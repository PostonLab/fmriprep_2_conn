%% ============================================================
%  batch_fmriprep_conn.m - WITH QC AND SCRUBBING (FIXED)
% ============================================================

clear batch;

% === User paths ===
fmriprep_dir = '/home/Neuron/data/derivatives/HIVPD23/derivatives_6subj';
project_dir  = '/home/Neuron/Documents/test_conn_6_subj';
participant_file = 'participant_ids.csv';  % Or specify full path
TR           = 2.0;

% Read subject IDs from CSV
if ~exist(participant_file, 'file')
    error('participant_ids.csv not found at: %s', participant_file);
end

% Read the CSV file
opts = detectImportOptions(participant_file);
participant_table = readtable(participant_file, opts);

% Extract subject IDs (assumes column is named 'participant_id' or 'subject_id')
if ismember('participant_id', participant_table.Properties.VariableNames)
    subjects = participant_table.participant_id;
elseif ismember('subject_id', participant_table.Properties.VariableNames)
    subjects = participant_table.subject_id;
else
    % If column name is different, use first column
    subjects = participant_table{:, 1};
end

% Convert to cell array of strings if needed
if ~iscell(subjects)
    subjects = cellstr(subjects);
end

% Remove any empty entries
subjects = subjects(~cellfun(@isempty, subjects));

fprintf('Loaded %d subjects from %s\n', length(subjects), participant_file);

% === QC Parameters ===
FD_threshold = 0.5;  % Framewise displacement threshold (mm)
DVARS_threshold = 1.5;  % DVARS threshold (standardized)

nsub = numel(subjects);

if ~exist(project_dir, 'dir')
    mkdir(project_dir);
end

batch.filename = fullfile(project_dir, 'conn_project.mat');
batch.Setup.isnew = 1;
batch.Setup.nsubjects = nsub;
batch.Setup.RT = TR;

%% ============================================================
%  Setup Anatomical ROIs
% ============================================================
batch.Setup.rois.names = {'Grey Matter', 'White Matter', 'CSF'};
batch.Setup.rois.files = cell(3, 1);

for i = 1:nsub
    subj_id = subjects{i};
    subj_dir = fullfile(fmriprep_dir, subj_id, 'anat');
    
    roi_paths = {
        fullfile(subj_dir, [subj_id '_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz']),
        fullfile(subj_dir, [subj_id '_space-MNI152NLin2009cAsym_label-WM_probseg.nii.gz']),
        fullfile(subj_dir, [subj_id '_space-MNI152NLin2009cAsym_label-CSF_probseg.nii.gz'])
    };
    
    for r = 1:3
        if exist(roi_paths{r}, 'file')
            batch.Setup.rois.files{r}{i} = roi_paths{r};
        else
            batch.Setup.rois.files{r}{i} = '';
        end
    end
end

%% ============================================================
%  Main Processing - WITH SCRUBBING
% ============================================================
fprintf('\n=== STARTING DATA COLLECTION ===\n');

% Define covariates: realignment + scrubbing
batch.Setup.covariates.names = {'realignment', 'scrubbing'};
batch.Setup.covariates.files = cell(2, 1);  % 2 covariates (each will contain nsub cells)
for cov = 1:2
    batch.Setup.covariates.files{cov} = cell(nsub, 1);
end
all_tasks = {};

% Storage for QC metrics
qc_summary = struct();

for i = 1:nsub
    subj_id = subjects{i};
    subj_dir = fullfile(fmriprep_dir, subj_id);
    fprintf('\n--- Subject: %s ---\n', subj_id);
    
    % T1w Anatomical
    anat_dir = fullfile(subj_dir, 'anat');
    anat_patterns = {
        [subj_id '_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz'],
        [subj_id '_*space-MNI*T1w.nii.gz']
    };
    
    anat_file = [];
    for p = 1:length(anat_patterns)
        anat_file = dir(fullfile(anat_dir, anat_patterns{p}));
        if ~isempty(anat_file)
            break;
        end
    end
    
    if isempty(anat_file)
        error('T1w not found for %s', subj_id);
    end
    
    batch.Setup.structurals{i} = fullfile(anat_file(1).folder, anat_file(1).name);
    fprintf('T1w: %s\n', anat_file(1).name);
    
    % Sessions
    session_dirs = dir(fullfile(subj_dir, 'ses-*'));
    if isempty(session_dirs)
        session_dirs = struct('name', {''}, 'folder', {subj_dir});
    end
    nses = numel(session_dirs);
    batch.Setup.nsessions(i) = nses;
    fprintf('Number of sessions: %d\n', nses);
    
    batch.Setup.functionals{i} = cell(nses, 1);
    batch.Setup.covariates.files{1}{i} = cell(nses, 1);  % realignment
    batch.Setup.covariates.files{2}{i} = cell(nses, 1);  % scrubbing
    
    for s = 1:nses
        if isempty(session_dirs(s).name)
            sess_label = '';
            sess_dir = fullfile(subj_dir, 'func');
        else
            sess_label = session_dirs(s).name;
            sess_dir = fullfile(subj_dir, sess_label, 'func');
        end
        
        fprintf('\n  Session %d: %s\n', s, sess_label);
        
        if ~exist(sess_dir, 'dir')
            fprintf('  WARNING: Directory does not exist\n');
            batch.Setup.functionals{i}{s} = {''};
            batch.Setup.covariates.files{1}{i}{s} = '';
            batch.Setup.covariates.files{2}{i}{s} = '';
            continue;
        end
        
        % Functional files
        if isempty(sess_label)
            func_pattern = [subj_id '_task-*space-MNI*bold.nii.gz'];
        else
            func_pattern = [subj_id '_' sess_label '_task-*space-MNI*bold.nii.gz'];
        end
        
        func_files = dir(fullfile(sess_dir, func_pattern));
        fprintf('  Found %d functional files\n', length(func_files));
        
        if ~isempty(func_files)
            n_runs = numel(func_files);
            func_paths = cell(n_runs, 1);
            
            % Storage for concatenated confounds
            session_motion_all = [];
            run_scrub_data = {};  % Store scrubbing data per run
            session_qc = struct('runs', [], 'mean_FD', [], 'mean_DVARS', [], ...
                              'n_scrubbed', [], 'pct_scrubbed', []);
            
            for r = 1:n_runs
                func_name = func_files(r).name;
                func_paths{r} = fullfile(func_files(r).folder, func_name);
                fprintf('    Run %d: %s\n', r, func_name);
                
                % Extract task name
                task_tokens = regexp(func_name, 'task-([^_]+)', 'tokens');
                task_name = 'rest';
                if ~isempty(task_tokens)
                    task_name = task_tokens{1}{1};
                    if ~ismember(task_name, all_tasks)
                        all_tasks{end+1} = task_name;
                    end
                end
                
                % Extract run number
                run_tokens = regexp(func_name, 'run-(\d+)', 'tokens');
                run_num = r;
                if ~isempty(run_tokens)
                    run_num = str2double(run_tokens{1}{1});
                end
                
                % Find confounds TSV
                if isempty(sess_label)
                    conf_pattern = sprintf('%s_task-%s_run-%02d*confounds*.tsv', subj_id, task_name, run_num);
                else
                    conf_pattern = sprintf('%s_%s_task-%s_run-%02d*confounds*.tsv', subj_id, sess_label, task_name, run_num);
                end
                
                conf_files = dir(fullfile(sess_dir, conf_pattern));
                
                if ~isempty(conf_files)
                    original_tsv = fullfile(conf_files(1).folder, conf_files(1).name);
                    fprintf('      Found confounds: %s\n', conf_files(1).name);
                    
                    % Extract motion parameters and scrubbing regressors
                    [motion_params, scrub_regs, qc_metrics] = extract_confounds_with_scrubbing(...
                        original_tsv, FD_threshold, DVARS_threshold);
                    
                    if ~isempty(motion_params)
                        session_motion_all = [session_motion_all; motion_params];
                        
                        % Store scrubbing data (will concatenate later with proper alignment)
                        run_scrub_data{r} = scrub_regs;
                        
                        % Store QC metrics
                        session_qc.runs(r) = run_num;
                        session_qc.mean_FD(r) = qc_metrics.mean_FD;
                        session_qc.mean_DVARS(r) = qc_metrics.mean_DVARS;
                        session_qc.n_scrubbed(r) = qc_metrics.n_scrubbed;
                        session_qc.pct_scrubbed(r) = qc_metrics.pct_scrubbed;
                        
                        fprintf('      QC: FD=%.3f mm, DVARS=%.3f, Scrubbed=%d (%.1f%%)\n', ...
                            qc_metrics.mean_FD, qc_metrics.mean_DVARS, ...
                            qc_metrics.n_scrubbed, qc_metrics.pct_scrubbed);
                    end
                else
                    fprintf('      WARNING: No confounds TSV found\n');
                    run_scrub_data{r} = [];
                end
            end
            
            % Save concatenated motion parameters
            if ~isempty(session_motion_all)
                session_motion_file = sprintf('rp_%s_session%02d.txt', subj_id, s);
                session_motion_path = fullfile(sess_dir, session_motion_file);
                dlmwrite(session_motion_path, session_motion_all, 'delimiter', '\t', 'precision', 8);
                batch.Setup.covariates.files{1}{i}{s} = session_motion_path;
                fprintf('  Created motion file: %s (%d timepoints)\n', ...
                    session_motion_file, size(session_motion_all, 1));
            else
                batch.Setup.covariates.files{1}{i}{s} = '';
            end
            
            % Concatenate scrubbing regressors with proper alignment
            if ~isempty(run_scrub_data) && any(~cellfun(@isempty, run_scrub_data))
                session_scrub_all = concatenate_scrubbing_regressors(run_scrub_data);
                
                session_scrub_file = sprintf('scrub_%s_session%02d.txt', subj_id, s);
                session_scrub_path = fullfile(sess_dir, session_scrub_file);
                dlmwrite(session_scrub_path, session_scrub_all, 'delimiter', '\t', 'precision', 0);
                batch.Setup.covariates.files{2}{i}{s} = session_scrub_path;
                fprintf('  Created scrubbing file: %s (%d timepoints, %d regressors)\n', ...
                    session_scrub_file, size(session_scrub_all, 1), size(session_scrub_all, 2));
            else
                batch.Setup.covariates.files{2}{i}{s} = '';
            end
            
            % Store QC summary
            qc_summary.(sprintf('sub%d_ses%d', i, s)) = session_qc;
            
            batch.Setup.functionals{i}{s} = func_paths;
            
        else
            fprintf('  No functional files found\n');
            batch.Setup.functionals{i}{s} = {''};
            batch.Setup.covariates.files{1}{i}{s} = '';
            batch.Setup.covariates.files{2}{i}{s} = '';
        end
    end
end

%% ============================================================
%  Save QC Summary
% ============================================================
qc_summary_file = fullfile(project_dir, 'qc_summary.mat');
save(qc_summary_file, 'qc_summary');
fprintf('\nQC summary saved to: %s\n', qc_summary_file);

% Print overall QC summary
fprintf('\n=== QC SUMMARY ===\n');
all_fds = [];
all_dvars = [];
all_pct_scrubbed = [];
for fn = fieldnames(qc_summary)'
    sess_qc = qc_summary.(fn{1});
    all_fds = [all_fds, sess_qc.mean_FD];
    all_dvars = [all_dvars, sess_qc.mean_DVARS];
    all_pct_scrubbed = [all_pct_scrubbed, sess_qc.pct_scrubbed];
end
fprintf('Mean FD across all runs: %.3f ± %.3f mm\n', mean(all_fds), std(all_fds));
fprintf('Mean DVARS across all runs: %.3f ± %.3f\n', mean(all_dvars), std(all_dvars));
fprintf('Mean %% scrubbed: %.1f%% ± %.1f%%\n', mean(all_pct_scrubbed), std(all_pct_scrubbed));

%% ============================================================
%  Setup Conditions
% ============================================================
fprintf('\n=== SETTING UP CONDITIONS ===\n');
fprintf('Tasks found: %s\n', strjoin(all_tasks, ', '));

if isempty(all_tasks)
    all_tasks = {'rest'};
end

batch.Setup.conditions.names = all_tasks;

for i = 1:nsub
    for s = 1:batch.Setup.nsessions(i)
        for c = 1:length(all_tasks)
            batch.Setup.conditions.onsets{c}{i}{s} = 0;
            batch.Setup.conditions.durations{c}{i}{s} = inf;
        end
    end
end

batch.Setup.conditions.missingdata = 1;

%% ============================================================
%  SETUP COMPLETION
% ============================================================
batch.Setup.done = 1;
batch.Setup.overwrite = 'Yes';

%% ============================================================
%  DENOISING SETUP
% ============================================================
batch.Denoising.filter = [0.01, 0.1];  % Band-pass filter
batch.Denoising.detrending = 1;         % Linear detrending
batch.Denoising.despiking = 0;          % Already done in fMRIPrep
batch.Denoising.regbp = 1;              % Bandpass before regression

% Confound regressors
batch.Denoising.confounds.names = {'realignment', 'scrubbing', 'White Matter', 'CSF'};
batch.Denoising.confounds.deriv = {0, 0, 0, 0};  % No derivatives

batch.Denoising.done = 1;
batch.Denoising.overwrite = 'Yes';

%% ============================================================
%  Save and Run
% ============================================================
fprintf('\n=== SAVING BATCH ===\n');
save(fullfile(project_dir, 'batch_prerun.mat'), 'batch');

fprintf('\n=== RUNNING CONN BATCH ===\n');
try
    conn_batch(batch);
    fprintf('\n✓ SUCCESS! Project created: %s\n', batch.filename);
    fprintf('✓ QC metrics saved to: %s\n', qc_summary_file);
catch ME
    fprintf('\n✗ FAILED: %s\n', ME.message);
    fprintf('\nFull error:\n');
    disp(getReport(ME, 'extended'));
    
    save(fullfile(project_dir, 'debug_batch_error.mat'), 'batch', 'ME');
    fprintf('\nError info saved\n');
end

%% ============================================================
%  HELPER FUNCTION - CONCATENATE SCRUBBING REGRESSORS
% ============================================================

function concatenated = concatenate_scrubbing_regressors(run_scrub_data)
    % Concatenate scrubbing regressors from multiple runs
    % Each run may have different number of outliers, so we need to
    % align them properly in time
    
    n_runs = length(run_scrub_data);
    
    % Get dimensions
    run_lengths = zeros(n_runs, 1);
    run_n_outliers = zeros(n_runs, 1);
    
    for r = 1:n_runs
        if ~isempty(run_scrub_data{r})
            run_lengths(r) = size(run_scrub_data{r}, 1);
            run_n_outliers(r) = size(run_scrub_data{r}, 2);
        end
    end
    
    total_timepoints = sum(run_lengths);
    total_outliers = sum(run_n_outliers);
    
    % Create empty matrix
    concatenated = zeros(total_timepoints, total_outliers);
    
    % Fill in scrubbing regressors, properly aligned
    current_timepoint = 1;
    current_outlier = 1;
    
    for r = 1:n_runs
        if ~isempty(run_scrub_data{r})
            n_tp = run_lengths(r);
            n_out = run_n_outliers(r);
            
            % Place this run's scrubbing regressors in the correct position
            concatenated(current_timepoint:(current_timepoint + n_tp - 1), ...
                        current_outlier:(current_outlier + n_out - 1)) = run_scrub_data{r};
            
            current_timepoint = current_timepoint + n_tp;
            current_outlier = current_outlier + n_out;
        end
    end
    
    % Remove any all-zero columns (shouldn't happen, but just in case)
    concatenated = concatenated(:, any(concatenated, 1));
    
    % If no scrubbing regressors at all, create a dummy column
    if isempty(concatenated) || size(concatenated, 2) == 0
        concatenated = zeros(total_timepoints, 1);
    end
end

%% ============================================================
%  HELPER FUNCTION - EXTRACT CONFOUNDS WITH SCRUBBING
% ============================================================

function [motion_params, scrub_regressors, qc_metrics] = extract_confounds_with_scrubbing(tsv_path, fd_thresh, dvars_thresh)
    % Initialize outputs
    motion_params = [];
    scrub_regressors = [];
    qc_metrics = struct('mean_FD', NaN, 'mean_DVARS', NaN, 'n_scrubbed', 0, 'pct_scrubbed', 0);
    
    try
        % Read confounds TSV
        opts = detectImportOptions(tsv_path, 'FileType', 'text', 'Delimiter', '\t');
        opts.VariableNamingRule = 'preserve';
        data_table = readtable(tsv_path, opts);
        
        n_timepoints = height(data_table);
        
        %% Extract 6 motion parameters
        motion_param_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
        motion_params = zeros(n_timepoints, 6);
        
        for p = 1:6
            param = motion_param_names{p};
            if ismember(param, data_table.Properties.VariableNames)
                col_data = data_table.(param);
                col_numeric = convert_to_numeric(col_data);
                col_numeric(isnan(col_numeric) | isinf(col_numeric)) = 0;
                motion_params(:, p) = col_numeric;
            end
        end
        
        %% Extract or compute FD and DVARS
        if ismember('framewise_displacement', data_table.Properties.VariableNames)
            FD = convert_to_numeric(data_table.framewise_displacement);
            FD(isnan(FD)) = 0;
        else
            FD = compute_FD(motion_params);
        end
        
        if ismember('dvars', data_table.Properties.VariableNames)
            DVARS = convert_to_numeric(data_table.dvars);
            DVARS(isnan(DVARS)) = 0;
        else
            DVARS = zeros(n_timepoints, 1);
        end
        
        %% Extract motion outliers from fMRIPrep
        outlier_cols = find(startsWith(data_table.Properties.VariableNames, 'motion_outlier'));
        
        if ~isempty(outlier_cols)
            % fMRIPrep already identified outliers
            scrub_regressors = zeros(n_timepoints, length(outlier_cols));
            for c = 1:length(outlier_cols)
                col_data = data_table{:, outlier_cols(c)};
                scrub_regressors(:, c) = convert_to_numeric(col_data);
            end
            n_scrubbed = length(outlier_cols);  % Number of outlier timepoints
            fprintf('        Using %d fMRIPrep motion_outlier regressors\n', length(outlier_cols));
        else
            % Create scrubbing regressors based on FD and DVARS thresholds
            scrub_idx = (FD > fd_thresh) | (DVARS > dvars_thresh);
            n_scrubbed = sum(scrub_idx);
            
            if n_scrubbed > 0
                scrub_regressors = zeros(n_timepoints, n_scrubbed);
                outlier_timepoints = find(scrub_idx);
                for i = 1:n_scrubbed
                    scrub_regressors(outlier_timepoints(i), i) = 1;
                end
                fprintf('        Created %d scrubbing regressors (FD>%.2f or DVARS>%.2f)\n', ...
                    n_scrubbed, fd_thresh, dvars_thresh);
            else
                scrub_regressors = zeros(n_timepoints, 1);  % Dummy column
                fprintf('        No timepoints exceed thresholds\n');
            end
        end
        
        %% Compute QC metrics
        qc_metrics.mean_FD = mean(FD(FD > 0));
        qc_metrics.mean_DVARS = mean(DVARS(DVARS > 0));
        qc_metrics.n_scrubbed = n_scrubbed;
        qc_metrics.pct_scrubbed = 100 * n_scrubbed / n_timepoints;
        
    catch ME
        warning('Failed to extract confounds from %s: %s', tsv_path, ME.message);
    end
end

%% ============================================================
%  HELPER FUNCTION - CONVERT TO NUMERIC
% ============================================================

function numeric_data = convert_to_numeric(col_data)
    if iscell(col_data)
        numeric_data = zeros(length(col_data), 1);
        for i = 1:length(col_data)
            if ischar(col_data{i}) || isstring(col_data{i})
                numeric_data(i) = str2double(col_data{i});
            else
                numeric_data(i) = double(col_data{i});
            end
        end
    else
        numeric_data = double(col_data);
    end
end

%% ============================================================
%  HELPER FUNCTION - COMPUTE FD
% ============================================================

function FD = compute_FD(motion_params)
    % Compute framewise displacement (Power et al. 2012)
    
    n_timepoints = size(motion_params, 1);
    FD = zeros(n_timepoints, 1);
    
    if n_timepoints > 1
        trans_diff = abs(diff(motion_params(:, 1:3), 1, 1));
        rot_diff = abs(diff(motion_params(:, 4:6), 1, 1)) * 50;
        FD(2:end) = sum([trans_diff, rot_diff], 2);
    end
end
