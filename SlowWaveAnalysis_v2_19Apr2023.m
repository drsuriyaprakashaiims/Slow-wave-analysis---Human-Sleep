%% Rereferencing and filtering data
% Rereferencing to mastoids
EEG = pop_reref( EEG, [57 100] ,'keepref','on');

% Filter specifications
fs = 250; % Sampling frequency (Hz)
fpass = [0.5 4]; % Passband frequency range (Hz)
fstop = [0.1 10]; % Stopband frequency range (Hz)
apass = 3; % Maximum passband ripple (dB)
astop = 25; % Minimum stopband attenuation (dB)

% Design the filter using cheby2
[n, Wn] = cheb2ord(fpass/(fs/2), fstop/(fs/2), apass, astop);
[b, a] = cheby2(n, astop, Wn);

% Apply the filter to the signal
EEG.data = filter(b, a, EEG.data);
EEG = pop_saveset( EEG, 'savemode','resave');
EEG = eeg_checkset( EEG );

% Identification of slow waves
% Define the threshold value of at least 80 microvolts
y_filtered = -( EEG.data);
threshold = 80;

% Initialize the neg_peak vector
neg_peak = [];

% Loop through each row of the data matrix
for i = 1:size(y_filtered, 1)
    % Find the peaks in the current row
    [pks, locs] = findpeaks(y_filtered(i,:), 'MinPeakHeight', threshold);
    
    % Append the peak locations to the neg_peak vector
    neg_peak = [neg_peak; locs(:)];
end

% Keep only unique values in the neg_peak vector
neg_peak = unique(neg_peak);

% Remove succeeding elements if value difference is less than 25 (100 ms) with 
% preceding value in neg_peak vector.
for i = length(neg_peak):-1:2
    if (neg_peak(i) - neg_peak(i-1)) < 25
        neg_peak(i) = [];
    end
end

% Initialize first_zero and second_zero vectors
first_zero = [];
second_zero = [];

% Loop through each row of the data matrix
for i = 1:size(y_filtered, 1)
    % Loop through each peak location in the neg_peak vector
    for j = 1:length(neg_peak)
        % Find index of peak location in current row
        idx = neg_peak(j);
        
        % Find index of zero crossing immediately before peak location in current row
        k = idx - 1;
        while k > 0 && sign(y_filtered(i,k)) == sign(y_filtered(i,idx))
            k = k - 1;
        end
        
        % Append index of zero crossing to first_zero vector if it exists
        if k > 0 && sign(y_filtered(i,k)) ~= sign(y_filtered(i,idx))
            first_zero(end+1) = k;
        end
        
        % Find index of zero crossing immediately after peak location in current row
        k = idx + 1;
        while k <= size(y_filtered,2) && sign(y_filtered(i,k)) == sign(y_filtered(i,idx))
            k = k + 1;
        end
        
        % Append index of zero crossing to second_zero vector if it exists
        if k <= size(y_filtered,2) && sign(y_filtered(i,k)) ~= sign(y_filtered(i,idx))
            second_zero(end+1) = k;
        end
    end
end

% Keep only unique values in the first_zero and second_zero vectors
first_zero = unique(first_zero);
second_zero = unique(second_zero);

% Remove succeeding elements if value difference is less than 75(=300 ms) with 
% preceding value in first and second zero crossings vectors separately.
for i = length(first_zero):-1:2
    if (first_zero(i) - first_zero(i-1)) < 75
        first_zero(i) = [];
    end
end

% Check if the first_zero array is empty
if isempty(first_zero)
    % If it is, print a message to the console
    disp('The first_zero array is empty. No slow waves. Stopping script.')
    % Use the return statement to stop the script and return to the command line
    return
end

% extraction of second_zero of slow wave based on the number of elements in
% first zero
Z = linkage(second_zero', 'ward');
num_clusters = length(first_zero);
cluster_idx = cluster(Z, 'maxclust', num_clusters);
cluster_medians = zeros(1, num_clusters);
for i = 1:num_clusters
    cluster_medians(i) = median(second_zero(cluster_idx == i));
end
cluster_medians = sort(cluster_medians);
second_zero = cluster_medians;
second_zero = int32(second_zero);

% Initialize last_zero vector for zero crossings subsequent to after zero crossings.
last_zero = [];

% Loop through each row of the data matrix
for i = 1:size(y_filtered, 1)
    % Loop through each after zero crossing location in the second_zero vector
    for j = 1:length(second_zero)
        % Find index of after zero crossing location in current row
        idx = second_zero(j);
        
        % Find index of zero crossing immediately after after zero crossing 
        % location in current row
        k = idx + 1;
        while k <= size(y_filtered,2) && sign(y_filtered(i,k)) == sign(y_filtered(i,idx))
            k = k + 1;
        end
        
        % Append index of zero crossing to last_zero vector if it exists
        if k <= size(y_filtered,2) && sign(y_filtered(i,k)) ~= sign(y_filtered(i,idx))
            last_zero(end+1) = k;
        end
    end
end

% Keep only unique values in the last_zero vectors
last_zero = unique(last_zero);

% extraction of last_zero of slow wave based on the number of elements in
% first zero
Z = linkage(last_zero', 'ward');
num_clusters = length(first_zero);
cluster_idx = cluster(Z, 'maxclust', num_clusters);
cluster_medians = zeros(1, num_clusters);
for i = 1:num_clusters
    cluster_medians(i) = median(last_zero(cluster_idx == i));
end
cluster_medians = sort(cluster_medians);
last_zero = cluster_medians;
last_zero = int32(last_zero);

% initialize the result cell array
result = cell(1, length(first_zero));
% loop through each element in first_zero
for i = 1:length(first_zero)
% find the difference between the first_zero and second_zero (i.e. negative peak
% duration)
diff = second_zero - first_zero(i);

% negative peak duration criteria (300 ms to 1000 ms)
idx = find(diff > 75 & diff < 250);

% check if any indices were found
if isempty(idx)
% if no indices were found, store an empty cell
result{i} = {};
else
% if indices were found, store them in the result cell array
result{i} = idx;
end
end

result = cellfun(@(x) num2str(x), result, 'UniformOutput', false);
[~, idx] = unique(result, 'stable');
result(setdiff(1:numel(result), idx)) = {[]};

result = cellfun(@str2double, result);
first_zero(isnan(result)) = [];
result = result(~isnan(result));
second_zero = second_zero(result);
last_zero = last_zero(result);
final_result = [first_zero; second_zero; last_zero];
clear diff
final_result = final_result(:, all(diff(final_result) >= 0));
final_result = double(final_result);
% positive peak duration criteria (less than 1000 ms)
fiff = final_result(3,:) - final_result(2,:);
final_result(:,fiff > 250) = [];

% convert time frames to time in milliseconds
final_result = final_result* 4;

% create event array for EEGLAB with type, latency, and duration
[rows, columns] = size(final_result);
i = 1:columns;
total_dur = final_result(3,:) - final_result(1,:);
events = [i; final_result(1,:); total_dur];
events = events';
for i = 1:size(events,1)
    if mod(events(i,2),4) ~= 0
        events(i,2) = round(events(i,2)/4)*4;
    end
end

% Extraction of epochs
if isempty(events)
    disp('No slow waves detected')
    return
else
    EEG = pop_importevent( EEG, 'event',events, ...
    'fields',{'type','latency','duration'},'timeunit',1E-3,'optimalign','on');
    % Extract EEG epochs of slow waves of 2 s duration (-0.5 s to 1.5 s)
    event_type = num2cell(events(:,1));
    EEG = pop_epoch( EEG, event_type, [-0.5         1.5], 'newname', ...
    'Data_Sleep epochs', 'epochinfo', 'yes');
    EEG = pop_rmbase( EEG, [-500 1496] ,[]);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'savemode','resave');
    EEG = eeg_checkset( EEG );
end

% Find and reject the epochs with amplitude value more than +75 and
% less than -250 microvolts, abnormal trend >150 microvolts/second
pop_rejmenu(EEG,1)




     