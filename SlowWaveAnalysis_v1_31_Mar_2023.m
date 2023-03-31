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


%% Identification of slow waves and extraction of epochs
% Define the threshold value
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

% Remove succeeding elements if value difference is less than 25 with 
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

% Remove succeeding elements if value difference is less than 25 with 
% preceding value in first and second zero crossings vectors separately.
for i = length(first_zero):-1:2
    if (first_zero(i) - first_zero(i-1)) < 25
        first_zero(i) = [];
    end
end

% Find large differences between consecutive elements
large_diff = diff(second_zero) > 50;
% Keep values at large difference indices and last element
second_zero = [second_zero(large_diff) second_zero(end)];

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
% Find large differences between consecutive elements
large_diff = diff(last_zero) > 50;
% Keep values at large difference indices and last element
last_zero = [last_zero(large_diff) last_zero(end)];

% convert time frames to time in seconds
fs = 250;
first_zero = first_zero/fs;
second_zero = second_zero/fs;
last_zero = last_zero/fs;
% reject based on the duration of negative peak - if the duration is less
% than 300 ms, more than 1000 ms
ngpk_dur = second_zero - first_zero;
% Find the indices of elements that meet the condition
indices = find(ngpk_dur < 0.3 | ngpk_dur > 1.0);
% Remove the elements from each vector
first_zero(indices) = [];
second_zero(indices) = [];
last_zero(indices) = [];
% create event array for EEGLAB with type, latency, and duration
i = 1:length(first_zero);
Total_dur = last_zero - first_zero;
first_zero = [i; first_zero; Total_dur];
first_zero = first_zero';
%dlmwrite('first_zero.txt', first_zero, 'delimiter', '\t');
EEG = pop_importevent( EEG, 'event','C:\\Toolbox\\eeglab2023.0\\first_zero.txt', ...
    'fields',{'type','latency','duration'},'timeunit',1,'optimalign','off');
% Extract EEG epochs of slow waves of 2 s duration (-0.5 s to 1.5 s)
EEG = pop_epoch( EEG, {  '1'  '2'  '3'  '4'  }, [-0.5         1.5], 'newname', ...
    'Data_Sleep4 epochs', 'epochinfo', 'yes');
EEG = pop_rmbase( EEG, [-500 1496] ,[]);
EEG = eeg_checkset( EEG );

% Find and reject the epochs with amplitude value more than +75
rej_epoch = find(max(max(EEG.data,[],1),[],2) > 75);
for i = 1:length(rej_epoch)
    epo_len(rej_epoch(i)) = 1;
end
EEG = pop_rejepoch( EEG, epo_len,0 );
EEG = pop_saveset( EEG, 'savemode','resave');

     