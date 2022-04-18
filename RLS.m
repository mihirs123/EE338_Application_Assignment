
%% Data Acquisition 

%Import Data 
data1 = load('Contaminated_Data.mat');   %Fz|Cz|O2|HEOG
data2 = load('HEOG.mat');
data3 = load('VEOG.mat');


%Data Extraction
eog_HEOG = data2.heog_3; %HEOG artifacts
eog_HEOG = eog_HEOG.';
eog_VEOG = data3.veog_3;
eog_VEOG = eog_VEOG.';
eeg_Fz = data1.sim3_con(1,:);   %Raw at Fz
eeg_Fz = eeg_Fz.';
size(eeg_Fz)

correct_fz = load('Pure_Data.mat');
%% RLS ALGORITHM EXECUTION FOR Denoising at Fz is a easy, stable and fast convergence method suitable for online removal of EOG artifacts 

% Execution for Fz electrode
% Variable Initialization
sample_no = size(eeg_Fz);    % No of samples/time points
order = 3;                   % Order of the Adaptive Filter (User Tunable)
sigma = 0.01;                % Initializing variable
lambda = 0.9999;             % Forgetting Factor for RLS Algorithm (User Tunable)
H = zeros(order,1);          % Initial filter Coefficients
R = sigma*eye(order,order);  % Initial value for Reference combination

correct_Fz = zeros(sample_no(1),1); 

correct_z = correct_fz.sim3_resampled(1,:);

% RLS Algorithm for Adaptive Denoising
for n = 1:sample_no          % Loop to simulate reality situation
   
    s = eeg_Fz(n,1);         % eeg @ Fz at that time point
    if n>=order               
            
        j = 1;               % calculation of last "order" reference vector
        for i = n:-1:(n+1-order)  
        
            r(j,1) = eog_HEOG(i,1) + eog_VEOG(i,1);
            j = j+1;
        end
        
        
        % Calculation of Filter Coefficients
        K = ((inv(R))*r)*(inv(lambda+r'*(inv(R))*r)); 
        e = s - (r'*H);
        H = H + (K*e);
        R = inv((inv(lambda)*inv(R)) - ((inv(lambda))*(K)*((r)')*(inv(R))));
        % Calculation Of RLS algorithm estimated signal
        correct_Fz(n,1) = s - (r'*H);
       

        
    end
end




% Fz Denoising plot

              %HEOG plot
figure(1)
subplot(3,2,1)
plot(eog_HEOG)
title('Horizontal EOG Artifact')
xlabel('sample number')
ylabel('HEOG mu volts')
ylim([-200 200])

subplot(3,2,2)
plot(eog_VEOG)
title('Vertical EOG Artifact')
xlabel('sample number')
ylabel('VEOG mu volts')
ylim([-200 200])


subplot(3,2,[3,4])

plot(eeg_Fz)
title('Raw EEG at Fz')
xlabel('sample number')
ylabel('EEG mu volts')
ylim([-200 200])

subplot(3,2,[5,6])
plot(correct_z)
title('Corrected EEG at Fz')
xlabel('sample number')
ylabel('EEG mu volts')
ylim([-200 200])



