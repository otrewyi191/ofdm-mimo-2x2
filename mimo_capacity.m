%MIMO信道容量的Matlab仿真程序

%% Simulation of Capacity of MIMO channel
% Author: alpswy
% Email: alpswy@gmail.com
% Revise History:
%   2009-08-22 Created

%% Parameters for normal simulation
ChanType = 'Rayleigh'; % AWGN, Rayleigh, CorrRayleigh, Rice
LOOP = 10000; % number of channel tap samples
const = 1/sqrt(2); % constant number
NT = 4; % number of receive antennas
NR = 4; % number of transmit antennas
SNRdB = 10; % SNR in dB
SNRdec = 10^(0.1*SNRdB); % SNR in decibel

%% Parameters for displaying
numGroup = 100; % variable used to count samples
Range_MIN = 0; % minimum value of figure
Range_MAX = 6; % maximum value of figure

%% Parameters for simulation of correlation of channel
rho_t = 0;
rho_r = 0;

%% Parameters for simulation of Rice channel
K = 1; % Rice factor
D = sqrt(K);

switch ChanType
    case 'AWGN'
        H = ones(NR, NT);
        C = log2(det(eye(NR) + SNRdec/NT*H*H'));
        C_1 = [C, C];
        plot(C_1, 0:1);
        axis([Range_MIN, Range_MAX, 0, 1]);
    case 'Rayleigh'
        C = zeros(1, LOOP);
        C_2 = zeros(1, numGroup);
        for idx_LOOP = 1:LOOP
            H = const*(randn(NR, NT) + 1j*randn(NR, NT));
            C(idx_LOOP) = log2(det(eye(NR) + SNRdec/NT*H*H'));
        end
        C_min = min(C);
        C_max = max(C);
        C_1 = hist(C, numGroup);
        C_2(1) = C_1(1);
        for idx = 2:numGroup
            C_2(idx) = C_2(idx-1) + C_1(idx);
        end
        C_2 = C_2 / LOOP;
        Range_MIN = C_min;
        Range_MAX = C_max;
        X_axe = Range_MIN + [0:(numGroup-1)]/numGroup*(Range_MAX-Range_MIN);
        plot(X_axe, C_2);
    case 'CorrRayleigh'
        C = zeros(1, LOOP);
        C_2 = zeros(1, numGroup);
        R_tx = [1, rho_t; rho_t, 1];% for 2*2 only
        R_rx = [1, rho_r; rho_r, 1];
        for idx_LOOP = 1:LOOP
            H = const*(randn(NR, NT) + 1j*randn(NR, NT));
            H_Corr = R_rx * H * R_tx;
            C(idx_LOOP) = log2(det(eye(NR) + SNRdec/NT*H_Corr*H_Corr'));
        end
        C_min = min(C);
        C_max = max(C);
        C_1 = hist(C, numGroup);
        C_2(1) = C_1(1);
        for idx = 2:numGroup
            C_2(idx) = C_2(idx-1) + C_1(idx);
        end
        C_2 = C_2 / LOOP;
        Range_MIN = C_min;
        Range_MAX = C_max;
        X_axe = Range_MIN + [0:(numGroup-1)]/numGroup*(Range_MAX-Range_MIN);
        plot(X_axe, C_2);
    case 'Rice'
        C = zeros(1, LOOP);
        C_2 = zeros(1, numGroup);
        for idx_LOOP = 1:LOOP
            H = const*(randn(NR, NT) + 1j*randn(NR, NT)) + D * ones(NR, NT);
            C(idx_LOOP) = log2(det(eye(NR) + SNRdec/NT*H*H'));
        end
        C_min = min(C);
        C_max = max(C);
        C_1 = hist(C, numGroup);
        C_2(1) = C_1(1);
        for idx = 2:numGroup
            C_2(idx) = C_2(idx-1) + C_1(idx);
        end
        C_2 = C_2 / LOOP;
        Range_MIN = C_min;
        Range_MAX = C_max;
        X_axe = Range_MIN + [0:(numGroup-1)]/numGroup*(Range_MAX-Range_MIN);
        plot(X_axe, C_2);
    otherwise
        disp('Other Channel Type is still under construction. Only AWGN and Rayleigh fading is available.')
        keyboard
end