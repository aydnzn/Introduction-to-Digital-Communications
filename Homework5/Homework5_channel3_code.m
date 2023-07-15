%%
% Aydın Uzun 2015401210
%HW_5
% channel 3
%% define variables
clear
N  = 10^5;         % number of bits
Eb_N0_dB = [0:1:20]; %  Eb/N0 values

K = 7; % constant to arrange MMSE taps

ref = [0 0 ;  0 1; 1 0  ; 1 1 ]; % states
ipLUT = [ 0   0   0   0;...
    0   0   0   0;...
    1   1   0   0;...
    0   0   1   1 ]; % input given current state and reference state
%%
for nEN = 1:length(Eb_N0_dB)
    
    %% Viterbi Part
    
    % Transmitter
    bits = rand(1,N)>0.5; % generating 0,1 with equal probability
    binary_data_symbol = 2*bits-1; % BPSK modulation 0  -1; 1 -> 1
    
    % convolutional coding
    cip1 = mod(conv(bits,[1 1 1 ]),2);
    cip2 = mod(conv(bits,[1 0 1 ]),2);
    cip = [cip1;cip2];
    cip = cip(:).';
    
    binary_data_symbol_viterbi = 2*cip-1; % BPSK modulation 0 -> -1; 1 -> 1
    
    n_viterbi = 1/sqrt(2)*[randn(size(cip)) + j*randn(size(cip))]; % create wgn
    
    % add wgn
    out_viterbi = binary_data_symbol_viterbi + 10^(-Eb_N0_dB(nEN)/20)*n_viterbi; % additive white gaussian noise
    % receiver - hard decision decoding
    cipHat = real(out_viterbi)>0;
    % Viterbi decoding
    pathMetric  = zeros(4,1);  % path metric
    survivorPath_v  = zeros(4,length(out_viterbi)/2); % survivor path
    
    for jj = 1:length(out_viterbi)/2
        r = cipHat(2*jj-1:2*jj); % taking 2 coded bits
        
        % computing the Hamming distance between ip coded sequence with [00;01;10;11]
        rv = kron(ones(4,1),r);
        hammingDist = sum(xor(rv,ref),2);
        if (jj == 1) || (jj == 2)
            
            % branch metric and path metric for state 0
            bm1 = pathMetric(1,1) + hammingDist(1);
            pathMetric_n(1,1)  = bm1;
            survivorPath(1,1)  = 1;
            % branch metric and path metric for state 1
            bm1 = pathMetric(3,1) + hammingDist(3);
            pathMetric_n(2,1) = bm1;
            survivorPath(2,1)  = 3;
            
            
            % branch metric and path metric for state 2
            bm1 = pathMetric(1,1) + hammingDist(4);
            pathMetric_n(3,1) = bm1;
            survivorPath(3,1)  = 1;
            
            % branch metric and path metric for state 3
            bm1 = pathMetric(3,1) + hammingDist(2);
            pathMetric_n(4,1) = bm1;
            survivorPath(4,1)  = 3;
            
        else
            % branch metric and path metric for state 0
            bm1 = pathMetric(1,1) + hammingDist(1);
            bm2 = pathMetric(2,1) + hammingDist(4);
            [pathMetric_n(1,1) idx] = min([bm1,bm2]);
            survivorPath(1,1)  = idx;
            
            % branch metric and path metric for state 1
            bm1 = pathMetric(3,1) + hammingDist(3);
            bm2 = pathMetric(4,1) + hammingDist(2);
            [pathMetric_n(2,1) idx] = min([bm1,bm2]);
            survivorPath(2,1)  = idx+2;
            
            % branch metric and path metric for state 2
            bm1 = pathMetric(1,1) + hammingDist(4);
            bm2 = pathMetric(2,1) + hammingDist(1);
            [pathMetric_n(3,1) idx] = min([bm1,bm2]);
            survivorPath(3,1)  = idx;
            
            % branch metric and path metric for state 3
            bm1 = pathMetric(3,1) + hammingDist(2);
            bm2 = pathMetric(4,1) + hammingDist(3);
            [pathMetric_n(4,1) idx] = min([bm1,bm2]);
            survivorPath(4,1)  = idx+2;
            
        end
        
        pathMetric = pathMetric_n;
        survivorPath_v(:,jj) = survivorPath;
        
    end
    % trace back unit
    currState = 1;
    ipHat_v = zeros(1,length(out_viterbi)/2);
    for kk = length(out_viterbi)/2:-1:1
        prevState = survivorPath_v(currState,kk);
        ipHat_v(kk) = ipLUT(currState,prevState);
        currState = prevState;
    end
    
    % counting the errors
    nErrViterbi(nEN) = size(find([bits- ipHat_v(1:N)]),2);
    
    %% 15-Tap-MMSE
    % Channel model
 ht = [0.802 0.536 0.265]; 
 L  = length(ht);
    
    channel_out = conv(binary_data_symbol,ht);
    n = 1/sqrt(2)*[randn(1,N+length(ht)-1) + j*randn(1,N+length(ht)-1)]; % create wgn
    
    % add wgn
    out_mmse = channel_out + 10^(-Eb_N0_dB(nEN)/20)*n; % additive white gaussian noise
    % mmse equalization
    hAutoCorr = conv(ht,fliplr(ht));
    hM = toeplitz([hAutoCorr([3:end]) zeros(1,2*K+1-L)], [ hAutoCorr([3:end]) zeros(1,2*K+1-L) ]);
    hM = hM + 1/2*10^(-Eb_N0_dB(nEN)/10)*eye(2*K+1);
    d  = zeros(1,2*K+1);
    d([-1:1]+K+1) = fliplr(ht); 
    c_mmse  = [inv(hM)*d.'].';
    yFilt_mmse = conv(out_mmse,c_mmse);
    yFilt_mmse = yFilt_mmse(K+2:end);
    yFilt_mmse = conv(yFilt_mmse,ones(1,1)); % convolution
    ySamp_mmse = yFilt_mmse(1:1:N);  % sampling 
    
    % receiver - hard decision decoding
    ipHat_mmse = real(ySamp_mmse)>0;
    
    % counting the errors
    nErrMMSE(1,nEN) = size(find([bits- ipHat_mmse]),2);
    
    
%%    
end
%%
simBer_mmse = nErrMMSE/N; % simulated ber - MMSE
simBer_Viterbi = nErrViterbi/N; % simulated ber - Viterbi 
%%
% plots
close all
figure
semilogy(Eb_N0_dB,simBer_mmse(1,:),'-o','Linewidth',1);
hold on
semilogy(Eb_N0_dB,simBer_Viterbi,'-x','LineWidth',1);
axis square;
grid on
set(gca,'FontSize',14);
legend( '15-Tap-MMSE','Viterbi');
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (Channel 3) ');
%% References
% my reference to : Krishna Pillai, http://www.dsplog.com
% http://www.dsplog.com/2009/01/04/viterbi/
% my reference to : Krishna Sankar, http://www.dsplog.com
% http://www.dsplog.com/2010/01/24/ber-bpsk-isi-channel-mmse-equalization/
%%

