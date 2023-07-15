clear all;
warning off;
%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%
Infty=50;
number_of_bits=4000;
N=(number_of_bits);
max_nframe =2000;
ferlim = 100;
snr_db=0:2:20;
%%%%%%%%%%%%%added
ES_N0_in_lin = 10.^(snr_db./10);
Eb_N0_in_dB = snr_db - 10*log10(2);
Eb_N0_in_lin = 10.^(Eb_N0_in_dB./10);
%%%%%%%%%%%%%added
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb4=[0 0;0 1; 1 0; 1 1];
tb8=[zeros(4,1) tb4; ones(4,1) tb4];
tb16=[zeros(8,1) tb8; ones(8,1) tb8];
%%%%%%%%%%%%%CORRELATION MATRICES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%ERROR MATRICES%%%%%%%%%%%%%%%%%%%%
errs=zeros(length(snr_db), 1);
%%%%%%%%%%%%%%%%%%%%%%%% added
% initialize symbol error vector
symbol_errs=zeros(length(snr_db), 1);
%%%%%%%%%%%%%%%%%%%%%%%% added
nframes=zeros(length(snr_db), 1);
ferrs=errs;
%%%%%%%%%%%%%%%%%%%%%%%% added
% initialize symbol error vector
ferrs_symbol = symbol_errs;
%%%%%%%%%%%%%%%%%%%%%%%% added
m_bit=1; n_bit=1;
M=2^(m_bit+n_bit);
bits=tb4;
code_symbols=exp(1j*2*pi*[0:(M-1)]/M);
tot_bits=m_bit+n_bit;
Nsy=N/tot_bits;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    snr_p=snr_db(nEN);
    en = 10^(snr_p/10); % convert SNR from unit db to normal numbers
    sigma = 1/sqrt(en); % standard deviation of AWGN noise
    nframe = 0;
    while (nframe<max_nframe) && (ferrs(nEN)<ferlim)
        err_count=0;
        nframe = nframe + 1;
        info_bits=round(rand(1,number_of_bits));
        info_part=reshape(info_bits, tot_bits, Nsy);
        info_matrix=info_part';
        sym_vec=ones(Nsy, 1);
        for v=1:tot_bits
            sym_vec=sym_vec+info_matrix(:,v).*2^(m_bit+n_bit-v);
        end
        sym_seq=code_symbols(sym_vec);
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        noise=1/sqrt(2)*[randn(1, Nsy) + 1j*randn(1,Nsy)];
        det_seq=zeros(1,N);
        rec_sig=sym_seq+sigma*noise;
        %%%%DETECTOR %%%%%%%%%%%%
        for k=1:Nsy
            min_metric=10^6; dm=zeros(1,M);
            for r=1:M
                x_r=code_symbols(r);
                dm(r)=norm(rec_sig(k)-x_r);
            end
            [rowmin, sym_ind]=min(dm);
            detected_bits=bits(sym_ind, :);
            err = length(find(info_part(:,k)~=detected_bits'));
            errs(nEN)=errs(nEN)+err;
            err_count=err_count+err;
        end
        if err_count~=0
            ferrs(nEN)=ferrs(nEN)+1;
        end
      
    end % End of while loop
    nframes(nEN)=nframe;
    sim_res=[errs nframes]
    %%%%%%%%%%%%%%%%%%%%%%%% added
    % show error vector step by step
    %%%%%%%%%%%%%%%%%%%%%%%% added
end %end for (SNR points)
sim_res=[errs nframes]
%%%%%%%%%%%%%%%%%%%%%%%% added
%save error vector
%%%%%%%%%%%%%%%%%%%%%%%% added
figure(1);
semilogy(snr_db, errs./nframes/number_of_bits, '-x'); %BER in Es/No
%%%%%%%%%%%%%%%%%%%%%%%% added
hold on;
grid on;

%%

%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%
Infty=50;
number_of_bits=4000;
N=(number_of_bits);
max_nframe =2000;
ferlim = 100;
snr_db=0:2:20;
%%%%%%%%%%%%%added
ES_N0_in_lin = 10.^(snr_db./10);
Eb_N0_in_dB = snr_db - 10*log10(2);
Eb_N0_in_lin = 10.^(Eb_N0_in_dB./10);
%%%%%%%%%%%%%added
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb4=[0 0;0 1; 1 0; 1 1];
tb8=[zeros(4,1) tb4; ones(4,1) tb4];
tb16=[zeros(8,1) tb8; ones(8,1) tb8];
%%%%%%%%%%%%%CORRELATION MATRICES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%ERROR MATRICES%%%%%%%%%%%%%%%%%%%%
errs=zeros(length(snr_db), 1);
nframes=zeros(length(snr_db), 1);
ferrs=errs;
m_bit=0; n_bit=1; % change
M=2^(m_bit+n_bit);
bits=[0;1]; % change
code_symbols=exp(1j*2*pi*[0:(M-1)]/M);
tot_bits=m_bit+n_bit;
Nsy=N/tot_bits;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    snr_p=snr_db(nEN);
    en = 10^(snr_p/10); % convert SNR from unit db to normal numbers
    sigma = 1/sqrt(en); % standard deviation of AWGN noise
    nframe = 0;
    while (nframe<max_nframe) && (ferrs(nEN)<ferlim)
        err_count=0;
        nframe = nframe + 1;
        info_bits=round(rand(1,number_of_bits));
        info_part=reshape(info_bits, tot_bits, Nsy);
        info_matrix=info_part';
        sym_vec=ones(Nsy, 1);
        for v=1:tot_bits
            sym_vec=sym_vec+info_matrix(:,v).*2^(m_bit+n_bit-v);
        end
        sym_seq=code_symbols(sym_vec);
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% added
                %  awgn : Add white Gaussian noise to a signal
        rec_sig = awgn(sym_seq,snr_db(nEN),'measured');
                % awgn measures the signal power before adding noise.
        %%%%%%%%%%%%%%%%%%%%%%%% added
        %%%%DETECTOR %%%%%%%%%%%%
        for k=1:Nsy
            min_metric=10^6; dm=zeros(1,M);
            for r=1:M
                x_r=code_symbols(r);
                dm(r)=norm(rec_sig(k)-x_r);
            end
            [rowmin, sym_ind]=min(dm);
            detected_bits=bits(sym_ind, :);
            err = length(find(info_part(:,k)~=detected_bits'));
            errs(nEN)=errs(nEN)+err;
            err_count=err_count+err;
        end
        if err_count~=0
            ferrs(nEN)=ferrs(nEN)+1;
        end
    end % End of while loop
    nframes(nEN)=nframe;
    sim_res=[errs nframes]
end %end for (SNR points)
sim_res=[errs nframes]

semilogy(snr_db, errs./nframes/number_of_bits, '-x'); %BER in Es/No


%%
%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%
Infty=50;
number_of_bits=4000;
N=(number_of_bits);
max_nframe =2000;
ferlim = 100;
snr_db=0:2:20;
%%%%%%%%%%%%%added
Es_N0_in_lin = 10.^(snr_db./10);
%%%%%%%%%%%%%added
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb4=[0 0;0 1; 1 0; 1 1];
tb8=[zeros(4,1) tb4; ones(4,1) tb4];
tb16=[zeros(8,1) tb8; ones(8,1) tb8];
%%%%%%%%%%%%%CORRELATION MATRICES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%ERROR MATRICES%%%%%%%%%%%%%%%%%%%%
errs=zeros(length(snr_db), 1);
nframes=zeros(length(snr_db), 1);
ferrs=errs;
m_bit=1; n_bit=1;
M=2^(m_bit+n_bit);
bits=tb4; % change
%%%%%%%%%%%%%%%%%%%%%%%% added
% Same Modulator from the previous homework
%The PAMModulator object modulates using M-ary pulse amplitude modulation.
%The output is a baseband representation of the modulated signal.
%The M-ary number parameter, M, represents the number of points in the
%signal constellation and requires an even integer.
% binary symbopl mapping is used and the average power is normalized to
% unity.
modulator_obj_4pam = comm.PAMModulator(M,'SymbolMapping', 'Binary', 'NormalizationMethod','Average Power');
code_symbols = step(modulator_obj_4pam,[0;1;2;3]);
code_symbols = code_symbols';
%%%%%%%%%%%%%%%%%%%%%%%% added
tot_bits=m_bit+n_bit;
Nsy=N/tot_bits;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    snr_p=snr_db(nEN);
    en = 10^(snr_p/10); % convert SNR from unit db to normal numbers
    sigma = 1/sqrt(en); % standard deviation of AWGN noise
    nframe = 0;
    while (nframe<max_nframe) && (ferrs(nEN)<ferlim)
        err_count=0;
        nframe = nframe + 1;
        info_bits=round(rand(1,number_of_bits));
        info_part=reshape(info_bits, tot_bits, Nsy);
        info_matrix=info_part';
        sym_vec=ones(Nsy, 1);
        for v=1:tot_bits
            sym_vec=sym_vec+info_matrix(:,v).*2^(m_bit+n_bit-v);
        end
        sym_seq=code_symbols(sym_vec);
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% added
        %  awgn : Add white Gaussian noise to a signal
        rec_sig = awgn(sym_seq,snr_db(nEN),'measured');
        % awgn measures the signal power before adding noise.
        %%%%%%%%%%%%%%%%%%%%%%%% added
        %%%%DETECTOR %%%%%%%%%%%%
        for k=1:Nsy
            min_metric=10^6; dm=zeros(1,M);
            for r=1:M
                x_r=code_symbols(r);
                dm(r)=norm(rec_sig(k)-x_r);
            end
            [rowmin, sym_ind]=min(dm);
            detected_bits=bits(sym_ind, :);
            err = length(find(info_part(:,k)~=detected_bits'));
            errs(nEN)=errs(nEN)+err;
            err_count=err_count+err;
        end
        if err_count~=0
            ferrs(nEN)=ferrs(nEN)+1;
        end
    end % End of while loop
    nframes(nEN)=nframe;
    sim_res=[errs nframes]
end %end for (SNR points)
sim_res=[errs nframes]

semilogy(snr_db , errs./nframes/number_of_bits, '-x'); %BER in Es/No(in dB)

%%
Infty=50;
number_of_bits=4000;
N=(number_of_bits);
max_nframe =2000;
ferlim = 100;
snr_db=0:2:20;
Es_N0_in_lin = 10.^(snr_db./10);
Eb_N0_in_dB = snr_db - 10*log10(1);
Eb_N0_in_lin = 10.^(Eb_N0_in_dB./10);
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb4=[0 0;0 1; 1 0; 1 1];
tb8=[zeros(4,1) tb4; ones(4,1) tb4];
tb16=[zeros(8,1) tb8; ones(8,1) tb8];
%%%%%%%%%%%%%CORRELATION MATRICES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%ERROR MATRICES%%%%%%%%%%%%%%%%%%%%
errs=zeros(length(snr_db), 1);
nframes=zeros(length(snr_db), 1);
ferrs=errs;
m_bit=0; n_bit=1; % change
M=2^(m_bit+n_bit);
bits=[0;1]; %change
%%%%%%%%%%%%%%%%%%%%%%%% added
% Same Modulator from the previous homework
% The FSKModulator object modulates using the M-ary frequency shift keying method.
% The output is a baseband representation of the modulated signal.
% Binary Symbol mapping is used.
modulator_obj_bfsk = comm.FSKModulator(M,50,'SymbolMapping', 'binary');
modulator_obj_bfsk.SamplesPerSymbol = 1;
data = [1;0];
modSignal = modulator_obj_bfsk(data);
code_symbols = modSignal';
%%%%%%%%%%%%%%%%%%%%%%%% added
tot_bits=m_bit+n_bit;
Nsy=N/tot_bits;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    snr_p=snr_db(nEN);
    en = 10^(snr_p/10); % convert SNR from unit db to normal numbers
    sigma = 1/sqrt(en); % standard deviation of AWGN noise
    nframe = 0;
    while (nframe<max_nframe) && (ferrs(nEN)<ferlim)
        err_count=0;
        nframe = nframe + 1;
        info_bits=round(rand(1,number_of_bits));
        info_part=reshape(info_bits, tot_bits, Nsy);
        info_matrix=info_part';
        sym_vec=ones(Nsy, 1);
        for v=1:tot_bits
            sym_vec=sym_vec+info_matrix(:,v).*2^(m_bit+n_bit-v);
        end
        sym_seq=code_symbols(sym_vec);
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% added
        %  awgn : Add white Gaussian noise to a signal
        rec_sig = awgn(sym_seq,snr_db(nEN),'measured');
        % awgn measures the signal power before adding noise.
        %%%%%%%%%%%%%%%%%%%%%%%% added
        %%%%DETECTOR %%%%%%%%%%%%
        for k=1:Nsy
            min_metric=10^6; dm=zeros(1,M);
            for r=1:M
                x_r=code_symbols(r);
                dm(r)=norm(rec_sig(k)-x_r);
            end
            [rowmin, sym_ind]=min(dm);
            detected_bits=bits(sym_ind, :);
            err = length(find(info_part(:,k)~=detected_bits'));
            errs(nEN)=errs(nEN)+err;
            err_count=err_count+err;
        end
        if err_count~=0
            ferrs(nEN)=ferrs(nEN)+1;
        end
    end % End of while loop
    nframes(nEN)=nframe;
    sim_res=[errs nframes]
end %end for (SNR points)
sim_res=[errs nframes]

semilogy(snr_db , errs./nframes/number_of_bits, '-x'); %BER in Es/No(in dB)
%%

%AYDIN UZUN
%2015401210
%EE 477 HW#3
%Please install Communications Toolbox to run this code. Because this code
%has some toolbox specific functions and classes.
%%
%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%
Infty=50;
number_of_bits=4000;
N=(number_of_bits);
max_nframe =2000;
ferlim = 100;
snr_db=0:2:20;
Es_N0_in_lin = 10.^(snr_db./10);
Eb_N0_in_dB = snr_db - 10*log10(4);
Eb_N0_in_lin = 10.^(Eb_N0_in_dB./10);
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb4=[0 0;0 1; 1 0; 1 1];
tb8=[zeros(4,1) tb4; ones(4,1) tb4];
tb16=[zeros(8,1) tb8; ones(8,1) tb8];
%%%%%%%%%%%%%CORRELATION MATRICES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%ERROR MATRICES%%%%%%%%%%%%%%%%%%%%
errs=zeros(length(snr_db), 1);
nframes=zeros(length(snr_db), 1);
ferrs=errs;
m_bit=2; n_bit=2; % change
M=2^(m_bit+n_bit);
bits=tb16; %change
%%%%%%%%%%%%%%%%%%%%%%%% added
% Same Modulator from the first homework
modulator_obj_16qam = comm.RectangularQAMModulator(M,'SymbolMapping', 'binary', 'NormalizationMethod','Average Power');
baseband_mod_out_symbol_16qam = step(modulator_obj_16qam,[0:15]');
code_symbols = baseband_mod_out_symbol_16qam';
%%%%%%%%%%%%%%%%%%%%%%%% added
tot_bits=m_bit+n_bit;
Nsy=N/tot_bits;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    snr_p=snr_db(nEN);
    en = 10^(snr_p/10); % convert SNR from unit db to normal numbers
    sigma = 1/sqrt(en); % standard deviation of AWGN noise
    nframe = 0;
    while (nframe<max_nframe) && (ferrs(nEN)<ferlim)
        err_count=0;
        nframe = nframe + 1;
        info_bits=round(rand(1,number_of_bits));
        info_part=reshape(info_bits, tot_bits, Nsy);
        info_matrix=info_part';
        sym_vec=ones(Nsy, 1);
        for v=1:tot_bits
            sym_vec=sym_vec+info_matrix(:,v).*2^(m_bit+n_bit-v);
        end
        sym_seq=code_symbols(sym_vec);
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% added
        %  awgn : Add white Gaussian noise to a signal
        rec_sig = awgn(sym_seq,snr_db(nEN),'measured');
        % awgn measures the signal power before adding noise.
        %%%%%%%%%%%%%%%%%%%%%%%% added
        %%%%DETECTOR %%%%%%%%%%%%
        for k=1:Nsy
            min_metric=10^6; dm=zeros(1,M);
            for r=1:M
                x_r=code_symbols(r);
                dm(r)=norm(rec_sig(k)-x_r);
            end
            [rowmin, sym_ind]=min(dm);
            detected_bits=bits(sym_ind, :);
            err = length(find(info_part(:,k)~=detected_bits'));
            errs(nEN)=errs(nEN)+err;
            err_count=err_count+err;
        end
        if err_count~=0
            ferrs(nEN)=ferrs(nEN)+1;
        end
    end % End of while loop
    nframes(nEN)=nframe;
    sim_res=[errs nframes]
end %end for (SNR points)
sim_res=[errs nframes]

semilogy(snr_db , errs./nframes/number_of_bits, '-x'); %BER in Es/No(in dB)

%%
clear all
warning off
%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%
Infty=50;
number_of_bits=4002;
N=(number_of_bits);
max_nframe =2000;
ferlim = 100;
snr_db=0:2:20;
Es_N0_in_lin = 10.^(snr_db./10);
Eb_N0_in_dB = snr_db - 10*log10(4);
Eb_N0_in_lin = 10.^(Eb_N0_in_dB./10);
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb4=[0 0;0 1; 1 0; 1 1];
tb8=[zeros(4,1) tb4; ones(4,1) tb4];
tb16=[zeros(8,1) tb8; ones(8,1) tb8];
%%%%%%%%%%%%%CORRELATION MATRICES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%ERROR MATRICES%%%%%%%%%%%%%%%%%%%%
errs=zeros(length(snr_db), 1);
nframes=zeros(length(snr_db), 1);
ferrs=errs;
m_bit=1; n_bit=2; % change
M=2^(m_bit+n_bit);
bits=tb8; %change
%%%%%%%%%%%%%%%%%%%%%%%% added
% Same Modulator from the first homework
modulator_obj_8psk = comm.PSKModulator(8, 0,'SymbolMapping', 'binary');
baseband_mod_out_symbol_8psk = step(modulator_obj_8psk,[0;1;2;3;4;5;6;7]);
code_symbols = baseband_mod_out_symbol_8psk';
%%%%%%%%%%%%%%%%%%%%%%%% added
tot_bits=m_bit+n_bit;
Nsy=N/tot_bits;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    snr_p=snr_db(nEN);
    en = 10^(snr_p/10); % convert SNR from unit db to normal numbers
    sigma = 1/sqrt(en); % standard deviation of AWGN noise
    nframe = 0;
    while (nframe<max_nframe) && (ferrs(nEN)<ferlim)
        err_count=0;
        nframe = nframe + 1;
        info_bits=round(rand(1,number_of_bits));
        info_part=reshape(info_bits, tot_bits, Nsy);
        info_matrix=info_part';
        sym_vec=ones(Nsy, 1);
        for v=1:tot_bits
            sym_vec=sym_vec+info_matrix(:,v).*2^(m_bit+n_bit-v);
        end
        sym_seq=code_symbols(sym_vec);
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%% added
        %  awgn : Add white Gaussian noise to a signal
        rec_sig = awgn(sym_seq,snr_db(nEN),'measured');
        % awgn measures the signal power before adding noise.
        %%%%%%%%%%%%%%%%%%%%%%%% added
        %%%%DETECTOR %%%%%%%%%%%%
        for k=1:Nsy
            min_metric=10^6; dm=zeros(1,M);
            for r=1:M
                x_r=code_symbols(r);
                dm(r)=norm(rec_sig(k)-x_r);
            end
            [rowmin, sym_ind]=min(dm);
            detected_bits=bits(sym_ind, :);
            err = length(find(info_part(:,k)~=detected_bits'));
            errs(nEN)=errs(nEN)+err;
            err_count=err_count+err;
        end
        if err_count~=0
            ferrs(nEN)=ferrs(nEN)+1;
        end
    end % End of while loop
    nframes(nEN)=nframe;
    sim_res=[errs nframes]
end %end for (SNR points)
sim_res=[errs nframes]

semilogy(snr_db , errs./nframes/number_of_bits, '-x'); %BER in Es/No(in dB)


%%




zoom on;
xlabel('Es/No(in dB)');
ylabel('BER');
legend('QPSK', 'BPSK', '4PAM', 'BFSK','16QAM','8PSK');
%%%%%%%%%%%%%%%%%%%%%%%% added

