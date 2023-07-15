%AYDIN UZUN
%2015401210
%EE 477 HW#2
%Please install Communications Toolbox to run this code. Because this code
%has some toolbox specific functions and classes.
%%
clear all
warning off
%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%
Infty=50;
number_of_bits=4000;
N=(number_of_bits);
max_nframe =2000;
ferlim = 100;
snr_db=0:2:20;
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
save BFSK_uniform_demo.mat sim_res

figure(1);
semilogy(snr_db , errs./nframes/number_of_bits, '-x'); %BER in Es/No(in dB)
xlabel('Es/No (in dB)');
ylabel('BER');

hold on;
%%
% calculation for gray mapping
% I expect to get the same result as for uniform mapping
% Because there is only 1 bit information, Binary 0 or binary 1
Infty=50;
number_of_bits=4000;
N=(number_of_bits);
max_nframe =2000;
ferlim = 100;
snr_db=0:2:20;
%%%%% SYMBOL CONSTELLATION %%%%%%%%%
tb4=[0 0;0 1; 1 0; 1 1];
tb8=[zeros(4,1) tb4; ones(4,1) tb4];
tb16=[zeros(8,1) tb8; ones(8,1) tb8];
%%%%%%%%%%%%%CORRELATION MATRICES %%%%%%%%%%%%%%
%%%%%%%%%%%%%%ERROR MATRICES%%%%%%%%%%%%%%%%%%%%
errs=zeros(length(snr_db), 1);
nframes=zeros(length(snr_db), 1);
ferrs=errs;
m_bit=0; n_bit=1;
M=2^(m_bit+n_bit);
bits=[0;1];
%%%%%%%%%%%%%%%%%%%%%%%% added
% Same Modulator from the previous homework
% The FSKModulator object modulates using the M-ary frequency shift keying method.
% The output is a baseband representation of the modulated signal.
% Binary Symbol mapping is used.
modulator_obj_bfsk = comm.FSKModulator(M,50,'SymbolMapping', 'Gray');
modulator_obj_bfsk.SamplesPerSymbol = 1;
modSignal = modulator_obj_bfsk([1;0]);
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
save BFSK_gray_demo.mat sim_res
semilogy(snr_db, errs./nframes/number_of_bits, '-x'); %BER in Es/No
xlabel('Es/No (in dB)');
ylabel('BER');
legend('BFSK BER versus Es/No UNIFORM mapping', 'BFSK BER versus Es/No GRAY mapping');

