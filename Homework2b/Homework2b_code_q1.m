%AYDIN UZUN
%2015401210
%EE 477 HW#2b
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
        %%%%%%%%%%%%%%%%%%%%%%%% added
        % symbol error counter
        err_sym_count=0;
        %%%%%%%%%%%%%%%%%%%%%%%% added
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
            %%%%%%%%%%%%%%%%%%%%%%%% added
            % check if there's an symbol error at this instant
            if err~=0
                symbol_err=1 ;
            else
                symbol_err=0 ;
            end
            % create symbol error vectors
            % completely analogous to bit error rate calculation
            % for example between '01' and '10' there's only one symbol
            % error but there's two bit errors.
            symbol_errs(nEN) = symbol_errs(nEN) + symbol_err;
            err_sym_count = err_sym_count + symbol_err;
            %%%%%%%%%%%%%%%%%%%%%%%% added
            errs(nEN)=errs(nEN)+err;
            err_count=err_count+err;
        end
        if err_count~=0
            ferrs(nEN)=ferrs(nEN)+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%% added
        % completely analogous to bit error rate calculation
        if err_sym_count~=0
            ferrs_symbol(nEN)=ferrs_symbol(nEN)+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%% added
    end % End of while loop
    nframes(nEN)=nframe;
    sim_res=[errs nframes]
    %%%%%%%%%%%%%%%%%%%%%%%% added
    % show error vector step by step
    sim_res_sym = [symbol_errs nframes]
    %%%%%%%%%%%%%%%%%%%%%%%% added
end %end for (SNR points)
sim_res=[errs nframes]
save QPSK_demo.mat sim_res
%%%%%%%%%%%%%%%%%%%%%%%% added
%save error vector
sim_res_sym = [symbol_errs nframes]
save QPSK_demo_sym.mat sim_res_sym
%%%%%%%%%%%%%%%%%%%%%%%% added
figure(1);
semilogy(snr_db, errs./nframes/number_of_bits, '-x'); %BER in Es/No
%%%%%%%%%%%%%%%%%%%%%%%% added
hold on;
grid on;
% BER performance versus Eb/No
% Es/N0 = Eb/N0 + 10log(k) where k = log2(M) % all of them in dB
% take Eb/N0 in dB
semilogy(snr_db - 10*log10(2), errs./nframes/number_of_bits, '-x'); % BER in Eb/No
semilogy(snr_db, symbol_errs./nframes/number_of_bits, '-x'); % SER in Es/No(dB)
semilogy(snr_db - 10*log10(2), symbol_errs./nframes/number_of_bits, '-x'); %SER  in Eb/No(dB)
%%%%%%%%%%%%%%%%added
% theoretical BER versus Es_N0
a = sqrt(2.*Eb_N0_in_lin);
y= qfunc(a);
x = snr_db;
plot(x,y);
%theoretical BER versus Eb_N0
x = Eb_N0_in_dB;
plot(x,y);
% theoretical SER versus Es_No
a = sqrt(2.*Eb_N0_in_lin);
y = 2.*qfunc(a).*(1-0.5.*qfunc(a));
x = snr_db;
plot(x,y);
% theoretical SER versus Eb_N0
x = Eb_N0_in_dB;
plot(x,y);
zoom on;
xlabel('Es/No or Eb/No (in dB)');
ylabel('SER or BER');
legend('BER versus Es/No', 'BER versus Eb/No', 'SER versus Es/No', 'SER versus Eb/No','theoretical BER versus Es_N0','theoretical BER versus Eb_N0', 'theoretical SER versus Es_N0','theoretical SER versus Eb_No');
%%%%%%%%%%%%%%%%%%%%%%%% added

