%AYDIN UZUN
%2015401210
%EE 477 HW#3
%Please install Communications Toolbox to run this code. Because this code
%has some toolbox specific functions and classes.
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
save 8PSK_uniform_demo.mat sim_res

figure(1);
semilogy(snr_db , errs./nframes/number_of_bits, '-x'); %BER in Es/No(in dB)
hold on;
grid on;

a = sqrt(Es_N0_in_lin).*sin(pi/8);
y = (1/3).*erfc(a);
x = snr_db;
plot(x,y);

zoom on;
xlabel('Es/No (in dB)');
ylabel('BER');
legend('8PSK BER versus Es/No UNIFORM', 'theoretical 8PSK BER versus Es/N0');

