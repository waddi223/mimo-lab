% C4 Lab1 : Capacity and Outage
% Author : Waddi Abdrahamane

% Parameters

rng("default")
param.SNR = 0:1:30; % SNR range
param.SNRlin = 10.^(0.1*param.SNR);
param.Ptot = 1;
param.A = 1;
param.P = 1;
param.nr = 1;
param.nt = 1;
param.N = 10000; % Channel realizations
param.mod = 'QPSK';


%% 1 Single-antenna design

%P.1.1

avgFadingCapacity = zeros(1, length(param.SNRlin));
avgAWGNCapacity = zeros(1, length(param.SNRlin));

h = sqrt(0.5)*(randn(param.N, 1) + sqrt(-1)*randn(param.N, 1)); % Channel realizations

for ii=1:length(param.SNRlin)

    fadingCapacity = log2( 1 + param.SNRlin(ii).*(abs(h)).^2 ); % Capacity with fading
    awgnCapacity = log2( 1 + param.SNRlin(ii) ); % Capacity with AWGN only

    avgFadingCapacity(ii) = mean(fadingCapacity);
    avgAWGNCapacity(ii) = mean(awgnCapacity);

end


figure('Name', 'Average Rates')
plot(param.SNR, avgFadingCapacity)
hold on
plot(param.SNR, avgAWGNCapacity)
grid on
xlabel('SNR (dB)')
ylabel('Average Rates')
legend('Capacity w/ fading', 'Capacity w/o fading')
title('Average Rates vs SNR')

%P1.2
% Pout given in the course notes. It's P( |h|^2 < (2^R -1)/snr )
% h has a Rayleigh distribution, which implies that |h|^2 has an
% exponantial distribution

param.R = 6; % rate

Pout = 1 - exp(- (2^param.R - 1)./param.SNRlin);

figure('Name', 'Outage Probability')
semilogy(param.SNR, Pout)
grid on
xlabel('SNR (dB)')
ylabel('Outage Probability')
title('Outage Probability vs SNR')


% T1.3 The average realized data rate is derived as follow:
% When there's no outage Pout=0, the average data rate is R, considered
% that R is the predefined rate of our transmission. 
% avgRealizedRate = (1-Pout).R with Pout=0, which gives avgRealizedRate = R
% In the presence of outage, Pout~=0, the average realized is as defined
% above : avgRealizedRate = (1-Pout).R
% By using the expression of Pout, the average realized rate becomes
% avgRealizedRate = R.[ 1 - (1-exp((2.^R -1)SNR)) ] = R.exp((2.^R -1)SNR)

% P1.3 For a given average SNR, say, 15 dB, plot the average realized rate versus the 
% encoding rate R by letting the encoding rate R range from 0 to a large number4.
% Find out numerically the optimal encoding rate R* and the corresponding 
% average realized rate

% From P1.1 results, for a SNR of 15 dB the Gaussian capacity is approximatively 
% 5 bit/s/Hz, so the double is 10 bit/s/Hz. The rate R will hence range from
% 0 to 10 bit/s/Hz

param.Rrange = 0:.1:10;

avgRealizedRate = param.Rrange.*exp( -(2.^param.Rrange - 1)./15 ); % 15 dB of SNR

figure('Name', 'Average Realized Rate')
plot(param.Rrange, avgRealizedRate)
grid on
xlabel('Rate')
ylabel('Average Realized Rate')
title('Average Realized Rate vs Rate')

% The optimal encoding rate and the corresponding average rate

[optAvgRealizedRate, optAvgRealizedRateIdx] = max(avgRealizedRate);
optRate = param.Rrange(optAvgRealizedRateIdx);

fprintf('The optimal rate is %.2f\n', optRate)

% P1.4 Plot the optimized average realized rate versus the SNR in dB

optAvgRealizedRateVec = optRate.*exp( -(2.^optRate - 1)./param.SNRlin ); % 15 dB of SNR

figure('Name', 'Average Rates 2')
plot(param.SNR, optAvgRealizedRateVec, '-o')
grid on
hold on
plot(param.SNR, avgFadingCapacity, '-^')
hold on
plot(param.SNR, avgAWGNCapacity, '-*')
grid on
xlabel('SNR (dB)')
ylabel('Average Rates')
legend('Optimal average realized rate', 'Capacity w/ fading', 'Capacity w/o fading')
title('Average Rates vs SNR')

% Observation: it's observed that in the full CSI knowledge case the capacity 
% is higher than the optimal average realized case, which is normal since in 
% the previously mentioned case the lost packet are discarded

%% Multiple antennas

param.nr = 4;
param.nt = 4;

% P2.1


function [P_opt] = waterfilling(lambda,P_tot) % function provided by Prof Sheng, need to be verified
lambda = reshape(lambda,1,[]); %reshape into a row vector
[lambda, sort_idx] = sort(lambda,'descend');
L = length(lambda);
P_opt = zeros(1,L);
th = (1:L-1)./lambda(2:L) - cumsum(1./lambda(1:L-1)); %threshold
L_opt = sum(P_tot>th)+1; %optimal number of levels with non-zero power
mu_inv = (P_tot + sum(1./lambda(1:L_opt)))./L_opt; %find out 1/mu
P_opt(1:L_opt) = mu_inv - 1./lambda(1:L_opt); %optimal power allocation
P_opt(sort_idx) = P_opt;
end


avgRate = zeros(1, length(param.SNRlin));

for ii=1:length(param.SNRlin)
    rateMIMO = 0;
    for jj=1:param.N
        H = sqrt(0.5).*( randn(param.nr, param.nt) + sqrt(-1)*randn(param.nr, param.nt) );
        [U , S, Vh] = svd(H);
        lambda = diag(S); % to store the singular values of the channel
        Popt = waterfilling(lambda, param.Ptot);
        rateMIMO = rateMIMO + sum( log2(1 + param.SNRlin(ii).*(lambda.').*Popt) );
    end
    avgRate(ii) = rateMIMO./param.N;
end

figure('Name', 'Average Rate Multi Antenna')
plot(param.SNR, avgRate)
grid on
% hold on
% plot(param.SNR, avgFadingCapacity)
xlabel('SNR (dB)')
ylabel('Average Rate')
legend('Average Rate 4 Tx and 4 Rx')
title('Average Rate vs SNR')


% P2.2 Ratio between average rate Tx 1 * Rx 1 and average rate Tx4 Rx4

figure('Name', 'Ratio of Average Rates')
plot(param.SNR, avgRate./avgFadingCapacity)
grid on
xlabel('SNR (dB)')
ylabel('Ratio')
title('Ratio of Average Rates')

% Two times better rate for multi antenna case

% P2.3 The outage probability


Qx = (1/param.nt)*eye(param.nt);
outageOutMIMO = zeros(1, length(param.SNRlin));
outageOutMISO = zeros(1, length(param.SNRlin));

for ii=1:length(param.SNRlin)
    outageMIMO = 0;
    outageMISO = 0;
    for jj=1:param.N

        H = sqrt(0.5).*( randn(param.nr, param.nt) + sqrt(-1)*randn(param.nr, param.nt) );
        h = sqrt(0.5).*( randn(1, param.nt) + sqrt(-1)*randn(1, param.nt) );

        rateMIMO = log2( det( eye(param.nt) + param.SNRlin(ii)*H*Qx*H' ) );
        rateMISO = log2( 1 + param.SNRlin(ii)*(norm(h)).^2 );

        outageMIMO = outageMIMO + (rateMIMO < param.R);
        outageMISO = outageMISO + (rateMISO < param.R);

    end
    outageOutMIMO(ii) = outageMIMO./param.N;
    outageOutMISO(ii) = outageMISO./param.N;
end

figure('Name', 'Outage Probabilities')
semilogy(param.SNR, outageOutMIMO, '-o')
hold on
semilogy(param.SNR, outageOutMISO, '-*')
semilogy(param.SNR, Pout, '-^')
grid on
legend('MIMO', 'MISO', 'SISO')
xlabel('SNR (dB)')
ylabel('Outage probabilities')
title('Outage probabilities vs SNR')

%%
% P2.4 Average realized rates


% avgRealizedRateMIMO = (1-outageOutMIMO(param.SNR == 15)).*param.Rrange;
% avgRealizedRateMISO = (1-outageOutMISO(param.SNR == 15)).*param.Rrange;

% 
% outageOutMIMOvarR = zeros(1, length(param.Rrange));
% outageOutMISOvarR = zeros(1, length(param.Rrange));
% 
% for ii=1:length(param.Rrange)
%     outageMIMOvarR = 0;
%     outageMISOvarR = 0;
%     for jj=1:param.N
% 
%         H = sqrt(0.5).*( randn(param.nr, param.nt) + sqrt(-1)*randn(param.nr, param.nt) );
%         h = sqrt(0.5).*( randn(1, param.nt) + sqrt(-1)*randn(1, param.nt) );
% 
%         rateMIMO = log2( det( eye(param.nt) + 15*H*Qx*H' ) ); % 15 dB for the SNR
%         rateMISO = log2( 1 + 15*(norm(h)).^2 ); % 15 dB for the SNR
%         outageMIMOvarR = outageMIMOvarR + (rateMIMO < param.Rrange(ii));
%         outageMISOvarR = outageMISOvarR + (rateMISO < param.Rrange(ii));
% 
%     end
%     outageOutMIMOvarR(ii) = outageMIMOvarR./param.N;
%     outageOutMISOvarR(ii) = outageMISOvarR./param.N;
% end
% 
% figure(20)
% plot(param.Rrange, outageOutMIMOvarR)
% hold on
% plot(param.Rrange, outageOutMISOvarR)

% figure('Name', 'Average Realized Rates')
% plot(param.Rrange, avgRealizedRateMIMO)
% hold on
% plot(param.Rrange, avgRealizedRateMISO)