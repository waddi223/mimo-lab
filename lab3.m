% C4 Lab3 Receiver
% Author : Waddi Abdrahamane

% We consider a block fading channel with coherence time T, namely, in each block b,
% b = 1, . . . , B,
% Ynr×T[b] = Hnr×nt[b]Xnt×T[b] + Znr×T[b]

% 1. Signal generation

% The matrix X is considered as a packet of length T and is composed of two parts: 
% the training (pilot) matrix Xp with nt columns and data matrix Xd with T′ = T − nt columns
% X = [X_p X_d];
% we assume that the pilot matrix Xp is unitary (e.g., identity matrix).
% The data matrices are generated randomly with independent QAMsymbols. We assume that nt = 2, nr = 8,
% T = 10. The average transmit power per channel use is P

% The size M of the constellation : 128 bits/(2{nr}*8{T-nr}) = 8 bits = 256QAM


% P1.1

% Parameters
param.B = 100;
param.SNR = 0:2:40; % SNR range
param.SNRlin = 10.^(0.1*param.SNR);
param.modOrder = 256; % 256QAM
param.mod = "256QAM";
param.nt = 2;
param.nr = 8;
param.T = 10;


% for ii=1:length(param.SNRlin)
% 
%     for b=1:param.B
% 
%         % Tx
%         msg = randi([0 param.mod-1], param.nt, param.T-param.nt);
%         Xd = qammod(msg, param.mod, 'UnitAveragePower', true); % Replace and add normalization factor
%        
%         
%         Xp = eye(param.nt);
%         X = [Xp, Xd];
% 
%         % Channel
%         H = sqrt(0.5) * ( randn(param.nr,param.nt) + 1i * randn(param.nr,param.nt) );
%         Z = sqrt(0.5) * ( randn(param.nr,param.T) + 1i * randn(param.nr,param.T) );
% 
%         % Rx
%         Y = ( sqrt(param.SNRlin(ii)) .* H * X ) + Z;
% 
%     end
% end

% PAMsize = 4;
% 
% const = qammod(0:PAMsize-1, PAMsize, "UnitAveragePower",true);
% bitTable = int2bit((0:PAMsize-1), log2(PAMsize)).';
% spheredecoder = comm.SphereDecoder(const.', bitTable, DecisionType="Hard");
% Need to implement my own sphere decoder algorithm
% % 2.Perfect CSIR
% P2.1
% e = zeros(param.B,1);
% per = zeros(length(param.SNR), 1);
% for ii=1:length(param.SNRlin)
% for k = 1:param.B
% H = sqrt(0.5*param.SNRlin(ii)).*randn(param.nr,param.nt)+1i*randn(param.nr,param.nt);
% x = ones(param.nt,1)*(3+1i);
% y = H*x + (randn(param.nr,1)+1i*randn(param.nr,1))*0.1;
% x_hat = spheredecoder(y.', reshape(H.', [1, param.nt,param.nr]));
% e(k) = isequal(x_hat,x);
% end
% p = 1-mean(e);
% per(ii) = p;
% end

% P2.2

% Below we can observe two normalization factors
% The first is done using the whole constellation, each symbol has an
% average power of P. The result for this alpha is a reduced error rate
% since we're transmitting with high power, this can be observed in the ber
% plot. But the problem is the high power required or used

% The second normalization factor give us an average power P per channel
% use. Since with this configuration we use less power compared to the
% previous one, we observe a reduced ber performance compared to the
% previous configuration ber. The present configuration will be used

ber = zeros(2, length(param.SNR));
for ii=1:length(param.SNRlin)

    %alphaXd = sqrt( (param.SNRlin(ii)*param.modOrder) / (sum(abs(qammod((0:param.modOrder-1), param.modOrder)) .^2))); % Normalization
    % factor constellation

    err = zeros(1, param.B);    
    for b=1:param.B

        % Tx
        msg = randi([0 param.modOrder-1], param.nt, param.T-param.nt); % Use bits instead

        Xd = qammod(msg, param.modOrder);
        alphaXd = sqrt( (param.SNRlin(ii)) ./ ( sum(abs(Xd).^2) )); %Normalization such that average transmit power per channel use is P
        Xd = alphaXd.*Xd;


        % Channel
        H = sqrt(0.5) * ( randn(param.nr,param.nt) + 1i * randn(param.nr,param.nt) );
        Z = sqrt(0.5) * ( randn(param.nr,param.T-param.nt) + 1i * randn(param.nr,param.T-param.nt) );

        % Rx
        Y = H * Xd  + Z;
        XHat = (pinv(H)*Y)./alphaXd;
        msgHat = qamdemod(XHat, param.modOrder);

        err(b) = sum(sum( int2bit(msg, log2(param.modOrder)) ~= int2bit(msgHat, log2(param.modOrder)) ))/numel(int2bit(msg, log2(param.modOrder)));
    end
    ber(1, ii) = mean(err);

end


% 3 Imperfect CSIR

for ii=1:length(param.SNRlin)
    %alphaXd = sqrt( (param.SNRlin(ii)*param.modOrder) / (sum(abs(qammod((0:param.modOrder-1), param.modOrder)) .^2)));

    err = zeros(1, param.B);    
    for b=1:param.B

        % Tx
        msg = randi([0 param.modOrder-1], param.nt, param.T-param.nt); % Use bits instead
        
        Xd = qammod(msg, param.modOrder);
        alphaXd = sqrt( (param.SNRlin(ii)) ./ ( sum(abs(Xd).^2) )); %Normalization such that average transmit power per channel use is P
        Xd = alphaXd.*Xd;

        alphaXp = param.SNRlin(ii);
        Xp = alphaXp.*eye(param.nt);

        X = [Xp Xd];

        % Channel
        H = sqrt(0.5) * ( randn(param.nr,param.nt) + 1i * randn(param.nr,param.nt) );
        Z = sqrt(0.5) * ( randn(param.nr,param.T) + 1i * randn(param.nr,param.T) );
        
        Y = H * X  + Z;


        % Rx
        Yp = Y( : , 1 : param.nt);
        Yd = Y( : , param.nt + 1 : end);

        HHat = (Yp*Xp'/(Xp*Xp'));

        XdHat = (pinv(HHat)*Yd)./alphaXd;
        msgHat = qamdemod(XdHat, param.modOrder);

        err(b) = sum(sum( int2bit(msg, log2(param.modOrder)) ~= int2bit(msgHat, log2(param.modOrder)) ))/numel(int2bit(msg, log2(param.modOrder)));
    end
    ber(2, ii) = mean(err);
end


for ii=1:length(param.SNRlin)
    %alphaXd = sqrt( (param.SNRlin(ii)*param.modOrder) / (sum(abs(qammod((0:param.modOrder-1), param.modOrder)) .^2)));

    err = zeros(1, param.B);    
    for b=1:param.B

        % Tx
        msg = randi([0 param.modOrder-1], param.nt, param.T-param.nt); % Use bits instead

        Xd = qammod(msg, param.modOrder);
        alphaXd = sqrt( (param.SNRlin(ii)) ./ ( sum(abs(Xd).^2) )); %Normalization such that average transmit power per channel use is P
        Xd = alphaXd.*Xd;

        alphaXp = param.SNRlin(ii);
        Xp = alphaXp.*eye(param.nt);
        Xp = [Xp Xp];

        X = [Xp Xd];

        % Channel
        H = sqrt(0.5) * ( randn(param.nr,param.nt) + 1i * randn(param.nr,param.nt) );
        Z = sqrt(0.5) * ( randn(param.nr,param.T + param.nt) + 1i * randn(param.nr,param.T + param.nt) );

        Y = H * X  + Z;


        % Rx
        Yp = Y( : , 1 : 2*param.nt);
        Yd = Y( : , 2*param.nt + 1 : end);

        HHat = Yp*Xp'/(Xp*Xp');

        XdHat = (pinv(HHat)*Yd)./alphaXd;
        msgHat = qamdemod(XdHat, param.modOrder);

        err(b) = sum(sum( int2bit(msg, log2(param.modOrder)) ~= int2bit(msgHat, log2(param.modOrder)) ))/numel(int2bit(msg, log2(param.modOrder)));
    end
    ber(3, ii) = mean(err);
end

figure('Name','BER vs SNR')
grid on
semilogy(param.SNR, ber(1, :), '-o')
grid on
hold on
semilogy(param.SNR, ber(2, :), '-^')
semilogy(param.SNR, ber(3, :), '-*')
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
legend('ZF Perfect CSIR', 'ZF Imperfect CSIR 1 Pilot mtx', 'ZF Imperfect CSIR 2 Pilot mtx')