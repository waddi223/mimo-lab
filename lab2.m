% C4 Lab2 : Space-time transmission
% Author : Waddi Abdrahamane

% Parameters
param.N = 10000;
param.SNR = 10:4:40; % SNR range
param.SNRlin = 10.^(0.1*param.SNR);


function param = setModParam(param)

    switch param.mod
        case "QPSK"
            param.modOrder = 4;
        case "16QAM"
            param.modOrder = 16;
        case "64QAM"
            param.modOrder = 64;
        case "256QAM"
            param.modOrder = 256;
        otherwise
            error('Unsupported modulation type.');
    end
    
    param.bitPerSymb = log2(param.modOrder);        
    param.modSymbSet = qammod(0:param.modOrder-1, param.modOrder, 'UnitAveragePower', true);

end

param.mod = "QPSK";
param = setModParam(param);

%% 1 Single-antenna channel

% P1.1 and P1.2 

ber = zeros(1, length(param.SNRlin));
per = zeros(1, length(param.SNRlin));

for ii=1:length(param.SNRlin)

    xHat = zeros(param.N, 1);

    msg = randi([0 1], param.N*param.bitPerSymb, 1);
    x = nrSymbolModulate(msg, param.mod); % directly normalized
    h = sqrt(0.5)*( randn(param.N, 1) + sqrt(-1)*randn(param.N, 1) );
    z = sqrt(0.5)*( randn(param.N, 1) + sqrt(-1)*randn(param.N, 1) );
    y = ( sqrt(param.SNRlin(ii)) .* h .* x ) + z;

    for jj=1:param.N
        eucNorm = (abs( y(jj) - h(jj).*param.modSymbSet)).^2;
        [~, minEucNormIdx] = min(eucNorm);
        xHat(jj) = param.modSymbSet(minEucNormIdx);
    end

    msgHat = nrSymbolDemodulate(xHat, param.mod, 'DecisionType','Hard');

    ber(ii) = sum(msg~=msgHat)./(param.N*param.bitPerSymb);
    per(ii) = sum(x~=xHat)./param.N;

end

figure('Name','Error Rates SISO')
semilogy(param.SNR, per, '-o')
grid on
hold on
semilogy(param.SNR, ber, '-*')
semilogy(param.SNR, (param.SNRlin).^-1, '-^')
legend('Packet Error Rate', 'Bit Error Rate', 'Theoritical Error Rate')
title('Error Rates vs SNR SISO')

% From the result, d is 1. And the diversity here is also 1.
% Explanation of the result : with higheer modulation order, symbols are
% more packed which leads to error even with slight noise. It can be said
% that the symbols are more sensitive to error and interference. But with a
% small modulation order, e.g QPSK the symbols are appart, making then more
% robust to noisy or fading channels compared to the ones using higher
% modulation order

%% 2*2 MIMO

% T2.1 Symbol rates
% Repetition code,  Symb Rate = 1/2
% Alamouti code,    Symb Rate = 1
% VBLAST,           Symb Rtae = 2
% Golden code,      Symb Rate = 2

param.SNR = 0:2:30; % SNR range
param.SNRlin = 10.^(0.1*param.SNR);
param.nt = 2;
param.nr = 2;
param.T = 2;
param.M = 256; % Message size

time.RP = 0;
time.AL = 0;
time.VB = 0;
time.GD = 0;


% param.nummc = 10; % Number of Monte Carlo

% Minimum modulation for RP is 256QAM since the message size is 256 and the same symbol is transmitted twice 
% Minimum modulation for AL is 16QAM = log2(256)/2, here two symbols transmitted twice
% Minimum modulation for VB is 4QAM = log2(256)/4, here four symbols transmitted once

% Below we'll use the minimum constellation for each scheme to respect the message size of 256

% Codebook inintialization
cb.RP = zeros(param.M, param.nr, param.nt);
cb.AL = zeros(param.M, param.nr, param.nt);
cb.VB = zeros(param.M, param.nr, param.nt);

param.modRP = "256QAM";
param.bitPerSymbRP = 8;
param.modOrderRP = 256;
param.modSymbSetRP = qammod(0:param.modOrderRP-1, param.modOrderRP, 'UnitAveragePower', true);

param.modAL = "16QAM";
param.bitPerSymbAL = 4;
param.modOrderAL = 16;
param.modSymbSetAL = qammod(0:param.modOrderAL-1, param.modOrderAL, 'UnitAveragePower', true);

param.modVB = "QPSK";
param.bitPerSymbVB = 2;
param.modOrderVB = 4;
param.modSymbSetVB = qammod(0:param.modOrderVB-1, param.modOrderVB, 'UnitAveragePower', true);


% Codebook Repetition Scheme
for m=1:param.M
    cb.RP(m, :, :) = param.modSymbSetRP(m).*eye(param.nr); % Codebook RP
end

% Codebook Alamouti
for i=1:param.modOrderAL
    s1 = param.modSymbSetAL(i);
    for j=1:param.modOrderAL
        s2 = param.modSymbSetAL(j);
        cb.AL((i-1) * param.modOrderAL + j, :, :) = [s1, -conj(s2); s2, conj(s1)]; 
    end
end

% Codebook V-BLAST
idxVB = 1;
for i=1:param.modOrderVB
    s1 = param.modSymbSetVB(i);
    for j=1:param.modOrderVB
        s2 = param.modSymbSetVB(j);
        for k=1:param.modOrderVB
            s3 = param.modSymbSetVB(k);
            for l=1:param.modOrderVB
                s4 = param.modSymbSetVB(l);
                cb.VB(idxVB, : , :) = [s1, s3; s2, s4];
                idxVB = idxVB + 1;
            end
        end
    end
end

per = zeros(4, length(param.SNRlin)); % packet error rate initialization

% RP simulation
tic
for ii=1:length(param.SNRlin)

    % for i=1:param.nummc

        msg = randi([0 param.modOrderRP-1], param.N, 1);
        x = qammod(msg, param.modOrderRP, 'UnitAveragePower', true);

        xHat = zeros(param.N, 1);

        for jj=1:param.N

            xRP = x(jj).*eye(param.T);

            h = sqrt(0.5)*( randn(param.nr, param.nt) + sqrt(-1)*randn(param.nr, param.nt) );
            z = sqrt(0.5)*( randn(param.nr, param.T) + sqrt(-1)*randn(param.nr, param.T) );


            yRP = ( sqrt(param.SNRlin(ii)) .* h * xRP ) + z;

            eucNorm = zeros(param.M, 1);
            for m=1:param.M
                % tempCBmtx : 2*2 matrix of the codebook that will be used for ML decoding. Previously use directly from cb.RP with squeeze, but squeeze is computationaly expensive
                tempCBmtx = cb.RP(m, :, :);
                tempCBmtx = reshape(tempCBmtx, [param.nr, param.nt]);
                eucNorm(m) = norm( yRP - sqrt(param.SNRlin(ii)).*h*tempCBmtx ).^2;
            end
            [~, minEucNormIdx] = min(eucNorm);
            xHat(jj) = cb.RP(minEucNormIdx, 1, 1);

        end
    per(1, ii) = per(1, ii) + sum(x~=xHat)./param.N;

    % end
    % per(ii) = per(ii)./param.nummc;
end
time.RP = toc;

% AL simulation
tic;
for ii=1:length(param.SNRlin)

        msg = randi([0 param.modOrderAL-1], param.N, 1);
        x = qammod(msg, param.modOrderAL, 'UnitAveragePower', true);

        xHat = zeros(param.N, 1);
        for jj=1:2:param.N
            xAL = [x(jj), -conj(x(jj+1)); x(jj+1), conj(x(jj))];

            h = sqrt(0.5)*( randn(param.nr, param.nt) + sqrt(-1)*randn(param.nr, param.nt) );
            z = sqrt(0.5)*( randn(param.nr, param.T) + sqrt(-1)*randn(param.nr, param.T) );

            yAL = ( sqrt(param.SNRlin(ii)) .* h * xAL ) + z;

            eucNormAL = zeros(param.M, 1);
            for m=1:param.M
                tempCBmtx = cb.AL(m, :, :);
                tempCBmtx = reshape(tempCBmtx, [param.nr, param.nt]);
                eucNormAL(m) = norm( yAL - sqrt(param.SNRlin(ii)).*h*tempCBmtx ).^2;
            end
            [~, minEucNormIdx] = min(eucNormAL);

            xHat(jj) = cb.AL(minEucNormIdx, 1, 1);
            xHat(jj+1) = cb.AL(minEucNormIdx, 2, 1);

        end
        per(2, ii) = per(2, ii) + sum(x~=xHat)./param.N;

end
time.AL = toc;

% VB simulation
tic;
for ii=1:length(param.SNRlin)

        msg = randi([0 param.modOrderVB-1], param.N, 1);
        x = qammod(msg, param.modOrderVB, 'UnitAveragePower', true);

        xHat = zeros(param.N, 1);
        for jj=1:4:param.N
            xVB = [x(jj), x(jj+2); x(jj+1), x(jj+3)];

            h = sqrt(0.5)*( randn(param.nr, param.nt) + sqrt(-1)*randn(param.nr, param.nt) );
            z = sqrt(0.5)*( randn(param.nr, param.T) + sqrt(-1)*randn(param.nr, param.T) );

            yVB = ( sqrt(param.SNRlin(ii)) .* h * xVB ) + z;

            eucNormVB = zeros(param.M, 1);
            for m=1:param.M
                tempCBmtx = cb.VB(m, :, :);
                tempCBmtx = reshape(tempCBmtx, [param.nr, param.nt]);
                eucNormVB(m) = norm( yVB - sqrt(param.SNRlin(ii)).*h*tempCBmtx ).^2;
            end
            [~, minEucNormIdx] = min(eucNormVB);

            xHat(jj:jj+3) = reshape(cb.VB(minEucNormIdx, :, :), [4, 1]);

        end
        per(3, ii) = per(3, ii) + sum(x~=xHat)./param.N;
end
time.VB = toc;

% GD parameters
cb.GD = zeros(param.M, param.nr, param.nt);
cb.GDsymbs = zeros(param.M, 4);
phi = (1+sqrt(5))/2;
phiHat = 1-phi;
alph = 1 + 1i*phi;
alphHat = 1 - 1i*phiHat;

% Use same modulation as VB
% Golden Code
idxGD = 1;
for i=1:param.modOrderVB
    s1 = param.modSymbSetVB(i);
    for j=1:param.modOrderVB
        s2 = param.modSymbSetVB(j);
        for k=1:param.modOrderVB
            s3 = param.modSymbSetVB(k);
            for l=1:param.modOrderVB
                s4 = param.modSymbSetVB(l);
                cb.GD(idxGD, : , :) = (1/sqrt(5)).*[alph.*(s1 + s2.*phi)         , alph.*(s3 + s4.*phi); ...
                                                    1i*alphHat.*(s3 + s4*phiHat) , alphHat.*(s1 + s2.*phiHat)];
                cb.GDsymbs(idxGD, : ) = [s1 s2 s3 s4];
                idxGD = idxGD + 1;
            end
        end
    end
end

tic
for ii=1:length(param.SNRlin)

        msg = randi([0 param.modOrderVB-1], param.N, 1);
        x = qammod(msg, param.modOrderVB, 'UnitAveragePower', true);

        xHat = zeros(param.N, 1);
        for jj=1:4:param.N

            s1 = x(jj); s2 = x(jj+1); s3 = x(jj+2); s4 = x(jj+3);
            xGD = (1/sqrt(5)).*[alph.*(s1 + s2.*phi)         , alph.*(s3 + s4.*phi); ...
                                                    1i*alphHat.*(s3 + s4*phiHat) , alphHat.*(s1 + s2.*phiHat)];

            h = sqrt(0.5)*( randn(param.nr, param.nt) + sqrt(-1)*randn(param.nr, param.nt) );
            z = sqrt(0.5)*( randn(param.nr, param.T) + sqrt(-1)*randn(param.nr, param.T) );

            yGD = ( sqrt(param.SNRlin(ii)) .* h * xGD ) + z;

            eucNormVB = zeros(param.M, 1);
            for m=1:param.M
                tempCBmtx = cb.GD(m, :, :);
                tempCBmtx = reshape(tempCBmtx, [param.nr, param.nt]);
                eucNormVB(m) = norm( yGD - sqrt(param.SNRlin(ii)).*h*tempCBmtx ).^2;
            end
            [~, minEucNormIdx] = min(eucNormVB);
            
            xHat(jj:jj+3) = cb.GDsymbs(minEucNormIdx, :);
            
        end
        per(4, ii) = per(4, ii) + sum(x~=xHat)./param.N;

end
time.GD = toc;

marker = {'-o', '-^','-+', '-d'};

figure('Name','Error Rates MIMO')

for p=1:4
semilogy(param.SNR, per(p, :), marker{p})
hold on
grid on
end
xlabel('SNR (dB)')
ylabel('Error Rates')
legend('RP', 'AL', 'VB', 'GD')
title('Error rates vs SNR MIMO')

% This simulation can be optimized after to reduce the complexity and
% running time