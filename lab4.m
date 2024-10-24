% C4 Lab4 Uplink and Downlink Communications
% Author : Waddi Abdrahamane

% Parameters
param.M = 8; % Number of antenna at BS
param.K = 6; % Number of user with single antenna
param.N = 1002; % Total time or number of channel realization
param.SNR = 0:4:40; % SNR range
param.SNRlin = 10.^(0.1*param.SNR);

% We consider TDD (timedivision duplex), i.e., uplink and downlink communications take place in the same
% band but successively. We are interested in the achievable uplink and downlink rates
% with linear precoding and linear detection.

% 1 Uplink
% 1.1 TDMA

% P1.1 & P1.2 The sum rate of the K users and average sum rate

% Each block here corresponds to one channel realization

% Time allocation
function param = setTimeAlloc(param, timeAlloc)
    switch timeAlloc
        case 'Random' % Allocation done randomly    
        % param.theta = rand(1, param.K);
        % param.theta = param.N*param.theta / sum(param.theta);
        % param.theta = round(param.theta);
        % 
        % while sum(param.theta) ~= param.N
        %     idx = randi(param.K);
        %     param.theta(idx) = param.theta(idx) + (param.N - sum(param.theta));
        % end
        % %disp(sum(param.theta))
        % 
        % param.theta = param.theta/param.N;
            param.timeAlloc = 'Random';
            param.theta = rand(1, param.K);
            param.theta = param.theta / sum(param.theta);

        case 'Uniform' % Uniform time allocation
            param.timeAlloc = 'Uniform';
            param.theta = ones(1, param.K)./param.K;
        otherwise
            param.timeAlloc = timeAlloc;
    end
end


function [avgSumRate] = genAvgSumRateTDMA(param, timeAlloc)

    param = setTimeAlloc(param, timeAlloc);
    avgSumRate = zeros(1, length(param.SNR));

    for ii=1:length(param.SNRlin)
        totalSumRate = 0;
        for jj=1:param.N

            H = sqrt(0.5) * (randn(param.M, param.K) + 1i*randn(param.M, param.K));
            if param.timeAlloc == "Optimal"
                param.theta = sum(abs(H).^2) ./ sum(sum(abs(H).^2));
            end
            sumRate = 0;
            for k=1:param.K
                sumRate = sumRate + param.theta(k) * log2(1 + (param.SNRlin(ii)* (norm(H(:, k)).^2))/param.theta(k) );
            end
            totalSumRate = totalSumRate + sumRate;
        end
        avgSumRate(ii) = totalSumRate/param.N;
    end
end

avgSumRateRandom1 = genAvgSumRateTDMA(param, "Random");
avgSumRateRandom2 = genAvgSumRateTDMA(param, "Random");
avgSumRateUniform = genAvgSumRateTDMA(param, "Uniform");
avgSumRateOptimal = genAvgSumRateTDMA(param, "Optimal");



% 1.2 SDMA

function [avgSumRate] = genAvgSumRateSDMA(param)
    avgSumRate = zeros(1, length(param.SNR));

    for ii=1:length(param.SNRlin)
        totalSumRate = 0;
        for jj=1:param.N

            H = sqrt(0.5) * (randn(param.M, param.K) + 1i*randn(param.M, param.K));
            sumRate = 0;
            for k=1:param.K
                orthoMtx = null(H(:, [1:k-1, k+1:end])');
                gk = orthoMtx(:, 1); % Just one vector, it can be any from the null space
                sumRate = sumRate + log2(1 + (param.SNRlin(ii)*abs(gk'*H(:, k)).^2)/norm(gk)^2 );
            end
            totalSumRate = totalSumRate + sumRate;
        end
        avgSumRate(ii) = totalSumRate/param.N;
    end
end

avgSumRateSDMA = genAvgSumRateSDMA(param);


function [avgSumCapacity] = genAvgSumCapacitySDMA(param)
    avgSumCapacity = zeros(1, length(param.SNR));
    
    for ii=1:length(param.SNRlin)
        totalSumCapacity = 0;
        for jj=1:param.N
    
            H = sqrt(0.5) * (randn(param.M, param.K) + 1i*randn(param.M, param.K));
            mtx = 0;
            for k=1:param.K
                mtx = mtx + H(:, k)*H(:, k)';
            end
            totalSumCapacity = totalSumCapacity + log2(det( eye(param.M) + param.SNRlin(ii)*mtx ));
        end
        avgSumCapacity(ii) = totalSumCapacity/param.N;
    end
end

avgSumCapacity = genAvgSumCapacitySDMA(param);


marker = {'-o', '-^','-+', '-d', '-*', '-s'};

figure("Name",'Average Sum Rates')
plot(param.SNR, avgSumRateRandom1, marker{1})
hold on 
grid on
plot(param.SNR, avgSumRateRandom2, marker{2})
plot(param.SNR, avgSumRateUniform, marker{3})
plot(param.SNR, avgSumRateOptimal, marker{4})
plot(param.SNR, avgSumRateSDMA, marker{5})
plot(param.SNR, avgSumCapacity, marker{6})
xlabel('SNR (dB)')
ylabel('Avg Sum Rate')
legend('TDMA Random time alloc 1', 'TDMA Random time alloc 2', 'TDMA Uniform time alloc', 'Optimal time alloc', 'SDMA Avg Sum Rate', 'SDMA Sum Capacity')
title('Average Sum Rates')