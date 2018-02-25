close all; clear all; clc

parD.U = 20; % number of UEs
parD.N = 256; % number of BS antennas
parD.trials = 128; % number of Monte-Carlo trials (transmissions)
parD.rHe = 0; % relative channel estimate error
parD.SNRdB_list = -15:3:15; % list of SNR [dB] values to be simulated 
parD.mod = '64QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
parD.precoder = {'MRT','MRTi', 'ZF','ZFi', 'WF', 'WFi'}; % 
%  , 'SQUID'
switch (parD.mod)
    case 'BPSK'
        parD.symbols = [ -1 1 ];
    case 'QPSK'
        parD.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
    case '16QAM'
        parD.symbols = [...
            -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM'
        parD.symbols = [...
            -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
end


parD.symbols = parD.symbols/sqrt(mean(abs(parD.symbols).^2));

% precompute bit labels
parD.card = length(parD.symbols); % cardinality
parD.bps = log2(parD.card); % number of bits per symbol
parD.bits = de2bi(0:parD.card-1,parD.bps,'left-msb'); % symbols-to-bits


% S = [];
% for t=1:parD.trials
%     
%         % generate random bit stream
%         b = randi([0 1],parD.U,parD.bps);
% 
%         % generate transmit symbols
%         idx = bi2de(b,'left-msb')+1;
%         s = parD.symbols(idx).';
%         S = [s S];
% end



t = 120; % iteration num
c = 20/128; % dimension ratio
dx = 0.05; 
sigma = 1;  
v = [];
for i = 1:t
        S = [];
        for t=1:parD.trials
    
        % generate random bit stream
                bb = randi([0 1],parD.U,parD.bps);

                % generate transmit symbols
                idx = bi2de(bb,'left-msb')+1;
                s = parD.symbols(idx).';
                S = [s S];
        end
%      X = randn(N,T);
     s = S*S';
     vv = eig(s)/128;
     v = [v vv];
end

v = sort(reshape(v,1,size(v,1)*size(v,2))); 
a = sigma*(1-sqrt(c))^2; b = sigma*(1+sqrt(c))^2;

[count, x] = hist(v, a:dx:b);
cla reset
bar(x, count/(t*20*dx), 'w'); % t*
hold on;
x = linspace (a, b);
plot(x, sqrt((x-a).*(b-x))./(2*pi*x*c*sigma),'--k', 'LineWidth', 2)
% axis([0 ceil(b) -0.1 1.5]);
xlabel('x');  ylabel('f(x)');
legend('EED', 'MP-law')
