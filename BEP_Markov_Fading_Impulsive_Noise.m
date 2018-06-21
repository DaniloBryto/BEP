% Cálculo da BEP e BER de um esquema de modulação M-QAM em um canal AWGN
% com desvanecimento modelado por meio de uma cadeia de Markov com até 3
% estados. Rotina contendo ruído impulsivo e possibilidade de modelar os
% combinações dos desvanecimentos Alpha/Eta/Kappa-mu

clc
close all
clear all

kappa = [0.3 2.7];

M = 1024; % Ordem da constelação M-QAM
N = 5e5; % Número de bits transmitidos
J = 3; % Número de estados da Cadeia (Por sempre 3 para simulação)
pa =@(x) [1e-5 x;1e-5 x;1e-5 x]; % Parâmetros das distribuições eta/kappa-mu
P = [0.70 0.15 0.15;
     0.10 0.90 0.00;
     0.10 0.00 0.90]; % Matriz de Probabilidades de transição
w = [1 1 1]; % Vetor de tempo de permanência de estado
w = w/sum(w);
dist = [4 4 4]; % Seletor de distribuição (1 - Eta-Mu, 2 - Kappa-Mu, 3 - Gaussiana, 4 -Possível Alpha-Mu)
sig_quad = [1 1 1]; % Potências dos Desvanecimento. POR SEMPRE 1!!!

% Parâmetros do Ruído Impulsivo
Ap = 0.10; % Alpha_p. No Modelo de Poisn Equivale ao Lambda
B = 1; % Beta. No Modelo de Poison Equivale ao t
p1 = 0.20; % No Modelo de Poison Equivale ao p
p2 = 1; % No Modelo de Poison Equivale ao q
SNI_dB = 20; % SNIs(dB) observáveis

SNR_dB = linspace(0,40,3); % SNRs(dB) observáveis


tic
ber = zeros(length(SNR_dB),length(M),length(SNI_dB),length(Ap),length(kappa));
ne = zeros(1,length(SNR_dB));

for kk = 1:length(kappa)
    kk
    for pi = 1:length(Ap);
        pi;
        for v = 1:length(SNI_dB)
            v;
            for mk = 1:length(M)
                mk;
                k = log2(M(mk)); % Número de bits por símbolo
                if mod(N,log2(M(mk))) ~= 0 % Correção no número de bits
                    N = N+log2(M(mk))-mod(N,log2(M(mk)));
                end

                Eav = 2*(M(mk)-1)/3; % Potência do Sinal
                Eb = Eav/log2(M(mk)); % Energia por bit
                N0 = Eb.*10.^(-0.1*SNR_dB); % SNR adimensional
                Ni = Eb.*10.^(-0.1*SNI_dB(v)); %SNI adimensional

                for i=1:length(SNR_dB)
                    i
                    dados = randi([0 1],N,1); % Geração da sequência binária
                    dados_s2p = reshape(dados,N/k,k); % Conversão Série-Paralelo
                    dados_dec = bi2de(dados_s2p); % Conversão binario-decimal para futura
                                                  % correlação coms os símbolos da
                                                  % constelação

                    s = qammod(0:1:M(mk)-1,M(mk),0,'gray'); % Geração dos símbolos da constelação
                                                            % utilizando
                                                            % mapeamento Gray imp_noise_poison(p,q,lam,t,N0,Ni,N)

                    n = imp_noise(Ap(pi),B,p1,p2,N0(i)/2,Ni/2,length(dados_dec))+...
                        1j*imp_noise(Ap(pi),B,p1,p2,N0(i)/2,Ni/2,length(dados_dec));

                    [g,p] = var_Markov_3State(P,w,dist,pa(kappa(kk)),sig_quad,length(dados_dec)); % Desvanecimento complexo
                                                          % modelado por meio de uma
                                                          % cadeia de Markov de 3 Estados

                    r = s(dados_dec+1) + n./(g(1:length(dados_dec)).*exp(1j.*p(1:length(dados_dec)))); % Sinal recebido

                    demod = qamdemod(r,M(mk),0,'gray'); % Sinal demodulado

                    dados_demod = de2bi(demod,k); % Conversão Decimal-Binário
                    dados_p2s = dados_demod(:); % Conversão Paralelo-Serial
                    [ne(i),ber(i,mk,v,pi,kk)] = biterr(dados,dados_p2s); % BER

                end
            end
        end
    end
end
toc


%%

dist = [1 1 1];
SNR = linspace(0,40,1e1);
Be =@(x,y,Ap)bep_mquam_imp_markov(J,M,pa(y),P,w,dist,sig_quad,x,SNI_dB,Ap,B,p1,p2);

figure(1)
semilogy(SNR,Be(SNR,kappa(1),Ap),'m',...
         SNR,Be(SNR,kappa(2),Ap),'r',...
         SNR_dB,ber(:,1,1,1,1),'mx',...
         SNR_dB,ber(:,1,1,1,2),'rx',...
         'linewidth',1.2)
axis([min(SNR) max(SNR) 1e-4 10])
grid on
legend('L=1','L=2','L=3','L=4')

%%
figure(2)
histnorm(g,1e2)