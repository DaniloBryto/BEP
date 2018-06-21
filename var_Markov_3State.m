function [G,Ph] = var_Markov_3State(P,w,dist,pa,sig_quad,N)
% Função para geração das variáveis aleatórias do desvanecimento complexo
% modelado segundo uma FSMC. Capaz de modelar as distribuições eta-mu,
% kappa-mu e normal.
% p0 representa o vetor de probabilidades iniciais
% P representa a matriz de probabilidade de transição
% w representa o tempo médio de permanência no estado
% dist determina qual a distribuição a ser utilizada
% pa representa os parâmetros das distribuições
% sig_quad representa a potência do desvanecimento

pe = [1 0 0]*(P^50);

lengthChain = ceil(N/sum(pe.*w)); % Número de Transições

aux = 1;
te = zeros(1,3);
while sum(te) < N % Assegurar VAs suficientes para os símbolos
    [channel] = markovChain3States(P,ceil(lengthChain/aux)); % Estados da Cadeia
    for i = 1:length(channel)
        if channel(i) == 1
            te(1) = te(1) + exprnd(w(1));% exprnd(w(1)); % Tempo Total de Permanência no Estado 1
        elseif channel(i) == 2
            te(2) = te(2) + exprnd(w(2)); % Tempo Total de Permanência no Estado 2
        else
            te(3) = te(3) + exprnd(w(3)); % Tempo Total de Permanência no Estado 3
        end
    end
    aux = aux + 1;
end
t_show = te/sum(te);

G = [];
Ph = [];
for i=1:3
    switch dist(i)
    case 1 % Modelo eta-mu
        g =@(x) var_envelope_eta_mu(pa(i,1),pa(i,2),sig_quad(i),x);
        p =@(x) var_fase_eta_mu(pa(i,1),pa(i,2),x);
    case 2 % Modelo kappa-mu
        g =@(x) var_envelope_kappa_mu(pa(i,1),pa(i,2),sig_quad(i),x);
        p =@(x) var_fase_kappa_mu(pa(i,1),pa(i,2),pi/13,x);
    case 3 % Modelo Gaussiano
        g =@(x) normrnd(pa(i,1),pa(i,2),[1 x]);
        p =@(x) normrnd(pa(i,1),pa(i,2),[1 x]);
    case 4
        g =@(x) randraw('Nakagami',pa(i,2),[1 x]);
        p =@(x) pi*(1-2*rand(1,x));
    otherwise
        fprintf('Distribuição não encontrada\n')
        fprintf('Vai dar errado!!!\n')
        g =@(x) zeros(1,x);
        p =@(x) zeros(1,x);
    end
    G = cat(2,G,g(round(te(i))));
    Ph = cat(2,Ph,p(round(te(i))));
end

end