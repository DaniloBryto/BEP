function g = var_envelope_kappa_mu(k,m,sig_quad,N)
% Função para geração das variáveis aleatórias da envoltória
% no modelo de desvanecimento Kappa-Mu

y = zeros(1,10001);
g = zeros(1,N);

% FDP da Envoltória no Modelo Kappa-Mu
C = (2*m*(1+k)^(0.5*m+0.5))./...
    ((k^(0.5*m-0.5))*exp(k*m)*sqrt(sig_quad)^(m+1));

pdf =@(x) C*(abs(x).^m).*exp((-m*(1+k)*abs(x).^2)/sig_quad).*...
          besseli(m-1,2*m*sqrt(k*(1+k)/sig_quad)*abs(x));

% Identificação do Máximo da FDP para Utilização
% no Limite Superior do Método de Monte Carlo
i=1;
for x = 0:0.001:5
   y(i) = pdf(x);
   i = i+1;
end

mx = max(y);

% Geração das Variáveis Aleatórias da Envoltória
% no Modelo Kappa-Mu
i=1;
while i <= N
    x1 = random('unif',0,5);
    y1 = random('unif',0,mx);
    if y1 <= pdf(x1)
        g(i)= x1;
        i=i+1;
    end
end

if N == 0
    g = [];
end

end