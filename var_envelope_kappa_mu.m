function g = var_envelope_kappa_mu(k,m,sig_quad,N)
% Fun��o para gera��o das vari�veis aleat�rias da envolt�ria
% no modelo de desvanecimento Kappa-Mu

y = zeros(1,10001);
g = zeros(1,N);

% FDP da Envolt�ria no Modelo Kappa-Mu
C = (2*m*(1+k)^(0.5*m+0.5))./...
    ((k^(0.5*m-0.5))*exp(k*m)*sqrt(sig_quad)^(m+1));

pdf =@(x) C*(abs(x).^m).*exp((-m*(1+k)*abs(x).^2)/sig_quad).*...
          besseli(m-1,2*m*sqrt(k*(1+k)/sig_quad)*abs(x));

% Identifica��o do M�ximo da FDP para Utiliza��o
% no Limite Superior do M�todo de Monte Carlo
i=1;
for x = 0:0.001:5
   y(i) = pdf(x);
   i = i+1;
end

mx = max(y);

% Gera��o das Vari�veis Aleat�rias da Envolt�ria
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