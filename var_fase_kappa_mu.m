function z = var_fase_kappa_mu(k,m,phi,N)
% Fun��o para Gera��o das Vari�veis Aleat�rias 
% da Fase no Modelo de Desvanecimento Kappa-Mu

z = zeros(1,N);

% FDP da Fase no Modelo Kappa-Mu
Gx =@(theta)(k^(1-0.5*m))*...
    abs(sin(2*theta)).^(0.5*m).*...
    exp(2*m*sqrt(k*(1+k))*cos(theta-phi)).*...
    besseli(0.5*m-1,2*m*sqrt(k*(1+k))*abs(cos(theta)*cos(phi))).*...
    besseli(0.5*m-1,2*m*sqrt(k*(1+k))*abs(sin(theta)*sin(phi))).*...
    sech(2*m*sqrt(k*(1+k)).*cos(theta)*cos(phi)).*...
    sech(2*m*sqrt(k*(1+k)).*sin(theta)*sin(phi));

S = integral(Gx,-pi,pi);

Fo =@(theta) (1/S)*abs(sin(2*theta)).^(0.5*m).*...
    exp(2*m*sqrt(k*(1+k))*cos(theta-phi)).*...
    besseli(0.5*m-1,2*m*sqrt(k*(1+k))*abs(cos(theta)*cos(phi))).*...
    besseli(0.5*m-1,2*m*sqrt(k*(1+k))*abs(sin(theta)*sin(phi))).*...
    sech(2*m*sqrt(k*(1+k)).*cos(theta)*cos(phi)).*...
    sech(2*m*sqrt(k*(1+k)).*sin(theta)*sin(phi));


% Identifica��o do M�ximo da FDP para Utiliza��o
% no Limite Superior do M�todo de Monte Carlo
i = 1;
for x = -pi:0.001:pi
   y(i) = Fo(x);
   i = i+1;
end

mx = max(y);

% Gera��o das Vari�veis Aleat�rias da Envolt�ria
% no Modelo Kappa-Mu
i = 1;
while i <= N
    x1 = random('unif',-pi,pi);
    y = random('unif',0,mx);
    if y <= Fo(x1)
        z(i)= x1;
        i=i+1;
    end
end

if N == 0
    z = [];
end

end
