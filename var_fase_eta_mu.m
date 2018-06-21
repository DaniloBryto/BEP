function z = var_fase_eta_mu(e,m,N)
% Função para geração das variáveis aleatórias da fase
% no modelo de desvanecimento eta-mu

h = (2+e+e^(-1))/4;
H = (e^(-1)-e)/4;

C = (h.^2-H.^2).*gamma(2*m)./((4.^m).*gamma(m).^2);

fdp_fase_eta_mu =@(theta) C*abs(sin(2*theta)).^(2*m-1)./...
                          ((h+H*cos(2*theta)).^(2*m));

i = 1;

for x = -pi:0.0001:pi
   y(i) = fdp_fase_eta_mu(x);
   i = i+1;
end

mx = max(y);

i = 1;
while i <= N
    x1 = random('unif',-pi,pi);
    y = random('unif',0,mx);
    if y <= fdp_fase_eta_mu(x1)
        z(i)= x1;
        i=i+1;
    end
end

if N == 0
    z = [];
end

end
