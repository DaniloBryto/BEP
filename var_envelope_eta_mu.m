function g = var_envelope_eta_mu(e,m,sig_quad,N)
% Função para geração das variáveis aleatórias da envoltória
% no modelo de desvanecimento Eta-Mu

h = (2+e+e^(-1))/4;
H = (e^(-1)-e)/4;

C = (4*sqrt(pi)*(h^(m))*(m^(m+0.5)))/...
    (gamma(m)*(H^(m-0.5))*(sig_quad^(m+0.5)));

fdp_envelope_eta_mu =@(x) C*(x.^(2*m)).*...
                          exp(-2*m*h*(x.^2)./sig_quad).*...
                          besseli(m-0.5,2*m*H*(x.^2)./sig_quad);

i=1;
for x = 0:0.0001:10
   y(i) = fdp_envelope_eta_mu(x);
   i = i+1;
end

mx = max(y);

i=1;

while i <= N
    x1 = random('unif',0,2*pi);
    y = random('unif',0,mx);
    if y <= fdp_envelope_eta_mu(x1)
        g(i)= x1;
        i=i+1;
    end
end

if N == 0
    g = [];
end

end