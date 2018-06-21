function [Pb] = bep_mquam_imp_markov(N,M,pa,P,t,dist,Omega,SNR_dB,SNI_dB,Ap,B,p1,p2)
% Fun��o para gera��o da BEP de um esquema de modula��o M-QAM com
% mapeamento Gray, em um canal AWGN com desvanecimento modelado segundo uma
% FSMC com N estados.
% M representa a ordem da constela��o M-QAM.
% SNR_dB representa a rela��o sinal ru�do por bit em dB.
% pa � uma matriz com N linhas, cujos elementos das linhas representam 
% respectivamento os par�metros eta e mu ou kappa e mu das distribui��es
% eta-mu e kappa-mu.
% P representa a matriz de probabilidades de transi��o.
% t � um vetor com N elementos que representa o tempo de perman�ncia no
% n-�simo estado da cadeia de Markov.
% dist � um vetor que determina qual distribui��o ser� utilizada.
% Omega � um vetor com N elementos que representa a pot�ncia do
% desvanecimento do n-�simo estado.

col = cat(2,1,zeros(1,N-1));
p = col*(P^50);
fac = sum(p.*t);
coef = (p.*t./fac);

w =@(i,k,M) (-1).^floor((i.*2.^(k-1))./sqrt(M)).*...
            (-floor(0.5+(i.*2.^(k-1))./sqrt(M))+2.^(k-1));
a =@(i,M)   3*((2*i+1).^2).*log2(M)./(M-1);

ga = 10.^(0.1*SNR_dB);
gi = 10.^(0.1*SNI_dB);
Pb = zeros(1,length(ga));
    
for j = 1:N
    switch dist(j)
    case 1
        h = (2+pa(j,1)+pa(j,1).^-1)/4;
        H = (1-pa(j,1)^2)./(4*pa(j,1));
        Mx =@(x) ((4*h*pa(j,2)^2)./...
            ((2*(h-H)*pa(j,2)+x*Omega(j)).*(2*(h+H)*pa(j,2)+x*Omega(j)))).^pa(j,2);
        f1 =@(a,ga) integral(@(o) Mx(0.5*a*ga*gi./((ga+gi).*sin(o).^2)),0,pi/2);
        f2 =@(a,ga) integral(@(o) Mx(0.5*a*ga./(sin(o).^2)),0,pi/2);
    case 2
        Mx =@(x) exp(-pa(j,1)*pa(j,2)+(pa(j,1)*(1+pa(j,1))*pa(j,2)^2)./(pa(j,2)*(1+pa(j,1))+x*Omega(j))).*...
                 (pa(j,2)*(1+pa(j,1))./(pa(j,2)*(1+pa(j,1))+x*Omega(j))).^(pa(j,2));
        f1 =@(a,ga) integral(@(o) Mx(0.5*a*ga*gi./((ga+gi).*sin(o).^2)),0,pi/2);
        f2 =@(a,ga) integral(@(o) Mx(0.5*a*ga./(sin(o).^2)),0,pi/2);
    otherwise
        disp('No defined Distribution - Wrong Value Expected')
        f1 =@(a,ga) 0;
        f2 =@(a,ga) 0;
    end

    for v = 1:length(ga)
        for k = 1:log2(sqrt(M))
            Pbk = 0;
            for i = 0:(1-2^-k)*sqrt(M)-1
                Pbk = Pbk + coef(j)*w(i,k,M)*...
                    (Ap*B*p1*p2*f1(a(i,M),ga(v))+(1-Ap*B*p1*p2)*f2(a(i,M),ga(v)));
            end
            Pb(v) = Pb(v) + Pbk;
        end
    end
end

Pb = 2*Pb./(pi*sqrt(M)*log2(sqrt(M)));

end