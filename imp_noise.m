function n = imp_noise(Ap,B,p1,p2,N0,Ni,N)

p = round(Ap*B*p1*p2*N);

n = cat(2,normrnd(0,sqrt(N0),[1 N-p]),normrnd(0,sqrt(Ni+N0),[1 p]));

end
