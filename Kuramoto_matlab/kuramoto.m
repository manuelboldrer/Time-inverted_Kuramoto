function f=kuramoto(x,K,N,Omega,Ad,Dd,thetaref)
f = ((Omega))+ (K/N).*sum(sin(x.*Ad-Ad'.*x'),2) ;
end
