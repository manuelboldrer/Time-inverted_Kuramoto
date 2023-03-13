function f=kuramoto(x,K,N,Omega,Ad,Dd,thetaref)
% % f=Omega+(K/N).*Ad*sum(sin(x*ones(1,N)-(ones(N,1)*x')))';
 f = ((Omega))+ (K/N).*sum(sin(x.*Ad-Ad'.*x'),2) ;
%   f=(Omega) + (K/N).*sum(sin(x.*Ad-Ad'.*x'),2) ;
% %  f=Omega+ sum(exp(x.*Ad-Ad.*x')-ones(N))';
% f = (Omega)+(thetaref-x);
% f =  Omega + K/N.*Dd*cos(Dd'*x);
end