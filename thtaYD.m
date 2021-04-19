%% follow A3N_Kmethod_5Dec2014
function p=thtaYD(x,x1,gx1,w,lafa,beta)
global cal_SS 
S=x1;Y=gx1; 
% u=calCov_K(x,S,lafa,beta)*(calCov_K(S,S,lafa,beta)\(Y-u1))+u1;
% D=calCov_K(x,x,lafa,beta)-calCov_K(x,S,lafa,beta)*inv(calCov_K(S,S,lafa,beta))*calCov_K(S,x,lafa,beta);
cal_xS=calCov_K(x,S,lafa,beta);
u=cal_xS*(cal_SS*(Y-Reg(w,S)))+Reg(w,x);
D=calCov_K(x,x,lafa,beta)-cal_xS*(cal_SS*cal_xS');% cal_dS*(cal_SS\cal_dS')
p=exp(-.5*sum(sum(x.*x)))...
    *(2*pi)^(-size(x,2)/2)*det(D)^(-0.5)*exp(-.5*(-u)'*(D\(-u)));
