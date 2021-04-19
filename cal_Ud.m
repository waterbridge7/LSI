%---------new calculation of Ud-------
function Ud=cal_Ud(xOut,xIn,d,x1,gx1)
global w lafa beta
Ud=zeros(length(xOut),1);n=size(d,1);
yIn=zeros(n,length(xIn));
panduan=floor(max(size(gx1))/2)>=1;
% *********预计算**********
cal_dd=calCov_K(d,d,lafa,beta);

if panduan==0
    cal_S_1_S_1=[];
    cal_S_1_d=[];
else
    cal_S_1_S_1=calCov_K(x1,x1,lafa,beta);
    cal_S_1_d=calCov_K(x1,d,lafa,beta);
end

parfor j=1:length(xIn)

if panduan==0
    S=xIn(j,:);Y=0;
else
    S=[x1;xIn(j,:)];Y=[gx1;0];
end

    cal_dS=[cal_S_1_d;calCov_K(S(end,:),d,lafa,beta)]';% cal_dS=calCov_K(S,d,lafa,beta)';
    cal_S_1_end=calCov_K(S(1:end-1,:),S(end,:),lafa,beta);cal_e_e=calCov_K(S(end,:),S(end,:),lafa,beta);
    cal_SS=[cal_S_1_S_1 cal_S_1_end;cal_S_1_end' cal_e_e]\eye(size(Y,1));
%     cal_SS = inv([cal_S_1_S_1 cal_S_1_end;cal_S_1_end' cal_e_e]);

    % [calCov_K(S(1:end-1,:),S(1:end-1,:),lafa,beta) calCov_K(S(1:end-1,:),S(end,:),lafa,beta);calCov_K(S(end,:),S(1:end-1,:),lafa,beta) calCov_K(S(end,:),S(end,:),lafa,beta)]
%     uIn=cal_dS*(cal_SS*(Y-u1))+u1;
    uIn=cal_dS(1,:)*(cal_SS*(Y-Reg(w,S)))+Reg(w,d);
%    原 DIn(:,:,j)=calCov_K(d,d,lafa,beta)-calCov_K(d,S,lafa,beta)*inv(calCov_K(S,S,lafa,beta))*calCov_K(S,d,lafa,beta);
    DIn=cal_dd-cal_dS*(cal_SS*cal_dS');
    DDIn=(DIn+DIn')/2;
    yIn(:,j)=mvnrnd (uIn,DDIn,1)';
end
% options = statset('Display','final','Maxiter',1000);
options = statset('Maxiter',500);
obj = gmdistribution.fit(yIn',1,'Options',options);

parfor i=1:length(xOut)

if panduan==0
    S=xOut(i,:);Y=0;
else
    S=[x1;xOut(i,:)];Y=[gx1;0];
end
cal_dS=[cal_S_1_d;calCov_K(S(end,:),d,lafa,beta)]';% cal_dS=calCov_K(S,d,lafa,beta)';
cal_S_1_end=calCov_K(S(1:end-1,:),S(end,:),lafa,beta);cal_e_e=calCov_K(S(end,:),S(end,:),lafa,beta);
cal_SS=[cal_S_1_S_1 cal_S_1_end;cal_S_1_end' cal_e_e]\eye(size(Y,1));
% cal_SS = inv([cal_S_1_S_1 cal_S_1_end;cal_S_1_end' cal_e_e]);

u=cal_dS(1,:)*(cal_SS*(Y-Reg(w,S)))+Reg(w,d);
D=cal_dd-cal_dS*(cal_SS*cal_dS');
%  原 D=calCov_K(d,d,lafa,beta)-calCov_K(d,S,lafa,beta)*inv(calCov_K(S,S,lafa,beta))*calCov_K(S,d,lafa,beta);
%%%%%%%%%%%%%%%%----------注意 由于计算精度问题经常导致不对称或不正定。下式目的就是将D变为对称
DD=(D+D')/2;
y=mvnrnd (u,DD,1)';
% ud1=-0.5*log(det(D))+(-.5*(y-u)'*((D)\(y-u)));%ln[p(y_i|ehta_i,d)]  % +log(C)
ud1=-.5*log(det(D));
obJ=pdf(obj,y');
ud2=log(obJ);
% ud2= -(y-mu_y)'*(D_y\(y-mu_y))/2-log(det(D_y))/2;
Ud(i)=ud1-ud2;
end
Ud=sum(real(Ud))/length(xOut);
% matlabpool close

