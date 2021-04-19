clear;clc;
% initial 10 data
% rand('seed',12);
xlati=16*(lhsdesign(10,2)-0.5);
x1=xlati;
gx1=g(x1(:,1),x1(:,2));
save x1_gx1 x1 gx1

% initial hyper-parameters of GP
global lafa beta cal_SS 
% lafa=1 ;beta=2;%% 0.8 0.2
global dimo par u0 D0 w
par=4; 
m=50; % for speed, we only use 50 iterations. 
dimo=2;u0=0;D0=1;
burnin = 5000;
% loop 10 iterations

% ---------------------initial-----------------------------------------
w0=rand(6,1); w=nlinfit(x1,gx1,@Reg,w0);
gprMdl = fitrgp(x1,gx1- Reg(w,x1),'KernelFunction','squaredexponential','KernelParameters',[1,2]);% ,'Sigma',1e-2
lafa = gprMdl.KernelInformation.KernelParameters(2,1);
beta = gprMdl.KernelInformation.KernelParameters(1,1);
%%--------------------MCMC for posterior_i (z)-------------------------
disp(' sample posterior p_i(z)')

cal_SS=calCov_K(x1,x1,lafa,beta)\eye(size(x1,1));
delta = .5;
pdf1 = @(xx) thtaYD(xx,x1,gx1,w,lafa,beta);
proprnd = @(xx) xx+rand(size(xx))*2*delta-delta ;   
nsamples =1.5*10^4;% 尽量得到2W左右的采样点
xx = mhsample([1,1],nsamples,'pdf',pdf1,'proprnd',proprnd,'symmetric',1);
pdf1 = @(x) thtaYD(x,x1,gx1,w,lafa,beta);
proprnd = @(x) x+rand(size(x))*2*delta-delta ;   
x = mhsample([-1,-1],nsamples,'pdf',pdf1,'proprnd',proprnd,'symmetric',1);
x=[x(burnin:end,:);xx(burnin:end,:)];
x = unique(x,'rows');

for i=1:15
    
    %  ---------------experimental design---------------
    disp(strcat('loop:',num2str(i)))
    x = x(randperm(length(x)),:);
    n_x=size(x,1);
    randindex=randperm(n_x); 
    if (mod(n_x,2)==1)
        xOut=x(randindex(1:(n_x+1)/2),:);xIn=x(randindex((n_x+1)/2:end),:);
    else
        xOut=x(randindex(1:n_x/2),:);xIn=x(randindex(n_x/2:end),:);
    end
    d=x(round([10, n_x*2/5, n_x*3/5, n_x-10]),:)+0.1; %随机取初值 
    % d0=[3 3;3 -3;-3 3;-3 -3];
    record_d=zeros(m,numel(d));record_Ud=zeros(m,1);
    save record record_d record_Ud
    k=1;
    disp(' optimal experimal design')
    for n=1:m
        Ud=cal_Ud(xOut,xIn,d,x1,gx1);
        record_Ud(n)=Ud;
        record_d(n,:)=reshape(d',1,dimo*par);
        disp([record_Ud(n) - record_Ud(1), n, record_d(n,:) - record_d(1,:)]);
        %-------------------d*=argmaxUd------------------------------->>
        a=.8*5 ;c=.1*10;A=100;alfa=0.602;gama=0.101; 
        ak=a/(A+k+1)^alfa;Ck=c/(k+1)^gama;
        % delta=2*ceil(rand(size(d,1),2)-0.5)-1;
        delta=round((rand(size(d,1),dimo)));
        d1=d+Ck*delta;d2=d-Ck*delta;
        Ud1=cal_Ud(xOut,xIn,d1,x1,gx1);
        Ud2=cal_Ud(xOut,xIn,d2,x1,gx1);
        % gk=(Ud1-Ud2)./(2*Ck*delta)*1;
        gk=(Ud1-Ud2)./(2*Ck)*delta;
        if Ud1>Ud || Ud2>Ud
            dnew=d+ak*gk;
        else
            dnew=d;
        end
        d=dnew;
        k=k+1;
    end
    
    %%-------------------evaluate simulation; get data------------------------
    d;gd=g(d(:,1),d(:,2));
    x1=[x1;d];gx1=[gx1;gd];

    % optimal hyper-parameters
    w0=rand(6,1); w=nlinfit(x1,gx1,@Reg,w0);
    gprMdl = fitrgp(x1,gx1- Reg(w,x1),'Basis','none','KernelFunction','squaredexponential','KernelParameters',[1,2],'Sigma',1e-2,'ConstantSigma',true);% 
    lafa = gprMdl.KernelInformation.KernelParameters(2,1);
    beta = gprMdl.KernelInformation.KernelParameters(1,1);
    %%--------------------MCMC for posterior_i (z)-------------------------
    disp(' sample posterior p_i(z)')

    cal_SS=calCov_K(x1,x1,lafa,beta)\eye(size(x1,1));
    delta = .5;
    pdf1 = @(xx) thtaYD(xx,x1,gx1,w,lafa,beta);
    proprnd = @(xx) xx+rand(size(xx))*2*delta-delta ;   
    nsamples =1.5*10^4;% 尽量得到2W左右的采样点
    xx = mhsample(x(round(n_x/5),:),nsamples,'pdf',pdf1,'proprnd',proprnd,'symmetric',1);
    pdf1 = @(x) thtaYD(x,x1,gx1,w,lafa,beta);
    proprnd = @(x) x+rand(size(x))*2*delta-delta ;   
    x = mhsample(x(round(n_x*4/5),:),nsamples,'pdf',pdf1,'proprnd',proprnd,'symmetric',1);
    x=[x(burnin:end,:);xx(burnin:end,:)];
    x = unique(x,'rows');

figure
[X,Y]=meshgrid(-8:.5:8);
Z=g(X,Y);
contour(X,Y,Z,[0,0]);
hold on

x_contour = [reshape(X,33*33,1),reshape(Y,33*33,1)];
cal_xS=calCov_K(x_contour,x1,lafa,beta);
U=cal_xS*(cal_SS*(gx1-Reg(w,x1)))+Reg(w,x_contour);
contour(X,Y,reshape(U,33,33),[0,0],'r');

hold on
plot(x(:,1),x(:,2),'.','color',[.5 .5 .5],'MarkerSize',4)
hold on
plot(x1(1:end,1),x1(1:end,2),'*','Color','b');hold on 
hold on
plot(x1(end-3:end,1),x1(end-3:end,2),'*','Color','r');
saveas(gcf,strcat(strcat('poster_',num2str(i)),'.png'))
end






