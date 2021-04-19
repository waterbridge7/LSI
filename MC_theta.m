%% MCMC
clear;clc 
global lafa beta w
%% 为下一轮xIn xOut 采样
% load record
% d=[record_d(end,1:2);record_d(end,3:4);record_d(end,5:6);record_d(end,7:8)];
% d;gd=g(d(:,1),d(:,2));x1=[x1;d];gx1=[gx1;gd];save x1_gx1 x1 gx1
load x1_gx1
x1=x1(1:10,:);gx1=gx1(1:10,:);
w0=rand(6,1); w=nlinfit(x1,gx1,@Reg,w0);

global cal_SS 
cal_SS=calCov_K(x1,x1,lafa,beta)\eye(size(x1,1));
delta = .5;
pdf1 = @(xx) thtaYD(xx,x1,gx1,w,lafa,beta);
proppdf = @(xx,y) 1*normrnd(xx,10);
proprnd = @(xx) xx+rand(size(xx))*2*delta-delta ;   
% proprnd = @(xx) normrnd(0,1,[1,2]) ;  % 貌似表现更好
nsamples =3*10^4;% 尽量得到2W左右的采样点
xx = mhsample([1,1],nsamples,'pdf',pdf1,'proprnd',proprnd,'symmetric',1);

pdf1 = @(x) thtaYD(x,x1,gx1,w,lafa,beta);
proppdf = @(x,y) 1*normrnd(x,10);
proprnd = @(x) x+rand(size(x))*2*delta-delta ;   
x = mhsample([-2,-2],nsamples,'pdf',pdf1,'proprnd',proprnd,'symmetric',1);
x=[x;xx];
x = unique(x,'rows');
save x x
figure(3)
plot(x(:,1),x(:,2),'.','color',[.5 .5 .5],'MarkerSize',4)
hold on
plot(x1(:,1),x1(:,2),'*','Color','red');
% hold on 
% plot(x1(1:4,1),x1(1:4,2),'*','Color','red');
% hold on
% plot(x1(5:8,1),x1(5:8,2),'x','Color','r');
% hold on
% plot(x1(9:12,1),x1(9:12,2),'+','Color','r');
% hold on
% plot(x1(13:16,1),x1(13:16,2),'rs','Color','r');
hold on



% %%  画等高线
[X,Y]=meshgrid(-8:.1:8);
Z=g(X,Y);
contour3(X,Y,Z,[0,0]);
