    % calculate covraiance
function X=calCov_K(X1,X2,lafa,beta)
% % % 要求X1 and X2 的格式分别为 X1=[1 1;2 2]; X2=[3 3;2 2;1 1];
for i=1:size(X2,1)
    x=X1-[X2(i,1).*ones(size(X1,1),1),X2(i,2).*ones(size(X1,1),1)];
    X(:,i)=lafa*exp(-(x(:,1).^2+x(:,2).^2)/beta);
end

% % % 要求X1 and X2 的格式分别为 X1=[1 1;2 2]; X2=[3 3;2 2;1 1];
% Xone=ones(size(X1,1),1);
% X=zeros(size(X1,1),size(X2,1));
% for i=1:size(X2,1)
%     x=X1-Xone*X2(i,:);
% % % % x=X1-[X2(i,1).*Xone,X2(i,2).*Xone,X2(i,3).*Xone,X2(i,4).*Xone,X2(i,5).*Xone,X2(i,6).*Xone,X2(i,7).*Xone,X2(i,8).*Xone];
% %     ABS_X=abs(x);
% %     X(:,i)=lafa*exp(-(ABS_X(:,1)/beta(1)+ABS_X(:,2)/beta(2))); %.^2
% 
%     X(:,i)=lafa*exp(-sum(abs(x),2)/beta(1)); %.^2
% end