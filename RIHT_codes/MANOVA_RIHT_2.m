function [RIHT] = MANOVA_RIHT_2(X1,X2,lambda)
%UNTITLED2 
[p,n1]=size(X1);
[p,n2]=size(X2);
n=n1+n2;
Ip=eye(p);
gamma=p/n;

Meanue1=mean(X1');
Snue1=cov(X1');
Meanue2=mean(X2');
Snue2=cov(X2');
Scom=((n1-1)*Snue1+(n2-1)*Snue2)/(n1+n2-2);
X=[X1,X2];
Mean_sum=mean(X');

RIHT1=n1*(Meanue1-Mean_sum)*(Scom+lambda*Ip)^(-1)*(Meanue1-Mean_sum)';
RIHT2=n2*(Meanue2-Mean_sum)*(Scom+lambda*Ip)^(-1)*(Meanue2-Mean_sum)';
RIHT0=RIHT1+RIHT2;

RIHT_mFnp=trace((Scom+lambda*Ip)^(-1))/p;
Theta1=(1-lambda*RIHT_mFnp)/(1-gamma*(1-lambda*RIHT_mFnp));
FIHT_mF_deriva=trace((Scom+lambda*Ip)^(-2))/p;
Theta21=(1-lambda*RIHT_mFnp)/((1-gamma+gamma*lambda*RIHT_mFnp)^3);
Theta22=lambda*(RIHT_mFnp-lambda*FIHT_mF_deriva)/((1-gamma+gamma*lambda*RIHT_mFnp)^4);
Theta2=Theta21-Theta22;

RIHT=sqrt(p)*(RIHT0/p-Theta1)/sqrt(2*Theta2);
end

