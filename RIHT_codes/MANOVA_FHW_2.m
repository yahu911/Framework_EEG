function [FHW] = MANOVA_FHW_2(X1,X2)
%UNTITLED4 
[p,n1]=size(X1);
[p,n2]=size(X2);

n=n1+n2;
Meanue1=mean(X1');
Snue1=cov(X1');
Meanue2=mean(X2');
Snue2=cov(X2');
X=[X1,X2];
Mean_sum=mean(X');

G=2;
e1=n-G;
e2=G-1;

bbH=n1*(Meanue1-Mean_sum)'*(Meanue1-Mean_sum)+n2*(Meanue2-Mean_sum)'*(Meanue2-Mean_sum);
bbB=(n1-1)*Snue1+(n2-1)*Snue2;

tilde_T_D=sqrt(p)*(e1*(trace(bbH))/(trace(bbB))-e2)*(trace(bbB))/(e1*p);
tilde_sigma_D=sqrt(2*e2*((trace(bbB^2))/(e1^2)-(trace(bbB))^2/(e1^3))/p);

FHW=tilde_T_D/tilde_sigma_D;
end