function [SK] = MANOVA_SK_2(X1,X2)
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

bbH=n1*(Meanue1-Mean_sum)'*(Meanue1-Mean_sum)+n2*(Meanue2-Mean_sum)'*(Meanue2-Mean_sum);
bbB=(n1-1)*Snue1+(n2-1)*Snue2;

bf_Sn=bbB/e1;
Ds=diag(diag(bf_Sn));
bf_Rn=Ds^(-0.5)*bf_Sn*Ds^(-0.5);
c_pn=1+((trace(bf_Rn^2))/(p^(1.5)));

SK=(trace(bbH*Ds^(-1))-e1*p*(e1-2)^(-1))/sqrt(2*c_pn*(trace(bf_Rn^2)-p^2/e1));
end