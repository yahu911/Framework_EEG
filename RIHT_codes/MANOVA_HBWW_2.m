function [HBWW] = MANOVA_HBWW_2(X1,X2)
%UNTITLED4 
[p,n1]=size(X1);
[p,n2]=size(X2);

n=n1+n2;
X=[X1,X2];
Mean_X1=mean(X1');
Mean_X2=mean(X2');

Tnk1=(Mean_X1-Mean_X2)*(Mean_X1-Mean_X2)';

Sn1=cov(X1');
Sn2=cov(X2');

Tnk2=(trace(Sn1)/n1+trace(Sn2)/n2);

Tnk=Tnk1-Tnk2;

Var_Tnk1_1=(n1-1)*(trace(Sn1^2)-(trace(Sn1))^2/(n1-1))/(n1*(n1+1)*(n1-2));
Var_Tnk1_2=(n2-1)*(trace(Sn2^2)-(trace(Sn2))^2/(n2-1))/(n2*(n2+1)*(n2-2));
Var_Tnk1=2*(Var_Tnk1_1+Var_Tnk1_2);

Var_Tnk2_1=trace(Sn1*Sn2);
Var_Tnk2=4*Var_Tnk2_1/(n1*n2);

HBWW=Tnk/sqrt(Var_Tnk1+Var_Tnk2);
end