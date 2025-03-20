function [Schott] = MANOVA_Schott_2(X1,X2)
%UNTITLED4 
[p,n1]=size(X1);
[p,n2]=size(X2);

n=n1+n2;
X=[X1,X2];
Mean_X1=mean(X1');
Mean_X2=mean(X2');
Mean_s=mean(X');
H=n1*(Mean_X1-Mean_s)'*(Mean_X1-Mean_s)+n2*(Mean_X2-Mean_s)'*(Mean_X2-Mean_s);

Sn1=cov(X1');
Sn2=cov(X2');
E=(n1-1)*Sn1+(n2-1)*Sn2;

e=n-2;

tnp=(trace(H)-trace(E)/e)/(sqrt(n-1));

a=(trace(E^2)-(trace(E))^2/e)/((e+2)*(e-1));
Schott=tnp/sqrt(2*a/e);
end
