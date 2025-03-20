clear all
tic
n1=200;
n2=200;
n=n1+n2;
p =150;

Theta_01=xlsread('Theta_generate_p150_low_sparsity.csv');
Theta_02=Theta_01;

theta1_ori=xlsread('Opt_thetahat_X1_p150_nk200_sparsity1.0.csv');
theta2_ori=xlsread('Opt_thetahat_X2_p150_nk200_sparsity1.0.csv');
X1_Ori=xlsread('Opt_sample_X1_p150_nk200_sparsity1.0.csv');
X1_Ori=X1_Ori(:,2:151);
X2_Ori=xlsread('Opt_sample_X2_p150_nk200_sparsity1.0.csv');
X2_Ori=X2_Ori(:,2:151);

N=500;
k1=1;
k2=1;
T_FGL=0;
R_FGL=zeros(1,N);
for i=1:N
  i
    X1=X1_Ori((i-1)*200+1:i*200,:);
    Sigma_hat1 = (X1.' * X1) / n1;
    X2=X2_Ori((i-1)*200+1:i*200,:);
    Sigma_hat2 = (X2.' * X2) / n2;
        
    Theta1_hat=theta1_ori((i-1)*150+1:i*150,:);
    Theta2_hat=theta2_ori((i-1)*150+1:i*150,:);
    
    Var=Theta1_hat(k1,k1)*Theta1_hat(k2,k2)+Theta1_hat(k1,k2)^2+Theta2_hat(k1,k1)*Theta2_hat(k2,k2)+Theta2_hat(k1,k2)^2;
    T_Theta=2*Theta1_hat-Theta1_hat*Sigma_hat1*Theta1_hat-(2*Theta2_hat-Theta2_hat*Sigma_hat2*Theta2_hat)-(Theta_01-Theta_02);
    
    R_FGL(1,i)=T_Theta(k1,k2)*sqrt(n/2)/sqrt(Var);
    if abs(R_FGL(1,i))>1.96;
        T_FGL=T_FGL+1;
    end
end

size=T_FGL/N

[fi,xi] = ecdf(R_FGL); 
ecdfhist(fi,xi,10);
axis([-4 4 0 0.5])
hold on
x = [-4:0.001:4];
y=normpdf(x,0,1);
plot(x,y,'r','linewidth',3)
hold on
xlabel('$T_{(1,1)}/\hat{\sigma}_{(1,1)}$','interpreter','latex')
ylabel('Density') 

% qqplot(R_FGL)

csvwrite('Figure_Opt_p150_nk200_sparsity_1_1.csv',R_FGL);

t=toc 