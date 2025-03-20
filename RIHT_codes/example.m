clear all % Empirical Power
tic
n1=100;
% n1=75;
% n1=100;
n2=20;
n=n1+n2;
p=63;

mu=zeros(p,n);
mu1=mu(:,1:n1);

Size_All_In=zeros(5,21);
Size_All_St=zeros(5,21);
hat_j=(0.003)*(0:20);
for kk=1:21
    kk
    mu_in=(-1).^((1:p)').*(sqrt(hat_j(kk))*randn(p,1));
    mu2=repmat(mu_in,1,n2);
    
    % hat_j=0.2;
    % Inter_5=round(p*0.05);
    % mu2=repmat([hat_j*ones(Inter_5,1);zeros(p-Inter_5,1)],1,n2);
    % mu3=repmat([zeros(Inter_5,1);-hat_j*ones(Inter_5,1);zeros(p-2*Inter_5,1)],1,n3);
    
    Sigma_In=eye(p);
    I=eye(p);
    Sigma_St=0.3*I+0.7*ones(p,p);
    
    N=5000;
    R_RIHT_In=zeros(1,N);
    T_RIHT_In=0;
    R_FHW_In=zeros(1,N);
    T_FHW_In=0;
    R_S_In=zeros(1,N);
    T_S_In=0;
    R_SK_In=zeros(1,N);
    T_SK_In=0;
    R_HBWW_In=zeros(1,N);
    T_HBWW_In=0;
    
    R_RIHT_St=zeros(1,N);
    T_RIHT_St=0;
    R_FHW_St=zeros(1,N);
    T_FHW_St=0;
    R_S_St=zeros(1,N);
    T_S_St=0;
    R_SK_St=zeros(1,N);
    T_SK_St=0;
    R_HBWW_St=zeros(1,N);
    T_HBWW_St=0;
    for i=1:N
        Asign1=sign(rand(p,n1)-0.5);
        Aexp1=0.5*exprnd(sqrt(2),p,n1);
        x_unif1=unifrnd(-sqrt(3),sqrt(3),p,n1);
        x1=(0.7827*x_unif1+0.6224*Asign1.*Aexp1);
        Asign2=sign(rand(p,n2)-0.5);
        Aexp2=0.5*exprnd(sqrt(2),p,n2);
        x_unif2=unifrnd(-sqrt(3),sqrt(3),p,n2);
        x2=(0.7827*x_unif2+0.6224*Asign2.*Aexp2);
        
        X1_In=(Sigma_In)^(0.5)*x1+mu1;
        X2_In=(Sigma_In)^(0.5)*x2+mu2;
        
        X1_St=(Sigma_St)^(0.5)*x1+mu1;
        X2_St=(Sigma_St)^(0.5)*x2+mu2;
        
        R_RIHT_In(1,i)=MANOVA_RIHT_2(X1_In,X2_In,1);
        if abs(R_RIHT_In(1,i))>1.96;
            T_RIHT_In=T_RIHT_In+1;
        end
        R_RIHT_St(1,i)=MANOVA_RIHT_2(X1_St,X2_St,1);
        if abs(R_RIHT_St(1,i))>1.96;
            T_RIHT_St=T_RIHT_St+1;
        end
        
        R_FHW_In(1,i)=MANOVA_FHW_2(X1_In,X2_In);
        if abs(R_FHW_In(1,i))>1.96;
            T_FHW_In=T_FHW_In+1;
        end
        R_FHW_St(1,i)=MANOVA_FHW_2(X1_St,X2_St);
        if abs(R_FHW_St(1,i))>1.96;
            T_FHW_St=T_FHW_St+1;
        end
        
        R_S_In(1,i)=MANOVA_Schott_2(X1_In,X2_In);
        if abs(R_S_In(1,i))>1.96;
            T_S_In=T_S_In+1;
        end
        R_S_St(1,i)=MANOVA_Schott_2(X1_St,X2_St);
        if abs(R_S_St(1,i))>1.96;
            T_S_St=T_S_St+1;
        end
        
        R_SK_In(1,i)=MANOVA_SK_2(X1_In,X2_In);
        if abs(R_SK_In(1,i))>1.96;
            T_SK_In=T_SK_In+1;
        end
        R_SK_St(1,i)=MANOVA_SK_2(X1_St,X2_St);
        if abs(R_SK_St(1,i))>1.96;
            T_SK_St=T_SK_St+1;
        end
        
        R_HBWW_In(1,i)=MANOVA_HBWW_2(X1_In,X2_In);
        if abs(R_HBWW_In(1,i))>1.96;
            T_HBWW_In=T_HBWW_In+1;
        end
        R_HBWW_St(1,i)=MANOVA_HBWW_2(X1_St,X2_St);
        if abs(R_HBWW_St(1,i))>1.96;
            T_HBWW_St=T_HBWW_St+1;
        end
    end
    
    Size_All_In(1,kk)=T_RIHT_In/N;
    Size_All_In(2,kk)=T_FHW_In/N;
    Size_All_In(3,kk)=T_S_In/N;
    Size_All_In(4,kk)=T_SK_In/N;
    Size_All_In(5,kk)=T_HBWW_In/N;
    
    Size_All_St(1,kk)=T_RIHT_St/N;
    Size_All_St(2,kk)=T_FHW_St/N;
    Size_All_St(3,kk)=T_S_St/N;
    Size_All_St(4,kk)=T_SK_St/N;
    Size_All_St(5,kk)=T_HBWW_St/N;
end

csvwrite('DenseAl_120_In.csv',Size_All_In);
csvwrite('DenseAl_120_St.csv',Size_All_St);

t=toc