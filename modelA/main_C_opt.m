clear all
clc
format long

%epsilon_vec = linspace(0,0.5,5); %epsilon is the probability of packet drop
epsilon_vec = [0 0.3 0.7];



S = 30; %Truncation size. I have denoted this by 'n' in my writeup
p = 0.3;
betaVec=[0.9];

err=10^(-4);

lambda = linspace(1,200,100);
k_opt = ones(length(epsilon_vec),length(lambda));
C_opt = zeros(length(epsilon_vec),length(lambda));
%alpha = linspace(0.001,0.99,50);
 d= zeros(2*S+1,1);
 ell= zeros(2*S+1,1);
K=1;



l_K_vec = [];
m_K_vec = [];

%m_opt = ones(length(epsilon_vec),length(alpha));
%lambda_cr = zeros(S,1);

%color = ['b','m','g','k','c'];
color = ['b', 'm', 'k'];

for b=1:length(betaVec)
    beta=betaVec(b);
    for ep=1:length(epsilon_vec)
        epsilon = epsilon_vec(ep);
        
        for k=1:S %thresholds
            
            %transition probability matrix of size 2*K+1 for birth-death markov chain
            P=zeros(2*S+1);
            P(1,:)=[1-2*p p zeros(1,2*S-1)];
            P(2*S+1,:)=[zeros(1,2*S-1) p 1-2*p];
            
            for J=2:2*S
                P(J,J-1)=p;
                P(J,J)=1-2*p;
                P(J,J+1)=p;
            end
            
            h=zeros(2*S+1,1);
            
            h(-k+S+2:k+S) = ones(2*k-1,1);
            h(1:-k+S+1)= epsilon*ones(S-k+1,1);
            h(k+S+1:2*S+1) = h(1:-k+S+1);
            
            P_had=hadamard_prod(h,P,k,S);
            
            %  d= zeros(2*S+1,1);
            
            for j=1:2*S+1
                if j<k+S+1 && j>-k+S+1
                    d(j) = abs(j-S-1); % per-step distortion func
                else d(j) = epsilon*abs(j-S-1);
                end
            end
            
            for j=1:2*S+1
                if j<k+S+1 && j>-k+S+1
                    ell(j) = 0; % per-step distortion func
                else ell(j) = 1;
                end
            end
            
            L_vec = (eye(2*S+1) - beta*P_had)^(-1)*d;
            M_vec = (eye(2*S+1) - beta*P_had)^(-1)*h;
            K_ell_vec = (eye(2*S+1) - beta*P_had)^(-1)*ell;
            
            l_add(k) = L_vec(S+1);
            m_add(k) = M_vec(S+1);
            K_ell_add(k) = K_ell_vec(S+1);
            
            trun_L(ep,k) = l_add(k);
            trun_M(ep,k) = m_add(k);
            trun_K_ell(ep,k) = K_ell_add(k) ;
            
            
            trun_D(ep,k) = l_add(k)/m_add(k);
            trun_N(ep,k) = K_ell_add(k)/m_add(k);
            
            
            %             for K=k:S
            %
            %
            %                 %transition probability matrix of size 2*K+1 for birth-death markov chain
            %                 P=zeros(2*K+1);
            %                 P(1,:)=[1-2*p p zeros(1,2*K-1)];
            %                 P(2*K+1,:)=[zeros(1,2*K-1) p 1-2*p];
            %
            %                 for J=2:2*K
            %                     P(J,J-1)=p;
            %                     P(J,J)=1-2*p;
            %                     P(J,J+1)=p;
            %                 end
            %
            %                 h=zeros(2*K+1,1);
            %
            %                 h(-k+K+2:k+K) = ones(2*k-1,1);
            %                 h(1:-k+K+1)= epsilon*ones(K-k+1,1);
            %                 h(k+K+1:2*K+1) = h(1:-k+K+1);
            %
            %                 P_had=hadamard_prod(h,P,k,K);
            %
            %                 d= zeros(2*K+1,1);
            %
            %               for j=1:2*K+1
            %                   if j<k+K+1 && j>-k+K+1
            %                       d(j) = abs(j-K-1); % per-step distortion func
            %                   else d(j) = epsilon*abs(j-K-1);
            %                   end
            %               end
            %               L_vec = (eye(2*K+1) - beta*P_had)^(-1)*d;
            %               M_vec = (eye(2*K+1) - beta*P_had)^(-1)*h;
            %
            %               l_add(k,K) = L_vec(K+1);
            %               m_add(k,K) = M_vec(K+1);
            %
            %             end
            %               trun_L(ep,k) = l_add(k,K);
            %               trun_M(ep,k) = m_add(k,K);
            %               trun_D(ep,k) = l_add(k,K)/m_add(k,K);
            %               trun_N(ep,k) = 1/m_add(k,K) - (1-beta);
            %
            %costly performance
            
            if k>=2
                lambda_cr(ep,k-1) = (trun_D(ep,k) - trun_D(ep,k-1))/(trun_N(ep,k-1) - trun_N(ep,k));
                %C(ep,k-1) = trun_D(ep,k) + lambda_cr(ep,k-1)* trun_N(ep,k-1);
            end
        end
    end
    %csvwrite(strcat('N_est_beta_', num2str(beta), '_trunSize_', num2str(S), '.csv'), trun_D);
end

for ep = 1:length(epsilon_vec)
    for ind = 1: length(lambda)
        for k = 1:S-2
            if k==1
                C_opt(ep,ind) = trun_D(ep,k) + lambda(ind)* trun_N(ep,k);
            end
            if lambda(ind) > lambda_cr(ep,k) && lambda(ind) <= lambda_cr(ep,k+1)
                k_opt(ep,ind) = k+1;
                C_opt(ep,ind) = trun_D(ep,k_opt(ep,ind)) + lambda(ind)* trun_N(ep,k_opt(ep,ind));
            end
        end
    end
end

for i=1:length(epsilon_vec)
    plot(lambda,C_opt(i,1:length(lambda))',color(i));
    hold on;
end


save('C_opt.mat','C_opt');
save('k_opt.mat','k_opt');
