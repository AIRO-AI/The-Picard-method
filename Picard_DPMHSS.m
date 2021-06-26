clear; clc;
% Ut-(alpha1+ibeta1)(Uxx+Uyy)+qU=(alpha2+ibeta2)U*exp(U)+sin(sqrt(1+Ux+Uy));
%% compare parameters
%eta=0.1;
%N=32;q=0.1;alpha=0.7;belta=1.8;
%N=32;q=1;alpha=0.7;belta=1.7;
%N=32;q=10;alpha=0.7;belta=1.7; 

%N=64;q=0.1;alpha=0.7;belta=1.7;
%N=64;q=1;alpha=0.7;belta=1.7;
%N=64;q=10;alpha=0.7;belta=1.7;

eta=0.01;
N=32;q=0.01;alpha=0.7;belta=1.7; 
%N=32;q=1;alpha=0.825;belta=0.827;
%N=32;q=10;alpha=0.808;belta=0.811;

%N=64;q=0.1;alpha=0.785;belta=0.785;
%N=64;q=1;alpha=0.784;belta=0.784;
%N=64;q=10;alpha=0.777;belta=0.778;

  %eta=0.001;
%N=32;q=0.1;alpha=0.827;belta=0.829;
%N=32;q=1;alpha=0.825;belta=0.827;
%N=32;q=10;alpha=0.808;belta=0.811;

%N=64;q=0.1;alpha=0.785;belta=0.785;
%N=64;q=1;alpha=0.784;belta=0.784;
%N=64;q=10;alpha=0.777;belta=0.778;

n=N^2; h=1/(N+1);dt=h;alpha1=1;beta1=1;alpha2=0.5;beta2=0.5;
%% problem describe
M=diag(2*ones(N,1)) + diag(-1*ones(N-1,1),1) + diag(-1*ones(N-1,1),-1);
IN=eye(N);In=eye(n);
A=kron(M,IN)+kron(IN,M);
W=h*(1+q*dt)*In+alpha1*dt/h*A;
T=beta1*A; 
A=W+1i*T;
V=W; aV=alpha*V;bV=belta*V;
[L1,U1]=lu(aV+W);
[L,U]=lu(bV+T);
%% give the initial value
X0=zeros(n,1);
F0=A*X0-phi(X0,alpha2,beta2,N,h,dt);

%% some coefficients
tmin=1000;
%   for belta=1.7
for m=1:1
%% initial steps
    outerror=1; l=0; k=0;
    X=X0; F=F0;
    tic
    %% iterative process
    while outerror>1e-6
        %===========================================
        Phi=phi(X,alpha2,beta2,N,h,dt);
        DF=A-d_phi(X,alpha2,beta2,N,h,dt);
        %===========================================
        X1=X;
        inerror=1;
        while inerror>eta
           b1=(aV-1i*T)*X+Phi;         
           X=U1\(L1\b1);
           b2=(bV+1i*W)*X-1i*Phi;
           X=U\(L\b2);
           s=X-X1;
           F1=DF*s+F;
           inerror1=norm(F1,2)/norm(F,2);
           if inerror1==inerror
           break;
           end
           inerror=inerror1;
           l=l+1;
        end  
        nl(k+1)=l;l=0;
        %  update the values of F and X_k 
        F=A*X-phi(X,alpha2,beta2,N,h,dt);
        outerror=norm(F,2)/norm(F0,2);
        k=k+1;
    end
    t=toc;
    if t<tmin
        tmin=t;
        sum=0;
       for i=1:k
           sum=nl(i)+sum;
       end       
    end   
end    

outerror
outIT=k
inIT=sum/k
IT=sum
tmin


