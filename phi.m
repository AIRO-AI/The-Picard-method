%phi(x)=(akpha2+ibeta2)U*exp(U)+sin(sqrt(1+Ux+Uy))
function [f]=phi(X,alpha2,beta2,N,h,dt)
f=(alpha2+1i*beta2)*X.*exp(X);
n=N^2;
%将X的n维列向量分解为(N+2)*(N+2)矩阵
for i=1:N+2
    for j=1:N+2
        if i==1 ||i==N+2||j==1 ||j==N+2
            y(i,j)=0;          
        else
            y(i,j)=X((i-2)*N+(j-1));
        end
    end
end
%计算sin(sqrt(1+Ux+Uy))
for k=1:n
    j=mod(k,N); 
    i=fix(k/N)+1;
    if j==0
        j=N;
        i=i-1;
    end
    %%a=y(i+1,j+1);
    b=y(i+2,j+1)-y(i,j+1);    
    c=y(i+1,j+2)-y(i+1,j);
    f(k)=f(k)+sin(sqrt(1+(b^2+c^2)/(4*h^2)));
end   
f=f*h*dt;

