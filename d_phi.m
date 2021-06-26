%phi(x)=(akpha2+ibeta2)U*exp(U)+sin(sqrt(1+Ux+Uy))
function [f]=d_phi(x,alpha2,beta2,N,h,dt)
f=diag((alpha2+1i*beta2)*(x+1).*exp(x));
n=N^2;
for i=1:N+2
    for j=1:N+2
        if i==1 ||i==N+2||j==1 ||j==N+2
            y(i,j)=0;
        else
            y(i,j)=x((i-2)*N+(j-1));
        end
    end
end
for i=1:n
        b1=mod(i,N);
        a1=fix(i/N)+1;
        if b1==0
           b1=N;
           a1=a1-1;
        end
        a=y(a1+1,b1+2)-y(a1+1,b1);
        b=y(a1+2,b1+1)-y(a1,b1+1);  
        c=2*h*sqrt(4*h^2+a^2+b^2);
        d=cos(sqrt(4*h^2+a^2+b^2)/(2*h));
        if i<n
            f(i,i+1)=d*a/c; 
            if mod(i,N)==0
            f(i,i+1)=0;
            end
        end        
        if i>1
            f(i,i-1)=-d*a/c;
            if mod(i,N)==1
            f(i,i-1)=0;
            end
        end
        if i<=n-N
            f(i,i+N)=d*b/c;
        end
        if i>N
            f(i,i-N)=-d*b/c;
        end
        
end
f=f*h*dt;
