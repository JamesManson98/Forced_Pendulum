function [lambda,Rdiag,x]=LyapunovQR(M,xini,N)
    xi=xini;    
    n=length(xi);   %we find the length of xi=xinitial
    Q=eye(n);       
    Rdiag=zeros(n,N);
    x=zeros(n,N);       %we preallocate Rdiag, x and lambda
    lambda=zeros(n,1);
    for i=1:N
        A=MyJacobian(M,xi,1e-6);    %we find A=partial(M(xi))
        [Q,Rnew]=qr(A*Q);           %we use QR decomposition
        r=diag(Rnew);               %we set r as the diagonal of Rnew
        for k=1:n
            if r(k)<0               %if r(k) is negative, we flip the sign of r(k) and Q(:,k)
                r(k)=-r(k);
                Q(:,k)=-Q(:,k);
            end
        end
        Rdiag(:,i)=r;   %we set Rdiag(:,i)=r;
        x(:,i)=xi;      %we input xi into the vector x for output
        xi=M(xi);       %we find the next xi by putting it into the map
    end
    for j=1:n
        lambda(j)=(1/N)*sum(log(Rdiag(j,:)));   %we then find all the lambdas 
    end
end