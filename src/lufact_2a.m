% Clear memory and screen
clear
clc

N = 8;
a = -2;
b = 1;
c = 1;
A = diag(a*ones(1,N)) + diag(b*ones(1,N-1),1) + diag(c*ones(1,N-1),-1);
B = [1 1 1 1 1 1 1 1]';

disp(A)
disp(B)
disp('The augmented matrix : ')
D = horzcat(A,B);
disp(D);
lu(A,B);

function lu(A,b)
% Solve the system Ax=b using the LU decomposition.
n=length(b);
y=zeros(n,1);
x=zeros(n,1);
fprintf('\n');
for i=1:n
    U(i,i)=1;
end
L(1,1)=A(1,1)/U(1,1);
for j=2:n
    L(j,1)=A(j,1)/U(1,1);
    U(1,j)=A(1,j)/L(1,1);
end
for i=2:n-1
    S=0;
    for k=1:i-1
        S=S+U(k,i)*L(i,k);
    end
    L(i,i)=(A(i,i)-S)/U(i,i);
    for j=i+1:n
        S=0;
        for k=1:i-1
            S=S+U(k,i)*L(j,k);
        end
        L(j,i)=(A(j,i)-S)/U(i,i);
        S=0;
        for k=1:i-1
            S=S+U(k,j)*L(i,k);
        end
        U(i,j)=(A(i,j)-S)/L(i,i);
    end
end
S=0;
for k=1:n-1
    S=S+U(k,n)*L(n,k);
end
L(n,n)=(A(n,n)-S)/U(n,n);
% Perform the forward substitution.
y(1)=b(1)/L(1,1);
for i=2:n
    S=b(i);
    for j=1:i-1
        S=S-L(i,j)*y(j);
    end
    y(i)=S/L(i,i);
end
% Perform the back substitution.
x(n)=y(n)/U(n,n);
for i=n-1:-1:1
    S=y(i);
    for j=i+1:n
        S=S-U(i,j)*x(j);
    end
    x(i)=S/U(i,i);
end
disp('The L matrix :')
disp(L)
disp('The forward substitution gives :')
disp(y)
disp('The U matrix :')
disp(U)
disp('The vector solution is :')
disp(x)
end