clear
clc

x=[1 -2 0];
y=[0 21 1];

% Create matrix from equation
xLength = length(x);
A = zeros(xLength);
for i=1:xLength
    A(i, 1)=x(i)^2;
    A(i, 2)=x(i);
    A(i, 3)=1;
end

b=y;
disp(A)
disp(b')
solution = ngaussel(A,b');

% Print parabola equation %
disp('The parabola equation solution is')
for i=1:xLength
    fprintf('a%d = %d\n', i, solution(i));
end

% Create parabola from input given
p=polyfit(x,y,(xLength - 1));
t = linspace(x(1),x(xLength));
yhat = polyval(p,t);
a1 = trapz(t,yhat);

figure
plot(x,y,'o');
hold on
plot(t,yhat)
legend('x,y','Equation')
hold off

function x = ngaussel(A,b)
    % Solve the system Ax=b using naive gaussian elimination
    n=length(b);
    x=zeros(n,1);
    fprintf('\n');
    disp('The augmented matrix is')
    augm=[A b];
    disp(augm)
    for k=1:n-1
        for i=k+1:n
            m=A(i,k)/A(k,k);
            for j=k+1:n
                A(i,j)=A(i,j)-m*A(k,j);
            end
            A(i,k)=m;
            b(i)=b(i)-m*b(k);
        end
    end
    x(n)=b(n)/A(n,n);
    for i=n-1:-1:1
        S=b(i);
        for j=i+1:n
            S=S-A(i,j)*x(j);
        end
        x(i)=S/A(i,i);
    end
    % Print the results
    fprintf('\n');
    disp('The transformed upper triangular augmented matrix C is =')
    fprintf('\n');
    for i=1:n
        for j=1:n
            if (j<i)
                A(i,j)=0; 
            end
        end
    end
    C=[A b];
    disp(C)
    fprintf('\n');
    disp('The vector solution is =')
    disp(x)
end