function[y, cf]=sso(B, A, x, ci);

cf=ci;
y=zeros(size(x));

    for i=1:length(x)

       N1=x(i)-A(2)*cf(1) -A(3)*cf(2); % A=x(i)-a1*M1 - a2*M2
       N2=N1*B(1)+B(2)*cf(1)+B(3)*cf(2);   % B=b0*A + b1*M1 + b2*M2
       %salida
       y(i)=N2;       %y(i)=B
       %actualización memoria
       cf(2)=cf(1);  %M2=M1;
       cf(1)=N1;     %M1=A;

    end

end
