function yp = IMEXRK_Matrix(y,N,S,d,locdt,A,Ahat,b,bhat,r)
Y = zeros(d,r);
Y(:,1)=y;
for i = 2:r
    temp=y;
    for j=1:i-1
        temp = temp + locdt *(A(i,j)*N + Ahat(i,j)*S)*Y(:,j); %%%
    end
    Y(:,i)=(speye(d) - locdt *Ahat(i,i)*S)\temp;
end
yp=y;
for i=1:r
    yp = yp+ locdt * (b(i)*N + bhat(i)*S)*Y(:,i); %%%
end
end