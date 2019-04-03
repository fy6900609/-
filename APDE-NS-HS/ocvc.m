function  [ocvqs,ocvhs,E] = ocvc( A,B,C ) %CSq,CSh,CS
%计算违反约束程度的函数
n=length(C(1,:));%约束数
l1=length(A(:,1));%种群中个体数
ocvqs=zeros(l1,1);
ocvhs=zeros(l1,1);
M=[A;B;C];
C=max(M);
E=C;
for i=1:l1
a=0;
b=0;
c=0;
d=0;
F=zeros(1,n);
M=find(C~=0);
     for s=M
     F(1,s)=1/C(1,s);
     b=b+F(1,s); a=a+F(1,s)*A(i,s);
     end
ocvqs(i)= a/b;

R=find(C~=0);
        for r=R
            F(1,r)=1/C(1,r);
            c=c+F(1,r); d=d+F(1,r)*B(i,r);
        end
    ocvhs(i)=d/c;
end
