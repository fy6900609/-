function y = computediv( A )
D=length(A(1,:));
ps=length(A(:,1));
B=mean(A);
div1=0;
for i=1:ps
    for j=1:D
        div1=div1+(A(i,j)-B(j)).^2;
    end

end
y=sqrt(div1/ps);
end

