function [ res ] = roulet( cutpoints )
%������Ϊ�����������������
h=length(cutpoints);
cp(1)= cutpoints(1);
for i=2:h
    cp(i)=cp(i-1)+cutpoints(i);
end
res= 1+ fix(sum(cp< rand(1)));
end

