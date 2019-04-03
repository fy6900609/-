function [ res ] = roulet( cutpoints )
%按概率为各个个体分配变异策略
h=length(cutpoints);
cp(1)= cutpoints(1);
for i=2:h
    cp(i)=cp(i-1)+cutpoints(i);
end
res= 1+ fix(sum(cp< rand(1)));
end

