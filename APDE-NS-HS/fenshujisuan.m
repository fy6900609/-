
Dimension=[9,16,18,10,10,10,10,10,10,10];
resave=10;
score=zeros(10,20);
for func_num=1:10
    D=Dimension(func_num);
eval(['load ' num2str(func_num) '_' num2str(D) '.mat finalvalue'])
    score1=sort(finalvalue);
for i=1:resave
        if score1(i)<1.000000001
            score(func_num,i)=10;
        elseif score1(i)<1.00000001
            score(func_num,i)=9;
        elseif score1(i)<1.0000001
            score(func_num,i)=8;
        elseif score1(i)<1.000001
            score(func_num,i)=7;
        elseif score1(i)<1.00001
            score(func_num,i)=6;
        elseif score1(i)<1.0001
            score(func_num,i)=5;
        elseif score1(i)<1.001
            score(func_num,i)=4;
        elseif score1(i)<1.01
            score(func_num,i)=3;
        elseif score1(i)<1.1
            score(func_num,i)=2;
        elseif score1(i)<2
            score(func_num,i)=1;
        else
            score(func_num,i)=0;
        end
end
end
scorefin1=sum(score,2)./resave;
scorefin=sum(scorefin1);