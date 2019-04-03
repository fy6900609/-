%  de_current_to_best binomial
function y=currenttopbestbintest(AR,velarch,PO,hodf,hodviol,p,F,CR,expt,a,b)
N=length(PO(:,1));        %个体数
d=length(PO(1,:));        %自变量维数 

pom=zeros(N,d+2);         %亲本个体信息矩阵
pom(:,1:d)=PO;            %基本信息
pom(:,d+1)=hodf;          %亲本目标函数值信息
pom(:,d+2)=hodviol;       %亲本约束违反值信息

pom=sortrows(pom,[d+2,d+1]);    %使各个体目标函数值升序排列


pbest=pom(1:p,1:d);       %提取前p个优秀个体基本信息
ktery=1+fix(p*rand(1));   %从p个优秀个体中随机选择一个个体
xpbest=pbest(ktery,:);
% prd1=size(PO)
% prd=expt(1)
xi=PO(expt(1),1:d);       %提取当前个体

vyb=nahvyb_expt(N,1,expt);%提取其他非当前个体的种群中其他个体
r1=PO(vyb,:);
expt=[expt,vyb];

vyb=nahvyb_expt(N+velarch,1,expt);%在当前种群和失败个体存档中提取第二个个体
sjed=[PO;AR];                     %总个体集
r2=sjed(vyb,:);

v=xi+F*(xpbest-xi)+F*(r1-r2);  %策略表达式
y=xi;
change=find(rand(1,d)<CR);     %交叉过程
if isempty(change) % at least one element is changed
    change=1+fix(d*rand(1));
end
y(change)=v(change);
y=zrcad_shade(y,xi,a,b);   %防止生成个体超越定义域
%
