%  de_current_to_best binomial
function y=currenttopbestbintest(AR,velarch,PO,hodf,hodviol,p,F,CR,expt,a,b)
N=length(PO(:,1));        %������
d=length(PO(1,:));        %�Ա���ά�� 

pom=zeros(N,d+2);         %�ױ�������Ϣ����
pom(:,1:d)=PO;            %������Ϣ
pom(:,d+1)=hodf;          %�ױ�Ŀ�꺯��ֵ��Ϣ
pom(:,d+2)=hodviol;       %�ױ�Լ��Υ��ֵ��Ϣ

pom=sortrows(pom,[d+2,d+1]);    %ʹ������Ŀ�꺯��ֵ��������


pbest=pom(1:p,1:d);       %��ȡǰp��������������Ϣ
ktery=1+fix(p*rand(1));   %��p��������������ѡ��һ������
xpbest=pbest(ktery,:);
% prd1=size(PO)
% prd=expt(1)
xi=PO(expt(1),1:d);       %��ȡ��ǰ����

vyb=nahvyb_expt(N,1,expt);%��ȡ�����ǵ�ǰ�������Ⱥ����������
r1=PO(vyb,:);
expt=[expt,vyb];

vyb=nahvyb_expt(N+velarch,1,expt);%�ڵ�ǰ��Ⱥ��ʧ�ܸ���浵����ȡ�ڶ�������
sjed=[PO;AR];                     %�ܸ��弯
r2=sjed(vyb,:);

v=xi+F*(xpbest-xi)+F*(r1-r2);  %���Ա��ʽ
y=xi;
change=find(rand(1,d)<CR);     %�������
if isempty(change) % at least one element is changed
    change=1+fix(d*rand(1));
end
y(change)=v(change);
y=zrcad_shade(y,xi,a,b);   %��ֹ���ɸ��峬Խ������
%
