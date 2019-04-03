%  de_rand binomial
%   RL random local, r1 by tournament, Kaelo, Ali 2007
function y=derand_RLe_cons(P,hodf,hodviol,F,CR,expt,rdiv)
cc=0.3;
n=3;
N=length(P(:,1));
d=length(P(1,:));
% prd1=size(P)
% prd=expt(1)
y=P(expt(1),1:d);
vyb=nahvyb_expt(N,n,expt);	% three random points without expt

r12345=P(vyb,:);
hodf12345=hodf(vyb);
hodviol12345=hodviol(vyb);

trivybrane=[r12345 hodf12345 hodviol12345];
trivybrane=sortrows(trivybrane,[d+2,d+1]);

r1=trivybrane(1,1:d);
if rand  < 0.5
    r2=trivybrane(2,1:d);
    r3=trivybrane(3,1:d);
%      r4=trivybrane(4,1:d);
%     r5=trivybrane(5,1:d);
else
    r2=trivybrane(3,1:d);
     r3=trivybrane(2,1:d);
%     r4=trivybrane(4,1:d);
%     r5=trivybrane(4,1:d);
end

v=r1+F*(r2-r3);%+cc*(rdiv)*(r1-r4);
change=find(rand(1,d)<CR);
if isempty(change) % at least one element is changed
    change=1+fix(d*rand(1));
end
y(change)=v(change);
%