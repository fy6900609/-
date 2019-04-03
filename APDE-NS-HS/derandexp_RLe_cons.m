%  de_rand exponential RL
%
function y=derandexp_RLe_cons(P,hodf,hodviol,F,CR,expt,rdiv,evals,maxiter)
cc=0.6;%0.5;%
n=5;
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
    r4=trivybrane(4,1:d);
    r5=trivybrane(5,1:d);
else
    r2=trivybrane(3,1:d);
    r3=trivybrane(2,1:d);
    r4=trivybrane(5,1:d);
    r5=trivybrane(4,1:d);
end
 if evals<0.2*maxiter
     v=r1+F*(r2-r3);
 else
 v=r1+F*(r2-r3)+cc*(1-rdiv)*(r4-r5);%r1+F*(r2-r3);%+cc*(1-rdiv)*(r1-r2);%
 end
%  change=find(rand(1,d)<CR);
% if isempty(change) % at least one element is changed
%     change=1+fix(d*rand(1));
% end
% y(change)=v(change);
L=1+fix(d*rand(1));  % starting position for crossover
change=L;
position=L;
while rand(1) < CR && length(change) < d
    position=position+1;
    if position <= d
        change(end+1)=position;
    else
        change(end+1)=mod(position,d);
    end
end
y(change)=v(change);
%