%  de_current_to_best binomial
function y=currenttopbestbin_izrc_cons(AR,velarch,PO,hodf,hodviol,p,F,CR,expt,a,b,rdiv)
N=length(PO(:,1));
d=length(PO(1,:));
cc=0.6;%0.3;%
pom=zeros(N,d+2);
pom(:,1:d)=PO;
pom(:,d+1)=hodf;
% pom(:,d+2)=hodviol;

pom=sortrows(pom,d+1);


pbest=pom(1:p,1:d);
ktery=1+fix(p*rand(1));
xpbest=pbest(ktery,:);
% prd1=size(PO)
% prd=expt(1)
xi=PO(expt(1),1:d);

vyb=nahvyb_expt(N,1,expt);
r1=PO(vyb,:);
expt=[expt,vyb];

% vyb=nahvyb_expt(N,1,expt);
% r3=PO(vyb,:);
% expt=[expt,vyb];

vyb=nahvyb_expt(N+velarch,1,expt);
sjed=[PO;AR];
r2=sjed(vyb,:);
expt=[expt,vyb];

% vyb=nahvyb_expt(N+velarch,1,expt);
% sjed=[PO;AR];
% r4=sjed(vyb,:);
% v=xi+F*(xpbest-xi)+F*(r1-r2);
v=xi+F*(xpbest-xi)+cc*(1-rdiv)*(r1-r2);
y=xi;
change=find(rand(1,d)<CR);
if isempty(change) % at least one element is changed
    change=1+fix(d*rand(1));
end
y(change)=v(change);
% L=1+fix(d*rand(1));  % starting position for crossover
% change=L;
% position=L;
% while rand(1) < CR && length(change) < d
%     position=position+1;
%     if position <= d
%         change(end+1)=position;
%     else
%         change(end+1)=mod(position,d);
%     end
% end
% y(change)=v(change);
y=zrcad_shade(y,xi,a,b);


%
