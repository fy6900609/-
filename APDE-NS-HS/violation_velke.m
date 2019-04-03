function [viol,H,G]=violation_velke(func,h,g,eps_viol,N)
% vraci sloupcovy vektor violation prislusny k maticim g,h, v nichz kazdy
% radek je vektorem g(resp.h) odpovidajim jednomu bodu
%
% h=h';
% g=g';

poctyhag=[1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;...
       1,1;0,2;0,3;1,1;1,1;1,1;1,1;1,2;0,2;0,2;...
       0,2;0,3;1,1;1,1;1,1;1,1;1,2;0,2];
H=zeros(N,6);
G=zeros(N,3);
mh=poctyhag(func,1);
mg=poctyhag(func,2);

sumH=zeros(N,1);
sumG=zeros(N,1);
if mh > 0
    for i=1:N
        h1=h(i,:);
        Hpom=abs(h1);
        Hpom(Hpom<=eps_viol)=0;
        H(i,1:mh)=Hpom;
        sumH(i,1)=sum(Hpom);
    end    
end    
if mg > 0
    for i=1:N
        g1=g(i,:);
        Gpom=g1;
        Gpom(Gpom<=0)=0;
        G(i,1:mg)=Gpom;
        sumG(i,1)=sum(Gpom);  
    end    
end    
viol=(sumH+sumG)/(mh+mg);
