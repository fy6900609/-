function [vystup,fminvec] = APDENS(Dim,Max_Gen,xmin,xmax,fnum,eps_viol)%
rand('state',sum(100*clock));
%binomialni a exponencialni krizeni soutezi ale jen v jedne generaci
global initial_flag

N_init=32*Dim;%26*Dim;%26*Dim;%24*Dim;
N_min=10;
ps=N_init;
p=0.11;
max_velikost_archivu=round(ps*2);
G=0; 
pesek1=0;pesek2=0;pesek3=0;pesek4=0;pesek6=0;pfeas=0;pesek7=0;
Gsearch=[];
Psearch=[];
rdiv=1;
 maxiter=Max_Gen*Dim;
 evals=ps;
 D=Dim;
 H=4;
CFpool=[0.9 0.9;0.5 0.5;0.9 0.3;0.6 0.8];
%   CFpool=[0.9 0.9;0.7 0.3;0.9 0.2;0.6 0.8];
 h=4;
  LP=150;
 sm = zeros(LP, h); fm = zeros(LP, h);     %学习时期的成功/失败记录
 smp1=zeros(LP,H); fmp1=zeros(LP,H);
 smp2=zeros(LP,H); fmp2=zeros(LP,H);
 smp3=zeros(LP,H); fmp3=zeros(LP,H);
 smp4=zeros(LP,H); fmp4=zeros(LP,H);
 epsil = 0.01;                             %一个微小余量（防止分母为0）
 puse = zeros(1, h)+ 1/h;                  %初始各策略的使用概率
 pusep=zeros(h, H)+ 1/H;
%  n0=2;
%  delta=1/(5*h);
%  ni=zeros(1,h)+n0;

 Xmin=repmat(xmin,1,D);
 Xmax=repmat(xmax,1,D);
 Xmin=repmat(Xmin,ps,1);
 Xmax=repmat(Xmax,ps,1);
 
 pos1=Xmin+(Xmax-Xmin).*rand(ps,D);

 P=zeros(ps,D+11);
 P(:,1:D)=pos1;
 initial_flag = 0;
 P(:,D+1)=cec19_func(pos1',fnum);
 %UMDAc algorithm
 
%  mu=5*Dim;
%  total=26*Dim;
%  e=P(:,D+1);
%  pos1=P(:,1:D);
%  weights=log(mu+1/2)-log(1:mu)';
%  weights=weights/sum(weights);
%  [a1,a2]=sort(e);
%  bestval=a1(1);
%  bestvec=pos1(a2(1),:);
%  pos=pos1(a2(1:mu),:);
%  meanval=mean(pos);
%  stdval=std(pos);
%  for k=1:total
%     pos(k,:)=meanval+stdval.*randn(1,D);
%  end
%  for k=1:total
%            for j=1:D
%              while pos(k,j)>xmax
%                  pos(k,j)=meanval(j)+stdval(j).*randn;
%              end
%              while pos(k,j)<xmin
%                  pos(k,j)=meanval(j)+stdval(j).*randn;
%              end
%            end
%  end
%  cc1=0;
%  for kk=1:3
%      e=cec19_func(pos',fnum);
%      [a1,a2]=sort(e);
%      if a1(1)<bestval
%                  bestval=a1(1);
%                  bestvec=pos(a2(1),:);
%      end
%      newpos=pos(a2(1:mu),:);
%      meanval=(newpos(:,1:D)'*weights)';
%      stdval=1*std(newpos);
%      FV(kk)=a1(1);
%      
%      if kk>30 
%                if mod(kk,20)==0
%                   [aa1,aa2]=min(FV);
%                   if aa2<kk-20
%                      cc1=1;
%                   end
%                end
%      end
%      for k=1:total
%                  if cc1==1      %kk>300
%                     a=0.96*randn(1,D);
%                  else
%                      a=randn(1,D);
%                  end
%                  pos(k,:)=meanval+stdval.*a;
%      end
%      for k=1:total
%                  for j=1:D
% %                     if pos(k,j)>xmax
% %                        pos(k,j)=mod(pos(k,j),xmax);
% %                     elseif pos(k,j)<xmin
% %                        pos(k,j)=mod(pos(k,j),xmin);
% %                     end
%                       while pos(k,j)>xmax
%                             pos(k,j)=meanval(j)+stdval(j).*randn;
%                       end
%                       while pos(k,j)<xmin
%                             pos(k,j)=meanval(j)+stdval(j).*randn;
%                       end
%                  end
%      end
%      
%      
%      
%      
%  end
%  
%  P(:,1:D)=pos;
%  P(:,D+1)=cec19_func(pos',fnum);
 
 
 
 
 
%  gf=zeros(ps,1);
%  hf=zeros(ps,1);
%  [P(:,D+2),P(:,D+3:D+8),P(:,D+9:D+11)]=violation_velke(fnum,hf,gf,eps_viol,ps);
%  CS=max(P(:,D+3:D+11));
 rdiv0=computediv(P(:,1:D));
%  [violmin, indexviolmin]=min(P(:,D+2)); %minimalizuju violation
%  pomocna=P(indexviolmin,:);
%  fmin=min(pomocna(:,D+1));

 %pro ulozeni prumernych hodnot v pameti (vsechny MCR a vsechny MF)
 %pro ulozeni Fmin-Fmax v 10 etapach
  
 
 velarchivu=0;
 A=[];

 while (evals<maxiter)
     G=G+1;
     Gsearch=[Gsearch,G];   
    strategie=zeros(1,ps);
    poolie=zeros(1,ps);
    poskon=zeros(ps,D);
 
    for i=1:ps  %VYTVORENI DALSI GENERACE
         [hh]= roulet(puse);                %为个体分配变异策略
%         [hh,p_min]=roulete(ni);
%         if p_min<delta
%             ni=zeros(1,h)+n0;
%         end  %reset%        r=1+fix(H*rand(1));
        switch hh
          case 1 %(CURRENTTORAND/BIN)
           strategie(1,i)=1;
           r=roulet(pusep(1,:)) ;
          CR=CFpool(r,1)+0.1*randn;
            if CR>1
                CR=1;
            elseif CR<0
                CR=0;
            end
            F=CFpool(r,2)+0.1*randn;
            while F<=0  
%                 F=rand*pi-pi/2; 
                F=0.1 * randn + CFpool(r,2);  
            end
            if F>1
                F=1;
            end
%             p = pmin+ (0.2-pmin) * rand;
            ppoc=round(p*ps);
            if ppoc==0 
                ppoc=1;
            end
           
            y=currenttopbestbin_izrc_cons(A,velarchivu,P(:,1:D),P(:,D+1),P(:,D+2),ppoc,F,CR,i,xmin,xmax,rdiv);
            poskon(i,:)=y;


          case 2  %(CURRENTTORAND/EXP)
            strategie(1,i)=2;
            r=roulet(pusep(2,:)) ;
           CR=CFpool(r,1)+0.1*randn;
           if CR>1
                CR=1;
           elseif CR<0
                CR=0;
           end
            F=CFpool(r,2)+0.1*randn;
            while F<=0 
%                 F=rand*pi-pi/2;
                  F=CFpool(r,2)+0.1*randn; 
            end
            if F>1
                F=1;
            end
%             p = pmin+ (0.2-pmin) * rand;
            ppoc=round(p*ps);
            if ppoc==0 
                ppoc=1;
            end
            y=currenttopbestexp_izrc_cons(A,velarchivu,P(:,1:D),P(:,D+1),P(:,D+2),ppoc,F,CR,i,xmin,xmax,rdiv);
            poskon(i,:)=y;

            case 3  %(RANDRL/BIN)
           strategie(1,i)=3;
           r=roulet(pusep(3,:)) ;
           CR=CFpool(r,1)+0.1*randn;
            if CR>1
                CR=1;
            elseif CR<0
                CR=0;
            end
             F=CFpool(r,2)+0.1*randn;
            while F<=0  
%                 F=rand*pi-pi/2; 
                  F=CFpool(r,2)+0.1*randn;
            end
            if F>1
                F=1;
            end
            ppoc=round(p*ps);
            if ppoc==0 
                ppoc=1;
            end
            y=currenttopbestbintest(A,velarchivu,P(:,1:D),P(:,D+1),P(:,D+2),ppoc,F,CR,i,xmin,xmax);%derand_RLe_cons(P(:,1:D),P(:,D+1),P(:,D+2),F,CR,i,rdiv);%
            poskon(i,:)=zrcad(y,xmin,xmax);
            
            case 4  %(RANDRL/EXP)
            strategie(1,i)=4;
            r=roulet(pusep(4,:)) ;
            CR=CFpool(r,1)+0.1*randn;
           if CR>1
                CR=1;
           elseif CR<0
                CR=0;
           end
            F=CFpool(r,2)+0.1*randn;
            while F<=0  
%                 F=rand*pi-pi/2; 
                 F=CFpool(r,2)+0.1*randn;
            end
            if F>1
                F=1;
            end
            y=derandexp_RLe_cons(P(:,1:D),P(:,D+1),P(:,D+2),F,CR,i,rdiv,evals,maxiter);%currenttopbestbintest(A,velarchivu,P(:,1:D),P(:,D+1),P(:,D+2),ppoc,F,CR,i,xmin,xmax);%
            poskon(i,:)=zrcad(y,xmin,xmax);
        end            
      poolie(1,i)=r;  
    end
   Q=zeros(ps,D+11);
   Q(:,1:D)=poskon;
   initial_flag = 0;
   Q(:,D+1)=cec19_func(poskon',fnum);
%    gf=zeros(ps,1);
%    hf=zeros(ps,1);
%    [Q(:,D+2),Q(:,D+3:D+8),Q(:,D+9:D+11)]=violation_velke(fnum,hf,gf,eps_viol,ps);
%    CSq=P(:,D+3:D+11); CSh=Q(:,D+3:D+11);
%    [qs,hs,E]=ocvc(CSq,CSh,CS);
%    CS=E;
% zjisteni, jak jsou na tom prvky Q
   jak=zeros(ps); 
   for i=1:ps
       
         if Q(i,D+1)<= P(i,D+1)
              % nahrad - y je uspesny
              jak(i)=1; 
         end
       
       
   end
    for i=1:ps
       if jak(i)==1 
            if velarchivu < max_velikost_archivu
                A=[A;P(i,1:D)];
                velarchivu=velarchivu+1;
            else
                ktere=nahvyb(velarchivu,1);
                A(ktere,:)=P(i,1:D);
            end
            P(i,:)=Q(i,:);
%             ni(strategie(1,i))=ni(strategie(1,i))+1;         % zmena prsti qi
           if mod(G, LP) ==0               %目前刚好为LP整数代
           sm(LP,strategie(1,i))= sm(LP, strategie(1,i))+ 1;       %当前变异策略成功次数加一
           else
           sm(mod(G, LP), strategie(1,i)) =sm(mod(G, LP), strategie(1,i))+ 1; %超过存储子代数LP则对存档进行更新
           end
           switch strategie(1,i)
               case 1
                    if mod(G, LP) ==0              
                    smp1(LP,poolie(1,i))= smp1(LP, poolie(1,i))+ 1;       
                    else
                    smp1(mod(G, LP), poolie(1,i)) =smp1(mod(G, LP), poolie(1,i))+ 1; 
                    end
               case 2
                    if mod(G, LP) ==0              
                    smp2(LP,poolie(1,i))= smp2(LP, poolie(1,i))+ 1;       
                    else
                    smp2(mod(G, LP), poolie(1,i)) =smp2(mod(G, LP), poolie(1,i))+ 1;
                    end
               case 3
                    if mod(G, LP) ==0              
                    smp3(LP,poolie(1,i))= smp3(LP, poolie(1,i))+ 1;       
                    else
                    smp3(mod(G, LP), poolie(1,i)) =smp3(mod(G, LP), poolie(1,i))+ 1; 
                    end
               case 4
                    if mod(G, LP) ==0              
                    smp4(LP,poolie(1,i))= smp4(LP, poolie(1,i))+ 1;       
                    else
                    smp4(mod(G, LP), poolie(1,i)) =smp4(mod(G, LP), poolie(1,i))+ 1; 
                    end
           end
       else
           if mod(G, LP)== 0
                fm(LP, hh) = fm(LP, hh) +1;
           else
                fm(mod(G, LP), hh)= fm(mod(G, LP), hh)+ 1;
           end
           switch strategie(1,i)
               case 1
                    if mod(G, LP) ==0              
                    fmp1(LP,poolie(1,i))= fmp1(LP, poolie(1,i))+ 1;       
                    else
                    fmp1(mod(G, LP), poolie(1,i)) =fmp1(mod(G, LP), poolie(1,i))+ 1; 
                    end
               case 2
                    if mod(G, LP) ==0              
                    fmp2(LP,poolie(1,i))= fmp2(LP, poolie(1,i))+ 1;       
                    else
                    fmp2(mod(G, LP), poolie(1,i)) =fmp2(mod(G, LP), poolie(1,i))+ 1;
                    end
               case 3
                    if mod(G, LP) ==0              
                    fmp3(LP,poolie(1,i))= fmp3(LP, poolie(1,i))+ 1;       
                    else
                    fmp3(mod(G, LP), poolie(1,i)) =fmp3(mod(G, LP), poolie(1,i))+ 1; 
                    end
               case 4
                    if mod(G, LP) ==0              
                    fmp4(LP,poolie(1,i))= fmp4(LP, poolie(1,i))+ 1;       
                    else
                    fmp4(mod(G, LP), poolie(1,i)) =fmp4(mod(G, LP), poolie(1,i))+ 1; 
                    end
           end
       end
    end
    rdivf=computediv(P(:,1:D));
    rdiv=rdivf/rdiv0;
%     Feasn=find(P(:,D+2)==0);
%     Lfeasn=length(Feasn);
%     pfeas=Lfeasn/ps;
    evals=evals+ps;
    if G>= LP                              %进行下一代使用各变异策略的概率进化
        ssm= sum(sm);                      %计算每种策略在当前代的总有益进化次数
        sfm= sum(fm);                      %计算每种策略在当前代的总失败进化次数
        ssmp1=sum(smp1);sfmp1=sum(fmp1);
        ssmp2=sum(smp2);sfmp2=sum(fmp2);
        ssmp3=sum(smp3);sfmp3=sum(fmp3);
        ssmp4=sum(smp4);sfmp4=sum(fmp4);
        ksucc= zeros(1, h);                %在当前代各策略的成功系数
        ksuccp1= zeros(1,H);
        ksuccp2= zeros(1,H);
        ksuccp3= zeros(1,H);
        ksuccp4= zeros(1,H);
        for j= 1:h
            ksucc(j) = ssm(j)/((ssm(j)+ sfm(j))+ epsil)+epsil;%定义成功系数（用来进化各变异策略概率）
        end
        for j=1:H
            ksuccp1(j) = ssmp1(j)/((ssmp1(j)+ sfmp1(j))+ epsil)+epsil;
            ksuccp2(j) = ssmp2(j)/((ssmp2(j)+ sfmp2(j))+ epsil)+epsil;
            ksuccp3(j) = ssmp3(j)/((ssmp3(j)+ sfmp3(j))+ epsil)+epsil;
            ksuccp4(j) = ssmp4(j)/((ssmp4(j)+ sfmp4(j))+ epsil)+epsil;
        end
        puse= ksucc/sum(ksucc);            %下一代使用第k种策略的概率
        pusep(1,:)=ksuccp1/sum(ksuccp1);
        pusep(2,:)=ksuccp2/sum(ksuccp2);
        pusep(3,:)=ksuccp3/sum(ksuccp3);
        pusep(4,:)=ksuccp4/sum(ksuccp4);
%         save puse
        sm(mod(G,LP)+1, :)= zeros(1, h);
        fm(mod(G,LP)+1, :)= zeros(1, h);   %将下一代要使用的成功记忆区置0
        smp1(mod(G,LP)+1, :)= zeros(1, h);
        fmp1(mod(G,LP)+1, :)= zeros(1, h);
        smp2(mod(G,LP)+1, :)= zeros(1, h);
        fmp2(mod(G,LP)+1, :)= zeros(1, h);
        smp3(mod(G,LP)+1, :)= zeros(1, h);
        fmp3(mod(G,LP)+1, :)= zeros(1, h);
        smp4(mod(G,LP)+1, :)= zeros(1, h);
        fmp4(mod(G,LP)+1, :)= zeros(1, h);
    end
    
    ps_minule=ps;
%      if evals<=0.5*maxiter
%         ps=round(((N_min-N_init)/(1.3*maxiter))*evals+N_init);
%     else
%         if pfeas>0.5&&pesek3==0
%         ps=round(((N_min-N_init)/(1.3*maxiter))*evals+N_init);
%         else
%             if pesek6==0
%                 psfix=ps_minule;
%                 pesek6=pesek6+1;
%                 evalsfix=evals;
%             end 
%             if pesek6==1&&pesek7==0&&pesek3==1
%                 psfix=ps_minule;
%                 pesek7=pesek7+1;
%                 evalsfix=evals;
%             end
%         if pfeas==0&&pesek3==1
%             ps=40;%round(N_init/6);
%         else
%             if 0<pfeas&&pfeas<0.5
%                 pfeas=0.5;
%             end
%             ps=round(psfix-((psfix-N_min)/(maxiter-evalsfix))*(evals-evalsfix)*pfeas);
%         end
%         end
%     end
%     if (evals> 0.9*maxiter)&&(pesek3==0)&&(ps>round(40))
%         ps=round(40);
%         pesek3=pesek3+1;
%     end
    ps=round(((N_min-N_init)/(maxiter))*evals+N_init);
    Psearch=[Psearch,ps];
    while ps>ps_minule
        ps=ps-1;    
    end
    if ps<ps_minule
        P=sortrows(P,[D+2,D+1]);
       %minimalizuju violation ale stejne tridim i podle fmin
        P=P(1:ps,:);
        max_velikost_archivu=round(ps*2);
        while velarchivu > max_velikost_archivu
            index_v_arch=nahvyb(velarchivu,1);
            A(index_v_arch,:)=[];
            velarchivu=velarchivu-1;
        end
    end

if (evals>=0.1*maxiter)&&(pesek1==0)
%     violmin1=min(P(:,D+2)); %minimalizuju violation
%     row=find(P(:,D+2)==violmin1);
%     pomocna=P(row,:);
    [fmin1,indexfmin]=min(P(:,D+1));
%     bodmin1=P(indexfmin,:);
    pesek1=1;
end
if (evals>=0.5*maxiter)&&(pesek2==0)
%     violmin2=min(P(:,D+2)); %minimalizuju violation
%     row=find(P(:,D+2)==violmin2);
%     pomocna=P(row,:);
    [fmin2,indexfmin]=min(P(:,D+1));
%     bodmin2=pomocna(indexfmin,:);
    pesek2=1;
end
      
 end
% violmin=min(P(:,D+2)); %minimalizuju violation
% row=find(P(:,D+2)==violmin);
% pomocna=P(row,:);
[fmin,indexfmin]=min(P(:,D+1));
fminvec=P(indexfmin,1:D);
% bodmin=pomocna(indexfmin,:);


% c1=spocitejc(bodmin1);
% c2=spocitejc(bodmin2);
% c=spocitejc(bodmin);
% por1=length(find(bodmin1(D+3:D+11)>0));
% por2=length(find(bodmin2(D+3:D+11)>0));
% por=length(find(bodmin(D+3:D+11)>0));
evals
% violmin 
fmin
vystup=[evals  fmin1 fmin2 fmin];


% plot(Gsearch,Psearch);
% xlabel('Iteration');
% ylabel('ps');
% drawnow;

