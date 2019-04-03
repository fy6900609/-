clear all;
global initial_flag
initial_flag = 0;
Dimension=[9,16,18,10,10,10,10,10,10,10];
Xmax = 100 + zeros(10,1);
% Xmax(4) = 10;
% Xmax(5) = 10;
% Xmax(6) = 20;
% Xmax(7) = 50;
% Xmax(9) = 10;
% Xmax(19) = 50;
% Xmax(28) = 50;
Xmax(1)=8192;
Xmax(2)=16384;
Xmax(3)=4;

Max_const_pro_FES=30000;

    tabulka=[];
    tabulka2=[];
    vystup2=[];
    vec=[];
    resave=25;
    score=zeros(10,resave);
    runs=50;
    strat=1;
    ftrx=[1 2 3 4 5 6 7 8 9 10];%[4];%[5];%
for func_num=1:10
    func_num1=ftrx(func_num);
    vystup2=[];
    func_num
    if func_num1==4
        Max_const_pro_FES=30000;
    end
    a=-Xmax(func_num1,1);
    b=Xmax(func_num1,1);
    eps_viol=0.0001;
    D=Dimension(func_num1);
     
    
    for j= 1:runs
        vec=[];
        [vystup,fvec] = LSHADE44const(D,Max_const_pro_FES,a,b,func_num1,eps_viol);
        vystup=[func_num j vystup];
        tabulka=[tabulka;vystup];
%         vec=[vec;fvec]

    end
        save vec1 vec
        clear vec
    tabulka1=sortrows(tabulka((func_num-strat)*runs+1:(func_num-strat)*runs+runs,:),6);
    tabulka1=tabulka1(1:resave,:);
    Best10=tabulka1(1,6);
    Worst10=tabulka1(end,6);
    Mean10=mean(tabulka1(1:end,6));%mean(tabulka1((func_num-strat)*runs+1:(func_num-strat)*runs+runs,6));
    Median10=median(tabulka1(1:end,6));
    Std10=std(tabulka1(1:end,6));%std(tabulka1((func_num-strat)*runs+1:(func_num-strat)*runs+runs,6));
%     Meanv=mean(tabulka((func_num-strat)*runs+1:(func_num-strat)*runs+runs,16));
%     feasn=find(tabulka((func_num-strat)*runs+1:(func_num-strat)*runs+runs,16)==0);
%     feasnum=length(feasn);
%     pfeas=feasnum/runs;
    score1=tabulka1(:,6);
    score1=score1';
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
    scorefin1=sum(score,2)./resave;
%     scorefin=sum(scorefin1);
    vystup2=[func_num Best10  Worst10 Mean10 Median10 Std10];
    tabulka2=[tabulka2;vystup2];
    save rua
end
scorefin=sum(scorefin1);
tabulka


    soubor=strcat('LSHADE44const_i_s_c_a_porusenim ',num2str(D),'_behu',num2str(runs),'.txt');
    fid = fopen (soubor, 'wt');
    fprintf(fid,' %14.11g   %14.11g   %14.11g  %14.11g   %14.11g   %14.11g \n', tabulka');
    fprintf(fid,' %14.11g   %14.11g   %14.11g  %14.11g   %14.11g   %14.11g \n',tabulka2');
    fprintf(fid,' %14.11g   %14.11g   %14.11g  %14.11g   %14.11g   %14.11g   %14.11g   %14.11g  %14.11g  %14.11g  %14.11g   %14.11g   %14.11g  %14.11g   %14.11g   %14.11g   %14.11g   %14.11g  %14.11g  %14.11g  %14.11g   %14.11g   %14.11g  %14.11g  %14.11g \n',score');
    fprintf(fid,' %14.11g   %14.11g   %14.11g  %14.11g   %14.11g   %14.11g   %14.11g   %14.11g  %14.11g  %14.11g \n',scorefin1');
    fprintf(fid,' %14.11g \n',scorefin');
    fclose(fid);

    
    clear all;

