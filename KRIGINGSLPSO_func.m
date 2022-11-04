function [gbest,gbestval,evacount,estcount,Itergbest,itereva,internum,convergegbest,convergeeva]= KRIGINGSLPSO_func(Dimension,Max_Gen,VRmin,VRmax,varargin)
%[gbest,gbestval,fitcount]= PSO_func('f8',3500,200000,30,30,-5.12,5.12)
rand('state',sum(100*clock));

maxEVA=Max_Gen;%z最大迭代次数
%ps=Particle_Number;%种群个数
D=Dimension;
estcount=0;%估计的总次数
func_num=problem;

%   alpha=0.5;%%样本选取范围
%%%%初始slpso参数
 M = 100;
    m = M + D/10;
    c3 = D/M*0.01;
    PL = zeros(m,1);

    for i = 1 : m
        PL(i) = (1 - (i - 1)/m)^log(sqrt(ceil(D/M)));%学习概率
    end
 ps=m;  
 fitest=zeros(ps,1);%每一代的对应的个体被估计时为1，否则为0
fiteva=zeros(ps,1);%每一代的对应的个体被实际计算时为1，否则为0 
%%%%%%%%%
%%%%%%%%%拉丁超立方体初始化
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=0.5*(VRmax-VRmin);%初始化
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
% Vmin=repmat(-mv,ps,1);
% Vmax=-Vmin;
 v = zeros(m,D);
% p=VRmin+(VRmax-VRmin).*rand(ps,D);  %用rand函数生成随机数
p=VRmin+(VRmax-VRmin).*lhsdesign(ps,D); %%%用拉丁超立方生成随机数


    e(:,1)=benchmark_func(p,func_num);%计算出的适应值是一个行向量
    fiteva(:,1)=1;
    fitest(:,1)=0;


posfit=e;

evacount=ps;%实际计算次数
% FITEVAALL=fiteva;%%%每一个的对应的个体被实际计算时为1，否则为0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以上为初始化
       [posfit, rank] = sort(posfit, 'descend');%降序排列
        p = p(rank,:);
        v = v(rank,:);
        besty = posfit(m);
        gbest = p(m, :);
        gbestval = besty;
        pastpos=p;%一列是一行数
        pastposfit=e;
        fprintf('Best fitness: %e\n',gbestval);
        Itergbest(1,1)=gbestval;
        itereva(1,1)=evacount; 
        convergegbest=gbestval;%%%%%%%%%%%%%%%%%%%%%第一次采样
        convergeeva=evacount;
    %population sorting
       
        center = ones(m,1)*mean(p);  %center position
        %random matrix
        %rand('seed', sum(100 * clock));
        randco1 = rand(m, D);
        %rand('seed', sum(100 * clock));
        randco2 = rand(m, D);
        %rand('seed', sum(100 * clock));
        randco3 = rand(m, D);
        winidxmask = repmat([1:m]', [1 D]);
        winidx = winidxmask + ceil(rand(m, D).*(m - winidxmask));%ceil函数的作用是朝正无穷方向取整
        %winidx = m - floor(0.5*rand(m, d).*(m - winidxmask));
        pwin = p;
        for j = 1:D
                pwin(:,j) = p(winidx(:,j),j);
        end
        
        %learning
         lpmask = repmat(rand(m,1) < PL, [1 D]);
         lpmask(m,:) = 0;
       v1 =  1*(randco1.*v + randco2.*(pwin - p) + c3*randco3.*(center - p));
% v1 =  1*(randco1.*v + randco2.*(pwin - p));
         p1 =  p + v1;   
         
         
         v = lpmask.*v1 + (~lpmask).*v;         
         p = lpmask.*p1 + (~lpmask).*p;
%          for i = 1:m
            p = max(p, VRmin);
            p = min(p, VRmax);
%         end
i2=2; 

e(:,1)=benchmark_func(p,func_num);%计算出的适应值是一个行向量
fiteva(:,1)=1;
fitest(:,1)=0;
posfit=e;

evacount=evacount+ps;%实际最优点是在变的
convergegbest=[convergegbest;gbestval];%%%%%%%%%%%%%%%%%%%%%第二次采样
convergeeva=[convergeeva;evacount];
 [posfit, rank] = sort(posfit, 'descend');%降序排列
        p = p(rank,:);
        v = v(rank,:);
        besty = posfit(m);
        gbest = p(m, :);
        gbestval = besty;
        fprintf('Best fitness: %e\n',gbestval);
        Itergbest(2,1)=gbestval;
        pastpos=[pastpos;p];%一列是一行数
        pastposfit=[pastposfit;e];
        itereva(2,1)=evacount;
   convergegbest=[convergegbest;gbestval];%%%%%%%%%%%%%%%%%%%%%第二次采样
convergeeva=[convergeeva;evacount];     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以上为第二代
    %population sorting

        center = ones(m,1)*mean(p);  %center position
        %random matrix
        %rand('seed', sum(100 * clock));
        randco1 = rand(m, D);
        %rand('seed', sum(100 * clock));
        randco2 = rand(m, D);
        %rand('seed', sum(100 * clock));
        randco3 = rand(m, D);
        winidxmask = repmat([1:m]', [1 D]);
        winidx = winidxmask + ceil(rand(m, D).*(m - winidxmask));%ceil函数的作用是朝正无穷方向取整
        %winidx = m - floor(0.5*rand(m, d).*(m - winidxmask));
        pwin = p;
        for j = 1:D
                pwin(:,j) = p(winidx(:,j),j);
        end
        
        %learning
         lpmask = repmat(rand(m,1) < PL, [1 D]);
         lpmask(m,:) = 0;
%          v1 =  1*(randco1.*v + randco2.*(pwin - p) + c3*randco3.*(center - p));
 v1 =  1*(randco1.*v + randco2.*(pwin - p));
         p1 =  p + v1;   
         
         
         v = lpmask.*v1 + (~lpmask).*v;         
         p = lpmask.*p1 + (~lpmask).*p;
%          for i = 1:m
            p = max(p, VRmin);
            p = min(p, VRmax);
%         end
i2=3; 
        GPfitness;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%从第三代开始用GP预测
% while(i2<=me&&evacount<=100*D)
while(evacount<=maxEVA)
%        GPbuilding=0;
%         [posfit, rank] = sort(posfit, 'descend');在GPfitness中已经排过序了
        p = p(gbestmid,:);
        v = v(gbestmid,:);
%         besty = posfit(m);
%         gbest = p(m, :);
%         gbestval = min(besty, gbestval);
    %population sorting
        center = ones(m,1)*mean(p);  %center position
        randco1 = rand(m, D);
        %rand('seed', sum(100 * clock));
        randco2 = rand(m, D);
        %rand('seed', sum(100 * clock));
        randco3 = rand(m, D);
        winidxmask = repmat([1:m]', [1 D]);
        winidx = winidxmask + ceil(rand(m, D).*(m - winidxmask));%ceil函数的作用是朝正无穷方向取整
        %winidx = m - floor(0.5*rand(m, d).*(m - winidxmask));
        pwin = p;
        for j = 1:D
                pwin(:,j) = p(winidx(:,j),j);
        end
        
        %learning
         lpmask = repmat(rand(m,1) < PL, [1 D]);
         lpmask(m,:) = 0;
         v1 =  1*(randco1.*v + randco2.*(pwin - p) + c3*randco3.*(center - p));
% v1 =  1*(randco1.*v + randco2.*(pwin - p));
         p1 =  p + v1;   
         
         
         v = lpmask.*v1 + (~lpmask).*v;         
         p = lpmask.*p1 + (~lpmask).*p;
%          for i = 1:m
            p = max(p, VRmin);
            p = min(p, VRmax);

            i2=i2+1;
GPfitness;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%用GP预测
%  Itergbest(i2,1)=gbestval;
%      itereva(i2,1)=evacount;  
end
internum=i2;
end