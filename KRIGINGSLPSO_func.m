function [gbest,gbestval,evacount,estcount,Itergbest,itereva,internum,convergegbest,convergeeva]= KRIGINGSLPSO_func(Dimension,Max_Gen,VRmin,VRmax,varargin)
%[gbest,gbestval,fitcount]= PSO_func('f8',3500,200000,30,30,-5.12,5.12)
rand('state',sum(100*clock));

maxEVA=Max_Gen;%z����������
%ps=Particle_Number;%��Ⱥ����
D=Dimension;
estcount=0;%���Ƶ��ܴ���
func_num=problem;

%   alpha=0.5;%%����ѡȡ��Χ
%%%%��ʼslpso����
 M = 100;
    m = M + D/10;
    c3 = D/M*0.01;
    PL = zeros(m,1);

    for i = 1 : m
        PL(i) = (1 - (i - 1)/m)^log(sqrt(ceil(D/M)));%ѧϰ����
    end
 ps=m;  
 fitest=zeros(ps,1);%ÿһ���Ķ�Ӧ�ĸ��屻����ʱΪ1������Ϊ0
fiteva=zeros(ps,1);%ÿһ���Ķ�Ӧ�ĸ��屻ʵ�ʼ���ʱΪ1������Ϊ0 
%%%%%%%%%
%%%%%%%%%�������������ʼ��
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=0.5*(VRmax-VRmin);%��ʼ��
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
% Vmin=repmat(-mv,ps,1);
% Vmax=-Vmin;
 v = zeros(m,D);
% p=VRmin+(VRmax-VRmin).*rand(ps,D);  %��rand�������������
p=VRmin+(VRmax-VRmin).*lhsdesign(ps,D); %%%���������������������


    e(:,1)=benchmark_func(p,func_num);%���������Ӧֵ��һ��������
    fiteva(:,1)=1;
    fitest(:,1)=0;


posfit=e;

evacount=ps;%ʵ�ʼ������
% FITEVAALL=fiteva;%%%ÿһ���Ķ�Ӧ�ĸ��屻ʵ�ʼ���ʱΪ1������Ϊ0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ϊ��ʼ��
       [posfit, rank] = sort(posfit, 'descend');%��������
        p = p(rank,:);
        v = v(rank,:);
        besty = posfit(m);
        gbest = p(m, :);
        gbestval = besty;
        pastpos=p;%һ����һ����
        pastposfit=e;
        fprintf('Best fitness: %e\n',gbestval);
        Itergbest(1,1)=gbestval;
        itereva(1,1)=evacount; 
        convergegbest=gbestval;%%%%%%%%%%%%%%%%%%%%%��һ�β���
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
        winidx = winidxmask + ceil(rand(m, D).*(m - winidxmask));%ceil�����������ǳ��������ȡ��
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

e(:,1)=benchmark_func(p,func_num);%���������Ӧֵ��һ��������
fiteva(:,1)=1;
fitest(:,1)=0;
posfit=e;

evacount=evacount+ps;%ʵ�����ŵ����ڱ��
convergegbest=[convergegbest;gbestval];%%%%%%%%%%%%%%%%%%%%%�ڶ��β���
convergeeva=[convergeeva;evacount];
 [posfit, rank] = sort(posfit, 'descend');%��������
        p = p(rank,:);
        v = v(rank,:);
        besty = posfit(m);
        gbest = p(m, :);
        gbestval = besty;
        fprintf('Best fitness: %e\n',gbestval);
        Itergbest(2,1)=gbestval;
        pastpos=[pastpos;p];%һ����һ����
        pastposfit=[pastposfit;e];
        itereva(2,1)=evacount;
   convergegbest=[convergegbest;gbestval];%%%%%%%%%%%%%%%%%%%%%�ڶ��β���
convergeeva=[convergeeva;evacount];     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ϊ�ڶ���
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
        winidx = winidxmask + ceil(rand(m, D).*(m - winidxmask));%ceil�����������ǳ��������ȡ��
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
        GPfitness;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ӵ�������ʼ��GPԤ��
% while(i2<=me&&evacount<=100*D)
while(evacount<=maxEVA)
%        GPbuilding=0;
%         [posfit, rank] = sort(posfit, 'descend');��GPfitness���Ѿ��Ź�����
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
        winidx = winidxmask + ceil(rand(m, D).*(m - winidxmask));%ceil�����������ǳ��������ȡ��
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
GPfitness;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��GPԤ��
%  Itergbest(i2,1)=gbestval;
%      itereva(i2,1)=evacount;  
end
internum=i2;
end