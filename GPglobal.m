%% ѵ����/���Լ�����
% if GPbuilding==1
%%%%�ⲿ�ֻ�ʹ����������210
%     trainnum=209;
%   [TrainL,TrainD]=size(pastpos); 
% IP_train = pastpos((TrainL-trainnum:TrainL),:);
% OT_train = pastposfit((TrainL-trainnum:TrainL),:);
% INOUT=[IP_train OT_train];
% [INOUT,~,~] = unique(INOUT,'rows');
% P_train =INOUT(:,(1:D));
% T_train =INOUT(:,D+1);
% T_test = p;
%%%%%%%%%%%%%%%
%%%%%%%%%ȷ������Ϊ210,����Ϊ��õ�210������
%  trainnum=2*m;
% %   [TrainL,TrainD]=size(pastpos); 
% IP_train = pastpos;
% OT_train = pastposfit;
% INput=[ OT_train IP_train];
% [INput,~,~] = unique(INput,'rows');
%  [PRETrainL,PRETrainD]=size(INput); 
%  if PRETrainL<2*m
%      PREIN=INput;
%  else
% PREIN=INput((1:210),:);
%  end
% P_train =PREIN(:,(2:D+1));
% T_train =PREIN(:,1);
% T_test = p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%ȷ������Ϊ210������Ϊ���ѡ���210������
%  trainnum=2*m;
% %   [TrainL,TrainD]=size(pastpos); 
% IP_train = pastpos;
% OT_train = pastposfit;
% INput=[ OT_train IP_train];
% [INput,~,~] = unique(INput,'rows');
%  [PRETrainL,PRETrainD]=size(INput); 
%  if PRETrainL<2*m
%      PREIN=INput;
%  else
%     randa=randperm(PRETrainL);
% PREIN=INput(randa(1:210),:);
%  end
% P_train =PREIN(:,(2:D+1));
% T_train =PREIN(:,1);
% T_test = p;
%%%%%%%%%%%%%%%%%%%%%ȷ������Ϊ210������Ϊ���µ�ʵ�ʼ���������
SampleL=4;
 trainnum=SampleL*m;
   [TrainL,TrainD]=size(pastpos); 
IP_train = pastpos;
OT_train = pastposfit;
INput=[ OT_train IP_train];
if TrainL<=trainnum
     [INput,~,~] = unique(INput,'rows');
     PREIN=INput;
      [PRETrainL,PRETrainD]=size(PREIN); 
 else
      [INput,~,~] = unique(INput,'rows');
      PREIN=INput;
      [PRETrainL,PRETrainD]=size(PREIN); 
      if PRETrainL>trainnum
                   
       PREIN=INput(((PRETrainL-trainnum):PRETrainL),:);%%�Է����ظ��ģ�����������Ϊ����
      [PREIN,~,~] = unique(PREIN,'rows');
      [PRETrainL,PRETrainD]=size(PREIN);
      end
 end

P_train =PREIN(:,(2:D+1));
T_train =PREIN(:,1);
T_test = p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% GP����/ѵ��
% gprMdl = fitrgp(Pn_train,Tn_train,'KernelFunction','ardsquaredexponential', 'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%     struct('AcquisitionFunctionName','expected-improvement-plus'));
%��ͬ�ĺ˺������� 'exponential''squaredexponential''matern32''matern52''rationalquadratic'
gprMdl = fitrgp(P_train,T_train,'KernelFunction','exponential','BasisFunction','pureQuadratic','Optimizer','fmincon');
%%%%%%%%%%%%%%%%�ⲿ������dace�Ĵ���
% lob=repmat(1e-3,1,D);
% upb=repmat(1e3,1,D);
% theta=std(pastpos);
% [dmodel,perf]=dacefit(P_train,T_train,@regpoly0,@corrgauss,theta,lob,upb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GP�������

[t_sim,mse]=predict(gprMdl,T_test);

posfit=t_sim;
posmse=mse;
PopObj=[t_sim mse];
% estcount=estcount+ps;%���ƴ���
% fiteva=zeros(ps,1);
% fitest=ones(ps,1);
minmse=min(posmse);

