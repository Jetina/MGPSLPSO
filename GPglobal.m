%% 训练集/测试集产生
% if GPbuilding==1
%%%%这部分会使得样本少于210
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
%%%%%%%%%确保样本为210,这里为最好的210个样本
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
%%%%%%%%%%%%%%确保样本为210，这里为随机选择的210个样本
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
%%%%%%%%%%%%%%%%%%%%%确保样本为210，这里为最新的实际计算后的样本
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
                   
       PREIN=INput(((PRETrainL-trainnum):PRETrainL),:);%%以防有重复的，导致索引变为负数
      [PREIN,~,~] = unique(PREIN,'rows');
      [PRETrainL,PRETrainD]=size(PREIN);
      end
 end

P_train =PREIN(:,(2:D+1));
T_train =PREIN(:,1);
T_test = p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% GP创建/训练
% gprMdl = fitrgp(Pn_train,Tn_train,'KernelFunction','ardsquaredexponential', 'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%     struct('AcquisitionFunctionName','expected-improvement-plus'));
%不同的核函数参数 'exponential''squaredexponential''matern32''matern52''rationalquadratic'
gprMdl = fitrgp(P_train,T_train,'KernelFunction','exponential','BasisFunction','pureQuadratic','Optimizer','fmincon');
%%%%%%%%%%%%%%%%这部分是用dace的代码
% lob=repmat(1e-3,1,D);
% upb=repmat(1e3,1,D);
% theta=std(pastpos);
% [dmodel,perf]=dacefit(P_train,T_train,@regpoly0,@corrgauss,theta,lob,upb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GP仿真测试

[t_sim,mse]=predict(gprMdl,T_test);

posfit=t_sim;
posmse=mse;
PopObj=[t_sim mse];
% estcount=estcount+ps;%估计次数
% fiteva=zeros(ps,1);
% fitest=ones(ps,1);
minmse=min(posmse);

