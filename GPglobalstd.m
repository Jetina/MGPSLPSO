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
 trainnum=2*m;
  [TrainL,TrainD]=size(pastpos); 
IP_train = pastpos;
OT_train = pastposfit;
INput=[ OT_train IP_train];

 if TrainL<=2*m
     [INput,~,~] = unique(INput,'rows');
     PREIN=INput;
      [PRETrainL,PRETrainD]=size(PREIN); 
 else
    
      PREIN=INput((TrainL-trainnum:TrainL),:);
      [PREIN,~,~] = unique(PREIN,'rows');
      [PRETrainL,PRETrainD]=size(PREIN); 
      T=0;
      while PRETrainL<=2*m
          T=T+10;
       PREIN=INput((TrainL-(trainnum+T):TrainL),:);
      [PREIN,~,~] = unique(PREIN,'rows');
      [PRETrainL,PRETrainD]=size(PREIN);
      end
 end
%%%%ԭ����ֱ������ѵ��
P_train =PREIN(:,(2:D+1));
T_train =PREIN(:,1);
T_test = p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%���ݱ�׼������̬����
%%�Ƚ��������ݺͲ������ݷ���һ���׼��
IP_train=[P_train;T_test];
ITP_train=mapstd(IP_train');%��mapstd�ǰ��б�׼��������Ҫ��ת��
SP_train=ITP_train';
SSP_train =SP_train((1:PRETrainL),:);
SST_test =SP_train((PRETrainL+1:end),:);
%% GP����/ѵ��


% gprMdl = fitrgp(P_train,T_train,'KernelFunction','matern32','BasisFunction','pureQuadratic','Optimizer','fmincon');
gprMdl = fitrgp(SSP_train,T_train,'KernelFunction','matern32','BasisFunction','pureQuadratic','Optimizer','fmincon');
%% GP�������

 [t_sim,mse]=predict(gprMdl,SST_test);
% [t_sim,mse]=predict(gprMdl,T_test);
posfit=t_sim;
posmse=mse;
PopObj=[t_sim mse];
estcount=estcount+ps;%���ƴ���
fiteva=zeros(ps,1);
fitest=ones(ps,1);
minmse=min(posmse);

