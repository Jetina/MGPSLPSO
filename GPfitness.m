      fitest=zeros(ps,1);%ÿһ���Ķ�Ӧ�ĸ��屻����ʱΪ1������Ϊ0
      fiteva=zeros(ps,1);%ÿһ���Ķ�Ӧ�ĸ��屻ʵ�ʼ���ʱΪ1������Ϊ0 
        minevafitpast=min(pastposfit);
%         changelessGBNUM=0;%������Ӧֵ����Ĵ���
    GPglobal;
    [msenum,~]=size(unique(mse));
    
        for mj=1:m
        k=ismember(p(mj,:),pastpos,'rows');
          if k==1
            pastid=find((ismember(p(mj,:),pastpos,'rows')));
            fiteva(mj,1)=1;%�˴�������ʵ�ʼ��������Ϊ�Ǹ�ֵ�����Բ�û�����Ӽ������
            posfit(mj,1)=pastposfit(pastid(1),1);
            posmse(mj,1)=minmse;
          end
        end

        %ѡ��ʵ�ʼ���ĸ���
        oneindex=find(fiteva==1);
        k1=length(oneindex); %%����1�ĸ���
        if k1==0 
        evap=[];
        evaposfit=[];
        else
        first1id=oneindex(1); %%��һ��ʵ�ʼ����Ԫ�ص�λ��
        last1id=oneindex(k1); %%�ڶ���ʵ�ʼ����Ԫ�ص�λ��
        evap=p(oneindex(1:k1),:);
        evaposfit=posfit(oneindex(1:k1),:);
        evafiteva=fiteva(oneindex(1:k1),1);
        evaposmse=posmse(oneindex(1:k1),1);
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ע��˴�������
        
        zero_index=find(fiteva==0) ;
        k0=length(zero_index); %%����0�ĸ���
        firstzeroid=zero_index(1); %%��һ��0Ԫ�ص�λ��
        lastzeroid=zero_index(k0); %%���һ��0Ԫ�ص�λ��
        
%         if k0~=0  һ���������0��������仰�������
        firstzeroid=zero_index(1); %%��һ��0Ԫ�ص�λ��
        lastzeroid=zero_index(k0); %%���һ��0Ԫ�ص�λ��
        estp=p(zero_index(1:k0),:);
        estposfit=posfit(zero_index(1:k0),:);
        estfiteva=fiteva(zero_index(1:k0),1);
        POPObj1=PopObj(zero_index(1:k0),:);
        estposmse=posmse(zero_index(1:k0),:);
%         end%
        
         [FrontNo,MaxFNo] = NDSort(POPObj1,k0);%ȥ��֮ǰ������ֵ����Ϊ֮ǰ������ֵ�Ѿ���ʵ�ʼ���
         difNo=unique(FrontNo);  %���ص��Ǻ�FrontNo��һ����ֵ������û���ظ�Ԫ�ء������Ľ����������������
         difNoL=length(difNo);%�ܹ����˶��ٸ���
         MinFNo=min(FrontNo);
% %          minmse=min(posmse);���֮ǰ����˲�������
         Fristpointsid=find (FrontNo==MinFNo);%�ҳ��ڵ�һ�����ϵ���±�
         Lastpointsid=find (FrontNo==MaxFNo);%�ҳ������һ�����ϵ���±�
%       [~, gbestsort] = sort(posfit, 'descend');
        FristL=length(Fristpointsid);
        LastL=length(Lastpointsid);
        if msenum==1
                estposfit(Fristpointsid,1)=benchmark_func(estp(Fristpointsid,:),problem);%�����һ�����ϵ�ֵ
                evacount=evacount+1;
                estfiteva(Fristpointsid,1)=1;
                pastpos=[pastpos;estp(Fristpointsid,:)];%%%������ʷ����
                pastposfit=[pastposfit;estposfit(Fristpointsid,1)];
                estposmse(Fristpointsid)=minmse;  
                minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
              
          
   
        elseif msenum~=1  
            if FristL==1%%��һ������Ϊһ�������
               estposfit(Fristpointsid,1)=benchmark_func(estp(Fristpointsid,:),problem);
                evacount=evacount+1;
                estfiteva(Fristpointsid,1)=1;
                pastpos=[pastpos;estp(Fristpointsid,:)];%%%������ʷ����
                pastposfit=[pastposfit;estposfit(Fristpointsid,1)];
                estposmse(Fristpointsid)=minmse;  
               minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
            elseif FristL~=1
         [~, SFsid1] = min(estposfit(Fristpointsid(1:FristL),1));
          [~, SFsid2] = max(estposmse(Fristpointsid(1:FristL),1));
          if SFsid1==SFsid2
         Fid1=Fristpointsid(SFsid1);
                estposfit(Fid1,1)=benchmark_func(estp(Fid1,:),problem);
                evacount=evacount+1;
                estfiteva(Fid1,1)=1;
                pastpos=[pastpos;estp(Fid1,:)];%%%������ʷ����
                pastposfit=[pastposfit;estposfit(Fid1,1)];
                estposmse(Fid1)=minmse;  
              minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
          
           
          else
%                fprintf('��һ�����ϵ�minfit��maxmse��Ĳ���һ��ֵ��');
              
            Fid1=Fristpointsid(SFsid1);
          estposfit(Fid1,1)=benchmark_func(estp(Fid1,:),problem);
                evacount=evacount+1;
                estfiteva(Fid1,1)=1;
                pastpos=[pastpos;estp(Fid1,:)];%%%������ʷ����
                pastposfit=[pastposfit;estposfit(Fid1,1)];
                estposmse(Fid1)=minmse;  
                minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
           
            Fid2=Fristpointsid(SFsid2);
           estposfit(Fid2,1)=benchmark_func(estp(Fid2,:),problem);
                evacount=evacount+1;
                estfiteva(Fid2,1)=1;
                pastpos=[pastpos;estp(Fid2,:)];%%%������ʷ����
                pastposfit=[pastposfit;estposfit(Fid2,1)];
                estposmse(Fid2)=minmse;  
               minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
          end
            end

         
          if LastL==1%%���һ������Ϊһ�������
                estposfit(Lastpointsid,1)=benchmark_func(estp(Lastpointsid,:),problem);
                evacount=evacount+1;
                estfiteva(Lastpointsid,1)=1;
                pastpos=[pastpos;estp(Lastpointsid,:)];%%%������ʷ����
                pastposfit=[pastposfit;estposfit(Lastpointsid,1)];
               estposmse(Lastpointsid)=minmse;  
               minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
           
  
           
           
        elseif LastL~=1
         [~, SLastid] = max(estposmse(Lastpointsid(1:LastL),1));
         Lid=Lastpointsid(SLastid);
          estposfit( Lid,1)=benchmark_func(estp( Lid,:),problem);
                evacount=evacount+1;
                estfiteva( Lid,1)=1;
                pastpos=[pastpos;estp( Lid,:)];%%%������ʷ����
                pastposfit=[pastposfit;estposfit( Lid,1)];
                estposmse( Lid)=minmse;  
              minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
          end

      end
           
       
        
       
       posfit=[evaposfit;estposfit];
       p=[evap;estp];
       posmse=[evaposmse;estposmse]; 
       fiteva=[evafiteva;estfiteva];
      [posfit, gbestmid] = sort(posfit, 'descend');%%gbestmidΪ��Ӧֵ�ɴ�С����
      p=p(gbestmid,:);

    gbestval =minevafitnow;
    fprintf('Best fitness: %e\n',gbestval);
     fprintf('evacount: %e\n',evacount);
