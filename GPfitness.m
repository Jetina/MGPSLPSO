      fitest=zeros(ps,1);%每一代的对应的个体被估计时为1，否则为0
      fiteva=zeros(ps,1);%每一代的对应的个体被实际计算时为1，否则为0 
        minevafitpast=min(pastposfit);
%         changelessGBNUM=0;%最优适应值不变的次数
    GPglobal;
    [msenum,~]=size(unique(mse));
    
        for mj=1:m
        k=ismember(p(mj,:),pastpos,'rows');
          if k==1
            pastid=find((ismember(p(mj,:),pastpos,'rows')));
            fiteva(mj,1)=1;%此处仅代表实际计算过，因为是赋值，所以并没有增加计算次数
            posfit(mj,1)=pastposfit(pastid(1),1);
            posmse(mj,1)=minmse;
          end
        end

        %选出实际计算的个体
        oneindex=find(fiteva==1);
        k1=length(oneindex); %%所有1的个数
        if k1==0 
        evap=[];
        evaposfit=[];
        else
        first1id=oneindex(1); %%第一个实际计算的元素的位置
        last1id=oneindex(k1); %%第二个实际计算的元素的位置
        evap=p(oneindex(1:k1),:);
        evaposfit=posfit(oneindex(1:k1),:);
        evafiteva=fiteva(oneindex(1:k1),1);
        evaposmse=posmse(oneindex(1:k1),1);
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%注意此处可利用
        
        zero_index=find(fiteva==0) ;
        k0=length(zero_index); %%所有0的个数
        firstzeroid=zero_index(1); %%第一个0元素的位置
        lastzeroid=zero_index(k0); %%最后一个0元素的位置
        
%         if k0~=0  一定不会等于0，所以这句话多余添加
        firstzeroid=zero_index(1); %%第一个0元素的位置
        lastzeroid=zero_index(k0); %%最后一个0元素的位置
        estp=p(zero_index(1:k0),:);
        estposfit=posfit(zero_index(1:k0),:);
        estfiteva=fiteva(zero_index(1:k0),1);
        POPObj1=PopObj(zero_index(1:k0),:);
        estposmse=posmse(zero_index(1:k0),:);
%         end%
        
         [FrontNo,MaxFNo] = NDSort(POPObj1,k0);%去掉之前的最优值，因为之前的最优值已经被实际计算
         difNo=unique(FrontNo);  %返回的是和FrontNo中一样的值，但是没有重复元素。产生的结果向量按升序排序。
         difNoL=length(difNo);%总共分了多少个面
         MinFNo=min(FrontNo);
% %          minmse=min(posmse);这个之前求过了不用再求
         Fristpointsid=find (FrontNo==MinFNo);%找出在第一个面上点的下标
         Lastpointsid=find (FrontNo==MaxFNo);%找出在最后一个面上点的下标
%       [~, gbestsort] = sort(posfit, 'descend');
        FristL=length(Fristpointsid);
        LastL=length(Lastpointsid);
        if msenum==1
                estposfit(Fristpointsid,1)=benchmark_func(estp(Fristpointsid,:),problem);%计算第一个面上的值
                evacount=evacount+1;
                estfiteva(Fristpointsid,1)=1;
                pastpos=[pastpos;estp(Fristpointsid,:)];%%%更新历史数据
                pastposfit=[pastposfit;estposfit(Fristpointsid,1)];
                estposmse(Fristpointsid)=minmse;  
                minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
              
          
   
        elseif msenum~=1  
            if FristL==1%%第一个面上为一个的情况
               estposfit(Fristpointsid,1)=benchmark_func(estp(Fristpointsid,:),problem);
                evacount=evacount+1;
                estfiteva(Fristpointsid,1)=1;
                pastpos=[pastpos;estp(Fristpointsid,:)];%%%更新历史数据
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
                pastpos=[pastpos;estp(Fid1,:)];%%%更新历史数据
                pastposfit=[pastposfit;estposfit(Fid1,1)];
                estposmse(Fid1)=minmse;  
              minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
          
           
          else
%                fprintf('第一个面上的minfit跟maxmse真的不是一个值啊');
              
            Fid1=Fristpointsid(SFsid1);
          estposfit(Fid1,1)=benchmark_func(estp(Fid1,:),problem);
                evacount=evacount+1;
                estfiteva(Fid1,1)=1;
                pastpos=[pastpos;estp(Fid1,:)];%%%更新历史数据
                pastposfit=[pastposfit;estposfit(Fid1,1)];
                estposmse(Fid1)=minmse;  
                minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
           
            Fid2=Fristpointsid(SFsid2);
           estposfit(Fid2,1)=benchmark_func(estp(Fid2,:),problem);
                evacount=evacount+1;
                estfiteva(Fid2,1)=1;
                pastpos=[pastpos;estp(Fid2,:)];%%%更新历史数据
                pastposfit=[pastposfit;estposfit(Fid2,1)];
                estposmse(Fid2)=minmse;  
               minevafitnow = min(pastposfit);
            everyEVA=[everyEVA;evacount];
            everyGBEST=[everyGBEST;minevafitnow ];
          end
            end

         
          if LastL==1%%最后一个面上为一个的情况
                estposfit(Lastpointsid,1)=benchmark_func(estp(Lastpointsid,:),problem);
                evacount=evacount+1;
                estfiteva(Lastpointsid,1)=1;
                pastpos=[pastpos;estp(Lastpointsid,:)];%%%更新历史数据
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
                pastpos=[pastpos;estp( Lid,:)];%%%更新历史数据
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
      [posfit, gbestmid] = sort(posfit, 'descend');%%gbestmid为适应值由大到小排序，
      p=p(gbestmid,:);

    gbestval =minevafitnow;
    fprintf('Best fitness: %e\n',gbestval);
     fprintf('evacount: %e\n',evacount);
