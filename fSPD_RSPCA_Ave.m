function fSPD_RSPCA_Ave(surface,ChiralInv,samplingNum,roundNum)

    if samplingNum<4 || roundNum<1
        disp('Input Parameter Wrong!');
        return;
    end
    timestart=cputime;
    disp('Result Calculating ...');

    %% Initialization
    MaxRound=roundNum;
    StepLength=1;
    StartNumofPoint=samplingNum;  
    MaxNumofPoint=samplingNum;     
    Dis_CI_temp=pdist2(ChiralInv',-ChiralInv','euclidean');% Distance Matrix
    conference_head=pca([surface.X,surface.Y,surface.Z]); nv_conf=conference_head(:,1); % normal vector reference
    total_num=size(surface.X,1);
    %% Test loop
    sym_p_num_result_n_size=1+floor((MaxNumofPoint-StartNumofPoint)/StepLength);
    %     sym_p_num_result_n_min=zeros(sym_p_num_result_n_size,3);
    %     sym_p_num_result_n_max=zeros(sym_p_num_result_n_size,3);
    sym_p_num_result_n_ave=zeros(sym_p_num_result_n_size,3);
    for sym_p_num=StartNumofPoint:StepLength:MaxNumofPoint
        round_result_n=zeros(MaxRound,3);
        for round_i=1:MaxRound
            disp(['P:',num2str(sym_p_num),'-R:',num2str(round_i)]);% P--the number of points  R: the number of repeated round
            %% Matching of sym_p_num random sampling points
            in_p_index=randi([1,total_num],1,sym_p_num);       
            Dis_CI=Dis_CI_temp(in_p_index,:);

            NumRes=zeros(sym_p_num,total_num);
    %         ValRes=zeros(sym_p_num,total_num);
            parfor i=1:sym_p_num
                [lamdnum,lamdpos]=sort(Dis_CI(i,:));
                NumRes(i,:)=lamdpos(1,:);
    %             ValRes(i,:)=lamdnum(1,:);
            end
            out_p_index=NumRes(:,1)';
            %% PCA->Normal vector
            p_in_num=[surface.X(in_p_index,1)';surface.Y(in_p_index,1)';surface.Z(in_p_index,1)'];
            p_out_num=[surface.X(out_p_index,1)';surface.Y(out_p_index,1)';surface.Z(out_p_index,1)'];
            x=(p_out_num-p_in_num)';
            [COEFF, SCORE]=pca(x);normal=COEFF(:,1)';

            round_result_n(round_i,:)=normal;
        end

        %计算在（1,0,0）方向投影后的最值
    %     nv_conf
    %         result_n_score=abs(round_result_n*nv_conf);
        result_n_i=1+(sym_p_num-StartNumofPoint)/StepLength;
    %         sym_p_num_result_n_min(result_n_i,:)=round_result_n((find(result_n_score==min(result_n_score))),:);
    %         sym_p_num_result_n_max(result_n_i,:)=round_result_n((find(result_n_score==max(result_n_score))),:);
        sym_p_num_result_n_ave(result_n_i,:)=sum(round_result_n,1)/MaxRound;
    end
    normal=sym_p_num_result_n_ave;
    %     res_name=['sym_p_num_result_n_',num2str(StartNumofPoint),'_',num2str(StepLength),'_',num2str(MaxNumofPoint),'_',num2str(MaxRound)];
    %     save(res_name,'StartNumofPoint','StepLength','MaxNumofPoint','MaxRound','sym_p_num_result_n_min','sym_p_num_result_n_max','sym_p_num_result_n_ave','nv_conf');

disp(['P:',num2str(samplingNum),',R:',num2str(roundNum),' SPD_RSPCA_ave time consumption: ',num2str(cputime-timestart),' s']);


%% Result Figure
x=[surface.X,surface.Y,surface.Z];
[COEFF, SCORE]=pca(x);
score=abs(normal*COEFF(:,1));
score=acos(score/180*pi);

figure('numbertitle','off','name',['SPD: ','P-',num2str(samplingNum),',R-',num2str(roundNum),' Average result of Random sampling + PCA',' Angle:',num2str(score)]);
surfaceDisplay(surface);
hold on;
% Symmetry plane
a=normal(1,1);b=normal(1,2);c=normal(1,3);
o=sum([surface.X,surface.Y,surface.Z])/size(surface.X,1);
if c~=0 
    z=min(surface.Z):0.1:max(surface.Z);
    y=min(surface.Y):0.1:max(surface.Y);
    [y,z]=meshgrid(y,z);
    x=-(b*(y-o(1,2))+c*(z-o(1,3)))/a+o(1,1);
    h=mesh(x,y,z);
    set(h,'EdgeColor','r','FaceColor','r','MarkerEdgecolor','k','MarkerFacecolor','w')
end


end



