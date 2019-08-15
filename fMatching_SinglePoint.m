function fMatching_SinglePoint(surface,ChiralInv)
try
    res_data=load('Matching_SinglePoing_Res.mat');
    disp('Result Loading ...');
    res_matching=res_data.res_matching;
    conference_xv=res_data.conference_xv;
    clear res_name res_data;    
catch
    disp('Result Calculating ...');
    %% Normal result of single point matching 
    Dis_CI_temp=pdist2(ChiralInv',-ChiralInv','euclidean');%Distance Matrix based on chiral invariants
    total_num=size(surface.X,1);
    NumRes_all=zeros(total_num,total_num);
    % ValRes_all=zeros(total_num,total_num);
    parfor i=1:total_num
        [lamdnum,lamdpos]=sort(Dis_CI_temp(i,:));
        NumRes_all(i,:)=lamdpos(1,:);
    %     ValRes_all(i,:)=lamdnum(1,:);
    end
    out_p_index=NumRes_all(:,1)';

    p_in_num=[surface.X';surface.Y';surface.Z'];
    p_out_num=[surface.X(out_p_index,1)';surface.Y(out_p_index,1)';surface.Z(out_p_index,1)'];
    x=(p_out_num-p_in_num)';% normal vector of symmetry plane of each point
    parfor i=1:size(x,1)
        x(i,:)=x(i,:)./norm(x(i,:));
    end   
    %% Calculate matching score 
    conference_head=pca(p_in_num');
    conference_xv=conference_head(:,1);
    res_matching=abs(x*conference_xv);%score
    save('Matching_SinglePoing_Res','res_matching','conference_xv');
end

%% Matching score visualization
res_matching=acos(res_matching)*180/pi;
surface.I=res_matching;
figure('numbertitle','off','name','Visualization of Single Point Matching Precision');
surfaceDisplay(surface);
%% Symmetry plane reference visualization
figure('numbertitle','off','name','Symmetry Plane Reference');
surface.I=surface.I*0;
surfaceDisplay(surface);
hold on;
x=[surface.X,surface.Y,surface.Z];
o=sum(x,1)/size(x,1);
normal=conference_xv';
a=normal(1,1);b=normal(1,2);c=normal(1,3);
if c~=0 
    z=min(surface.Z):0.1:max(surface.Z);
    y=min(surface.Y):0.1:max(surface.Y);
    [y,z]=meshgrid(y,z);
    x=-(b*(y-o(1,2))+c*(z-o(1,3)))/a+o(1,1);
    h=mesh(x,y,z);
    set(h,'EdgeColor','r','FaceColor','r','MarkerEdgecolor','k','MarkerFacecolor','w')
end

end