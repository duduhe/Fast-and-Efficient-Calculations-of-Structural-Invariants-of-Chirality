function fSPD_RSPCA(surface,ChiralInv,rsNum)
% rsNum--the number of points in random sampling
timestart=cputime;

%% Matching 

% Calculate the distance matrix of rsNum points      
sym_p_num=rsNum;total_num=size(surface.X,1);
in_p_index=randi([1,total_num],1,sym_p_num);
Dis_CI_temp=pdist2(ChiralInv',-ChiralInv','euclidean');% Distance Matrix 
Dis_CI=Dis_CI_temp(in_p_index,:);

% Distance ranking for matching
NumRes=zeros(sym_p_num,total_num);
% ValRes=zeros(sym_p_num,total_num);
parfor i=1:sym_p_num
    [lamdnum,lamdpos]=sort(Dis_CI(i,:));
    NumRes(i,:)=lamdpos(1,:);
%     ValRes(i,:)=lamdnum(1,:);
end
out_p_index=NumRes(:,1)';

%% PCA->Normal vector

p_in_num=[surface.X(in_p_index,1)';surface.Y(in_p_index,1)';surface.Z(in_p_index,1)'];
p_out_num=[surface.X(out_p_index,1)';surface.Y(out_p_index,1)';surface.Z(out_p_index,1)'];

x=(p_out_num-p_in_num)';
[COEFF, SCORE]=pca(x);
% parfor i=1:size(x,1)
%     SCORE(i,:)=SCORE(i,:)./norm(SCORE(i,:));
% end        
normal=COEFF(:,1)';% normal vector
disp([num2str(rsNum),' points',' SPD_RSPCA time consumption: ',num2str(cputime-timestart),' s']);
%% Result Figure
x=[surface.X,surface.Y,surface.Z];
[COEFF, SCORE]=pca(x);
score=abs(normal*COEFF(:,1));
score=acos(score)*180/pi;

figure('numbertitle','off','name',['SPD: ',num2str(rsNum),' Random sampling + PCA',' Angle:',num2str(score)]);
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

