function fMatching_AssignedPoint(surface,ChiralInv,p_ISN)
%% Searching Red Point
disp('Searching Start¡¤¡¤¡¤');

total_num=size(surface.Z,1);sym_p_num=1;in_p_index=p_ISN;

Dis_CI=zeros(sym_p_num,total_num);%Distance Matrix
for j=1:total_num
    Dis_CI(1,j)=pdist2(ChiralInv(:,in_p_index(1))',-ChiralInv(:,j)','euclidean');   
end     

NumRes=zeros(sym_p_num,total_num);
% ValRes=zeros(sym_p_num,total_num);
[lamdnum,lamdpos]=sort(Dis_CI(1,:));
NumRes(1,:)=lamdpos(1,:);
% ValRes(1,:)=lamdnum(1,:);

out_p_index=NumRes(:,1)';
disp('Searching Finished£¡');
%% Calculate normal vector

p_in_num=[surface.X(in_p_index,1)';surface.Y(in_p_index,1)';surface.Z(in_p_index,1)'];
p_out_num=[surface.X(out_p_index,1)';surface.Y(out_p_index,1)';surface.Z(out_p_index,1)'];
normal=(p_out_num(:,1)-p_in_num(:,1))./norm(p_out_num(:,1)-p_in_num(:,1));
%% Result Figure
x=[surface.X,surface.Y,surface.Z];
[COEFF, SCORE]=pca(x);
score=abs(normal'*COEFF(:,1));
score=acos(score)*180/pi;

figure('numbertitle','off','name',['SPD: Assigned White Point ',num2str(p_ISN),' Angle:',num2str(score)]);
surfaceDisplay(surface);
hold on;
% White point
plot3(p_in_num(1,1),p_in_num(2,1),p_in_num(3,1),'w.','markersize',20);   
hold on;
% Red poing
plot3(p_out_num(1,1),p_out_num(2,1),p_out_num(3,1),'r.','markersize',20);   
hold on;  
% Symmetry plane
a=normal(1,1);b=normal(2,1);c=normal(3,1);
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

