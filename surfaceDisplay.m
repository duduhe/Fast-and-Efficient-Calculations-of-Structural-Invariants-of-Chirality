function [  ] = surfaceDisplay( surface )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
h=trisurf(surface.TRIV,surface.X,surface.Y,surface.Z,surface.I);
% xlabel('X','Rotation',15,'fontsize',20)
xlabel('X','fontsize',10);
ylabel('Y','fontsize',10);
zlabel('Z','fontsize',10);
axis equal; colormap(cool);
set(h,'LineWidth',0.25)
set(h,'EdgeColor','k','FaceColor','interp','MarkerEdgecolor','k','MarkerFacecolor','w')
view(0,-90);
% view(0,-0);
% view(180,90);

% view(30,-80);
colorbar;
end

