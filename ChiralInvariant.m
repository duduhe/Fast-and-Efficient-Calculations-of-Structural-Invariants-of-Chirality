function [ ChiralInv ] = ChiralInvariant( surface )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
   
    ChiralInv=zeros(5,size(surface.X,1));%用于存储各点的手性不变量

% % 以上为采用临近点进行曲面拟合后，相关参数的存储--End   
    NumofBoundary=0;
    disp(['正在去除',num2str(NumofBoundary),'层边界点···']);
    [ IsBoundaryVertex ] = GetBoundaryVertex( surface,NumofBoundary );
    NumofOutVert=sum(IsBoundaryVertex);
    NumofBoundary=sum(IsBoundaryVertex);
    parfor index_vertex=1:size(surface.X,1)
        disp(['顶点处理进度：',num2str(index_vertex),'/',num2str(size(surface.X,1))]);
        CenterVertex_index=index_vertex;
        %边界点不参与运算
        if IsBoundaryVertex(1,index_vertex)==1
%             disp(['顶点：',num2str(index_vertex),'为边界点，不参与计算！']);
            continue;
        end
%以下为采用临近点进行曲面拟合后，进行相关参数的计算--Begin
        numofRing=2;
        [ VertexInRing,FaceInRingIndex]=GetVertexInRing(surface,CenterVertex_index,numofRing);
        %当前计算的物体
        obj=zeros(size(VertexInRing,2)+1,3);
        
        %中心点的坐标与f的取值
        fo_x=surface.X(CenterVertex_index,1);
        fo_y=surface.Y(CenterVertex_index,1);
        fo_z=surface.Z(CenterVertex_index,1);
        obj(size(VertexInRing,2)+1,1:3)=[fo_x,fo_y,fo_z];
        for index_VertexInRing=1:size(VertexInRing,2)
            verIsn=VertexInRing(1,index_VertexInRing);
            ver_x=surface.X(verIsn,1);
            ver_y=surface.Y(verIsn,1);
            ver_z=surface.Z(verIsn,1);
            
            obj(index_VertexInRing,1:3)=[ver_x,ver_y,ver_z];
        end
       %计算这些点构成物体的手性不变量
       Point=obj';
       Density=ones(1,size(VertexInRing,2)+1);
       
       Inv = CalInv( Point,Density );      
        
    %以上为采用临近点进行曲面拟合后，进行相关参数的计算--End
    ChiralInv(:,index_vertex)=Inv;
    end
    disp(['共有',num2str(NumofOutVert),'个点未参与计算，其中边界点有',num2str(NumofBoundary),'个']);
end

%函数1：计算当前顶点所在位置周边N环内的点,及相应环内的面片
%输出：VertexInRing，见注释；FaceInRingIndex，见注释。二者皆为行向量。
function [ VertexInRing,FaceInRingIndex]=GetVertexInRing(surface,CenterVertex_index,numofRing)

    VertexInRing=[CenterVertex_index];%查找结果是一个行向量
    ListAdd=zeros(size(surface.TRIV,1),1);%提高搜索效率，已经找过的面片不再进行查找
    %找的次数
    for indexRing=1:numofRing
        tempRingsVertexIndex=VertexInRing;%记录找到的结果的中间变量    
        %遍历所有的面片
        for indexFace=1:size(surface.TRIV,1)
            %判断当前面片是否已经添加到了搜索结果中
            if ListAdd(indexFace,1)==0
                %判断当前面片中的顶点是否在已经找到的结果中
                if sum(ismember(surface.TRIV(indexFace,:),VertexInRing))
                    tempRingsVertexIndex=union(tempRingsVertexIndex,surface.TRIV(indexFace,:));
                    ListAdd(indexFace,1)=1;
                end
            end
        end
        VertexInRing=tempRingsVertexIndex; %更新最新的查找结果    
    end    
    VertexInRing=setdiff(VertexInRing,CenterVertex_index);%去除掉中心点
    FaceInRingIndex=find(ListAdd==1)';%给出n环内的面片编号
end
%函数2：拟合曲面并计算相应的偏导数
function [DiffX,DiffY,DiffXX,DiffYY,DiffXY] = GetDiff(surface,CenterVertex_index,VertexInRing,UpOrder)
%需要在函数1中注释掉去除掉中心点的代码   倒数第二行  即将所有的点的坐标放到方程组中参与计算
    XY=zeros(size(VertexInRing,2),(UpOrder+1)^2);
    if size(XY,1)<size(XY,2)
        disp('当前进行最小二乘法拟合时没有解！')
    end
    Z=zeros(size(VertexInRing,2),1);
    for index_VertexInRing=1:size(VertexInRing,2)
        verIsn=VertexInRing(1,index_VertexInRing);
        ver_x=surface.X(verIsn,1);
        ver_y=surface.Y(verIsn,1);
        ver_z=surface.Z(verIsn,1);
        for i=1:(UpOrder+1)
            for j=1:(UpOrder+1)
                XY(index_VertexInRing,(i-1)*4+j)=(ver_x^(i-1))*(ver_y^(j-1));
            end
        end             
        Z(index_VertexInRing,1)=ver_z;
    end
%     P=inv(XY'*XY)*XY'*Z;
%     P=XY\Z;
    P=pinv(XY)*Z;
    ver_cen_x=surface.X(CenterVertex_index,1);
    ver_cen_y=surface.Y(CenterVertex_index,1);
    %计算Z对X方向上的偏导数
    tempDiffX=zeros(1,(UpOrder+1)^2);
    for i=1:(UpOrder+1)
        for j=1:(UpOrder+1)
            if i==1
%                 if j==1
%                     tempDiffX(1,(i-1)*4+j)=0;
%                 else
%                     %tempDiffX(1,(i-1)*4+j)=(ver_x^(i-1))*(ver_y^(j-1));
%                     tempDiffX(1,(i-1)*4+j)=(ver_cen_y^(j-1));
%                 end      
                tempDiffX(1,(i-1)*4+j)=0;
            else
                tempDiffX(1,(i-1)*4+j)=(i-1)*(ver_cen_x^(i-2))*(ver_cen_y^(j-1));
            end
        end
    end  
    DiffX=tempDiffX*P;    
    %计算Z对Y方向上的偏导数  
    tempDiffY=zeros(1,(UpOrder+1)^2);
    for i=1:(UpOrder+1)
        for j=1:(UpOrder+1)
            if j==1
%                 if i==1
%                     tempDiffY(1,(i-1)*4+j)=0;
%                 else
%                     tempDiffY(1,(i-1)*4+j)=(ver_cen_x^(i-1));
%                 end
                tempDiffY(1,(i-1)*4+j)=0;
            else
                tempDiffY(1,(i-1)*4+j)=(j-1)*(ver_cen_x^(i-1))*(ver_cen_y^(j-2));
            end
        end
    end  
    DiffY=tempDiffY*P;    
    %计算Z对XX的导数
    tempDiffXX=zeros(1,(UpOrder+1)^2);
    for i=1:(UpOrder+1)
        for j=1:(UpOrder+1)
            if i==1||i==2 
                tempDiffXX(1,(i-1)*4+j)=0;
            else
                tempDiffXX(1,(i-1)*4+j)=(i-1)*(i-2)*(ver_cen_x^(i-3))*(ver_cen_y^(j-1));
            end
        end
    end  
    DiffXX=tempDiffXX*P;    
    %计算Z对YY的导数
    tempDiffYY=zeros(1,(UpOrder+1)^2);
    for i=1:(UpOrder+1)
        for j=1:(UpOrder+1)
            if j==1||j==2 
                tempDiffYY(1,(i-1)*4+j)=0;
            else
                tempDiffYY(1,(i-1)*4+j)=(j-1)*(j-2)*(ver_cen_x^(i-1))*(ver_cen_y^(j-3));
            end
        end
    end  
    DiffYY=tempDiffYY*P;
    %计算在Z对XY的导数
    tempDiffXY=zeros(1,(UpOrder+1)^2);
    for i=1:(UpOrder+1)
        for j=1:(UpOrder+1)
            if i==1||j==1
                tempDiffXY(1,(i-1)*4+j)=0;
            else
                tempDiffXY(1,(i-1)*4+j)=(i-1)*(j-1)*(ver_cen_x^(i-2))*(ver_cen_y^(j-2));
            end
        end
    end  
    DiffXY=tempDiffXY*P;
    
    
end
%函数2-2：采用微分几何P179页泰勒展开式，拟合曲面并计算相应的偏导数
function [DiffX,DiffY,DiffXX,DiffXY,DiffYY] = GetDiff_2_Taylor(surface,CenterVertex_index,VertexInRing)
    %需要在函数1中 放开 去除中心点的代码   倒数第二行  即将中心点的坐标单独放出来作为输入
  
    XY=zeros(size(VertexInRing,2),5);
    if size(XY,1)<size(XY,2)
        disp('当前进行最小二乘法拟合时没有解！')
    end
    Z=zeros(size(VertexInRing,2),1);
    %中心点的坐标与f的取值
    fo_x=surface.X(CenterVertex_index,1);
    fo_y=surface.Y(CenterVertex_index,1);
    fo_z=surface.Z(CenterVertex_index,1);
    %构造方程组
    for index_VertexInRing=1:size(VertexInRing,2)
        verIsn=VertexInRing(1,index_VertexInRing);
        ver_x=surface.X(verIsn,1);
        ver_y=surface.Y(verIsn,1);
        ver_z=surface.Z(verIsn,1);
        
        XY(index_VertexInRing,1)=ver_x-fo_x;%x
        XY(index_VertexInRing,2)=ver_y-fo_y;%y
        XY(index_VertexInRing,3)=0.5*(ver_x-fo_x)^2;% 1/2*x^2
        XY(index_VertexInRing,4)=(ver_x-fo_x)*(ver_y-fo_y);%xy
        XY(index_VertexInRing,5)=0.5*(ver_y-fo_y)^2;%1/2*y^2
         
        Z(index_VertexInRing,1)=ver_z-fo_z;
    end
%     P=inv(XY'*XY)*XY'*Z;
%     P=XY\Z;
    P=pinv(XY)*Z;
    DiffX=P(1,1);
    DiffY=P(2,1);
    DiffXX=P(3,1);
    DiffXY=P(4,1);
    DiffYY=P(5,1);    
end



%函数3：计算当前顶点的1-ring面片上的dA--用三角网格的面积等价
function [dA] = GetVertexdA(surface,CenterVertex_index) 
    [RingsVertexIndex,FaceIndex]  = GetVertexInRing(surface,CenterVertex_index,1);
    %计算当前点与周边点的edge所跨越的面片的编号
    [ EdgeAtFaceIndex ] = EdgeAtFace( surface,CenterVertex_index,RingsVertexIndex,FaceIndex);
    %计算当前点与周边点的edge所跨越的两个三角形的对角和面积
    [EdgeBasedAngleIndex,EdgeBasedAreaIndex]=EdgeBasedAngle_Ares(surface,EdgeAtFaceIndex);
    dA=(sum(EdgeBasedAreaIndex(:,3))+sum(EdgeBasedAreaIndex(:,4)));
end
%函数4：找到指定点与周边1环内的点组成的边所跨越的面片的序号
%输出：EdgeAtFaceIndex，见注释
function [ EdgeAtFaceIndex ] = EdgeAtFace( surface,CenterVertex_index,RingsVertexIndex,FaceIndex)
    %第一列为CenterVertex_index，第二列为RingsVertexIndex，
    EdgeAtFaceIndex=zeros(size(RingsVertexIndex,2),4);
    EdgeAtFaceIndex(:,1)=ones(size(RingsVertexIndex,2),1)*CenterVertex_index;
    EdgeAtFaceIndex(:,2)=RingsVertexIndex';
    %第3、4列为当前边所跨越的面片的序号
    %遍历所有的边
    for indexCenterVertex_index=1:size(RingsVertexIndex,2)
        tempEdge=EdgeAtFaceIndex(indexCenterVertex_index,1:2);
        signofFaceNum=0;
        %遍历所有的面
        for indexFaceIndex=1:size(FaceIndex,2)
            tempFace=surface.TRIV(FaceIndex(1,indexFaceIndex),:);
            if sum(ismember(tempFace,tempEdge))==2
                if signofFaceNum==0
                    EdgeAtFaceIndex(indexCenterVertex_index,3)=FaceIndex(1,indexFaceIndex);        
                    signofFaceNum=signofFaceNum+1;                    
                elseif signofFaceNum==1
                    EdgeAtFaceIndex(indexCenterVertex_index,4)=FaceIndex(1,indexFaceIndex);        
                    signofFaceNum=signofFaceNum+1;                     
                else
                    display('ERROR:EdgeAtFace function error')
                end                                
            end
        end
        if signofFaceNum==0
            display('ERROR:当前边无所在的面！')
        elseif signofFaceNum==1
%             display('边界edge！')
        end
    end  
end
%函数5：根据EdgeAtFaceIndex计算edge所跨越的两个三角形的对角和面积
%输出：EdgeBasedAngleIndex，见注释；EdgeBasedAreaIndex，见注释
function [EdgeBasedAngleIndex,EdgeBasedAreaIndex]=EdgeBasedAngle_Ares(surface,EdgeAtFaceIndex)
    %索引格式与EdgeAtFaceIndex相似，不同的是对应的3、4列分别为角度的cot和面积
    EdgeBasedAngleIndex=zeros(size(EdgeAtFaceIndex,1),size(EdgeAtFaceIndex,2));%边对应的角度索引
    EdgeBasedAreaIndex=zeros(size(EdgeAtFaceIndex,1),size(EdgeAtFaceIndex,2));%边对应的面积索引
    
    EdgeBasedAngleIndex(:,1:2)=EdgeAtFaceIndex(:,1:2);
    EdgeBasedAreaIndex(:,1:2)=EdgeAtFaceIndex(:,1:2);
    
    for indexEdgeAtFaceIndex=1:size(EdgeAtFaceIndex,1)
        edge=EdgeAtFaceIndex(indexEdgeAtFaceIndex,1:2);%得到一条边
        
        %得到一个三角形    
        faceindex1=EdgeAtFaceIndex(indexEdgeAtFaceIndex,3);  
        if faceindex1==0
%             display(['边  ',num2str(edge),' 的一边没有面片！']);
        else
            triangle1=surface.TRIV(faceindex1,:); 
            triangle1_1_2=find(ismember(triangle1,edge));%第1、2个点在三角形中的位置
            triangle1_3=find(ismember(triangle1,edge)==0);%第3个点在三角形中的位置
            %得到三角形点的坐标
            triangle1_vertex3=[surface.X(triangle1(1,triangle1_3),1);surface.Y(triangle1(1,triangle1_3),1);surface.Z(triangle1(1,triangle1_3),1)];
            triangle1_vertex1=[surface.X(triangle1(1,triangle1_1_2(1,1)),1);surface.Y(triangle1(1,triangle1_1_2(1,1)),1);surface.Z(triangle1(1,triangle1_1_2(1,1)),1)];    
            triangle1_vertex2=[surface.X(triangle1(1,triangle1_1_2(1,2)),1);surface.Y(triangle1(1,triangle1_1_2(1,2)),1);surface.Z(triangle1(1,triangle1_1_2(1,2)),1)];
            triangle1_edge12=triangle1_vertex2-triangle1_vertex1;%列向量
            triangle1_edge31=triangle1_vertex1-triangle1_vertex3;
            triangle1_edge32=triangle1_vertex2-triangle1_vertex3;
            %得到三角形的边长
            triangle1_edge12_length=sqrt(triangle1_edge12'*triangle1_edge12);
            triangle1_edge31_length=sqrt(triangle1_edge31'*triangle1_edge31);
            triangle1_edge32_length=sqrt(triangle1_edge32'*triangle1_edge32);
            %计算面积  使用半周长的方法来计算三角形面积  http://mathworld.wolfram.com/Semiperimeter.html
            triangle1_s_half=(triangle1_edge12_length+triangle1_edge31_length+triangle1_edge32_length)/2;
            triangle1_Area=sqrt(triangle1_s_half*(triangle1_s_half-triangle1_edge12_length)*(triangle1_s_half-triangle1_edge31_length)*(triangle1_s_half-triangle1_edge32_length)); 
            EdgeBasedAreaIndex(indexEdgeAtFaceIndex,3)=triangle1_Area;
            %计算角度
            %triangle1_Angle=asin(2*triangle1_Area/triangle1_edge31_length/triangle1_edge32_length);% *180/pi  这里角度的单位是弧度
            %半角定理
            temp_tri1_tan=sqrt((triangle1_s_half-triangle1_edge31_length)*(triangle1_s_half-triangle1_edge32_length)/triangle1_s_half/(triangle1_s_half-triangle1_edge12_length));
            triangle1_cot=1/(temp_tri1_tan*2/(1-temp_tri1_tan^2));
            EdgeBasedAngleIndex(indexEdgeAtFaceIndex,3)=triangle1_cot;            
        end       
        %得到另外一个三角形
        faceindex2=EdgeAtFaceIndex(indexEdgeAtFaceIndex,4);
        if faceindex2==0
%             display(['边  ',num2str(edge),' 的一边没有面片！']);
        else
            triangle2=surface.TRIV(faceindex2,:); 
            triangle2_1_2=find(ismember(triangle2,edge));%第1、2个点在三角形中的位置
            triangle2_3=find(ismember(triangle2,edge)==0);%第3个点在三角形中的位置
            %得到三角形点的坐标
            triangle2_vertex3=[surface.X(triangle2(1,triangle2_3),1);surface.Y(triangle2(1,triangle2_3),1);surface.Z(triangle2(1,triangle2_3),1)];
            triangle2_vertex1=[surface.X(triangle2(1,triangle2_1_2(1,1)),1);surface.Y(triangle2(1,triangle2_1_2(1,1)),1);surface.Z(triangle2(1,triangle2_1_2(1,1)),1)];    
            triangle2_vertex2=[surface.X(triangle2(1,triangle2_1_2(1,2)),1);surface.Y(triangle2(1,triangle2_1_2(1,2)),1);surface.Z(triangle2(1,triangle2_1_2(1,2)),1)];
            triangle2_edge12=triangle2_vertex2-triangle2_vertex1;%列向量
            triangle2_edge31=triangle2_vertex1-triangle2_vertex3;
            triangle2_edge32=triangle2_vertex2-triangle2_vertex3;
            %得到三角形的边长
            triangle2_edge12_length=sqrt(triangle2_edge12'*triangle2_edge12);
            triangle2_edge31_length=sqrt(triangle2_edge31'*triangle2_edge31);
            triangle2_edge32_length=sqrt(triangle2_edge32'*triangle2_edge32);
            %计算面积  使用半周长的方法来计算三角形面积  http://mathworld.wolfram.com/Semiperimeter.html
            triangle2_s_half=(triangle2_edge12_length+triangle2_edge31_length+triangle2_edge32_length)/2;
            triangle2_Area=sqrt(triangle2_s_half*(triangle2_s_half-triangle2_edge12_length)*(triangle2_s_half-triangle2_edge31_length)*(triangle2_s_half-triangle2_edge32_length)); 
            EdgeBasedAreaIndex(indexEdgeAtFaceIndex,4)=triangle2_Area;
            %计算角度
            %triangle2_Angle=asin(2*triangle2_Area/triangle2_edge31_length/triangle2_edge32_length);% *180/pi  这里角度的单位是弧度
            %半角定理
            temp_tri2_tan=sqrt((triangle2_s_half-triangle2_edge31_length)*(triangle2_s_half-triangle2_edge32_length)/triangle2_s_half/(triangle2_s_half-triangle2_edge12_length));
            triangle2_cot=1/(temp_tri2_tan*2/(1-temp_tri2_tan^2));            
            EdgeBasedAngleIndex(indexEdgeAtFaceIndex,4)=triangle2_cot;
        end        
        if faceindex1==0 && faceindex2==0
            display(['ERROR:  边  ',num2str(edge),' 没有所在的面片！']);
        end        
    end
end
%函数6：找到三角网格最外圈的边界点
function [ IsBoundaryVertex ] = GetBoundaryVertex( surface,NumofBoundary )
%UNTITLED2 此处显示有关此函数的摘要
%   NumofBoundary指去除的边界点的圈数
    M=zeros(size(surface.X,1),size(surface.X,1));
    for i=1:size(surface.TRIV,1)
        v1=surface.TRIV(i,1);
        v2=surface.TRIV(i,2);
        v3=surface.TRIV(i,3);
        M(v1,v2)=M(v1,v2)+1;M(v2,v1)=M(v2,v1)+1;
        M(v1,v3)=M(v1,v3)+1;M(v3,v1)=M(v3,v1)+1;
        M(v2,v3)=M(v2,v3)+1;M(v3,v2)=M(v3,v2)+1;         
    end

    tempBoundary=[];
    for i=1:size(surface.X,1)
        for j=i+1:size(surface.X,1)
            if M(i,j)==1
                [RingsVertexIndexI,FaceIndexI]  = GetVertexInRing(surface,i,NumofBoundary-1);
                [RingsVertexIndexJ,FaceIndexJ]  = GetVertexInRing(surface,j,NumofBoundary-1);
                
                tempBoundary=union(tempBoundary,RingsVertexIndexI);
                tempBoundary=union(tempBoundary,RingsVertexIndexJ);
              
            end
        end
    end
    IsBoundaryVertex=zeros(1,size(surface.X,1));
    for i=1:size(tempBoundary,1)
        IsBoundaryVertex(1,tempBoundary(i,1))=1;
    end    
end
%函数7：求解一个点附近的dA--采用Mixed Cell的方式进行计算
function [dA] = GetVertexdA_MixedCell(surface,CenterVertex_index,FaceInRingIndex,symbol)
    dA=0;
    if symbol==1
        surface.Z=surface.Z*0;%这里将网格投影到XY平面上，然后求相关的面积  针对paper中2D情形进行处理
    end    
    for i=1:size(FaceInRingIndex,2)
        dA=dA+GetVertexdA_MixedCell_TridA(surface,CenterVertex_index,FaceInRingIndex(1,i));
    end
end 
%函数8：根据三角形的边长求解对应的角度
function [deg1,deg2,deg3]=GetDegFromLength(surface,vert1_index,vert2_index,vert3_index)
    vert1=[surface.X(vert1_index,1);surface.Y(vert1_index,1);surface.Z(vert1_index,1)];
    vert2=[surface.X(vert2_index,1);surface.Y(vert2_index,1);surface.Z(vert2_index,1)];
    vert3=[surface.X(vert3_index,1);surface.Y(vert3_index,1);surface.Z(vert3_index,1)];
    L12=sqrt((vert2-vert1)'*(vert2-vert1));
    L13=sqrt((vert3-vert1)'*(vert3-vert1));
    L23=sqrt((vert3-vert2)'*(vert3-vert2));
    deg1=acos((L12^2+L13^2-L23^2)/(2*L12*L13))/pi*180;
    deg2=acos((L12^2+L23^2-L13^2)/(2*L12*L23))/pi*180;
    deg3=acos((L13^2+L23^2-L12^2)/(2*L13*L23))/pi*180;
end
%函数9：给定三角形的顶点序号，求三角形的外心,并求当前三角形分配给vert1_index的面积——锐角三角形的处理
function [dA]=GetCircumcenterAnddA_1(surface,vert1_index,vert2_index,vert3_index)
    vert1=[surface.X(vert1_index,1);surface.Y(vert1_index,1);surface.Z(vert1_index,1)];
    vert2=[surface.X(vert2_index,1);surface.Y(vert2_index,1);surface.Z(vert2_index,1)];
    vert3=[surface.X(vert3_index,1);surface.Y(vert3_index,1);surface.Z(vert3_index,1)];
    Vec12=vert2-vert1;
    Vec13=vert3-vert1;
    Vec23=vert3-vert2;
    
    l12=sqrt(Vec12'*Vec12);
    l13=sqrt(Vec13'*Vec13);
    l23=sqrt(Vec23'*Vec23); 
    
    s=0.5*(l12+l13+l23);
    S=sqrt(s*(s-l12)*(s-l13)*(s-l23));%http://mathworld.wolfram.com/Semiperimeter.html
    R=l12*l13*l23/4/S;%手册P46
    n=cross(Vec12,Vec13);
    Vertical_12=cross(Vec12,n);%边12的垂向量
    Vertical_13=cross(Vec13,n);%边12的垂向量
    
    C12=sqrt((R^2-(l12/2)^2)/(Vertical_12'*Vertical_12));%垂直平分线的向量系数
    C13=sqrt((R^2-(l13/2)^2)/(Vertical_13'*Vertical_13));%垂直平分线的向量系数
        
    l12_Ver=sqrt((C12*Vertical_12)'*(C12*Vertical_12));
    l13_Ver=sqrt((C13*Vertical_13)'*(C13*Vertical_13));

    dA=(l12/2)*l12_Ver/2+(l13/2)*l13_Ver/2;%计算分配的面积
    
%     Circum=vert1+0.5*Vec12+length12*Vertical_12;%计算外心坐标

end

%函数10：根据三角形的顶点序号，求三角形最长边的中点，并求当前三角形分配给vert1_index的面积——钝角三角形的处理
function [dA]=GetCircumcenterAnddA_2(surface,vert1_index,vert2_index,vert3_index)
    vert1=[surface.X(vert1_index,1);surface.Y(vert1_index,1);surface.Z(vert1_index,1)];
    vert2=[surface.X(vert2_index,1);surface.Y(vert2_index,1);surface.Z(vert2_index,1)];
    vert3=[surface.X(vert3_index,1);surface.Y(vert3_index,1);surface.Z(vert3_index,1)];
    Vec12=vert2-vert1;
    Vec13=vert3-vert1;
    Vec23=vert3-vert2;    
    l12=sqrt(Vec12'*Vec12);
    l13=sqrt(Vec13'*Vec13);
    l23=sqrt(Vec23'*Vec23);       
    s=0.5*(l12+l13+l23);
    S=sqrt(s*(s-l12)*(s-l13)*(s-l23));%http://mathworld.wolfram.com/Semiperimeter.html

    [deg1,deg2,deg3]=GetDegFromLength(surface,vert1_index,vert2_index,vert3_index);    
    
    if deg1>=90
        dA=0.5*S;
    else
        dA=0.25*S;
    end
end

%函数11：根据三角形的边长来确定当前三角形中的dA
function [dA]= GetVertexdA_MixedCell_TridA(surface,VertIndex,TriIndex)
    vert1_index=VertIndex;
    vert23_index=find(surface.TRIV(TriIndex,:)>VertIndex|surface.TRIV(TriIndex,:)<VertIndex);
    vert2_index=vert23_index(1,1);
    vert3_index=vert23_index(1,2);    
    deg123=GetDegFromLength(surface,vert1_index,vert2_index,vert3_index);
    maxdeg=max(deg123);
    if maxdeg>=90 && maxdeg<180%当前顶点的角度为钝角
        dA=GetCircumcenterAnddA_2(surface,vert1_index,vert2_index,vert3_index);
    elseif maxdeg>0 && maxdeg<90%当前顶点的角度为锐角
        dA=GetCircumcenterAnddA_1(surface,vert1_index,vert2_index,vert3_index);
    end
end
%函数12：计算当前顶点的法向量-PCA
function [Normal]= GetVertexNormal_PCA(surface,CenterVertex_index,VertexInRing)
    X=zeros(3,size(VertexInRing,2)+1);
    X(1,size(VertexInRing,2)+1)=surface.X(CenterVertex_index,1);
    X(2,size(VertexInRing,2)+1)=surface.Y(CenterVertex_index,1);
    X(3,size(VertexInRing,2)+1)=surface.Z(CenterVertex_index,1);
    for i=1:size(VertexInRing,2)
        X(1,i)=surface.X(VertexInRing(1,i),1);
        X(2,i)=surface.Y(VertexInRing(1,i),1);
        X(3,i)=surface.Z(VertexInRing(1,i),1);
    end
    X_Aver=sum(X,2)/(size(VertexInRing,2)+1);
    X=X-repmat(X_Aver,1,size(VertexInRing,2)+1);
    [V,D]=eig(X*X');%计算结果默认按照特征值从小到大进行排序
    Normal=V(:,1);   
end
%函数13：根据边所跨越的两个面片的索引信息，计算得到两个对角的大小
function [EdgeBasedAngleIndex]=GetEdgeBasedAngle(surface,EdgeAtFaceIndex)
 %索引格式与EdgeAtFaceIndex相似，不同的是对应的3、4列分别为角度的cot和面积
    EdgeBasedAngleIndex=zeros(size(EdgeAtFaceIndex,1),size(EdgeAtFaceIndex,2));%边对应的角度索引    
    EdgeBasedAngleIndex(:,1:2)=EdgeAtFaceIndex(:,1:2);
    for indexEdgeAtFaceIndex=1:size(EdgeAtFaceIndex,1)
        edge=EdgeAtFaceIndex(indexEdgeAtFaceIndex,1:2);%得到一条边       
        %得到一个三角形    
        faceindex1=EdgeAtFaceIndex(indexEdgeAtFaceIndex,3);  
        if faceindex1==0
%             display(['边  ',num2str(edge),' 的一边没有面片！']);
        else
            triangle1=surface.TRIV(faceindex1,:); 
            triangle1_1_2=find(ismember(triangle1,edge));%第1、2个点在三角形中的位置
            triangle1_3=find(ismember(triangle1,edge)==0);%第3个点在三角形中的位置
            %得到三角形点的坐标
            triangle1_vertex3=[surface.X(triangle1(1,triangle1_3),1);surface.Y(triangle1(1,triangle1_3),1);surface.Z(triangle1(1,triangle1_3),1)];
            triangle1_vertex1=[surface.X(triangle1(1,triangle1_1_2(1,1)),1);surface.Y(triangle1(1,triangle1_1_2(1,1)),1);surface.Z(triangle1(1,triangle1_1_2(1,1)),1)];    
            triangle1_vertex2=[surface.X(triangle1(1,triangle1_1_2(1,2)),1);surface.Y(triangle1(1,triangle1_1_2(1,2)),1);surface.Z(triangle1(1,triangle1_1_2(1,2)),1)];
            triangle1_edge12=triangle1_vertex2-triangle1_vertex1;%列向量
            triangle1_edge31=triangle1_vertex1-triangle1_vertex3;
            triangle1_edge32=triangle1_vertex2-triangle1_vertex3;
            %得到三角形的边长
            triangle1_edge12_length=sqrt(triangle1_edge12'*triangle1_edge12);
            triangle1_edge31_length=sqrt(triangle1_edge31'*triangle1_edge31);
            triangle1_edge32_length=sqrt(triangle1_edge32'*triangle1_edge32);
            %计算角度
            triangle1_Angle=acos((triangle1_edge31_length^2+triangle1_edge32_length^2-triangle1_edge12_length^2)/2/triangle1_edge31_length/triangle1_edge32_length)/pi*180;          
            EdgeBasedAngleIndex(indexEdgeAtFaceIndex,3)=triangle1_Angle;%这里为角度的大小，单位为度          
        end       
        %得到另外一个三角形
        faceindex2=EdgeAtFaceIndex(indexEdgeAtFaceIndex,4);
        if faceindex2==0
%             display(['边  ',num2str(edge),' 的一边没有面片！']);
        else
            triangle2=surface.TRIV(faceindex2,:); 
            triangle2_1_2=find(ismember(triangle2,edge));%第1、2个点在三角形中的位置
            triangle2_3=find(ismember(triangle2,edge)==0);%第3个点在三角形中的位置
            %得到三角形点的坐标
            triangle2_vertex3=[surface.X(triangle2(1,triangle2_3),1);surface.Y(triangle2(1,triangle2_3),1);surface.Z(triangle2(1,triangle2_3),1)];
            triangle2_vertex1=[surface.X(triangle2(1,triangle2_1_2(1,1)),1);surface.Y(triangle2(1,triangle2_1_2(1,1)),1);surface.Z(triangle2(1,triangle2_1_2(1,1)),1)];    
            triangle2_vertex2=[surface.X(triangle2(1,triangle2_1_2(1,2)),1);surface.Y(triangle2(1,triangle2_1_2(1,2)),1);surface.Z(triangle2(1,triangle2_1_2(1,2)),1)];
            triangle2_edge12=triangle2_vertex2-triangle2_vertex1;%列向量
            triangle2_edge31=triangle2_vertex1-triangle2_vertex3;
            triangle2_edge32=triangle2_vertex2-triangle2_vertex3;
            %得到三角形的边长
            triangle2_edge12_length=sqrt(triangle2_edge12'*triangle2_edge12);
            triangle2_edge31_length=sqrt(triangle2_edge31'*triangle2_edge31);
            triangle2_edge32_length=sqrt(triangle2_edge32'*triangle2_edge32);            
            %计算角度
            triangle2_Angle=acos((triangle2_edge31_length^2+triangle2_edge32_length^2-triangle2_edge12_length^2)/2/triangle2_edge31_length/triangle2_edge32_length)/pi*180;%这里角度的单位是度              
            EdgeBasedAngleIndex(indexEdgeAtFaceIndex,4)=triangle2_Angle;
        end        
        if faceindex1==0 && faceindex2==0
            display(['ERROR:  边  ',num2str(edge),' 没有所在的面片！']);
        end        
    end
end





