function [ ChiralInv ] = ChiralInvariant( surface )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
   
    ChiralInv=zeros(5,size(surface.X,1));%���ڴ洢��������Բ�����

% % ����Ϊ�����ٽ������������Ϻ���ز����Ĵ洢--End   
    NumofBoundary=0;
    disp(['����ȥ��',num2str(NumofBoundary),'��߽�㡤����']);
    [ IsBoundaryVertex ] = GetBoundaryVertex( surface,NumofBoundary );
    NumofOutVert=sum(IsBoundaryVertex);
    NumofBoundary=sum(IsBoundaryVertex);
    parfor index_vertex=1:size(surface.X,1)
        disp(['���㴦����ȣ�',num2str(index_vertex),'/',num2str(size(surface.X,1))]);
        CenterVertex_index=index_vertex;
        %�߽�㲻��������
        if IsBoundaryVertex(1,index_vertex)==1
%             disp(['���㣺',num2str(index_vertex),'Ϊ�߽�㣬��������㣡']);
            continue;
        end
%����Ϊ�����ٽ������������Ϻ󣬽�����ز����ļ���--Begin
        numofRing=2;
        [ VertexInRing,FaceInRingIndex]=GetVertexInRing(surface,CenterVertex_index,numofRing);
        %��ǰ���������
        obj=zeros(size(VertexInRing,2)+1,3);
        
        %���ĵ��������f��ȡֵ
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
       %������Щ�㹹����������Բ�����
       Point=obj';
       Density=ones(1,size(VertexInRing,2)+1);
       
       Inv = CalInv( Point,Density );      
        
    %����Ϊ�����ٽ������������Ϻ󣬽�����ز����ļ���--End
    ChiralInv(:,index_vertex)=Inv;
    end
    disp(['����',num2str(NumofOutVert),'����δ������㣬���б߽����',num2str(NumofBoundary),'��']);
end

%����1�����㵱ǰ��������λ���ܱ�N���ڵĵ�,����Ӧ���ڵ���Ƭ
%�����VertexInRing����ע�ͣ�FaceInRingIndex����ע�͡����߽�Ϊ��������
function [ VertexInRing,FaceInRingIndex]=GetVertexInRing(surface,CenterVertex_index,numofRing)

    VertexInRing=[CenterVertex_index];%���ҽ����һ��������
    ListAdd=zeros(size(surface.TRIV,1),1);%�������Ч�ʣ��Ѿ��ҹ�����Ƭ���ٽ��в���
    %�ҵĴ���
    for indexRing=1:numofRing
        tempRingsVertexIndex=VertexInRing;%��¼�ҵ��Ľ�����м����    
        %�������е���Ƭ
        for indexFace=1:size(surface.TRIV,1)
            %�жϵ�ǰ��Ƭ�Ƿ��Ѿ���ӵ������������
            if ListAdd(indexFace,1)==0
                %�жϵ�ǰ��Ƭ�еĶ����Ƿ����Ѿ��ҵ��Ľ����
                if sum(ismember(surface.TRIV(indexFace,:),VertexInRing))
                    tempRingsVertexIndex=union(tempRingsVertexIndex,surface.TRIV(indexFace,:));
                    ListAdd(indexFace,1)=1;
                end
            end
        end
        VertexInRing=tempRingsVertexIndex; %�������µĲ��ҽ��    
    end    
    VertexInRing=setdiff(VertexInRing,CenterVertex_index);%ȥ�������ĵ�
    FaceInRingIndex=find(ListAdd==1)';%����n���ڵ���Ƭ���
end
%����2��������沢������Ӧ��ƫ����
function [DiffX,DiffY,DiffXX,DiffYY,DiffXY] = GetDiff(surface,CenterVertex_index,VertexInRing,UpOrder)
%��Ҫ�ں���1��ע�͵�ȥ�������ĵ�Ĵ���   �����ڶ���  �������еĵ������ŵ��������в������
    XY=zeros(size(VertexInRing,2),(UpOrder+1)^2);
    if size(XY,1)<size(XY,2)
        disp('��ǰ������С���˷����ʱû�н⣡')
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
    %����Z��X�����ϵ�ƫ����
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
    %����Z��Y�����ϵ�ƫ����  
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
    %����Z��XX�ĵ���
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
    %����Z��YY�ĵ���
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
    %������Z��XY�ĵ���
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
%����2-2������΢�ּ���P179ҳ̩��չ��ʽ��������沢������Ӧ��ƫ����
function [DiffX,DiffY,DiffXX,DiffXY,DiffYY] = GetDiff_2_Taylor(surface,CenterVertex_index,VertexInRing)
    %��Ҫ�ں���1�� �ſ� ȥ�����ĵ�Ĵ���   �����ڶ���  �������ĵ�����굥���ų�����Ϊ����
  
    XY=zeros(size(VertexInRing,2),5);
    if size(XY,1)<size(XY,2)
        disp('��ǰ������С���˷����ʱû�н⣡')
    end
    Z=zeros(size(VertexInRing,2),1);
    %���ĵ��������f��ȡֵ
    fo_x=surface.X(CenterVertex_index,1);
    fo_y=surface.Y(CenterVertex_index,1);
    fo_z=surface.Z(CenterVertex_index,1);
    %���췽����
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



%����3�����㵱ǰ�����1-ring��Ƭ�ϵ�dA--���������������ȼ�
function [dA] = GetVertexdA(surface,CenterVertex_index) 
    [RingsVertexIndex,FaceIndex]  = GetVertexInRing(surface,CenterVertex_index,1);
    %���㵱ǰ�����ܱߵ��edge����Խ����Ƭ�ı��
    [ EdgeAtFaceIndex ] = EdgeAtFace( surface,CenterVertex_index,RingsVertexIndex,FaceIndex);
    %���㵱ǰ�����ܱߵ��edge����Խ�����������εĶԽǺ����
    [EdgeBasedAngleIndex,EdgeBasedAreaIndex]=EdgeBasedAngle_Ares(surface,EdgeAtFaceIndex);
    dA=(sum(EdgeBasedAreaIndex(:,3))+sum(EdgeBasedAreaIndex(:,4)));
end
%����4���ҵ�ָ�������ܱ�1���ڵĵ���ɵı�����Խ����Ƭ�����
%�����EdgeAtFaceIndex����ע��
function [ EdgeAtFaceIndex ] = EdgeAtFace( surface,CenterVertex_index,RingsVertexIndex,FaceIndex)
    %��һ��ΪCenterVertex_index���ڶ���ΪRingsVertexIndex��
    EdgeAtFaceIndex=zeros(size(RingsVertexIndex,2),4);
    EdgeAtFaceIndex(:,1)=ones(size(RingsVertexIndex,2),1)*CenterVertex_index;
    EdgeAtFaceIndex(:,2)=RingsVertexIndex';
    %��3��4��Ϊ��ǰ������Խ����Ƭ�����
    %�������еı�
    for indexCenterVertex_index=1:size(RingsVertexIndex,2)
        tempEdge=EdgeAtFaceIndex(indexCenterVertex_index,1:2);
        signofFaceNum=0;
        %�������е���
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
            display('ERROR:��ǰ�������ڵ��棡')
        elseif signofFaceNum==1
%             display('�߽�edge��')
        end
    end  
end
%����5������EdgeAtFaceIndex����edge����Խ�����������εĶԽǺ����
%�����EdgeBasedAngleIndex����ע�ͣ�EdgeBasedAreaIndex����ע��
function [EdgeBasedAngleIndex,EdgeBasedAreaIndex]=EdgeBasedAngle_Ares(surface,EdgeAtFaceIndex)
    %������ʽ��EdgeAtFaceIndex���ƣ���ͬ���Ƕ�Ӧ��3��4�зֱ�Ϊ�Ƕȵ�cot�����
    EdgeBasedAngleIndex=zeros(size(EdgeAtFaceIndex,1),size(EdgeAtFaceIndex,2));%�߶�Ӧ�ĽǶ�����
    EdgeBasedAreaIndex=zeros(size(EdgeAtFaceIndex,1),size(EdgeAtFaceIndex,2));%�߶�Ӧ���������
    
    EdgeBasedAngleIndex(:,1:2)=EdgeAtFaceIndex(:,1:2);
    EdgeBasedAreaIndex(:,1:2)=EdgeAtFaceIndex(:,1:2);
    
    for indexEdgeAtFaceIndex=1:size(EdgeAtFaceIndex,1)
        edge=EdgeAtFaceIndex(indexEdgeAtFaceIndex,1:2);%�õ�һ����
        
        %�õ�һ��������    
        faceindex1=EdgeAtFaceIndex(indexEdgeAtFaceIndex,3);  
        if faceindex1==0
%             display(['��  ',num2str(edge),' ��һ��û����Ƭ��']);
        else
            triangle1=surface.TRIV(faceindex1,:); 
            triangle1_1_2=find(ismember(triangle1,edge));%��1��2�������������е�λ��
            triangle1_3=find(ismember(triangle1,edge)==0);%��3�������������е�λ��
            %�õ������ε������
            triangle1_vertex3=[surface.X(triangle1(1,triangle1_3),1);surface.Y(triangle1(1,triangle1_3),1);surface.Z(triangle1(1,triangle1_3),1)];
            triangle1_vertex1=[surface.X(triangle1(1,triangle1_1_2(1,1)),1);surface.Y(triangle1(1,triangle1_1_2(1,1)),1);surface.Z(triangle1(1,triangle1_1_2(1,1)),1)];    
            triangle1_vertex2=[surface.X(triangle1(1,triangle1_1_2(1,2)),1);surface.Y(triangle1(1,triangle1_1_2(1,2)),1);surface.Z(triangle1(1,triangle1_1_2(1,2)),1)];
            triangle1_edge12=triangle1_vertex2-triangle1_vertex1;%������
            triangle1_edge31=triangle1_vertex1-triangle1_vertex3;
            triangle1_edge32=triangle1_vertex2-triangle1_vertex3;
            %�õ������εı߳�
            triangle1_edge12_length=sqrt(triangle1_edge12'*triangle1_edge12);
            triangle1_edge31_length=sqrt(triangle1_edge31'*triangle1_edge31);
            triangle1_edge32_length=sqrt(triangle1_edge32'*triangle1_edge32);
            %�������  ʹ�ð��ܳ��ķ������������������  http://mathworld.wolfram.com/Semiperimeter.html
            triangle1_s_half=(triangle1_edge12_length+triangle1_edge31_length+triangle1_edge32_length)/2;
            triangle1_Area=sqrt(triangle1_s_half*(triangle1_s_half-triangle1_edge12_length)*(triangle1_s_half-triangle1_edge31_length)*(triangle1_s_half-triangle1_edge32_length)); 
            EdgeBasedAreaIndex(indexEdgeAtFaceIndex,3)=triangle1_Area;
            %����Ƕ�
            %triangle1_Angle=asin(2*triangle1_Area/triangle1_edge31_length/triangle1_edge32_length);% *180/pi  ����Ƕȵĵ�λ�ǻ���
            %��Ƕ���
            temp_tri1_tan=sqrt((triangle1_s_half-triangle1_edge31_length)*(triangle1_s_half-triangle1_edge32_length)/triangle1_s_half/(triangle1_s_half-triangle1_edge12_length));
            triangle1_cot=1/(temp_tri1_tan*2/(1-temp_tri1_tan^2));
            EdgeBasedAngleIndex(indexEdgeAtFaceIndex,3)=triangle1_cot;            
        end       
        %�õ�����һ��������
        faceindex2=EdgeAtFaceIndex(indexEdgeAtFaceIndex,4);
        if faceindex2==0
%             display(['��  ',num2str(edge),' ��һ��û����Ƭ��']);
        else
            triangle2=surface.TRIV(faceindex2,:); 
            triangle2_1_2=find(ismember(triangle2,edge));%��1��2�������������е�λ��
            triangle2_3=find(ismember(triangle2,edge)==0);%��3�������������е�λ��
            %�õ������ε������
            triangle2_vertex3=[surface.X(triangle2(1,triangle2_3),1);surface.Y(triangle2(1,triangle2_3),1);surface.Z(triangle2(1,triangle2_3),1)];
            triangle2_vertex1=[surface.X(triangle2(1,triangle2_1_2(1,1)),1);surface.Y(triangle2(1,triangle2_1_2(1,1)),1);surface.Z(triangle2(1,triangle2_1_2(1,1)),1)];    
            triangle2_vertex2=[surface.X(triangle2(1,triangle2_1_2(1,2)),1);surface.Y(triangle2(1,triangle2_1_2(1,2)),1);surface.Z(triangle2(1,triangle2_1_2(1,2)),1)];
            triangle2_edge12=triangle2_vertex2-triangle2_vertex1;%������
            triangle2_edge31=triangle2_vertex1-triangle2_vertex3;
            triangle2_edge32=triangle2_vertex2-triangle2_vertex3;
            %�õ������εı߳�
            triangle2_edge12_length=sqrt(triangle2_edge12'*triangle2_edge12);
            triangle2_edge31_length=sqrt(triangle2_edge31'*triangle2_edge31);
            triangle2_edge32_length=sqrt(triangle2_edge32'*triangle2_edge32);
            %�������  ʹ�ð��ܳ��ķ������������������  http://mathworld.wolfram.com/Semiperimeter.html
            triangle2_s_half=(triangle2_edge12_length+triangle2_edge31_length+triangle2_edge32_length)/2;
            triangle2_Area=sqrt(triangle2_s_half*(triangle2_s_half-triangle2_edge12_length)*(triangle2_s_half-triangle2_edge31_length)*(triangle2_s_half-triangle2_edge32_length)); 
            EdgeBasedAreaIndex(indexEdgeAtFaceIndex,4)=triangle2_Area;
            %����Ƕ�
            %triangle2_Angle=asin(2*triangle2_Area/triangle2_edge31_length/triangle2_edge32_length);% *180/pi  ����Ƕȵĵ�λ�ǻ���
            %��Ƕ���
            temp_tri2_tan=sqrt((triangle2_s_half-triangle2_edge31_length)*(triangle2_s_half-triangle2_edge32_length)/triangle2_s_half/(triangle2_s_half-triangle2_edge12_length));
            triangle2_cot=1/(temp_tri2_tan*2/(1-temp_tri2_tan^2));            
            EdgeBasedAngleIndex(indexEdgeAtFaceIndex,4)=triangle2_cot;
        end        
        if faceindex1==0 && faceindex2==0
            display(['ERROR:  ��  ',num2str(edge),' û�����ڵ���Ƭ��']);
        end        
    end
end
%����6���ҵ�������������Ȧ�ı߽��
function [ IsBoundaryVertex ] = GetBoundaryVertex( surface,NumofBoundary )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   NumofBoundaryָȥ���ı߽���Ȧ��
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
%����7�����һ���㸽����dA--����Mixed Cell�ķ�ʽ���м���
function [dA] = GetVertexdA_MixedCell(surface,CenterVertex_index,FaceInRingIndex,symbol)
    dA=0;
    if symbol==1
        surface.Z=surface.Z*0;%���ｫ����ͶӰ��XYƽ���ϣ�Ȼ������ص����  ���paper��2D���ν��д���
    end    
    for i=1:size(FaceInRingIndex,2)
        dA=dA+GetVertexdA_MixedCell_TridA(surface,CenterVertex_index,FaceInRingIndex(1,i));
    end
end 
%����8�����������εı߳�����Ӧ�ĽǶ�
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
%����9�����������εĶ�����ţ��������ε�����,����ǰ�����η����vert1_index�����������������εĴ���
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
    R=l12*l13*l23/4/S;%�ֲ�P46
    n=cross(Vec12,Vec13);
    Vertical_12=cross(Vec12,n);%��12�Ĵ�����
    Vertical_13=cross(Vec13,n);%��12�Ĵ�����
    
    C12=sqrt((R^2-(l12/2)^2)/(Vertical_12'*Vertical_12));%��ֱƽ���ߵ�����ϵ��
    C13=sqrt((R^2-(l13/2)^2)/(Vertical_13'*Vertical_13));%��ֱƽ���ߵ�����ϵ��
        
    l12_Ver=sqrt((C12*Vertical_12)'*(C12*Vertical_12));
    l13_Ver=sqrt((C13*Vertical_13)'*(C13*Vertical_13));

    dA=(l12/2)*l12_Ver/2+(l13/2)*l13_Ver/2;%�����������
    
%     Circum=vert1+0.5*Vec12+length12*Vertical_12;%������������

end

%����10�����������εĶ�����ţ�����������ߵ��е㣬����ǰ�����η����vert1_index����������۽������εĴ���
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

%����11�����������εı߳���ȷ����ǰ�������е�dA
function [dA]= GetVertexdA_MixedCell_TridA(surface,VertIndex,TriIndex)
    vert1_index=VertIndex;
    vert23_index=find(surface.TRIV(TriIndex,:)>VertIndex|surface.TRIV(TriIndex,:)<VertIndex);
    vert2_index=vert23_index(1,1);
    vert3_index=vert23_index(1,2);    
    deg123=GetDegFromLength(surface,vert1_index,vert2_index,vert3_index);
    maxdeg=max(deg123);
    if maxdeg>=90 && maxdeg<180%��ǰ����ĽǶ�Ϊ�۽�
        dA=GetCircumcenterAnddA_2(surface,vert1_index,vert2_index,vert3_index);
    elseif maxdeg>0 && maxdeg<90%��ǰ����ĽǶ�Ϊ���
        dA=GetCircumcenterAnddA_1(surface,vert1_index,vert2_index,vert3_index);
    end
end
%����12�����㵱ǰ����ķ�����-PCA
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
    [V,D]=eig(X*X');%������Ĭ�ϰ�������ֵ��С�����������
    Normal=V(:,1);   
end
%����13�����ݱ�����Խ��������Ƭ��������Ϣ������õ������ԽǵĴ�С
function [EdgeBasedAngleIndex]=GetEdgeBasedAngle(surface,EdgeAtFaceIndex)
 %������ʽ��EdgeAtFaceIndex���ƣ���ͬ���Ƕ�Ӧ��3��4�зֱ�Ϊ�Ƕȵ�cot�����
    EdgeBasedAngleIndex=zeros(size(EdgeAtFaceIndex,1),size(EdgeAtFaceIndex,2));%�߶�Ӧ�ĽǶ�����    
    EdgeBasedAngleIndex(:,1:2)=EdgeAtFaceIndex(:,1:2);
    for indexEdgeAtFaceIndex=1:size(EdgeAtFaceIndex,1)
        edge=EdgeAtFaceIndex(indexEdgeAtFaceIndex,1:2);%�õ�һ����       
        %�õ�һ��������    
        faceindex1=EdgeAtFaceIndex(indexEdgeAtFaceIndex,3);  
        if faceindex1==0
%             display(['��  ',num2str(edge),' ��һ��û����Ƭ��']);
        else
            triangle1=surface.TRIV(faceindex1,:); 
            triangle1_1_2=find(ismember(triangle1,edge));%��1��2�������������е�λ��
            triangle1_3=find(ismember(triangle1,edge)==0);%��3�������������е�λ��
            %�õ������ε������
            triangle1_vertex3=[surface.X(triangle1(1,triangle1_3),1);surface.Y(triangle1(1,triangle1_3),1);surface.Z(triangle1(1,triangle1_3),1)];
            triangle1_vertex1=[surface.X(triangle1(1,triangle1_1_2(1,1)),1);surface.Y(triangle1(1,triangle1_1_2(1,1)),1);surface.Z(triangle1(1,triangle1_1_2(1,1)),1)];    
            triangle1_vertex2=[surface.X(triangle1(1,triangle1_1_2(1,2)),1);surface.Y(triangle1(1,triangle1_1_2(1,2)),1);surface.Z(triangle1(1,triangle1_1_2(1,2)),1)];
            triangle1_edge12=triangle1_vertex2-triangle1_vertex1;%������
            triangle1_edge31=triangle1_vertex1-triangle1_vertex3;
            triangle1_edge32=triangle1_vertex2-triangle1_vertex3;
            %�õ������εı߳�
            triangle1_edge12_length=sqrt(triangle1_edge12'*triangle1_edge12);
            triangle1_edge31_length=sqrt(triangle1_edge31'*triangle1_edge31);
            triangle1_edge32_length=sqrt(triangle1_edge32'*triangle1_edge32);
            %����Ƕ�
            triangle1_Angle=acos((triangle1_edge31_length^2+triangle1_edge32_length^2-triangle1_edge12_length^2)/2/triangle1_edge31_length/triangle1_edge32_length)/pi*180;          
            EdgeBasedAngleIndex(indexEdgeAtFaceIndex,3)=triangle1_Angle;%����Ϊ�ǶȵĴ�С����λΪ��          
        end       
        %�õ�����һ��������
        faceindex2=EdgeAtFaceIndex(indexEdgeAtFaceIndex,4);
        if faceindex2==0
%             display(['��  ',num2str(edge),' ��һ��û����Ƭ��']);
        else
            triangle2=surface.TRIV(faceindex2,:); 
            triangle2_1_2=find(ismember(triangle2,edge));%��1��2�������������е�λ��
            triangle2_3=find(ismember(triangle2,edge)==0);%��3�������������е�λ��
            %�õ������ε������
            triangle2_vertex3=[surface.X(triangle2(1,triangle2_3),1);surface.Y(triangle2(1,triangle2_3),1);surface.Z(triangle2(1,triangle2_3),1)];
            triangle2_vertex1=[surface.X(triangle2(1,triangle2_1_2(1,1)),1);surface.Y(triangle2(1,triangle2_1_2(1,1)),1);surface.Z(triangle2(1,triangle2_1_2(1,1)),1)];    
            triangle2_vertex2=[surface.X(triangle2(1,triangle2_1_2(1,2)),1);surface.Y(triangle2(1,triangle2_1_2(1,2)),1);surface.Z(triangle2(1,triangle2_1_2(1,2)),1)];
            triangle2_edge12=triangle2_vertex2-triangle2_vertex1;%������
            triangle2_edge31=triangle2_vertex1-triangle2_vertex3;
            triangle2_edge32=triangle2_vertex2-triangle2_vertex3;
            %�õ������εı߳�
            triangle2_edge12_length=sqrt(triangle2_edge12'*triangle2_edge12);
            triangle2_edge31_length=sqrt(triangle2_edge31'*triangle2_edge31);
            triangle2_edge32_length=sqrt(triangle2_edge32'*triangle2_edge32);            
            %����Ƕ�
            triangle2_Angle=acos((triangle2_edge31_length^2+triangle2_edge32_length^2-triangle2_edge12_length^2)/2/triangle2_edge31_length/triangle2_edge32_length)/pi*180;%����Ƕȵĵ�λ�Ƕ�              
            EdgeBasedAngleIndex(indexEdgeAtFaceIndex,4)=triangle2_Angle;
        end        
        if faceindex1==0 && faceindex2==0
            display(['ERROR:  ��  ',num2str(edge),' û�����ڵ���Ƭ��']);
        end        
    end
end





