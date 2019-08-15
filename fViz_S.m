function fViz_S(surface,ChiralInv)

for i=1:5
    temp=ChiralInv(i,:);
    % normalization of S
%     temp=temp-(max(temp)+min(temp))/2;
    temp=temp./max(abs(temp));
    surface.I=nthroot(temp,5);
    figure('numbertitle','off','name',['Visualization of S',num2str(i)]);
    surfaceDisplay(surface);
end

end

