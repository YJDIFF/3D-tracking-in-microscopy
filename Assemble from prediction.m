%% Section 1: Reading cropped predictions and features and generates full-image prediction & deep feature map
clear
load('F:\Mo\my3D_matlab\Tracking\colormap.mat','map')
for time=1:69
    disp(time)
    tt=num2str(time);
    addr=strcat('F:\Mo\my3D_matlab\Prediction\FC-DenseNet\',tt,'\');
    addr2=strcat('F:\Mo\my3D_matlab\Tracking\',tt,'\');
    if ~exist(addr2,'dir')
        mkdir(addr2);
    end
    Files1=dir(strcat(addr,'*.nii'));
    Fullsize=zeros(280,512,13);
    Fullsize_regression=zeros(280,512,13);
    Weights=zeros(280,512,13,20);
    %Weights=zeros(280,512,13,256);
    for it=1:3:length(Files1) % prediction class
        clo=32*floor((it-1)/3/8)+1;
        row=35*rem((it-1)/3,8)+1;
        V = niftiread(strcat(addr,Files1(it).name));
        V=1-V;
        V2= uint8(V);
        V3= imresize3(V2, [35 32 13],'linear');
        Fullsize(row:row+34,clo:clo+31,:)=V3;
    end
%     for it=3:3:length(Files1) % regression score
%         clo=32*floor((it-2)/3/8)+1;
%         row=35*rem((it-2)/3,8)+1;
%         V = niftiread(strcat(addr,Files1(it).name));
%         V2= double(V(:,:,:,1));
%         V3= imresize3(V2, [35 32 13],'linear');
%         Fullsize_regression(row:row+34,clo:clo+31,:)=V3;
%     end
    for it=2:3:length(Files1) % deep feature maps
        display(it/length(Files1))
        clo=32*floor((it-2)/3/8)+1;
        row=35*rem((it-2)/3,8)+1;
        V = niftiread(strcat(addr,Files1(it).name));
        for iy=1:64  %%all=256
            V2= double(V(:,:,:,iy));
            %V3= imresize3(V2, [35 32 13],'linear');
            Weights(row:row+34,clo:clo+31,:,iy)=V2;
        end
    end
    
    Fullsize2=logical(Fullsize);
    for it=1:13
        img=Fullsize2(:,:,it);
        [f,orgnum] = bwlabel(img);
        g = regionprops(f, 'Area');
        area_values=[g.Area];
        % idx = find ((25<= area_values & 800>= area_values));
        idx = find ((5<= area_values));
        h = ismember (f, idx);
        Fullsize2(:,:,it)=h;
    end
    Fullsize2=double(Fullsize2);
    
    mor_ref=zeros(280,512,13);
    addr3=strcat('F:\Mo\ECad_image_sequence\');
    Files3=dir(strcat(addr3,'*.tif'));
    for i1=13*(time-1)+1:13*(time-1)+13
        img=imread(strcat(addr3,Files3(i1).name));
        img=img*(255/45);
%         img=imgaussfilt(img,1);
        mor_ref(:,:,i1-13*(time-1))=img;
    end
            
    %% exclude object that z=1
    stack_after=Fullsize2;
    [y, x, z] = size(Fullsize);
  
    stack_after_BW=logical(stack_after);
    [stack_after_label,orgnum]=bwlabeln(Fullsize2);
    CC = bwconncomp(Fullsize2,6);
    stats = regionprops3(CC,'BoundingBox','VoxelList','ConvexHull','Centroid');
    j=height(stats);
    for i=1:j %delect single layer object
        if stats.BoundingBox(i,6)==1 || stats.BoundingBox(i,4)==1 || stats.BoundingBox(i,5)==1 
            stats.BoundingBox(i,:) = 0;
            %stack_after_BW(stack_after_label==i)=0;  
            %stack_after_label(stack_after_label==i)=0;  
            [a,b]=size(stats.VoxelList{i,1});
            for m=1:a
                stack_after(stats.VoxelList{i,1}(m,2),stats.VoxelList{i,1}(m,1),stats.VoxelList{i,1}(m,3))=0;
                Fullsize2(stats.VoxelList{i,1}(m,2),stats.VoxelList{i,1}(m,1),stats.VoxelList{i,1}(m,3))=0;
            end
            stats.VoxelList{i,1}= [0,0,0];
        end
        
    end
    [stack_after_label,orgnum]=bwlabeln(stack_after, 6);
    CC = bwconncomp(stack_after,6);
    stats1 = regionprops3(CC,'BoundingBox','VoxelList','ConvexHull','Centroid');
    %stack_after(stack_after==0)=nan;
    niftiwrite(stack_after_label,strcat(addr2,'Fullsize_label','_',tt,'.nii'));
    for i=1:height(stats1)
        b=stats1.VoxelList{i,1};
        [x,y]=size(b);
        for i1=1:x
            stack_after_label(stack_after_label(b(i1,2),b(i1,1),b(i1,3))>0)=i;
        end
    end
    
    % draw 3D segmentation and labelling figures
    %stack_after_label(stack_after_label~=57)=0;
    stack_after_label(stack_after_label==0)=nan;
    h=figure;
    [X,Y,Z] = ndgrid(1:size(stack_after_label,1), 1:size(stack_after_label,2), 1:size(stack_after_label,3));
    pointsize = 5;
%     scatter3(X(:), Y(:), Z(:), pointsize, stack_after_label(:),'filled');
    colormap(map);
    % zlim([0 100]);
    colorbar;
    hold on
    grid on
    
    
    %% fill object with new assigned ID
    %% generate 3D figure
    
    for i=1:height(stats1)

        b=stats1.VoxelList{i,1};

        k = boundary(b);
        hold on
        value=i;
        Registration(value,1)=value;
        Registration(value,2:4)=stats1.Centroid(i,:);
%         if isnan(stack_after_label(round(stats1.Centroid(i,2)),round(stats1.Centroid(i,1)),round(stats1.Centroid(i,3))))==0
% %             value=stack_after_label(round(stats1.Centroid(i,2)),round(stats1.Centroid(i,1)),round(stats1.Centroid(i,3)));
%             value=i;
%             Registration(value,1)=value;
%             Registration(value,2:4)=stats1.Centroid(i,:);
%         else           
%             value=0;
%             while(value==0)
%                 a=ceil(rand*length(b(:,1)));
%                 if isnan(stack_after_label(round(b(a,2)),round(b(a,1)),round(b(a,3))))==0
%                     value=stack_after_label(round(b(a,2)),round(b(a,1)),round(b(a,3)));
%                     Registration(value,1)=value;
%                     Registration(value,2:4)=stats1.Centroid(i,:);
%                 end
%             end
%                 
%         end
        trisurf(k,b(:,1),b(:,2),b(:,3),'Facecolor',map(value,1:3),'FaceAlpha',0.3,'Edgecolor',map(value,1:3),'EdgeAlpha',0.3)
        text(b(end,1),b(end,2),b(end,3), num2str(value), 'Rotation',+15)
    end
    
    
    hold off
    view([0 0 1]);
    set(gca, 'YDir','reverse')
    xlim([0 512]);
    ylim([0 280]);
    zlim([0 13]);
    niftiwrite(Fullsize2,strcat(addr2,'Fullsize','_',tt,'.nii'));
    niftiwrite(Registration,strcat(addr2,'Registration','_',tt,'.nii'));
%    niftiwrite(Fullsize_regression,strcat(addr2,'Fullsize_regression','_',tt,'.nii'));
    niftiwrite(Weights,strcat(addr2,'Weights','_',tt,'.nii'));
    %niftiwrite(Weights,strcat(addr2,'Weights','_',tt,'.nii'));
    savefig(h,strcat(addr2,tt,'_3Dconnection2.fig'));
    saveas(h,strcat(addr2,tt,'_3Dconnection2.png'))
    close(h);
end
disp('finish')
