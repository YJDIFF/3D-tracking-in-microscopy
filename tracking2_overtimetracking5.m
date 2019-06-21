%% write full segmentation file
clear
for time=1:69
    addr=strcat('F:\Mo\my3D_1\Tracking\',num2str(time),'\');
    F=niftiread(strcat(addr,'Fullsize','_',num2str(time),'.nii'));
    F2(:,:,:,time)=F;
end
niftiwrite(F2,strcat('F:\Mo\my3D_1\Tracking\','FullSegmentation','.nii'));

%% plot single time segmetnation result
clear
time=60;
load('F:\Mo\my3D_matlab\Tracking\colormap.mat','map');
addr1=strcat('F:\Mo\my3D_matlab\Tracking\',num2str(time),'\');
Fullsize_1 = niftiread(strcat(addr1,'Fullsize_label_',num2str(time),'.nii'));
if time<10
    tt=strcat('00',num2str(time));
else
    tt=strcat('0',num2str(time));
end

threeDimg = niftiread(strcat('F:\Mo\my3D_matlab\3Dimage\','imgz',tt,'.nii'));

stats = regionprops3(Fullsize_1,'BoundingBox','VoxelList','ConvexHull','Centroid');
Fullsize_1(Fullsize_1>1)=1;
Fullsize_1(Fullsize_1==0)=nan;
h=figure;
% [X,Y,Z] = ndgrid(1:size(Fullsize_1,1), 1:size(Fullsize_1,2), 1:size(Fullsize_1,3));
xlim([0 512]);%0-512
ylim([0 280]);%0-280
zlim([0 13]);%0-13
% hold on
% grid on
for z=1:1:13
  img = imrotate(threeDimg(:,:,z),0);
  g = hgtransform('Matrix',makehgtform('translate',[0 0 z]));
  ii=image(g,img);
  ii.AlphaData = 0.15;
  hold on
end
for i=1:height(stats)
    
    b=stats.VoxelList{i,1};
    
    k = boundary(b);
    hold on
%     trisurf(k,b(:,1),b(:,2),b(:,3),'Facecolor',map(i,1:3),'FaceAlpha',0.5,'Edgecolor',map(i,1:3),'EdgeAlpha',0.5)
    trisurf(k,b(:,1),b(:,2),b(:,3),'Facecolor','r','FaceAlpha',0.5,'Edgecolor','r','EdgeAlpha',0.5)
%     text(b(end,1),b(end,2),b(end,3),num2str(i),'Color','black', 'Rotation',+15);
end
set(gca, 'XDir','reverse')
% set(gca, 'ZDir','reverse')
view([0.5 2.0 8]);
grid
%%
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
    for it=1:3:length(Files1)
        clo=32*floor((it-1)/3/8)+1;
        row=35*rem((it-1)/3,8)+1;
        V = niftiread(strcat(addr,Files1(it).name));
        V=1-V;
        V2= uint8(V);
        V3= imresize3(V2, [35 32 13],'linear');
        Fullsize(row:row+34,clo:clo+31,:)=V3;
    end
%     for it=3:3:length(Files1)
%         clo=32*floor((it-2)/3/8)+1;
%         row=35*rem((it-2)/3,8)+1;
%         V = niftiread(strcat(addr,Files1(it).name));
%         V2= double(V(:,:,:,1));
%         V3= imresize3(V2, [35 32 13],'linear');
%         Fullsize_regression(row:row+34,clo:clo+31,:)=V3;
%     end
    for it=2:3:length(Files1)
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
%     %Threshold
%     maxIntensity=max(max(max(mor_ref)));
%     Thre=0.4*maxIntensity;
%     Fullsize2(mor_ref<Thre)=0;
        
    
    
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
%% calculating corelation 
clear
load('F:\Mo\my3D_1\Tracking\colormap.mat','map')
disp('start')
padding=20;
time=clock;
filename = strcat('F:\Mo\my3D_1\Tracking\TrackingID',num2str(time),'.xls');
xlswrite(filename,cellstr('time'),2, 'A1');
xlswrite(filename,cellstr('old'),2, 'B1');
xlswrite(filename,cellstr('new'),2, 'C1');
xlswrite(filename,cellstr('split'),2, 'D1');
xlswrite(filename,cellstr('fusion'),2, 'E1');
trackbackT=2;

% load('F:\Mo\my3D_1\Tracking\colormap.mat','map')
% filename = strcat('F:\Mo\my3D_1\Tracking\TrackingID2019              2             27             16             37          7.995.xls');
% load(strcat('F:\Mo\my3D_1\Tracking\','xlswriter1','.mat'),'xlswriter1');
% load(strcat('F:\Mo\my3D_1\Tracking\','xlswriter2','.mat'),'xlswriter2');
% load(strcat('F:\Mo\my3D_1\Tracking\','xlswriter3','.mat'),'xlswriter3');
% load(strcat('F:\Mo\my3D_1\Tracking\','xlswriter4','.mat'),'xlswriter4');
% load(strcat('F:\Mo\my3D_1\Tracking\','xlswriter5','.mat'),'xlswriter5');
% load(strcat('F:\Mo\my3D_1\Tracking\','xlswriter6','.mat'),'xlswriter6');
% load(strcat('F:\Mo\my3D_1\Tracking\','xlswriter7','.mat'),'xlswriter7');
% load(strcat('F:\Mo\my3D_1\Tracking\','xlswriter8','.mat'),'xlswriter8');
% load(strcat('F:\Mo\my3D_1\Tracking\','xlswriter9','.mat'),'xlswriter9');
depth=64;
spatial_extend_matrix=zeros(10,10,3,depth);
for i1=1:10
    for i2=1:10
        for i3=1:3
            spatial_extend_matrix(i1,i2,i3,:)=exp(((i1-5)+(i2-5)+(i3-2))/20);
        end
    end
end
for time=1:68
    tic;
    disp(strcat('time point: ', num2str(time)))
    t1=num2str(time);
    t2=num2str(time+1);
    xlswriter1(1,(time)*2-1)=cellstr(t1);
    xlswriter1(1,(time)*2)=cellstr(t2);
    xlswriter3(1,(time)*2)=cellstr(t2);
    xlswriter4(1,(time)*2)=cellstr(t2);
    xlswriter5(1,(time)*2)=cellstr(t2);
    xlswriter6(1,(time)*2)=cellstr(t2);
    xlswriter7(1,(time)*2)=cellstr(t2);
    xlswriter8(1,(time)*2)=cellstr(t2);
    addr1=strcat('F:\Mo\my3D_1\Tracking\',t1,'\');
    addr2=strcat('F:\Mo\my3D_1\Tracking\',t2,'\');
    Files1=dir(strcat(addr1,'*.nii'));
    Files2=dir(strcat(addr2,'*.nii'));
%     if time<trackbackT+1 %calculating correlation
%         for i1=1:time
%             Fullsize_2 = niftiread(strcat(addr2,'Fullsize_label_',t2,'.nii'));
%             Fullsize_regression_2 = niftiread(strcat(addr2,'Weights_',t2,'.nii'));
%             t1=num2str(time+1-i1);
%             addr1=strcat('F:\Mo\my3D_1\Tracking\',t1,'\');
%             if i1==time
%                 Fullsize_1 = niftiread(strcat(addr1,'Fullsize_label_',t1,'.nii'));
%                 Fullsize_regression_1 = niftiread(strcat(addr1,'Weights_',t1,'.nii'));
%             else
%                 Fullsize_1 = niftiread(strcat(addr1,'Fullsize_2_aftertracking_',t1,'.nii'));
%                 Fullsize_regression_1 = niftiread(strcat(addr1,'Weights_',t1,'.nii'));
%             end
%             correlation(Fullsize_1,Fullsize_2,Fullsize_regression_1,Fullsize_regression_2,t2,i1,spatial_extend_matrix);
%         end
%     else
%         for i1=1:trackbackT
%             Fullsize_2 = niftiread(strcat(addr2,'Fullsize_label_',t2,'.nii'));
%             Fullsize_regression_2 = niftiread(strcat(addr2,'Weights_',t2,'.nii'));
%             t1=num2str(time+1-i1);
%             addr1=strcat('F:\Mo\my3D_1\Tracking\',t1,'\');
%             Fullsize_1 = niftiread(strcat(addr1,'Fullsize_2_aftertracking_',t1,'.nii'));
%             Fullsize_regression_1 = niftiread(strcat(addr1,'Weights_',t1,'.nii'));
%             correlation(Fullsize_1,Fullsize_2,Fullsize_regression_1,Fullsize_regression_2,t2,i1,spatial_extend_matrix);
%         end
%     end
%         
%     
%     
% clear Fullsize_1 Fullsize_regression_1 Fullsize_2 Fullsize_regression_2  Fullsize_1_padding Fullsize_2_padding ...
%         Fullsize_regression_1_padding Fullsize_regression_2_padding Fullsize_1_label Fullsize_2_label
    
    % plot tracking
    t1=num2str(time);
    t2=num2str(time+1);

    correlation_map_padding_show1 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t2,'\correlation_map_padding_show_traceback1_',t2,'.nii'));
    correlation_map_padding_hide1 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t2,'\','correlation_map_padding_show_traceback1','_',t2,'.nii'));
    
    if time<trackbackT+1 && time>1
        for i1=2:time
            Registration1 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t1,'\Registration2_tracking_',t1,'.nii'));
            correlation_map_padding_show1_2 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t2,'\','correlation_map_padding_show_traceback',num2str(i1),'_',t2,'.nii'));
            correlation_map_padding_hide1_2 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t2,'\','correlation_map_padding_hide_traceback',num2str(i1),'_',t2,'.nii'));
            for i2=1:320
                for i3=1:552
                    for i4=1:17
                        if correlation_map_padding_hide1(i2,i3,i4)<correlation_map_padding_hide1_2(i2,i3,i4) && correlation_map_padding_show1_2(i2,i3,i4)~=0
                            correlation_map_padding_show1(i2,i3,i4)=correlation_map_padding_show1_2(i2,i3,i4);
                        end
                    end
                end
            end
%             correlation_map_padding_show1(correlation_map_padding_show1==0)=correlation_map_padding_show1_2(correlation_map_padding_show1==0);
        end
    elseif time >=trackbackT+1 && time>1
        for i1=2:trackbackT
            Registration1 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t1,'\Registration2_tracking_',t1,'.nii'));
            correlation_map_padding_show1_2 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t2,'\','correlation_map_padding_show_traceback',num2str(i1),'_',t2,'.nii'));
            correlation_map_padding_hide1_2 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t2,'\','correlation_map_padding_hide_traceback',num2str(i1),'_',t2,'.nii'));
            for i2=1:320
                for i3=1:552
                    for i4=1:17
                        if correlation_map_padding_hide1(i2,i3,i4)<correlation_map_padding_hide1_2(i2,i3,i4) && correlation_map_padding_show1_2(i2,i3,i4)~=0
                            correlation_map_padding_show1(i2,i3,i4)=correlation_map_padding_show1_2(i2,i3,i4);
                        end
                    end
                end
            end
%                         correlation_map_padding_show1(correlation_map_padding_show1==0)=correlation_map_padding_show1_2(correlation_map_padding_show1==0);
        end
    else
        Registration1 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t1,'\Registration_',t1,'.nii'));
    end
    
    Fullsize_2 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',t2,'\Fullsize_label_',t2,'.nii'));
    %Fullsize_1 = niftiread(strcat('F:\Mo\my3D_1\Tracking\',num2str(time-1),'\Fullsize_2_aftertracking_',num2str(time-1),'.nii'));
    Fullsize_2_2= zeros(size(Fullsize_2));
%     Fullsize_2_2(Fullsize_2_2>0)=0;
    correlation_map_padding_show2 =correlation_map_padding_show1(21:end-20,21:end-20,3:15);
    Fullsize_2_mark=correlation_map_padding_show2;
    if time>1
        correlation_map_padding_show2_2 =correlation_map_padding_show1_2(21:end-20,21:end-20,3:15);
        Fullsize_1 =correlation_map_padding_show2_2;
        Fullsize_1(Fullsize_1==0)=nan;
        detector_fusion_old=load(strcat('F:\Mo\my3D_1\Tracking\',t1,'\fusion_tracking_',t1,'.mat'),'detector3_fusion');
        for i1=2:2:size(detector_fusion_old.detector3_fusion,1)
            detector_fusion_old.detector3_fusion(i1,:)=0;
        end
    end
    Fullsize_2_mark(Fullsize_2==0)=0;

    clear correlation_map_padding_show1 correlation_map_padding_show1_2 correlation_map_padding_hide1 correlation_map_padding_hide1_2
% disp('largest num of corr %d/nlargest num of reg %d:')
% disp(max(max(max(correlation_map_padding_show2))));
% disp(max(Registration1(:,1)));
    
%draw figures starts  
% correlation_map_padding_show2(Fullsize_2==0)=0;
%     stack_after_BW=logical(correlation_map_padding_show2);
%     stats = regionprops3(stack_after_BW,'BoundingBox','VoxelList','ConvexHull');
%     disp('Draw overlap figure:')
% 
%     correlation_map_padding_show2(correlation_map_padding_show2==0)=nan;
%     [X,Y,Z] = ndgrid(1:size(correlation_map_padding_show2,1), 1:size(correlation_map_padding_show2,2), 1:size(correlation_map_padding_show2,3));
%     pointsize = 5;
%     h1=figure;
%     scatter3(X(:), Y(:), Z(:), pointsize, correlation_map_padding_show2(:),'filled');
%     colormap(map);
%     % zlim([0 100]);
%     colorbar;
%     xlim([0 280]);
%     ylim([0 512]);
%     zlim([0 13]);
%     hold on
% 
%     grid on
% 
%     for i=1:height(stats)
% 
%         b=stats.VoxelList{i,1};
%         %     for j=length(stats.ConvexHull{i,1})
%         %         plot3(a(:,2),a(:,1),a(:,3),'r');
%         k = boundary(b);
%         value=correlation_map_padding_show2(b(end,2),b(end,1),b(end,3));
% 
%         hold on
%         trisurf(k,b(:,2),b(:,1),b(:,3),'Facecolor',map(value,1:3),'FaceAlpha',0.1,'Edgecolor',map(value,1:3),'EdgeAlpha',0.1)
% 
% 
%         %text(b(end,2),b(end,1),b(end,3), num2str(value), 'Rotation',+15)
%     end
% 
% 
%     hold off
%     savefig(h1,strcat(addr2,t2,'_corroverlap.fig'));
%     close(h1);
%draw figures ends    

    Fullsize_2_mark_BW=Fullsize_2_mark;
    Fullsize_2_mark_BW(Fullsize_2_mark_BW>0)=1;
    Fullsize_2_mark_BW=logical(Fullsize_2_mark_BW);
    [Fullsize_2_mark_label,orgnum]=bwlabeln(Fullsize_2_mark);
    stats1 = regionprops3(Fullsize_2,'BoundingBox','VoxelList','ConvexHull','Centroid');
    
    sizelist=zeros;
    stats2=table;
    for i1=1:size(stats1,1)%sort the order in terms of size
        sizelist(i1,1)=size(stats1.VoxelList{i1,1},1);
    end
    [sizelistB,sizelistIndex] = sort(sizelist,'descend');
    for i1=1:size(stats1,1)
        stats2.Centroid(i1,1:3)=[stats1.Centroid(sizelistIndex(i1,1),1) stats1.Centroid(sizelistIndex(i1,1),2) stats1.Centroid(sizelistIndex(i1,1),3)];
        stats2.BoundingBox(i1,1:3)=[stats1.BoundingBox(sizelistIndex(i1,1),1) stats1.BoundingBox(sizelistIndex(i1,1),2) stats1.BoundingBox(sizelistIndex(i1,1),3)];
        stats2.VoxelList{i1,1}=stats1.VoxelList{sizelistIndex(i1,1),1};
        stats2.ConvexHull{i1,1}=stats1.ConvexHull{sizelistIndex(i1,1),1};
    end
        

%     for i1=2:2:length(detector_fusion)
%         for i2=1:length(detector_fusion(i1,:))
%             if detector_fusion(i1,i2)<1
%                 detector_fusion(i1,i2)=0;
%                 detector_fusion(i1-1,i2)=0;
%             end
%         end
%         if length(unique(detector_fusion(i1,:)))<3
%             detector_fusion(i1-1:i1,:)=0;
%         end
%     end
    
    
    stack_after_label(Fullsize_2_mark>0)=0;
    
    %[stack_after_label,orgnum]=bwlabeln(Fullsize2);
    %stats1 = regionprops3(stack_after_BW,'BoundingBox','VoxelList','ConvexHull');
    %stack_after(stack_after==0)=nan;
    
    %     stack_after_label(stack_after_label==0)=nan;
    %     h=figure;
    %     [X,Y,Z] = ndgrid(1:size(stack_after_label,1), 1:size(stack_after_label,2), 1:size(stack_after_label,3));
    %     pointsize = 5;
    %     scatter3(X(:), Y(:), Z(:), pointsize, stack_after_label(:),'filled');
    %     colormap([0 0 0]);
    %     % zlim([0 100]);
    %     %colorbar;
    %     hold on
    
    %stack_after_BW=logical(Fullsize_2_mark);
    %[stack_after_label,orgnum]=bwlabeln(Fullsize_2_mark);
    %stats = regionprops3(stack_after_BW,'BoundingBox','VoxelList','ConvexHull');
    
    %[stack_after_label,orgnum]=bwlabeln(Fullsize2);
    %stats1 = regionprops3(stack_after_BW,'BoundingBox','VoxelList','ConvexHull');
    %stack_after(stack_after==0)=nan;
    
    Fullsize_2_mark(Fullsize_2_mark==0)=nan;

    h2=figure;
%     [X,Y,Z] = ndgrid(1:size(Fullsize_2_mark,1), 1:size(Fullsize_2_mark,2), 1:size(Fullsize_2_mark,3));
% %     set(gca,'Color','k')
%     xlim([0 300]);
%     ylim([0 600]);
%     zlim([0 13]);
%     %scatter3(X(:), Y(:), Z(:), pointsize, Fullsize_2_mark(:),'filled');
%     colormap(map);
%     hold on
%     grid on
    
    
    newc=0;
    l=length(Registration1);
    Registration2=zeros(l,4);
    detector_old=zeros;
    detector_new=zeros;
    detector_split=zeros;
    detector3_fusion=zeros;
    detector_numbering=zeros;
    c1=1;
    c2=1;
    c3=1;
    c_numbering=0;
    cc=zeros;
    for i=1:size(stats2.VoxelList,1)
        %         for i1=1:length(Registration1)
        %             Dist(i1)=sqrt((stats.Centroid(i,1)-Registration1(i1,2))^2+(stats.Centroid(i,2)-Registration1(i1,3))^2+(stats.Centroid(i,3)-Registration1(i1,4))^2);
        %         end
        %         [D,index]=min(Dist);
%         average_object_intensity=0;
        max_object_intensity=0;
        b=stats2.VoxelList{i,1};
        if time+1<10
            add1='z00';
        elseif time+1<100
            add1='z0';
        end
        threeDimg=niftiread(strcat('F:\Mo\my3D_1\3Dimage\','img',add1,num2str(time+1),'.nii'));
        threeDimgPixellist=zeros;
        for i1=1:size(b,1)
            threeDimgPixellist(i1,1)=threeDimg(b(i1,2),b(i1,1),b(i1,3));
%             average_object_intensity=average_object_intensity+threeDimg(b(i1,2),b(i1,1),b(i1,3));
            if threeDimg(b(i1,2),b(i1,1),b(i1,3))>max_object_intensity
                max_object_intensity=threeDimg(b(i1,2),b(i1,1),b(i1,3));
            end
        end
        threeDimgPixellist = sort(threeDimgPixellist,'descend');
        % 80% top-valued pixel to average
        average_object_intensity=sum(threeDimgPixellist(threeDimgPixellist>threeDimgPixellist(round(size(threeDimgPixellist,1)*0.8),1)))/(round(size(threeDimgPixellist,1)*0.8)-1);
%         average_object_intensity=average_object_intensity/size(b,1);
        a=zeros;
        a_t_1=zeros;
        k = boundary(b);
        for i1=1:size(b,1)
            a(i1,1)=Fullsize_2_mark(b(i1,2),b(i1,1),b(i1,3));
        end
        [value,Value_f]=mode(a,'all');
        countnan=sum(isnan(a));
        if countnan>Value_f
            value=nan;
        end
        if time>1 % numbering overtracking
            for i1=1:size(b,1)
                a_t_1(i1,1)=Fullsize_1(b(i1,2),b(i1,1),b(i1,3));
            end
            [value_t_1,Value_f_t_1]=mode(a_t_1,'all');
            if ~isnan(value_t_1) && isempty(intersect(value_t_1,Registration1(:,1))) && ~isempty(intersect(value_t_1,detector_fusion_old.detector3_fusion)) %merge happed in last time point
                c_numbering=c_numbering+1;
                detector_numbering(c_numbering,1:2)=[value value_t_1];
                value=value_t_1;                          
                disp(value)
            end
        end
        if ~isempty(intersect(value,Registration2(:,1)))
%             disp(value)
            value2=setdiff(a,Registration2(:,1));
            if ~isempty(value2) && size(value2,2)>0 && ~isempty(intersect(value2,Registration1(:,1)))
                value=value2(1,ceil(rand*size(value2,2)));
            end
        end    
        if isnan(value)==1
            
            color=[0 0 0];
            newc=newc+1;
            Registration2(l+newc,1)=l+newc;
            Registration2(l+newc,2:4)=stats2.Centroid(i,:);
            value=l+newc;
            for i1=1:size(b,1)
                Fullsize_2_mark(b(i1,2),b(i1,1),b(i1,3))=value;
                Fullsize_2_2(b(i1,2),b(i1,1),b(i1,3))=value;
            end
            txt=strcat('NEW ',num2str(value));
            detector_new(c1,1)=value;
            c1=c1+1;
            xlswriter1(l+newc,(time)*2-1)=cellstr('new');
            xlswriter1(l+newc,(time)*2)=cellstr(num2str(l+newc));
            xlswriter3(l+newc,(time)*2)=cellstr(num2str(max_object_intensity));
            xlswriter4(l+newc,(time)*2)=cellstr(num2str(average_object_intensity));
            xlswriter5(l+newc,(time)*2)=cellstr(num2str(size(b,1)));
            xlswriter6(l+newc,(time)*2)=cellstr(num2str(stats2.Centroid(i,2)));
            xlswriter7(l+newc,(time)*2)=cellstr(num2str(stats2.Centroid(i,1)));
            xlswriter8(l+newc,(time)*2)=cellstr(num2str(stats2.Centroid(i,3)));
            draw_text(value)=text(b(end,1),b(end,2),b(end,3), txt, 'Rotation',+15);
        elseif isnan(value)==0 && value>0

            if isempty(intersect(value,Registration2(:,1)))
                
                %                 value=index;
                color=map(value,1:3);
%                 Registration2(value,:)=Registration1(value,:);
                Registration2(value,1)=value;
                Registration2(value,2:4)=stats2.Centroid(i,:);
                txt=strcat('OLD ',num2str(value));
                for i1=1:size(b,1)
                    Fullsize_2_2(b(i1,2),b(i1,1),b(i1,3))=value;
                end
                detector_old(c2,1)=value;
                c2=c2+1;
                xlswriter1(value,(time)*2-1)=cellstr(num2str(value));
                xlswriter1(value,(time)*2)=cellstr(num2str(value));
                xlswriter3(value,(time)*2)=cellstr(num2str(max_object_intensity));
                xlswriter4(value,(time)*2)=cellstr(num2str(average_object_intensity));
                xlswriter5(value,(time)*2)=cellstr(num2str(size(b,1)));
                xlswriter6(value,(time)*2)=cellstr(num2str(stats2.Centroid(i,2)));
                xlswriter7(value,(time)*2)=cellstr(num2str(stats2.Centroid(i,1)));
                xlswriter8(value,(time)*2)=cellstr(num2str(stats2.Centroid(i,3)));
                draw_forsure=0;                
                for i2=1:size(xlswriter1,2)
                    if iscellstr(xlswriter1(value,i2))
                        if string(xlswriter1(value,i2))~="new"
                            if str2num(string(xlswriter1(value,i2)))~=value
                                txt=strcat('OLD',string(xlswriter1(value,i2)),'(',num2str(value),')');
                                draw_text(value)=text(b(end,1),b(end,2),b(end,3), txt, 'Rotation',+15);
                                color=map(str2num(string(xlswriter1(value,i2))),1:3);
                                draw_forsure=1;
                                break
                            end
                            break
                        end
                    end
                end
                if ~draw_forsure
                    draw_text(value)=text(b(end,1),b(end,2),b(end,3), txt, 'Rotation',+15);
                end

                %draw_textt(value)=text(b(end,2),b(end,1),b(end,3), txt, 'Rotation',+15);
            else
                color=map(value,1:3);
                newc=newc+1;
                Registration2(l+newc,1)=l+newc;
                Registration2(l+newc,2:4)=stats2.Centroid(i,:);
                detector_split(c3,1)=value;
                detector_split(c3,2)=l+newc;
                xlswriter10(value,(time)*2)=cellstr(num2str(value));
                xlswriter10(l+newc,(time)*2)=cellstr(num2str(value));
                c3=c3+1; 
                for i1=1:size(b,1)
                    Fullsize_2_2(b(i1,2),b(i1,1),b(i1,3))=l+newc;
                    Fullsize_2_mark(b(i1,2),b(i1,1),b(i1,3))=l+newc;
                end
                xlswriter1(l+newc,1:(time-1)*2)=xlswriter1(value,1:(time-1)*2);
                xlswriter1(l+newc,(time)*2-1)=cellstr(num2str(value));
                xlswriter1(l+newc,(time)*2)=cellstr(num2str(l+newc));
                xlswriter3(l+newc,(time)*2)=cellstr(num2str(max_object_intensity));
                xlswriter4(l+newc,(time)*2)=cellstr(num2str(average_object_intensity));
                xlswriter5(l+newc,(time)*2)=cellstr(num2str(size(b,1)));
                xlswriter6(l+newc,(time)*2)=cellstr(num2str(stats2.Centroid(i,2)));
                xlswriter7(l+newc,(time)*2)=cellstr(num2str(stats2.Centroid(i,1)));
                xlswriter8(l+newc,(time)*2)=cellstr(num2str(stats2.Centroid(i,3)));
                for i2=1:size(xlswriter1,2)
                    if iscellstr(xlswriter1(value,i2))
                        if string(xlswriter1(value,i2))~="new"
                            value=str2num(string(xlswriter1(value,i2)));
                            break
                        end
                    end
                end
                txt=strcat('OLD',num2str(value),'=',num2str(l+newc));
                color=map(value,1:3);
                draw_text(l+newc)=text(b(end,1),b(end,2),b(end,3), txt, 'Rotation',+15);
                value=l+newc;
            end
            
        end
        trisurf(k,b(:,1),b(:,2),b(:,3),'Facecolor',color,'FaceAlpha',0.3,'Edgecolor',color,'EdgeAlpha',0.3);
        hold on
        if i==1
            draw_text(value)=text(b(end,1),b(end,2),b(end,3), txt, 'Rotation',+15);
        end
            
        
        
        
    end

    %pointsize = 5;
    %stack_after_label(stack_after_label==0)=nan;
    %scatter3(X(:), Y(:), Z(:), pointsize, stack_after_label(:),'filled');
    %scatter3(X(:), Y(:), Z(:), pointsize, Fullsize_2_mark(:),'filled');
    colormap(map);
    

    niftiwrite(Registration2,strcat(addr2,'Registration2_tracking','_',t2,'.nii'));
    niftiwrite(Fullsize_2_2,strcat(addr2,'Fullsize_2_aftertracking','_',t2,'.nii'));
    
        c=1;  %% fusion alarm part1
    for i1=1:size(stats2.VoxelList,1)
        b=stats2.VoxelList{i1,1};
        UNIQUECOUNT=zeros;
        for i2=1:size(b,1)
            UNIQUECOUNT(i2,1)=Fullsize_2_mark(b(i2,2),b(i2,1),b(i2,3));
            if isnan(UNIQUECOUNT(i2,1))
                UNIQUECOUNT(i2,1)=0;
            end
        end
        [C,ia,ic] = unique(UNIQUECOUNT);
        a_counts = accumarray(ic,1);
        value_counts = [C, a_counts];
        [x,y]=size(value_counts);
        
        if length(C)>1
%             if sum(C(:)==0)  %% value replacement for 0s within an object
%                 if length(find(C))==1
%                     replacement=C(find(C));
%                 elseif length(find(C))>1
%                     a=find(C);
%                     replacement=C(a(1));
%                 end
%                 for i2=1:size(b,1)
%                     if Fullsize_2_mark(b(i2,2),b(i2,1),b(i2,3))==0
%                         Fullsize_2_mark(b(i2,2),b(i2,1),b(i2,3))=replacement;
%                     end
%                 end
%             else
            detector_fusion(c:c+1,1:size(C,1))=value_counts';
            c=c+2;
%             end
        end
    end
%     for i1=1:2:size(detector_fusion,1) %% fusion size filter
%         for i2=1:1:size(detector_fusion,2)
%             if detector_fusion(i1+1,i2)<5
%                 detector_fusion(i1:i1+1,i2:end-1)=detector_fusion(i1:i1+1,i2+1:end);
%             end
%         end
%     end

    for i1=1:2:size(detector_fusion,1) %% fusion 0 filter
        if detector_fusion(i1,1)==0
            detector_fusion(i1:i1+1,1:end-1)=detector_fusion(i1:i1+1,2:end);
        end
    end
%     for i1=1:2:size(detector_fusion,1) %% fusion 0 filter
%         if detector_fusion(i1,1)==0
%             detector_fusion(i1:i1+1,1:end-1)=detector_fusion(i1:i1+1,2:end);
%         end
%     end
    detector2_fusion=detector_fusion;
%     for i1=1:length(Registration2(:,1)) %% fusion alarm part2
%         detector2_fusion(detector2_fusion==Registration2(i1,1))=0;
%     end
    for i1=1:2:size(detector2_fusion,1) %% fusion alarm part2
        for i2=1:size(detector2_fusion,2)        
            if intersect(detector2_fusion(i1,i2),Registration2(:,1))
                detector2_fusion(i1:i1+1,i2)=0;
            end
        end
        for i2=1:1:size(detector2_fusion,2)% fusion size filter
            if detector2_fusion(i1+1,i2)<5
                detector2_fusion(i1:i1+1,i2)=0;
            end
        end
    end

    c=1;
    for i1=1:2:length(detector2_fusion(:,1))
        if ~isempty(nonzeros(detector2_fusion(i1,:)))
            detector3_fusion(c:c+1,1:size(detector_fusion,2))=detector_fusion(i1:i1+1,:);
            c=c+2;
        end
    end
%     if exist('detector3_fusion','var')
%         detector3_fusion(:,1)=[];    %% fusion alarm part1
%     else
%         detector3_fusion=0;
%     end

    for i2=1:2:size(detector3_fusion,1)
        detector3_fusion_exist=0;
        for i1=1:size(detector3_fusion,2)
            if detector3_fusion(i2,i1)>0
                if  isa(draw_text(detector3_fusion(i2,i1)),'matlab.graphics.primitive.Text')
                    draw_text(detector3_fusion(i2,i1)).Color = 'red';
                    if detector3_fusion_exist~=0
                        if detector3_fusion(i2+1,i1)>c
                            detector3_fusion_exist=detector3_fusion(i2,i1);
                            c=detector3_fusion(i2+1,i1);
                        end
                    else
                        detector3_fusion_exist=detector3_fusion(i2,i1);
                        c=detector3_fusion(i2+1,i1);
                    end
                end
            end
        end
        if detector3_fusion_exist==0 %numbering for fusion debug (situation: continue merge that old number is recovered)
            for i1=1:size(detector3_fusion,2)
                if detector3_fusion(i2,i1)>0
                    if ~isempty(intersect(detector3_fusion(i2,i1),detector_numbering))
                        [ii,jj]=find(detector_numbering==detector3_fusion(i2,i1));
                        if size(jj,1)==1
                            value_numbering_recover=detector_numbering(ii,jj+1);
                            
                            txt=draw_text(value_numbering_recover).String;
                            txt=strcat(txt,'M','(',num2str(detector3_fusion(i2,i1)),')');
                            draw_text(value_numbering_recover).String=txt;
                            draw_text(value_numbering_recover).Color='red';
                            detector3_fusion_exist=value_numbering_recover;
                            disp('----')
                            disp(value_numbering_recover)
                        end
                    end
                end
            end
        end
        for i1=1:size(detector3_fusion,2)
            if detector3_fusion(i2,i1)>0 && exist('detector3_fusion_exist','var') 
                if  ~isa(draw_text(detector3_fusion(i2,i1)),'matlab.graphics.primitive.Text')
                    xlswriter9(detector3_fusion(i2,i1),time*2-1) = xlswriter1(detector3_fusion_exist,time*2-1);
                end
            end
        end
        
    end
    view([0 0 1]);
    set(gca, 'YDir','reverse')
    xlim([0 512]);
    ylim([0 280]);
    zlim([0 13]);
    hold off
    savefig(h2,strcat(addr2,t2,'_tracking.fig'));
    save(strcat(addr2,'fusion_tracking','_',t2,'.mat'),'detector3_fusion');
    save(strcat(addr2,'split_tracking','_',t2,'.mat'),'detector_split');
    save(strcat(addr2,'detector_numbering','_',t2,'.mat'),'detector_numbering');
    save(strcat(addr2,'draw_text','_',t2,'.mat'),'draw_text');
    close(h2);  
    xlswriter2(time+1,1)=cellstr(num2str(time));
    xlswriter2(time+1,2)=cellstr(num2str(size(detector_old,1)));
    xlswriter2(time+1,3)=cellstr(num2str(size(detector_new,1)));
    xlswriter2(time+1,4)=cellstr(num2str(size(detector_split,1)));
    xlswriter2(time+1,5)=cellstr(num2str(size(detector3_fusion,1)/2));
    
    timecount(time)=toc;
    disp(timecount(time))
    clear detector_old detector_new detector_split detector_fusion detector2_fusion detector3_fusion draw_text detector_numbering
end
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter1','.mat'),'xlswriter1');
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter2','.mat'),'xlswriter2');
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter3','.mat'),'xlswriter3');
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter4','.mat'),'xlswriter4');
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter5','.mat'),'xlswriter5');
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter6','.mat'),'xlswriter6');
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter7','.mat'),'xlswriter7');
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter8','.mat'),'xlswriter8');
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter9','.mat'),'xlswriter9');
save(strcat('F:\Mo\my3D_1\Tracking\','xlswriter10','.mat'),'xlswriter10');
xlswrite(filename,xlswriter1,1, 'A1');
xlswrite(filename,xlswriter2,2, 'A2');
xlswrite(filename,xlswriter3,3, 'A1');
xlswrite(filename,xlswriter4,4, 'A1');
xlswrite(filename,xlswriter5,5, 'A1');
xlswrite(filename,xlswriter6,6, 'A1');
xlswrite(filename,xlswriter7,7, 'A1');
xlswrite(filename,xlswriter8,8, 'A1');
xlswrite(filename,xlswriter9,9, 'A1');
xlswrite(filename,xlswriter10,10, 'A1');
disp('finished')
    %% plot the ground truth
    addr2=strcat('D:\Mo\gt\slices-refined\bww\');
    Files2=dir(strcat(addr2,'*.tif'));
    I_stack=zeros(280,512,13);
    for i1=1:length(Files2)
        I=imread(strcat(addr2,Files2(i1).name));
        I_stack(:,:,i1)=I;
    end
    stack_after_BW=logical(I_stack);
    stats = regionprops3(stack_after_BW,'BoundingBox','VoxelList','ConvexHull','Centroid');
    [stack_after_label,orgnum]=bwlabeln(I_stack);
    stack_after_label(stack_after_label==0)=nan;
    h=figure;
    [X,Y,Z] = ndgrid(1:size(stack_after_label,1), 1:size(stack_after_label,2), 1:size(stack_after_label,3));
    pointsize = 5;
    scatter3(X(:), Y(:), Z(:), pointsize, stack_after_label(:),'filled');
    colormap(map);
    % zlim([0 100]);
    colorbar;
    hold on
    grid on
    
    
    
    for i=1:height(stats)
        Registration(i,:)=stats.Centroid(i,:);
        b=stats.VoxelList{i,1};
        %     for j=length(stats.ConvexHull{i,1})
        %         plot3(a(:,2),a(:,1),a(:,3),'r');
        k = boundary(b);
        hold on
        trisurf(k,b(:,2),b(:,1),b(:,3),'Facecolor',map(i,1:3),'FaceAlpha',0.1,'Edgecolor',map(i,1:3),'EdgeAlpha',0.1)
        %         hold on
        %     end
        value=stack_after_label(b(end,2),b(end,1),b(end,3));
        text(b(end,2),b(end,1),b(end,3), num2str(i), 'Rotation',+15)
    end
    
    
    hold off
%%
for i=1:10000
    i1=rand();
    i2=rand();
    i3=rand();
    map(i,1:3)=[i1 i2 i3];
end
save('F:\Mo\my3D_1\Tracking\colormap.mat','map')
