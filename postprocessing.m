%%
clear
c1=1;
load('G:\Mo\my3D_matlab\colormap.mat','map')
addr='G:\Mo\my3D_matlab\Test\';
% "G:\Mo\my3D_matlab\Test\try\FC-DenseNet\",...
addr2=[
    "G:\Mo\my3D_matlab\Test\c01\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\c02\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\c03\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\c04\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\c05\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\w01\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\w02\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\w03\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\w04\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\w05\FC-DenseNet\",...
    "G:\Mo\my3D_matlab\Test\w06\FC-DenseNet\"];
sample=["cg2a_ ubiEcadGFP-01 okay strong z2-12 t20.nii",...
    "cg2a_ ubiEcadGFP-03 good long time from 191202 z2-12+t39.nii",...
    "cg2a_ ubiEcadGFP-04 okay weak z2-12+t27.nii",...
    "cg2a_ ubiEcadGFP-05 good from 191202 z2-12+t23.nii",...
    "cg14427_ ubiEcadGFP-01 okay from 191203 z2-12+t9.nii",...
    "sqhAx3_ ubiEcadGFP-01 good from 191203 z2-12+t20.nii",...
    "sqhAx3_ ubiEcadGFP-01 okay z3-11 t11.nii",...
    "sqhAx3_ ubiEcadGFP-02 okay from 191203 z3-12+t13.nii",...
    "sqhAx3_ ubiEcadGFP-03 okay z3-11 t11.nii",...
    "sqhAx3_ ubiEcadGFP-05 okay z2-3-12 t24.nii",...
    "sqhAx3_ubiEcad GFP CyO-02okay z3-16 t13.nii"];
timeclip_start=[4,25,17,13,1,4,1,1,1,9,1];
timeclip_end=[25,44,36,34,18,25,18,18,16,30,20];
z1=[2,2,2,2,2,2,2,3,2,2,2];
z2=[10,10,10,10,10,10,11,10,10,10,10];
Dsize=[512,256,13];

for ii=1:length(sample)
    disp(ii)
    c=1;
    V_gt=niftiread(strcat(addr,sample(ii)));
    Files1=dir(strcat(addr2(ii),'*.nii'));
    for i1=1:timeclip_start(ii)-timeclip_start(ii)+1:timeclip_end(ii)-timeclip_start(ii)+1
        disp(i1)
        V=zeros(Dsize);
        V_weight=zeros([Dsize,64]);
        V_gt_time=squeeze(V_gt(:,:,z1:z2,i1-1+timeclip_start(ii)));
        V_gt_time=imresize3(V_gt_time,Dsize);
        if i1<10
            tt=strcat('00',num2str(i1));
        elseif i1<100
            tt=strcat('0',num2str(i1));
        elseif i1<1000
            tt=strcat('',num2str(i1));
        end
        if ~exist(strcat(addr2(ii),'\',tt,'\'), 'dir')
            mkdir(strcat(addr2(ii),'\',tt,'\'))
        end
        if i1<10
            tt=strcat('00',num2str(i1));
        elseif i1<100
            tt=strcat('0',num2str(i1));
        elseif i1<1000
            tt=strcat('',num2str(i1));
        end
        if ~exist(addr2(ii), 'dir')
            mkdir(addr2(ii))
        end
        for m=1:16:Dsize(1)-15
            for n=1:16:Dsize(2)-15        
                Mat = niftiread(strcat(addr2(ii),Files1(c).name));
                Mat_weight = niftiread(strcat(addr2(ii),Files1(c+1).name));
                Mat_weight = Mat_weight(:,:,:,1:64); %all is 256
                Mat = imresize3(uint8(Mat), [16,16,13]);
                V(m:m+15, n:n+15,:)=1-Mat; 
                V_weight(m:m+15, n:n+15,:,:)=Mat_weight;
                c=c+2;
            end
        end
        Fullsize2=logical(V);
        for it=1:13
            img=Fullsize2(:,:,it);
            [f,orgnum] = bwlabel(img);
            g = regionprops(f, 'Area');
            area_values=[g.Area];
            % idx = find ((25<= area_values & 800>= area_values));
            idx = find ((3<= area_values));
            h = ismember (f, idx);
            Fullsize2(:,:,it)=h;
        end
        Fullsize2=double(Fullsize2);
               
        stack_after=Fullsize2;
        
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
        niftiwrite(stack_after_label,strcat(addr2(ii),'\',tt,'\','Fullsize_label','_',tt,'.nii'));
        for i=1:height(stats1)
            b=stats1.VoxelList{i,1};
            [x,y]=size(b);
            for i1=1:x
                stack_after_label(stack_after_label(b(i1,2),b(i1,1),b(i1,3))>0)=i;
            end
        end
        
        stack_after_label(stack_after_label==0)=nan;
        h=figure;
        [X,Y,Z] = ndgrid(1:size(stack_after_label,1), 1:size(stack_after_label,2), 1:size(stack_after_label,3));
        pointsize = 5;
        hold on
        grid on
                     
        
        for i=1:height(stats1)
            
            b=stats1.VoxelList{i,1};
            
            k = boundary(b);
            hold on
            value=i;
            Registration(value,1)=value;
            Registration(value,2:4)=stats1.Centroid(i,:);
            trisurf(k,b(:,1),b(:,2),b(:,3),'Facecolor',map(value,1:3),'FaceAlpha',0.3,'Edgecolor',map(value,1:3),'EdgeAlpha',0.3)
            
%             text(b(end,1),b(end,2),b(end,3), num2str(value), 'Rotation',+15, 'Color', c2)
        end
        
        
        hold off
        view([0 0 1]);
        set(gca, 'YDir','reverse')
%         set(gca, 'XDir','reverse')
        set(gca, 'ZDir','reverse')
        xlim([0 256]);
        ylim([0 512]);
        zlim([0 13]);

        niftiwrite(Fullsize2,strcat(addr2(ii),'\',tt,'\','Fullsize','_',tt,'.nii'));
        niftiwrite(V_gt_time,strcat(addr2(ii),'\',tt,'\','GT','_',tt,'.nii'));
        niftiwrite(Registration,strcat(addr2(ii),'\',tt,'\','Registration','_',tt,'.nii'));
        niftiwrite(V_weight,strcat(addr2(ii),'\',tt,'\','Weights','_',tt,'.nii'));
        savefig(h,strcat(addr2(ii),'\',tt,'\',tt,'_3Dconnection2.fig'));
        saveas(h,strcat(addr2(ii),'\',tt,'\',tt,'_3Dconnection2.png'))
        close(h);

    end
    
end
display('finish');