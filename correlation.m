function correlation(Fullsize_1,Fullsize_2,Fullsize_regression_1,Fullsize_regression_2,t2,time,spatial_extend_matrix)
depth=size(Fullsize_regression_1,4);
if ~exist('time','var')
     % third parameter does not exist, so default it to something
      time=1;
end
addr2=strcat('F:\Mo\my3D_1\Tracking\',t2,'\');
% for ix=1:280
%     for iy=1:512
%         for iz=1:13
%             if Fullsize_1(ix,iy,iz)==0
%                 Fullsize_regression_1(ix,iy,iz,:)=0;
%             end
%             if Fullsize_2(ix,iy,iz)==0
%                 Fullsize_regression_2(ix,iy,iz,:)=0;
%             end
%         end
%     end
% end

padding=20;
[x,y,z]=size(Fullsize_1);
Fullsize_1_padding=zeros(x+padding*2, y+padding*2, z+4);
Fullsize_2_padding=zeros(x+padding*2, y+padding*2, z+4);
Fullsize_regression_1_padding=zeros(x+padding*2, y+padding*2, z+4,depth);
Fullsize_regression_2_padding=zeros(x+padding*2, y+padding*2, z+4,depth);
Fullsize_1_padding(padding+1:x+padding,padding+1:y+padding,3:z+2)=Fullsize_1;
Fullsize_2_padding(padding+1:x+padding,padding+1:y+padding,3:z+2)=Fullsize_2;
Fullsize_regression_1_padding(padding+1:x+padding,padding+1:y+padding,3:z+2,:)=Fullsize_regression_1(:,:,:,1:depth);
Fullsize_regression_2_padding(padding+1:x+padding,padding+1:y+padding,3:z+2,:)=Fullsize_regression_2(:,:,:,1:depth);

%correlation_map_padding=zeros(x+padding*2, y+padding*2, z+4,max(max(max(Fullsize_1))));
correlation_map_padding_corr=zeros(x+padding*2, y+padding*2, z+4);
correlation_map_padding_show=zeros(x+padding*2, y+padding*2, z+4);
clear Fullsize_regression_1 Fullsize_regression_2 Fullsize_1 Fullsize_2

Fullsize_1_label=Fullsize_1_padding;
% Fullsize_2_label=Fullsize_2_padding;
disp(max(max(max(Fullsize_1_padding))))
stats1 = regionprops3(logical(Fullsize_1_padding),'BoundingBox','VoxelList','ConvexHull');
% stats2 = regionprops3(logical(Fullsize_2_padding),'BoundingBox','VoxelList','ConvexHull');
%filename = strcat(addr2,t2,'_tracking.xls');
disp('next')

for i=1:height(stats1)
    if rem(i,50)==0
        disp(i/height(stats1))
    end
    if size(stats1.VoxelList{i,1},1)<30
        stepsize=1;
    else
        stepsize=3;
    end
    for n1=1:stepsize:size(stats1.VoxelList{i,1},1)
        if stepsize==1
            index=n1;
        else
            index=ceil(rand()*size(stats1.VoxelList{i,1},1));
        end
        if stats1.VoxelList{i,1}(index,2)<=300 && stats1.VoxelList{i,1}(index,1)<=532 && stats1.VoxelList{i,1}(index,3)<=15
            Feature_map1=Fullsize_regression_1_padding(stats1.VoxelList{i,1}(index,2)-4:stats1.VoxelList{i,1}(index,2)+5,stats1.VoxelList{i,1}(index,1)-4:stats1.VoxelList{i,1}(index,1)+5,stats1.VoxelList{i,1}(index,3)-1:stats1.VoxelList{i,1}(index,3)+1,:);
            
            Feature_map1=Feature_map1.*spatial_extend_matrix;

            for m1=-1:1
                x=2*m1;
                for m2=-1:1
                    y=2*m2;
                    for m3=-1:1
                        z=m3;
                        Feature_map2=Fullsize_regression_2_padding(stats1.VoxelList{i,1}(index,2)+x-4:stats1.VoxelList{i,1}(index,2)+x+5,stats1.VoxelList{i,1}(index,1)+y-4:stats1.VoxelList{i,1}(index,1)+y+5,stats1.VoxelList{i,1}(index,3)+z-1:stats1.VoxelList{i,1}(index,3)+z+1,:);
                        Feature_map2=Feature_map2.*spatial_extend_matrix;
                        %Feature_map1=Feature_map1/mean2(Feature_map1);
                        %Feature_map2=Feature_map2/mean2(Feature_map2);
                        %corr=convn(Feature_map1,Feature_map2(end:-1:1,end:-1:1,end:-1:1));
                        
                        Feature_map1_flatten = Feature_map1(:);
%                                                 Feature_map1_flatten = Feature_map1_flatten+1;
                        Feature_map2_flatten = Feature_map2(:);
%                                                 Feature_map2_flatten = Feature_map2_flatten+1;
                        %corr = dot(Feature_map1_flatten,Feature_map2_flatten)/(norm(Feature_map1_flatten)*norm(Feature_map2_flatten));
                        corr = corr2(Feature_map1_flatten,Feature_map2_flatten);
                        
                        %                             ele=nnz(corr);
                        %                             if ele==0
                        %                                 ele=1;
                        %                             end
                        %                             corr=sum(sum(sum(corr)))/ele^2;
                        

%                         b=stats1.VoxelList{i,1}; %debug
%                         a=zeros;
%                         for i1=1:size(b,1)
%                             a(i1,1)=Fullsize_1_label(b(i1,2),b(i1,1),b(i1,3));
%                         end
%                         value=mode(a,'all');
%                         if value==3538
%                             disp(value)
%                             countzero=size(a(a==0),1);
%                             if countzero>value
%                                 value=0;
%                             end
%                             disp(value)
%                             disp(corr)
%                             disp('')
%                         end
                        if corr>0.2
                            %xlswrite(filename,cellstr(stats1.VoxelList{i,1}(index,2)+x-4),i, strcat('A',num2str(count)));
                            %xlswrite(filename,cellstr(VoxelList{i,1}(index,1)+y-4),i, strcat('B',num2str(count)));
                            %xlswrite(filename,cellstr(VoxelList{i,1}(index,3)+y-4),i, strcat('C',num2str(count)));
                            %count=count+1;
                            b=stats1.VoxelList{i,1};
                            a=zeros;
                            for i1=1:size(b,1)
                                a(i1,1)=Fullsize_1_label(b(i1,2),b(i1,1),b(i1,3));
                            end
                            value=mode(a,'all');
                            countzero=size(a(a==0),1);
                            if countzero>value
                                value=0;
                            end
                            %                             value=Fullsize_1_label(b(end,2),b(end,1),b(end,3));
                            correlation_map_padding_corr_local=correlation_map_padding_corr(stats1.VoxelList{i,1}(index,2)+x-4:stats1.VoxelList{i,1}(index,2)+x+5,stats1.VoxelList{i,1}(index,1)+y-4:stats1.VoxelList{i,1}(index,1)+y+5,stats1.VoxelList{i,1}(index,3)+z-1:stats1.VoxelList{i,1}(index,3)+z+1);
                            correlation_map_padding_show_local=correlation_map_padding_show(stats1.VoxelList{i,1}(index,2)+x-4:stats1.VoxelList{i,1}(index,2)+x+5,stats1.VoxelList{i,1}(index,1)+y-4:stats1.VoxelList{i,1}(index,1)+y+5,stats1.VoxelList{i,1}(index,3)+z-1:stats1.VoxelList{i,1}(index,3)+z+1);
                            %correlation_map_padding(stats1.VoxelList{i,1}(index,2)+x-4:stats1.VoxelList{i,1}(index,2)+x+5,stats1.VoxelList{i,1}(index,1)+y-4:stats1.VoxelList{i,1}(index,1)+y+5,stats1.VoxelList{i,1}(index,3)+z-1:stats1.VoxelList{i,1}(index,3)+z+1,value)=corr;
                            %correlation_map_padding_show(stats1.VoxelList{i,1}(index,2)+x-4:stats1.VoxelList{i,1}(index,2)+x+5,stats1.VoxelList{i,1}(index,1)+y-4:stats1.VoxelList{i,1}(index,1)+y+5,stats1.VoxelList{i,1}(index,3)+z-1:stats1.VoxelList{i,1}(index,3)+z+1)=value;
                            %                                 for i1=stats1.VoxelList{i,1}(index,2)+x-4:stats1.VoxelList{i,1}(index,2)+x+5
                            %                                     for i2=stats1.VoxelList{i,1}(index,1)+y-4:stats1.VoxelList{i,1}(index,1)+y+5
                            %                                         for i3=stats1.VoxelList{i,1}(index,3)+z-1:stats1.VoxelList{i,1}(index,3)+z+1
                            %                                             if correlation_map_padding_corr(i1,i2,i3)<corr
                            %
                            %                                                 correlation_map_padding_corr(i1,i2,i3)=corr;
                            %                                                 correlation_map_padding_show(i1,i2,i3)=value;
                            %                                             end
                            %                                         end
                            %                                     end
                            %                                 end
                            correlation_map_padding_show_local(correlation_map_padding_corr_local<corr)=value;
                            correlation_map_padding_corr_local(correlation_map_padding_corr_local<corr)=corr;
                            correlation_map_padding_corr(stats1.VoxelList{i,1}(index,2)+x-4:stats1.VoxelList{i,1}(index,2)+x+5,stats1.VoxelList{i,1}(index,1)+y-4:stats1.VoxelList{i,1}(index,1)+y+5,stats1.VoxelList{i,1}(index,3)+z-1:stats1.VoxelList{i,1}(index,3)+z+1)=correlation_map_padding_corr_local;
                            correlation_map_padding_show(stats1.VoxelList{i,1}(index,2)+x-4:stats1.VoxelList{i,1}(index,2)+x+5,stats1.VoxelList{i,1}(index,1)+y-4:stats1.VoxelList{i,1}(index,1)+y+5,stats1.VoxelList{i,1}(index,3)+z-1:stats1.VoxelList{i,1}(index,3)+z+1)=correlation_map_padding_show_local;
                            
                            
                        end
                    end
                end
            end
        end
    end
end
%niftiwrite(correlation_map_padding,strcat(addr2,'correlation_map_padding','_',t2,'.nii'));
disp(max(max(max(correlation_map_padding_show))))
niftiwrite(correlation_map_padding_show,strcat(addr2,'correlation_map_padding_show_traceback',num2str(time),'_',t2,'.nii'));
niftiwrite(correlation_map_padding_corr,strcat(addr2,'correlation_map_padding_hide_traceback',num2str(time),'_',t2,'.nii'));
end

