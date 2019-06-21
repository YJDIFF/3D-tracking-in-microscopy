# 3D-tracking-in-microscopy
2019 CVPR. Multi-Object Portion Tracking in 4D Fluorescence Microscopy Imagery with Deep Feature Maps. Paper link: http://openaccess.thecvf.com/content_CVPRW_2019/papers/CVMI/Jiao_Multi-Object_Portion_Tracking_in_4D_Fluorescence_Microscopy_Imagery_With_Deep_CVPRW_2019_paper.pdf

To use this project,
1. you need segmentation result and deep feature maps from the DNN. Because in this project, the 3D image is cropped into 32x35x13, I firstly run "Assemble.m" to assemble the prediction result into the original size 280x512x13. "Assemble.m" will also draw 3D labelling figures.

2. After assembling all information, run "tracking.m". The code will tracking each independent 3D object through 2 time points. It will correlate deep feature maps (By correlation.m).
