% compute the image limitation based on data collected
data = getData();
c=3e8;
Wr   = c / (2*(data.deltaF*1e9))
Waz  = ((c / (data.FK*1e9)) / (2*data.deltaAz))
delta_r = c / (2*(data.K-1)*data.deltaF*1e9)
delta_Az = data.lambdac / (2*(data.Np-1)*data.deltaAz)

% form image based on image limitations
[data.x_mat,data.y_mat,data.z_mat] = meshgrid(-5:0.02:5,-5:0.02:5,0);


img = backproject_data(data);
%img_mf = matchfilter_data(data);

figure; imagesc(abs(img));
%figure; imagesc(abs(img_mf)); 


