function data = getData()
    load('backhoe_el000_az350to100.mat');
    data = struct;
    data.deltaF = FGHz(2) - FGHz(1); 
    data.minF = repmat(FGHz(1),[size(HH,2),1]);
    
    
    
    % useful convenience computations
    data.deltaAz = deg2rad( AZ(2) - AZ(1) ); % radians
    data.FK = FGHz(end); 
    data.K = length(FGHz);
    data.Fc = FGHz(data.K/2)*1e9;
    data.lambdac = 3e8 / data.Fc;
    data.Np = length(AZ);
    
    
    %arbitrary range to scene center (everything is mocomp to scene center so
    %is relative to the angles anyway)
    rc = 10000;
    data.AntX = cos(deg2rad(AZ))*rc;
    data.AntY = sin(deg2rad(AZ))*rc;
    data.AntZ = zeros(size(AZ));
    data.RO = sqrt(data.AntX.^2 + data.AntY.^2);
    data.phdata = VV;
    data.Nfft = 10*data.K;
    data.AZ = deg2rad(AZ);
end