function img = matchedfilter_data(data)
    c = 2.99792458e8;
    data.K = size(data.phdata,1);
    data.Np = size(data.phdata,2);
    data.AntAz = (unwrap(atan2(data.AntY,data.AntX)));
    data.deltaAz = abs(mean(diff(data.AntAz)));

    % total azimuth an((c / (data.FK*1e9))gle of the aperture
    data.totalAz = max(data.AntAz) - min(data.AntAz);
    data.maxWr = c/(2*data.deltaF*1e9);
    data.maxWx = c/(2*data.deltaAz*mean(data.FK*1e9));

    % determine resolution of the image
    data.dr = c/(2*data.deltaF*1e9*data.K);
    data.dx = c/(2*data.totalAz*(data.Fc));

    % display maximum scene size and resolution
    fprintf('Maximum Scene Size: %.2f m range, %.2f m cross-range\n', data.maxWr, data.maxWx);
    fprintf('Resolution: %.3f m range, %.3f m cross-range\n', data.dr, data.dx);
    
    data.im_final = zeros(size(data.x_mat));
    
    % set up a vector to keep execution times for every pulse
    t = zeros(1,data.Np);
    
    % loop through every pulse
    for ii = 1:data.Np
        
        % Display status of the imaging process
        if (ii>1)
            t_sofar = sum(t(1:ii-1));
            t_est = (t_sofar*data.Np/(ii-1) - t_sofar)/60;
            fprintf('Pulse %d of %d, %.02f minutes remaining\n',ii,data.Np,t_est);
        else 
            fprintf('Pulse %d of %d\n',ii,data.Np);
        end 
        tic;
        % compute the differential range for every pixel in the image (m)
        dR = sqrt((data.AntX(ii) - data.x_mat).^2 + ...
             (data.AntY(ii) - data.y_mat).^2 + ...
             (data.AntZ(ii) - data.z_mat).^2) - data.RO(ii);
         
        % calculate the frequency of each sample in the pulse (Hz)
        freq = (data.minF(ii) + (0:(data.K-1))*data.deltaF)*1e9;
        
        % perform matched filter operation
        for jj = 1:data.K 
            data.im_final = data.im_final + data.phdata(jj,ii)*exp(1i*4*pi*freq(jj)/c*dR);
        end
        % determine the execution for this pulse
        t(ii) = toc;
    end 
    img = data.im_final;
       
        
end
