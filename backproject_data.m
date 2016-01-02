function img = backproject_data(data)
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
    
    % calculate the range to every bin in range profile (m)
    data.r_vec = linspace(-data.Nfft/2,data.Nfft/2-1,data.Nfft)*data.maxWr/data.Nfft;
    
    % initialize image with all zero values
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
        
        % Form range profile with zero padding added
        rc = fftshift(ifft(data.phdata(:,ii),data.Nfft));
        % compute differential range for each pixel in the image (m)
        dR = sqrt((data.AntX(ii) - data.x_mat).^2 + ...
            (data.AntY(ii) - data.y_mat).^2 + ...
            (data.AntZ(ii) - data.z_mat).^2) - data.RO(ii);
        
        % Calculate the phase correction for image
        phCorr = exp(1i*4*pi*data.minF(ii)/c*dR);
        
        % Determine which pixels fall within the range swath
        I = find(and(dR > min(data.r_vec), dR < max(data.r_vec)));
        
        % Update the image using linear interpolation
        data.im_final(I) = data.im_final(I) + interp1(data.r_vec,rc,dR(I),'linear') .* phCorr(I);
        
        % Determine the execution time for this pulse
        t(ii) = toc;
    end
    
    img = data.im_final;
    
   
    
end