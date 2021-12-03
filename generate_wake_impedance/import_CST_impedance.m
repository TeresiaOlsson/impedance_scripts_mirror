function struct = import_CST_impedance(file,sampling_points,betas)
%% Import impedance from CST output and resample to required sampling points

    %% Read in data
    
    % Read file
    fileID = fopen(file);
    data = textscan(fileID,'%f %f %s %s','CommentStyle','#');
    fclose(fileID);        
        
    % Split file for planes
    index = find(data{1} == 0);
    data_lon = cellfun(@(v)v (index),data,'UniformOutput',false);

    index = find(data{1} == 1);
    data_hor = cellfun(@(v)v (index),data,'UniformOutput',false);

    index = find(data{1} == 2);       
    data_ver = cellfun(@(v)v (index),data,'UniformOutput',false);    
    
    %% Create output wakes

    ImpedanceFreq = zeros(length(sampling_points),1);
    ImpedanceRealZ = zeros(length(sampling_points),1);
    ImpedanceImagZ = zeros(length(sampling_points),1);
    ImpedanceRealX = zeros(length(sampling_points),1);
    ImpedanceImagX = zeros(length(sampling_points),1);  
    ImpedanceRealY = zeros(length(sampling_points),1);
    ImpedanceImagY = zeros(length(sampling_points),1);     
    
    %% Longitudinal impedance
    
    fprintf('Reading CST impedance files\n');
    
    fprintf('Reading longitudinal files\n');
    for i = 1:size(data_lon{1},1)
        
        filename = char(data_lon{4}(i));
        
        if ~isempty(filename)
            fileID = fopen(filename);
            data = textscan(fileID,'%f %f %f','CommentStyle','#');
            fclose(fileID);          
            lon_impedance = cell2mat(data);  % Units GHz, Ohm

            % Change units to Hz
            freq = lon_impedance(:,1).*1e9;
            real_impedance = lon_impedance(:,2);
            imag_impedance = lon_impedance(:,3);

            % Resample impedance
            impedance_real = interp1(freq,real_impedance,sampling_points,'linear',0);
            impedance_imag = interp1(freq,imag_impedance,sampling_points,'linear',0);

            % Put imaginary DC component to zero since required by Elegant
            index = find(freq == 0);
    %        impedance_real(index) = 0; 
            impedance_imag(index) = 0;        

            ImpedanceRealZ = ImpedanceRealZ + impedance_real;
            ImpedanceImagZ = ImpedanceImagZ + impedance_imag;
        end

    end
    fprintf('Finished reading longitudinal files. Read %d lines.\n\n',i);    
    
    %% Horizontal impedance
           
    if ~isempty(betas)                
        % Calculate average horizontal beta over element lengths       
        element_length = data_hor{2};
        element_s = [0; cumsum(element_length)]';         
        average_betax = zeros(1,length(element_s)-1);         
        for i = 1:length(element_s)-1                   
            average_betax(i) = integrate(betas.betax,element_s(i+1),element_s(i))./element_length(i);            
        end
    else
        element_length = data_hor{2};
        element_s = [0; cumsum(element_length)]';
        average_betax = ones(1,length(element_s)-1);
    end 
     
    fprintf('Reading horizontal files\n');     
    for i = 1:size(data_hor{1},1)
        
        filename = char(data_hor{4}(i));
        
        if ~isempty(filename)
            fileID = fopen(filename);
            data = textscan(fileID,'%f %f %f','CommentStyle','#');
            fclose(fileID);          
            hor_impedance = cell2mat(data);  % Units GHz, Ohm/mm

            % Change units to Hz and Ohm?m
            freq = hor_impedance(:,1).*1e9;
            real_impedance = hor_impedance(:,2).*1e3;
            imag_impedance = hor_impedance(:,3).*1e3;

            % Resample impedance
            impedance_real_interp = interp1(freq,real_impedance,sampling_points,'linear',0);
            impedance_imag_interp = interp1(freq,imag_impedance,sampling_points,'linear',0);
            
            % Change parts and signs to match Elegant conventions
            impedance_real = -impedance_imag_interp;
            impedance_imag  = impedance_real_interp;
                    
            ImpedanceRealX = ImpedanceRealX + impedance_real.*average_betax(i);
            ImpedanceImagX = ImpedanceImagX + impedance_imag.*average_betax(i);
        end

    end
    fprintf('Finished reading horizontal files. Read %d lines.\n\n',i');    
             
    %% Vertical impedance
    
    if ~isempty(betas)                
        % Calculate average vertical beta over element lengths       
        element_length = data_ver{2};
        element_s = [0; cumsum(element_length)]';         
        average_betay = zeros(1,length(element_s)-1);         
        for i = 1:length(element_s)-1                   
            average_betay(i) = integrate(betas.betay,element_s(i+1),element_s(i))./element_length(i);            
        end
    else
        element_length = data_ver{2};
        element_s = [0; cumsum(element_length)]'; 
        average_betay = ones(1,length(element_s)-1); 
    end
    
    fprintf('Reading vertical files\n');     
    for i = 1:size(data_ver{1},1)
        
        filename = char(data_ver{4}(i));
        
        if ~isempty(filename)
            fileID = fopen(filename);
            data = textscan(fileID,'%f %f %f','CommentStyle','#');
            fclose(fileID);          
            ver_impedance = cell2mat(data);  % Units GHz, Ohm/mm

            % Change units to Hz and Ohm?m
            freq = ver_impedance(:,1).*1e9;
            real_impedance = ver_impedance(:,2).*1e3;
            imag_impedance = ver_impedance(:,3).*1e3;

            % Resample impedance
            impedance_real_interp = interp1(freq,real_impedance,sampling_points,'linear',0);
            impedance_imag_interp = interp1(freq,imag_impedance,sampling_points,'linear',0);

            % Change parts and signs to match Elegant conventions
            impedance_real = -impedance_imag_interp;
            impedance_imag  = impedance_real_interp;      

            ImpedanceRealY = ImpedanceRealY + impedance_real.*average_betay(i);
            ImpedanceImagY = ImpedanceImagY + impedance_imag.*average_betay(i);
        end

    end
    fprintf('Finished reading vertical files. Read %d lines.\n\n',i);
    
    %% Save local beta functions
        
    save('local_beta','element_s','average_betax','average_betay');     
    
    %% Create output struct

    struct.ImpedanceFreq = sampling_points;
    struct.ImpedanceRealZ = ImpedanceRealZ;
    struct.ImpedanceImagZ = ImpedanceImagZ;
    struct.ImpedanceRealX = ImpedanceRealX;
    struct.ImpedanceImagX = ImpedanceImagX;
    struct.ImpedanceRealY = ImpedanceRealY;
    struct.ImpedanceImagY = ImpedanceImagY;

end