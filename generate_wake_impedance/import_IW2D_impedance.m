function struct = import_IW2D_impedance(file,sampling_points,betas)
%% Process impedanceWake2D input

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
    
    fprintf('Reading IW2D impedance files\n');
    
    fprintf('Reading longitudinal files\n');
    for i = 1:size(data_lon{1},1)
        
        RW_length = data_lon{2}(i);
        filename = char(data_lon{4}(i));
        
        if ~isempty(filename)
            lon_impedance = importdata(filename);
            lon_impedance = lon_impedance.data;

            % Resample impedance
            impedance_real = interp1(lon_impedance(:,1),lon_impedance(:,2),sampling_points,'linear',0);
            impedance_imag = interp1(lon_impedance(:,1),lon_impedance(:,3),sampling_points,'linear',0);        

            % Multiply impedance with length
            impedance_real  = impedance_real.*RW_length;        
            impedance_imag  = impedance_imag.*RW_length;

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
        
        RW_length = data_hor{2}(i);
        filename = char(data_hor{4}(i));
        
        hor_impedance = importdata(filename);
        hor_impedance = hor_impedance.data;
        
        % Resample impedance
        impedance_real_interp = interp1(hor_impedance(:,1),hor_impedance(:,2),sampling_points,'linear',0);
        impedance_imag_interp = interp1(hor_impedance(:,1),hor_impedance(:,3),sampling_points,'linear',0);
        
        % Multiply impedance with length and change parts and signs to
        % match Elegant conventions
        % Elegant imaginary = - real IW2D
        % Elegant real = imaginary IW2D
        impedance_real  = impedance_imag_interp.*RW_length;        
        impedance_imag  = -impedance_real_interp.*RW_length;
               
        ImpedanceRealX = ImpedanceRealX + impedance_real.*average_betax(i);
        ImpedanceImagX = ImpedanceImagX + impedance_imag.*average_betax(i);        
                                 
    end
    fprintf('Finished reading horizontal files. Read %d lines.\n\n',i);
    
    %% Vertical impedance
    
    if ~isempty(betas)                
        % Calculate average horizontal beta over element lengths       
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
        
        RW_length = data_ver{2}(i);
        filename = char(data_ver{4}(i));
                    
        ver_impedance = importdata(filename);
        ver_impedance = ver_impedance.data;
        
        % Resample impedance
        impedance_real_interp = interp1(ver_impedance(:,1),ver_impedance(:,2),sampling_points,'linear',0);
        impedance_imag_interp = interp1(ver_impedance(:,1),ver_impedance(:,3),sampling_points,'linear',0);
        
        % Multiply impedance with length and change parts and signs to
        % match Elegant conventions
        % Elegant imaginary = - real IW2D
        % Elegant real = imaginary IW2D
        impedance_real  = impedance_imag_interp.*RW_length;         
        impedance_imag  = -impedance_real_interp.*RW_length;
                      
        ImpedanceRealY = ImpedanceRealY + impedance_real.*average_betay(i);
        ImpedanceImagY = ImpedanceImagY + impedance_imag.*average_betay(i);          
        
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