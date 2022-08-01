%% Calculate loss/kick factors for specified bunch length
%% The calculation uses the impedance files in the impedance lattice

%% Load lattice
impedance_lattice_openIDs_noNEG = impedance_M_H6BA_34_1_1(0,0);

%% For 3 mm bunch length

sigma_s = 3e-3;

output_3mm = loss_kick_factors(impedance_lattice_openIDs_noNEG,sigma_s);

fprintf('\nFor %10.5f mm bunch length:\n', sigma_s.*1e3)
fprintf('Total loss factor (V/pC):\n')
fprintf('RW: %10.5f, geometric: %10.5f\n',output_3mm.tot_RW_loss_factor.*1e-12,output_3mm.tot_geom_loss_factor.*1e-12)

fprintf('Total hor kick factor (V/pC/mm):\n')
fprintf('RW: %10.5f, geometric: %10.5f\n',output_3mm.tot_RW_hor_kick_factor.*1e-12.*1e-3,output_3mm.tot_geom_hor_kick_factor.*1e-12.*1e-3)

fprintf('Total ver kick factor (V/pC/mm):\n')
fprintf('RW: %10.5f, geometric: %10.5f\n',output_3mm.tot_RW_ver_kick_factor.*1e-12.*1e-3,output_3mm.tot_geom_ver_kick_factor.*1e-12.*1e-3)

%% For 12 mm bunch length

sigma_s = 12e-3;

output_12mm = loss_kick_factors(impedance_lattice_openIDs_noNEG,sigma_s);

fprintf('\nFor %10.5f mm bunch length:\n', sigma_s.*1e3)
fprintf('Total loss factor (V/pC):\n')
fprintf('RW: %10.5f, geometric: %10.5f\n',output_12mm.tot_RW_loss_factor.*1e-12,output_12mm.tot_geom_loss_factor.*1e-12)

fprintf('Total hor kick factor (V/pC/mm):\n')
fprintf('RW: %10.5f, geometric: %10.5f\n',output_12mm.tot_RW_hor_kick_factor.*1e-12.*1e-3,output_12mm.tot_geom_hor_kick_factor.*1e-12.*1e-3)

fprintf('Total ver kick factor (V/pC/mm):\n')
fprintf('RW: %10.5f, geometric: %10.5f\n',output_12mm.tot_RW_ver_kick_factor.*1e-12.*1e-3,output_12mm.tot_geom_ver_kick_factor.*1e-12.*1e-3)


%%
function output = loss_kick_factors(IMPRING,sigma_s)

        c = 299792458;
        sigma_t = sigma_s./c;

    % Find components and amount             
    components = unique(atgetfieldvalues(IMPRING,'FamName'));
    component_amount = zeros(length(components),1);
    for j = 1:length(components)
        elements = findcells(IMPRING,'FamName',components{j});
        component_amount(j) = length(elements);
    end
    
    loss_factors = zeros(length(components),2);
    hor_kick_factors = zeros(length(components),2);
    ver_kick_factors = zeros(length(components),2);
    
    
    for i = 1:length(components)
            
        % Get element in family
        element_index = findcells(IMPRING,'FamName',components{i});
        element_index = element_index(1);
            
        % Resistive-wall    
        if ~isempty(IMPRING{element_index}.RW_impedance_files)
                
            % Get length
            RW_length = IMPRING{element_index}.Length;
                
            % --- Hor kick factor ---
                
            % Using impedance
            hor_data = importdata(IMPRING{element_index}.RW_impedance_files{1}).data;
            % Switch to Elegant conventions
            freq = hor_data(:,1);
            real_hor_impedance = hor_data(:,3).*RW_length;
            imag_hor_impedance = -hor_data(:,2).*RW_length;               
            hor_kick_factors(i,1) = calculate_kick_factor('impedance',[freq,real_hor_impedance,imag_hor_impedance],sigma_t);   
                               
            % --- Ver kick factor ---

            % Using impedance
            ver_data = importdata(IMPRING{element_index}.RW_impedance_files{2}).data;
            % Switch to Elegant conventions
            freq = ver_data(:,1);
            real_ver_impedance = ver_data(:,3).*RW_length;
            imag_ver_impedance = -ver_data(:,2).*RW_length;               
            ver_kick_factors(i,1) = calculate_kick_factor('impedance',[freq,real_ver_impedance,imag_ver_impedance],sigma_t); 
                                 
            % --- Loss factor ---
                                
            % Using impedance
            lon_data = importdata(IMPRING{element_index}.RW_impedance_files{3}).data;
            freq = lon_data(:,1);
            real_lon_impedance = lon_data(:,2).*RW_length; % Multiply with length of element
            imag_lon_impedance = lon_data(:,3).*RW_length; % Multiply with length of element
            loss_factors(i,1) = calculate_loss_factor('impedance',[freq,real_lon_impedance,imag_lon_impedance],sigma_t);                        

        end            
                       
        % Geometric        
        if ~isempty(IMPRING{element_index}.Geom_impedance_files)
                
            % --- Hor kick factor ---
                
            % Using impedance
            fileID = fopen(IMPRING{element_index}.Geom_impedance_files{1});
            data = textscan(fileID,'%f %f %f','CommentStyle','#');
            fclose(fileID);
            hor_data = cell2mat(data);  % Units GHz, Ohm/mm
            % Change units to Hz and Ohm/m and switch to Elegant
            % conventions
            freq = hor_data(:,1).*1e9;
            real_hor_impedance = -hor_data(:,3).*1e3;
            imag_hor_impedance = hor_data(:,2).*1e3;               
            hor_kick_factors(i,2) = calculate_kick_factor('impedance',[freq,real_hor_impedance,imag_hor_impedance],sigma_t);                
                       
            % --- Ver kick factor ---
                                
            % Using impedance
            fileID = fopen(IMPRING{element_index}.Geom_impedance_files{2});
            data = textscan(fileID,'%f %f %f','CommentStyle','#');
            fclose(fileID);
            ver_data = cell2mat(data);  % Units GHz, Ohm/mm
            % Change units to Hz and Ohm/m and switch to Elegant
            % conventions
            freq = ver_data(:,1).*1e9;
            real_ver_impedance = -ver_data(:,3).*1e3;
            imag_ver_impedance = ver_data(:,2).*1e3;               
            ver_kick_factors(i,2) = calculate_kick_factor('impedance',[freq,real_ver_impedance,imag_ver_impedance],sigma_t);                
                         
            % --- Loss factor ---
                
            % Using impedance
            fileID = fopen(IMPRING{element_index}.Geom_impedance_files{3});
            data = textscan(fileID,'%f %f %f','CommentStyle','#');
            fclose(fileID);
            lon_data = cell2mat(data);  % Units GHz, Ohm
            % Change units to Hz
            freq = lon_data(:,1).*1e9;
            real_lon_impedance = lon_data(:,2);
            imag_lon_impedance = lon_data(:,3);
            loss_factors(i,2) = calculate_loss_factor('impedance',[freq,real_lon_impedance,imag_lon_impedance],sigma_t);
                                                  
        end
                                                             
    end
    
    % Total loss/kick factors per component
    tot_hor_kick_factors = component_amount.*hor_kick_factors;
    tot_ver_kick_factors = component_amount.*ver_kick_factors;        
    tot_loss_factors = component_amount.*loss_factors;
    
    % Total RW loss/kick factor 
    tot_RW_hor_kick_factor = sum(tot_hor_kick_factors(:,1));
    tot_RW_ver_kick_factor = sum(tot_ver_kick_factors(:,1));            
    tot_RW_loss_factor = sum(tot_loss_factors(:,1));    
    
    % Total geometric loss/kick factor 
    tot_geom_hor_kick_factor = sum(tot_hor_kick_factors(:,2));
    tot_geom_ver_kick_factor = sum(tot_ver_kick_factors(:,2));            
    tot_geom_loss_factor = sum(tot_loss_factors(:,2));      
    
    output.components = components;
    output.component_amount = component_amount;
    output.loss_factors = loss_factors;
    output.hor_kick_factors = hor_kick_factors;
    output.ver_kick_factors = ver_kick_factors;
    output.tot_RW_hor_kick_factor = tot_RW_hor_kick_factor;
    output.tot_RW_ver_kick_factor = tot_RW_ver_kick_factor;
    output.tot_RW_loss_factor = tot_RW_loss_factor;
    output.tot_geom_hor_kick_factor = tot_geom_hor_kick_factor;
    output.tot_geom_ver_kick_factor = tot_geom_ver_kick_factor;
    output.tot_geom_loss_factor = tot_geom_loss_factor;    
    
end