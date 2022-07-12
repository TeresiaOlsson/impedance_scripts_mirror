function impedance_database_summary(IMPRING,aperture_flag,pie_flag,impedance_flag)
%% Produce summary statistics for impedance database
% aperture_flag - plot apertures
% pie_flag - calculate loss/kick factoran and plot pie charts
% impedance_flag - plot impedance

    %% Resistive wall statistics

    % Get all materials in database
    values = atgetfieldvalues(IMPRING,'Material');
    
    % Handle several materials
    for i = 1:length(values)
        if length(values{i}) > 1
            str = '';
            for j = 1:length(values{i})
                if j < length(values{i})
                    new_str = [values{i}{j}, ' + '];
                else
                    new_str = values{i}{j};
                end
                str = append(str,new_str);
            end
            values{i} = str;
        else
            values{i} = char(values{i});
        end           
    end
    
    materials = values;
    materials_unique = unique(materials);

    % Get total length of each material
    material_lengths = zeros(length(materials_unique),1);

    for i = 1:length(materials_unique)        
        elements = find(strcmp(materials_unique(i),materials));
        material_lengths(i) = sum(atgetfieldvalues(IMPRING,elements,'Length'));
    end

    % Total length
    total_length = sum(material_lengths);

    % Print out
    fprintf('Resistive wall:\n')
    for i = 1:length(materials_unique)
        fprintf('%s: %.5f m (%.2f%%)\n',materials_unique{i},material_lengths(i),(material_lengths(i)./total_length).*100)
    end
    fprintf('Total length = %f\n',total_length)
    fprintf('\n')

    %% Geometric impedance statistics
    
    % Get all classes in database
    classes = unique(atgetfieldvalues(IMPRING,'Class'));
    
    % Print out
    fprintf('Components:\n')
    for i = 1:length(classes)
        
        % Get elements in class
        elements_in_class = findcells(IMPRING,'Class',classes{i});
        
        % Get components in class and amount of each
        components = unique(atgetfieldvalues(IMPRING(elements_in_class),'FamName'));
        component_amount = zeros(length(components),1);
        for j = 1:length(components)
            elements = findcells(IMPRING,'FamName',components{j});
            component_amount(j) = length(elements);
        end
                
        fprintf('%s:\n',classes{i})
        for j = 1:length(components)
            fprintf('%s: %d \n',components{j},component_amount(j))
        end
        fprintf('\n')        
    end    %% Apertures
    
    if aperture_flag

        color=get(gca,'ColorOrder');

        % Find position of entrance of each element
        s = findspos(IMPRING,1:length(IMPRING))';
        lengths = atgetfieldvalues(IMPRING,'Length');

        apertures = atgetfieldvalues(IMPRING,'Apertures');
        apertures  = cat(1,apertures{:});
        entrance_apertures  = apertures(1:2:end,:);
        exit_apertures  = apertures(2:2:end,:);

        figure(1)
        plot(s,entrance_apertures(:,1)*1e3,'.-','Color',color(1,:))
        hold on
        plot(s+lengths,exit_apertures(:,1)*1e3,'.-','Color',color(1,:))
        plot(s,entrance_apertures(:,2)*1e3,'.-','Color',color(1,:))
        plot(s+lengths,exit_apertures(:,2)*1e3,'.-','Color',color(1,:))
        hold off
        xlabel('Position [m]')
        ylabel('Horizontal aperture [mm]')
        xlim([0,total_length])
        %xlim([0,total_length./6])

        figure(2)
        plot(s,entrance_apertures(:,3)*1e3,'.-','Color',color(1,:))
        hold on
        plot(s+lengths,exit_apertures(:,3)*1e3,'.-','Color',color(1,:))
        plot(s,entrance_apertures(:,4)*1e3,'.-','Color',color(1,:))
        plot(s+lengths,exit_apertures(:,4)*1e3,'.-','Color',color(1,:))
        hold off
        xlabel('Position [m]')
        ylabel('Vertical aperture [mm]')
        xlim([0,total_length])
        %xlim([0,total_length./6])
    end
    
    %% Plot apertures
    
    if aperture_flag

        color=get(gca,'ColorOrder');

        % Find position of entrance of each element
        s = findspos(IMPRING,1:length(IMPRING))';
        lengths = atgetfieldvalues(IMPRING,'Length');

        apertures = atgetfieldvalues(IMPRING,'Apertures');
        apertures  = cat(1,apertures{:});
        upstream_apertures  = apertures(1:2:end,:);
        downstream_apertures  = apertures(2:2:end,:);

        figure(1)
        plot(s,upstream_apertures(:,1)*1e3,'.-','Color',color(1,:))
        hold on
        plot(s+lengths,downstream_apertures(:,1)*1e3,'.-','Color',color(2,:))
        plot(s,upstream_apertures(:,2)*1e3,'.-','Color',color(1,:))
        plot(s+lengths,downstream_apertures(:,2)*1e3,'.-','Color',color(2,:))
        hold off
        xlabel('Position [m]')
        ylabel('Horizontal aperture [mm]')
        xlim([0,total_length])
        legend('Upstream','Downstream')
        %xlim([0,total_length./6])
        
        datacursormode on
        dcm = datacursormode(gcf);         
        set(dcm,'UpdateFcn',@(t,e) myupdatefcn(t,e,IMPRING,s) );
        
        figure(2)
        plot(s,upstream_apertures(:,3)*1e3,'.-','Color',color(1,:))
        hold on
        plot(s+lengths,downstream_apertures(:,3)*1e3,'.-','Color',color(2,:))
        plot(s,upstream_apertures(:,4)*1e3,'.-','Color',color(1,:))
        plot(s+lengths,downstream_apertures(:,4)*1e3,'.-','Color',color(2,:))
        hold off
        xlabel('Position [m]')
        ylabel('Vertical aperture [mm]')
        xlim([0,total_length])
        %xlim([0,total_length./6])
        
        datacursormode on
        dcm = datacursormode(gcf);         
        set(dcm,'UpdateFcn',@(t,e) myupdatefcn(t,e,IMPRING,s) );        
    end    
                
    %% Calculate loss/kick factors and plot pie charts
    % IW2D input - loss/kick factors calculated using impedance since wake
    % already convoluted with unknown bunch length
    % CST input - loss/kick factors calculated using impedance but factors
    % using wake also calculated for comparison with CST. For the wake
    % calculations sigma_t must match the bunch length that was used in CST
    % or the result will be wrong.
              
    if pie_flag
                
        % Find components and amount             
        components = unique(atgetfieldvalues(IMPRING,'FamName'));
        component_amount = zeros(length(components),1);
        for j = 1:length(components)
            elements = findcells(IMPRING,'FamName',components{j});
            component_amount(j) = length(elements);
        end
        
        % Calculate loss/kick factors for each component
        sigma_s = 3e-3; % This must match bunch length used in CST calculations.
        c = 299792458;
        sigma_t = sigma_s./c;
        
        loss_factors = zeros(length(components),3);
        hor_kick_factors = zeros(length(components),3);
        ver_kick_factors = zeros(length(components),3);
                    
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
                             
                % Using wake
                fileID = fopen(IMPRING{element_index}.Geom_wake_files{1});
                data = textscan(fileID,'%f %f','CommentStyle','#');
                fclose(fileID);          
                hor_data = cell2mat(data);  % Units mm, V/pC/mm
                % Change units to m and V/C/m and sign to AT convention
                hor_s = hor_data(:,1).*1e-3;
                hor_wake = -hor_data(:,2).*1e12.*1e3;
                hor_kick_factors(i,3) = calculate_kick_factor('wake',[hor_s,hor_wake],sigma_t);
                       
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
                             
                % Using wake
                fileID = fopen(IMPRING{element_index}.Geom_wake_files{2});
                data = textscan(fileID,'%f %f','CommentStyle','#');
                fclose(fileID);          
                ver_data = cell2mat(data);  % Units mm, V/pC/mm
                % Change units to m and V/C/m and sign to AT convention
                ver_s = ver_data(:,1).*1e-3;
                ver_wake = -ver_data(:,2).*1e12.*1e3;
                ver_kick_factors(i,3) = calculate_kick_factor('wake',[ver_s,ver_wake],sigma_t);
                                          
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
                                      
                % Using wake
                fileID = fopen(IMPRING{element_index}.Geom_wake_files{3});
                data = textscan(fileID,'%f %f','CommentStyle','#');
                fclose(fileID);  
                lon_data = cell2mat(data);  % Units mm, V/pC
                % Change units to m and V/C and sign of wake to match AT
                % conventions
                lon_s = lon_data(:,1).*1e-3;
                lon_wake = -lon_data(:,2).*1e12;
                loss_factors(i,3) = calculate_loss_factor('wake',[lon_s,lon_wake],sigma_t);                

            end
                                                             
        end
        
        % Total loss/kick factors per component
        tot_hor_kick_factors = component_amount.*hor_kick_factors;
        tot_ver_kick_factors = component_amount.*ver_kick_factors;        
        tot_loss_factors = component_amount.*loss_factors;
                
        % Total geometric loss/kick factor per class
        classes = unique(atgetfieldvalues(IMPRING,'Class'));
        tot_geom_hor_kick_factor = zeros(length(classes),1);
        tot_geom_ver_kick_factor = zeros(length(classes),1);        
        tot_geom_loss_factor = zeros(length(classes),1);
        
        for i = 1:length(classes)
            
            % Find families in class
             elements = findcells(IMPRING,'Class',classes{i});
             families = unique(atgetfieldvalues(IMPRING,elements,'FamName'));     
             
             [~,index] = ismember(families,components);
             tot_geom_hor_kick_factor(i) = sum(tot_hor_kick_factors(index,2));
             tot_geom_ver_kick_factor(i) = sum(tot_ver_kick_factors(index,2));                 
             tot_geom_loss_factor(i) = sum(tot_loss_factors(index,2));
            
        end
        
        % Total RW loss/kick factor 
        tot_RW_hor_kick_factor = sum(tot_hor_kick_factors(:,1));
        tot_RW_ver_kick_factor = sum(tot_ver_kick_factors(:,1));            
        tot_RW_loss_factor = sum(tot_loss_factors(:,1));
                      
        % Plot pie charts
        
        % Find four largest classes and combine the rest into others
        [~,max_index_hor] = maxk(abs(tot_geom_hor_kick_factor),4);
        [~,max_index_ver] = maxk(abs(tot_geom_ver_kick_factor),4);
        [~,max_index_lon] = maxk(tot_geom_loss_factor,4);
        
        hor_kick_factor_plot = abs([tot_RW_hor_kick_factor;tot_geom_hor_kick_factor(max_index_hor);sum(tot_geom_hor_kick_factor(setdiff(1:end, max_index_hor)))]);
        ver_kick_factor_plot = abs([tot_RW_ver_kick_factor;tot_geom_ver_kick_factor(max_index_ver);sum(tot_geom_ver_kick_factor(setdiff(1:end, max_index_ver)))]);
        loss_factor_plot = [tot_RW_loss_factor;tot_geom_loss_factor(max_index_lon);sum(tot_geom_loss_factor(setdiff(1:end, max_index_lon)))];
                

        figure(3)
        pie(hor_kick_factor_plot)
        legend(['RW';classes(max_index_hor);'OTHER'],'Location','NorthWestOutside')
        title('Horizontal kick factor')

        figure(4)
        pie(ver_kick_factor_plot)
        legend(['RW';classes(max_index_ver);'OTHER'],'Location','NorthWestOutside')
        title('Vertical kick factor')

        figure(5)
        pie(loss_factor_plot)
        legend(['RW';classes(max_index_lon);'OTHER'],'Location','NorthWestOutside')
        title('Loss factor')
        
        % Write kick/loss factors
        fileID = fopen('kick-loss-factors_resistive_wall.txt','w');
        fprintf(fileID,'Component Number in ring Loss factor impedance (V/pC) Hor. kick factor (V/pC/mm) Ver. kick factor (V/pC/mm)\n');
        for i = 1:length(components)
            fprintf(fileID,'%-15s %10d %10.5f %10.5f %10.5f\n',components{i},component_amount(i),loss_factors(i,1).*1e-12,hor_kick_factors(i,1).*1e-12.*1e-3,ver_kick_factors(i,1).*1e-12.*1e-3);
        end
        fprintf(fileID,'%-15s %10d %10.5f %10.5f %10.5f\n','','',tot_RW_loss_factor.*1e-12,tot_RW_hor_kick_factor.*1e-12.*1e-3,tot_RW_ver_kick_factor.*1e-12.*1e-3);
        fclose(fileID);
                
        fileID = fopen('kick-loss-factors_geometric.txt','w');
        fprintf(fileID,'Component Number in ring Loss factor impedance (V/pC) Loss factor wake (V/pC) Hor. kick factor impedance (V/pC/mm) Hor. kick factor wake (V/pC/mm) Ver. kick factor impedance (V/pC/mm) Ver. kick factor wake (V/pC/mm)\n');
        for i = 1:length(components)
            fprintf(fileID,'%-15s %10d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',components{i},component_amount(i),loss_factors(i,2).*1e-12,loss_factors(i,3).*1e-12,hor_kick_factors(i,2).*1e-12.*1e-3,hor_kick_factors(i,3).*1e-12.*1e-3,ver_kick_factors(i,2).*1e-12.*1e-3,ver_kick_factors(i,3).*1e-12.*1e-3);
        end
        fprintf(fileID,'%-15s %10d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n','','',sum(tot_loss_factors(:,2)).*1e-12,sum(tot_loss_factors(:,3)).*1e-12,sum(tot_hor_kick_factors(:,2)).*1e-12.*1e-3,sum(tot_hor_kick_factors(:,3)).*1e-12.*1e-3,sum(tot_ver_kick_factors(:,2).*1e-12.*1e-3),sum(tot_ver_kick_factors(:,3).*1e-12.*1e-3));
        fclose(fileID);
        
    end

    %% Plot cumulative impedance
    
    if impedance_flag

        % Interpolation parameters
        sampling_freq = linspace(0,35,1e4)'; %Units are in GHz
    
        % Find classes
        classes = unique(atgetfieldvalues(IMPRING,'Class'));
        labels = classes;

        % Get impedance per class
        hor_real_impedance = zeros(length(sampling_freq),length(classes));
        hor_imag_impedance = zeros(length(sampling_freq),length(classes));
        ver_real_impedance = zeros(length(sampling_freq),length(classes));
        ver_imag_impedance = zeros(length(sampling_freq),length(classes));
        lon_real_impedance = zeros(length(sampling_freq),length(classes));
        lon_imag_impedance = zeros(length(sampling_freq),length(classes));

        for i = 1:length(classes)

            elements = findcells(IMPRING,'Class',classes{i});

            for j = 1:length(elements)

                if ~isempty(IMPRING{elements(j)}.Geom_impedance_files)

                    % Import horizontal impedance data
                    hor_data = importdata(IMPRING{elements(j)}.Geom_impedance_files{1}).data;

                    real_impedance = interp1(hor_data(:,1),hor_data(:,2),sampling_freq,'linear',0);
                    imag_impedance = interp1(hor_data(:,1),hor_data(:,3),sampling_freq,'linear',0);

                    hor_real_impedance(:,i) = hor_real_impedance(:,i) + real_impedance;
                    hor_imag_impedance(:,i) = hor_imag_impedance(:,i) + imag_impedance;

                    % Import vertical impedance data
                    ver_data = importdata(IMPRING{elements(j)}.Geom_impedance_files{2}).data;

                    real_impedance = interp1(ver_data(:,1),ver_data(:,2),sampling_freq,'linear',0);
                    imag_impedance = interp1(ver_data(:,1),ver_data(:,3),sampling_freq,'linear',0);

                    ver_real_impedance(:,i) = ver_real_impedance(:,i) + real_impedance;
                    ver_imag_impedance(:,i) = ver_imag_impedance(:,i) + imag_impedance; 

                    % Import longitudinal impedance data
                    lon_data = importdata(IMPRING{elements(j)}.Geom_impedance_files{3}).data;

                    real_impedance = interp1(lon_data(:,1),lon_data(:,2),sampling_freq,'linear',0);
                    imag_impedance = interp1(lon_data(:,1),lon_data(:,3),sampling_freq,'linear',0);

                    lon_real_impedance(:,i) = lon_real_impedance(:,i) + real_impedance;
                    lon_imag_impedance(:,i) = lon_imag_impedance(:,i) + imag_impedance;               

                end

            end

        end

        figure(6)
        subplot(211)
        area(sampling_freq,-hor_real_impedance.*1e3.*1e-3,'LineStyle','None')
        ylabel('Re Z_x [k\Omega/m]')
        xlabel('Frequency [GHz]')
        legend(labels)

        subplot(212)
        area(sampling_freq,-hor_imag_impedance.*1e3.*1e-3,'LineStyle','None')
        ylabel('Im Z_x [k\Omega/m]')
        xlabel('Frequency [GHz]')
        legend(labels)

        figure(7)
        subplot(211)
        area(sampling_freq,-ver_real_impedance.*1e3.*1e-3,'LineStyle','None')
        ylabel('Re Z_y [k\Omega/m]')
        xlabel('Frequency [GHz]')
        legend(labels)

        subplot(212)
        area(sampling_freq,-ver_imag_impedance.*1e3.*1e-3,'LineStyle','None')
        ylabel('Im Z_y [k\Omega/m]')
        xlabel('Frequency [GHz]')
        legend(labels)

        figure(8)
        subplot(211)
        area(sampling_freq,lon_real_impedance.*1e-3,'LineStyle','None')
        ylabel('Re Z_L [k\Omega]')
        xlabel('Frequency [GHz]')
        legend(labels)

        subplot(212)
        area(sampling_freq,lon_imag_impedance.*1e-3,'LineStyle','None')
        ylabel('Im Z_L [k\Omega]')
        xlabel('Frequency [GHz]')
        legend(labels)
               
    end


end

function txt = myupdatefcn(~,event,ring,s)
    pos = get(event,'Position');
    i = find(s==pos(1));
    element = cell2mat(atgetfieldvalues(ring,i,'FamName'));
    txt = {[element], ...
        ['x: ',num2str(pos(1))],...
        ['y: ',num2str(pos(2))]};        
end


    
