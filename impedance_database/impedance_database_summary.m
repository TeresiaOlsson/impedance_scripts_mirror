function impedance_database_summary(IMPRING,aperture_flag,pie_flag,impedance_flag)
%% Produce summary statistics for impedance database
% aperture_flag - plot apertures
% pie_flag - plot pie charts
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
                
    %% Plot pie charts
    
    if pie_flag

        % Find classes
        classes = unique(atgetfieldvalues(IMPRING,'Class'));

        % Get total k factors per class
        tot_k_factors = zeros(length(classes),3);

        for i = 1:length(classes)

            elements = findcells(IMPRING,'Class',classes{i});    
            k_factors = cell2mat(atgetfieldvalues(IMPRING,elements,'Geom_k_factors'));    
            tot_k_factors(i,:) = tot_k_factors(i,:) + sum(abs(k_factors));

        end
        
        % Plot 4 largest
        
        labels = classes;
        
        [~,max_index_hor] = maxk(tot_k_factors(:,1),4);
        [~,max_index_ver] = maxk(tot_k_factors(:,2),4);
        [~,max_index_lon] = maxk(tot_k_factors(:,3),4);        

        figure(3)
        pie(tot_k_factors(max_index_hor,1))
        legend(labels(max_index_hor),'Location','NorthWestOutside')
        title('Horizontal kick factor')
        colormap('lines')

        figure(4)
        pie(tot_k_factors(max_index_ver,2))
        legend(labels(max_index_ver),'Location','NorthWestOutside')
        title('Vertical kick factor')
        colormap('lines')

        figure(5)
        pie(tot_k_factors(max_index_lon,3))
        legend(labels(max_index_lon),'Location','NorthWestOutside')
        title('Loss factor')
        colormap('lines')
        
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


    
