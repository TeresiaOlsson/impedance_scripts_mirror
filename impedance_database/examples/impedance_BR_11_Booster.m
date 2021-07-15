function IMPRING = impedance_BR_11(ID_config,neg_model)
%% Impedance database for lattice M-H6BA-34-1-1

%% ID configurations
% ID_config = 0 - All IDs open
% ID_config = 1 - All IDs closed

switch ID_config
    
    case 0
        ID_ver_aper = [-0.015 0.015];        
    case 1
        ID_ver_aper = [-0.002 0.002];
        
end
       
%% NEG models
% neg_model = 0 - no NEG coating
% neg_model = 1 - conductivity 1e6 S/m, thickness 1.0 um
% neg_model = 2 - conductivity 1e5 S/m, thickness 1.0 um
% neg_model = 3 - conductivity 1e6 S/m, thickness 1.0 um

switch neg_model
    
    case 0
        coating_material = ['Copper'];
        coating_conductivity = [5.96e7];
        coating_thickness = [Inf];
        
        RW_impedance_circular = {'/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-range0.1m-step1um/ZxdipWDiamond-II_1layers10.00mm_round_copper_1layer.dat'...
                           '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-range0.1m-step1um/ZydipWDiamond-II_1layers10.00mm_round_copper_1layer.dat'...
                           '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-range0.1m-step1um/ZlongWDiamond-II_1layers10.00mm_round_copper_1layer.dat'};
        
        RW_wake_circular = {'/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-range0.1m-step1um/WxdipWDiamond-II_1layers10.00mm_round_copper_1layer.dat'...
                            '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-range0.1m-step1um/WydipWDiamond-II_1layers10.00mm_round_copper_1layer.dat'...
                            '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-range0.1m-step1um/WlongWDiamond-II_1layers10.00mm_round_copper_1layer.dat'};
                        
        RW_impedance_rectangular = {'/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-range0.1m-step1um/ZxdipWDiamond-II_1layers15.00mm_rectangular_copper_1layer.dat'...
                              '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-range0.1m-step1um/ZydipWDiamond-II_1layers15.00mm_rectangular_copper_1layer.dat'...
                              '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-range0.1m-step1um/ZlongWDiamond-II_1layers15.00mm_rectangular_copper_1layer.dat'};
        
        RW_wake_rectangular = {'/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-range0.1m-step1um/WxdipWDiamond-II_1layers15.00mm_rectangular_copper_1layer.dat'...
                               '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-range0.1m-step1um/WyipWDiamond-II_1layers15.00mm_rectangular_copper_1layer.dat'...
                               '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-range0.1m-step1um/WlongWDiamond-II_1layers15.00mm_rectangular_copper_1layer.dat'};                                
    case 1
        coating_material = ['NEG','Copper'];
        coating_conductivity = [1e6,5.96e7];
        coating_thickness = [1e-6,Inf];
        
        RW_impedance_circular = {'/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e6-1.0um-range0.1m-step1um/ZxdipWDiamond-II_2layers10.00mm_round_copper_neg.dat'...
                           '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e6-1.0um-range0.1m-step1um/ZydipWDiamond-II_2layers10.00mm_round_copper_neg.dat'...
                           '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e6-1.0um-range0.1m-step1um/ZlongWDiamond-II_2layers10.00mm_round_copper_neg.dat'};
        
        RW_wake_circular = {'/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e6-1.0um-range0.1m-step1um/WxdipWDiamond-II_2layers10.00mm_round_copper_neg.dat'...
                            '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e6-1.0um-range0.1m-step1um/WydipWDiamond-II_2layers10.00mm_round_copper_neg.dat'...
                            '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e6-1.0um-range0.1m-step1um/WlongWDiamond-II_2layers10.00mm_round_copper_neg.dat'};
                        
        RW_impedance_rectangular = {'/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e6-1.0um-range0.1m-step1um/ZxdipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                              '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e6-1.0um-range0.1m-step1um/ZydipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                              '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e6-1.0um-range0.1m-step1um/ZlongWDiamond-II_2layers15.00mm_round_copper_neg.dat'};
        
        RW_wake_rectangular = {'/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e6-1.0um-range0.1m-step1um/WxdipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                               '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e6-1.0um-range0.1m-step1um/WydipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                               '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e6-1.0um-range0.1m-step1um/WlongWDiamond-II_2layers15.00mm_round_copper_neg.dat'};         
    case 2
        coating_material = ['NEG','Copper'];
        coating_conductivity = [1e5,5.96e7];
        coating_thickness = [1e-6,Inf];
        
        RW_impedance_circular = {'/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e5-1.0um-range0.1m-step1um/ZxdipWDiamond-II_2layers10.00mm_round_copper_neg.dat'...
                           '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e5-1.0um-range0.1m-step1um/ZydipWDiamond-II_2layers10.00mm_round_copper_neg.dat'...
                           '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e5-1.0um-range0.1m-step1um/ZlongWDiamond-II_2layers10.00mm_round_copper_neg.dat'};
        
        RW_wake_circular = {'/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e5-1.0um-range0.1m-step1um/WxdipWDiamond-II_2layers10.00mm_round_copper_neg.dat'...
                            '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e5-1.0um-range0.1m-step1um/WydipWDiamond-II_2layers10.00mm_round_copper_neg.dat'...
                            '/dls/physics/Impedance/IW2D/copper-circular-10mm-inf-neg-1e5-1.0um-range0.1m-step1um/WlongWDiamond-II_2layers10.00mm_round_copper_neg.dat'};
                        
        RW_impedance_rectangular = {'/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/ZxdipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                              '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/ZydipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                              '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/ZlongWDiamond-II_2layers15.00mm_round_copper_neg.dat'};
        
        RW_wake_rectangular = {'/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/WxdipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                               '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/WydipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                               '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/WlongWDiamond-II_2layers15.00mm_round_copper_neg.dat'};                 
        
    case 3
        coating_material = ['Steel316'];
        coating_conductivity = [1.3514e6];
        coating_thickness = [Inf];
        
        RW_impedance_circular = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Resistive wall\\rw_1m_12.5mm_stainless_impx.txt'...
                           'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Resistive wall\\rw_1m_12.5mm_stainless_impx.txt'...
                           'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Resistive wall\\rw_1m_12.5mm_stainless_impz.txt'};
        
        RW_wake_circular = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Resistive wall\\rw_1m_12.5mm_stainless_wakex.txt'...
                            'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Resistive wall\\rw_1m_12.5mm_stainless_wakex.txt'...
                            'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Resistive wall\\rw_1m_12.5mm_stainless_wakez.txt'};
                        
        RW_impedance_rectangular = {'/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/ZxdipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                              '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/ZydipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                              '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/ZlongWDiamond-II_2layers15.00mm_round_copper_neg.dat'};
        
        RW_wake_rectangular = {'/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/WxdipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                               '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/WydipWDiamond-II_2layers15.00mm_round_copper_neg.dat'...
                               '/dls/physics/Impedance/IW2D/copper-rectangular-36mm-15mm-inf-neg-1e5-1.0um-range0.1m-step1um/WlongWDiamond-II_2layers15.00mm_round_copper_neg.dat'};                 
        
end

%% Drifts

% --- Additional drifts ---  
name = 'DRIFT1';
class = 'DRIFT';
length = 0.3;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT1 = {impedance(name, class, length, material, conductivity, thickness, apertures,...
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

name = 'DRIFT2';
class = 'DRIFT';
length = 0.2;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT2 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

name = 'DRIFT3';
class = 'DRIFT';
length = 0.015;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT3 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

name = 'DRIFT4';
class = 'DRIFT';
length = 1.0;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT4 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

name = 'DRIFT5';
class = 'DRIFT';
length = 0.15;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT5 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

name = 'DRIFT6';
class = 'DRIFT';
length = 0.1;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT6 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

name = 'DRIFT7';
class = 'DRIFT';
length = 0.65;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT7 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

name = 'DRIFT8';
class = 'DRIFT';
length = 0.07875;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT8 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT9';
class = 'DRIFT';
length = 0.185;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT9 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT10';
class = 'DRIFT';
length = 0.04;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT10 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT11';
class = 'DRIFT';
length = 0.025;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT11 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT12';
class = 'DRIFT';
length = 0.0529;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT12 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT13';
class = 'DRIFT';
length = 0.0925;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT13 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT14';
class = 'DRIFT';
length = 0.085;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT14 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT15';
class = 'DRIFT';
length = 0.0031;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT15 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT16';
class = 'DRIFT';
length = 0.0409;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT16 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT17';
class = 'DRIFT';
length = 0.0341;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT17 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT18';
class = 'DRIFT';
length = 0.01;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT18 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT19';
class = 'DRIFT';
length = 0.0619;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT19 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT20';
class = 'DRIFT';
length = 0.03;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT20 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT21';
class = 'DRIFT';
length = 0.0775;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT21 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT22';
class = 'DRIFT';
length = 0.05;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT22 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DRIFT24';
class = 'DRIFT';
length = 0.0216;
material = coating_material;
conductivity = coating_conductivity;
thickness = coating_thickness;
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = RW_impedance_circular;
RW_wake_files = RW_wake_circular;
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [0 0 0];
geom_fit = struct([]);

DRIFT24 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

%% Dipole vessel
% Includes BPM3, keyhole flange, pump

name = 'DIPOLE1';
class = 'DIPOLE';
length = 1.3;
material = ['Aluminium'];
conductivity = [3.77e7];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [-0.195,-0.0301,0.117];
geom_fit = struct([]);

DIPOLE1 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'DIPOLE2';
class = 'DIPOLE';
length = 1.25;
material = ['Aluminium'];
conductivity = [3.77e7];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [-0.195,-0.0301,0.117];
geom_fit = struct([]);

DIPOLE2 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
%% Septum vessel

name = 'INJSEPT';
class = 'SETPUM';
length = 0.3;
material = ['Aluminium'];
conductivity = [3.77e7];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [];
geom_fit = struct([]);

INJSEPT = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'EXSEPT';
class = 'SETPUM';
length = 1.67;
material = ['Aluminium'];
conductivity = [3.77e7];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [];
geom_fit = struct([]);

EXSEPT = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
%% Kicker

name = 'HCM';
class = 'KICKER';
length = 0.075;
material = ['Aluminium'];
conductivity = [3.77e7];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [];
geom_fit = struct([]);

HCM = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
name = 'VCM';
class = 'KICKER';
length = 0.075;
material = ['Aluminium'];
conductivity = [3.77e7];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [];
geom_fit = struct([]);

VCM = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
%% Quad

name = 'QUAD';
class = 'QUADRUPOLE';
length = 0.45;
material = ['Aluminium'];
conductivity = [3.77e7];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [];
geom_fit = struct([]);

QUAD = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
%% Sextupole

name = 'SEXT';
class = 'SEXTUPOLE';
length = 0.05;
material = ['Aluminium'];
conductivity = [3.77e7];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [];
geom_fit = struct([]);

SEXT = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

%% BPM assemblies

% --- BPM 1 ---
name = 'BPM1';
class = 'BPM';
length = 0.05;
material = ['Steel316'];
conductivity = [1.3514e6];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\BPM\\boosterII_bpm_10.9mm_pec_impx.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\BPM\\boosterII_bpm_10.9mm_pec_impy.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\BPM\\boosterII_bpm_10.9mm_pec_impz.txt'};
geom_wake_files = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\BPM\\boosterII_bpm_10.9mm_pec_wx.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\BPM\\boosterII_bpm_10.9mm_pec_wy.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\BPM\\boosterII_bpm_10.9mm_pec_wz.txt'};
geom_k_factors = [-0.0298,-0.0127,0.0806];
geom_fit = struct([]);

BPM1 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
%% Cavities

% --- Main cavity ---

name = 'MAINCAV';
class = 'CAVITY';
length = 1.8;
material = ['Copper'];
conductivity = [5.96e7];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\NCC\\NBCWF_ImpX_ReIm.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\NCC\\NBCWF_ImpY_ReIm.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\NCC\\NBCWF_ImpZ_ReIm.txt'};
geom_wake_files = {};
geom_k_factors = [];
geom_fit = struct([]);

MAINCAV = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};


%% Valves

name = 'VALVE1';
class = 'VALVE';
length = 0.1;
material = ['Steel316'];
conductivity = [1.3514e6];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {};
geom_wake_files = {};
geom_k_factors = [];
geom_fit = struct([]);

VALVE1 = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

%% Vacuum Pump

name = 'VPCHAMBER';
class = 'PUMP';
length = 0.15;
material = ['Steel316'];
conductivity = [1.3514e6];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Pump\\boosterII_pumpsmall_v2_pec_impx.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Pump\\boosterII_pumpsmall_v2_pec_impy.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Pump\\boosterII_pumpsmall_v2_pec_impz.txt'};
geom_wake_files = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Pump\\boosterII_pumpsmall_v2_pec_wx.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Pump\\boosterII_pumpsmall_v2_pec_wy.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Pump\\boosterII_pumpsmall_v2_pec_wz.txt'};
geom_k_factors = [];
geom_fit = struct([]);

VPCHAMBER = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};  

%% Diagnostics

name = 'SCREEN';
class = 'SCREEN';
length = 0.2;
material = ['Steel316'];
conductivity = [1.3514e6];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Screen\\boosterII_screen_v2_impx_pec.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Screen\\boosterII_screen_v2_impy_pec.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Screen\\boosterII_screen_v2_impz_pec.txt'};
geom_wake_files = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Screen\\boosterII_screen_v2_wx_pec.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Screen\\boosterII_screen_v2_wy_pec.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Screen\\boosterII_screen_v2_wz_pec.txt'};
geom_k_factors = [-0.00564,-0.0126,0.0821];
geom_fit = struct([]);

SCREEN = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};
  
%% Ceramic Break

name = 'CBREAK';
class = 'CBREAK';
length = 0.0381;
material = ['Steel316'];
conductivity = [1.3514e6];
thickness = [Inf];
apertures = [-0.01 0.01 -0.01 0.01; -0.01 0.01 -0.01 0.01];
RW_impedance_files = {};
RW_wake_files = {};
RW_k_factors = [0 0 0];
geom_impedance_files = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Ceramic\\boosterII_ceramic_impx.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Ceramic\\boosterII_ceramic_impx.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Ceramic\\boosterII_ceramic_impz.txt'};
geom_wake_files = {'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Ceramic\\boosterII_ceramic_wx.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Ceramic\\boosterII_ceramic_wx.txt',...
                        'E:\\DIAMOND\\Booster\\Booster_Models\\Impedance_Database\\Impedance\\Ceramic\\boosterII_ceramic_wz.txt'};
geom_k_factors = [-0.00564,-0.0126,0.0821];
geom_fit = struct([]);

CBREAK = {impedance(name, class, length, material, conductivity, thickness, apertures,....
      'RW_impedance_files',RW_impedance_files,'RW_wake_files',RW_wake_files,'RW_k_factors',RW_k_factors,...
      'Geom_impedance_files',geom_impedance_files,'Geom_wake_files',geom_wake_files,'Geom_k_factors',geom_k_factors,'Geom_fit',geom_fit)};

%% Cells
S1 = [DRIFT11; SEXT; DRIFT11];
S2 = [DRIFT11; SEXT; DRIFT11];
VACPUMP = [DRIFT11; VPCHAMBER];
VACPUMPREV = [VPCHAMBER; DRIFT11];
VACGATEVALVE = [VALVE1; DRIFT19];
CBVESSEL = [DRIFT19; CBREAK; DRIFT19];
CORRECTORLINE = [DRIFT16; VCM; DRIFT17; BPM1; DRIFT17; HCM; DRIFT16];
RFCAV = [DRIFT8; VALVE1; DRIFT8; MAINCAV; DRIFT8; VALVE1; DRIFT8];

GIRDER1 = [DRIFT18; S2; DRIFT19; DIPOLE1; CORRECTORLINE; DIPOLE1; VACPUMP; DRIFT15];
GIRDER13 = [GIRDER1; GIRDER1; GIRDER1];
GIRDER14 = [GIRDER1; GIRDER1; GIRDER1; GIRDER1];
GIRDER2 = [DRIFT18; VACGATEVALVE; DIPOLE1; CORRECTORLINE; DIPOLE1; VACPUMP; DRIFT15];
GIRDER3 = [DRIFT18; CBVESSEL; DIPOLE1; CORRECTORLINE; DIPOLE1; VACPUMP; DRIFT15];
GIRDER4 = [DRIFT18; S2; DRIFT19; DIPOLE1; VACPUMP; DRIFT15]; 
GIRDER51 = [DRIFT6; DRIFT12; BPM1; DRIFT21; VCM; DRIFT14; HCM; DRIFT13; S1; DRIFT6; QUAD; DRIFT6; S1; DRIFT6; S1; DRIFT13; HCM; DRIFT14; VCM; VACPUMP; DRIFT21];
GIRDER52 = [VACPUMP; DRIFT21; VCM; DRIFT14; HCM; DRIFT13; S1; DRIFT6; S1; DRIFT6; QUAD; DRIFT6; S1; DRIFT13; HCM; DRIFT14; VCM; DRIFT10; BPM1; DRIFT24; VACPUMP];
GIRDER6 = [VACPUMPREV; DIPOLE2; DRIFT6; BPM1; DRIFT22; QUAD];
GIRDER7 = [QUAD; DRIFT6; BPM1; DRIFT22; DIPOLE2; VACPUMP];
GIRDER8 = [BPM1; DRIFT22; QUAD; DRIFT2; QUAD; DRIFT22; VCM; DRIFT22; HCM; DRIFT22; QUAD; DRIFT18; VACPUMP];
GIRDER8REV = [VACPUMP; DRIFT18; QUAD; DRIFT22; HCM; DRIFT22; VCM; DRIFT22; QUAD; DRIFT2; QUAD; DRIFT22; BPM1];
LATTICEGIRDER1 = [GIRDER6; GIRDER52; GIRDER14; GIRDER2; GIRDER14; GIRDER3; GIRDER14; GIRDER2; GIRDER13; GIRDER4; GIRDER51; GIRDER7];

INJSTRAIGHT1 = [DRIFT7; SCREEN; DRIFT22; SCREEN; DRIFT10; GIRDER8REV; DRIFT2; RFCAV; DRIFT9];
INJSTRAIGHT2 = [DRIFT1; DRIFT1; DRIFT2; INJSEPT; DRIFT22; SCREEN; DRIFT22; VACPUMP; DRIFT1; DRIFT11; DRIFT1; INJSEPT; DRIFT2; GIRDER8; DRIFT4; DRIFT6; DRIFT10];
EXTSTRAIGHT = [DRIFT5; DRIFT11; INJSEPT; DRIFT2; INJSEPT; DRIFT3; DRIFT5; GIRDER8REV; DRIFT6; DRIFT11; DRIFT11; INJSEPT; DRIFT11; DRIFT11; RFCAV; DRIFT2; DRIFT22; SCREEN; DRIFT22; EXSEPT; DRIFT20; DRIFT2; DRIFT9; GIRDER8; DRIFT10; DRIFT6; DRIFT4];

LATTICE = [INJSTRAIGHT2; LATTICEGIRDER1; EXTSTRAIGHT; LATTICEGIRDER1; INJSTRAIGHT1];

IMPRING = LATTICE;

