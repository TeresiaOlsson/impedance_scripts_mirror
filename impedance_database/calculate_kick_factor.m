function kick_factor = calculate_kick_factor(type,input,sigma_t)
% Wake should be wake function convoluted by sigma_t (e.g. the wake
% potential) according to AT conventions (W initially negative as t > 0)
% Impedance should be according to Elegant conventions

switch type
    
    case 'wake'
        s = input(:,1);
        wake = input(:,2);
                
        c = 299792458;        
        sigma_s = c.*sigma_t;              
        lambda = 1./(sqrt(2.*pi).*sigma_s).*exp(-s.^2./(2.*sigma_s.^2));
        
        y = -lambda.*wake;        
        kick_factor = trapz(s,y);
        
    case 'impedance'
        
        freq = input(:,1);
        real = input(:,2);
        imag = input(:,3);

        omega = 2.*pi.*freq;

        % For transverse impedance imaginary part cancel using Elegant's
        % impedance convention
        y = real.*h(omega,sigma_t);
        kick_factor = 1./(2.*pi).*trapz(omega,y).*2;
end

end

function h = h(x,sigma_t)
    h = exp(-x.^2*sigma_t.^2);
end
