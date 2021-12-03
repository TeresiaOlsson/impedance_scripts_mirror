function save_elegant_impedance_file(elegant_impedance)


    fileID = fopen('elegant_real_lon_impedance.txt','w');
    for i = 1:length(elegant_impedance.ImpedanceFreq)
        fprintf(fileID,'%10e %10e\n',elegant_impedance.ImpedanceFreq(i),elegant_impedance.ImpedanceRealZ(i));
    end
    fclose(fileID);
    system('plaindata2sdds elegant_real_lon_impedance.txt impedance_longitudinal_real_elegant.sdds -inputMode=ascii -outputMode=ascii -separator=" " -noRowCount -column=Freq,double,units=Hz  -column=RealZ,double,units=Ohm');
    delete('elegant_real_lon_impedance.txt');

    fileID = fopen('elegant_imag_lon_impedance.txt','w');
    for i = 1:length(elegant_impedance.ImpedanceFreq)
        fprintf(fileID,'%10e %10e\n',elegant_impedance.ImpedanceFreq(i),elegant_impedance.ImpedanceImagZ(i));
    end
    fclose(fileID);
    system('plaindata2sdds elegant_imag_lon_impedance.txt impedance_longitudinal_imag_elegant.sdds -inputMode=ascii -outputMode=ascii -separator=" " -noRowCount -column=Freq,double,units=Hz  -column=ImagZ,double,units=Ohm');
    delete('elegant_imag_lon_impedance.txt');

    fileID = fopen('elegant_transverse_impedance.txt','w');
    for i = 1:length(elegant_impedance.ImpedanceFreq)
        fprintf(fileID,'%10e %10e %10e %10e %10e\n',elegant_impedance.ImpedanceFreq(i),elegant_impedance.ImpedanceRealX(i),elegant_impedance.ImpedanceImagX(i),elegant_impedance.ImpedanceRealY(i),elegant_impedance.ImpedanceImagY(i));
    end
    fclose(fileID);
    system('plaindata2sdds elegant_transverse_impedance.txt impedance_transverse_elegant.sdds -inputMode=ascii -outputMode=ascii -separator=" " -noRowCount -column=Freq,double,units=Hz  -column=RealX,double,units=Ohm/m -column=ImagX,double,units=Ohm/m -column=RealY,double,units=Ohm/m -column=ImagY,double,units=Ohm/m');
    delete('elegant_transverse_impedance.txt');

end