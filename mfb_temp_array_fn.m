function [temp_array,linear_centre_freq]=mfb_temp_array_fn(fs,N)

    number_channels=40;                                                     % number of filter channels
    linear_freq=0:1:fs/2;                                                   % Range of the linear frequency
    mel_scale_freq=1127*log(1+linear_freq/700);                             % corresponding mel scale frequency
    band_width=max(mel_scale_freq)/(number_channels+1);                     % bandwidth of each channel
  
    freqz=0:max(mel_scale_freq)/(number_channels+1):max(mel_scale_freq);    % Band edge frequencies

    mel_Low_band_edges=freqz(1:number_channels);                            % Low band edge frequency vector (mel scale) 
    mel_upp_band_edges=freqz(3:number_channels+2);                          % high band edge frequency vector (mel scale)
    mel_centre_freq=freqz(2:number_channels+1);                             % centre frequencies (mel scale)
    
    linear_Low_band_edges=(exp(mel_Low_band_edges./1127)-1).*700;           % Low band edge frequency vector (linear scale)  
    linear_upp_band_edges=(exp(mel_upp_band_edges./1127)-1).*700;           % high band edge frequency vector (linear scale)
    linear_centre_freq=(exp(mel_centre_freq./1127)-1).*700;                 % centre frequencies (mel scale) (linear scale)
  
    % bin frequencies corresponding to the fft points
    linear_fft_bin_freq=0:fs/N:fs-1;                                        % OR freq_hz=0:fs/N:((fs/2)*(N-1)/N);
    
    linear_Low_band_bins=(linear_Low_band_edges.*N)./fs;
    linear_upp_band_bins=(linear_upp_band_edges.*N)./fs;
    % llbb=round(linear_Low_band_bins);
    % lubb=round(linear_upp_band_bins);
    
    temp_array=zeros(40,N/2);
    temp=1;
    
    for j=1:number_channels % covers the total no of channels
    count=1;                                 
    clear z;        

        for i=temp:N/2 % covers the fft range         
        bin_freq=linear_fft_bin_freq(i); % take the bin frequency corresponding to fft index 'i'
               
            % Group those values of FFT which falls within the channel band        
            if(bin_freq<=linear_upp_band_edges(j))
            
                % this loop is for the rising edge of the triangle 
                if(bin_freq<linear_centre_freq(j))
                    % calculate the slope value at bin freq index i
                    slope=1/(linear_centre_freq(j)-linear_Low_band_edges(j));
                    offset=(linear_Low_band_edges(j)/(linear_Low_band_edges(j)-linear_centre_freq(j)));
                    y=slope*bin_freq+offset;
                    z(i)=y;
                    % multiply the slope value and bin magnitude
                    mod_bin_mag=y;               
                    if j>1
                        temp_array(j,count+temp-1)=mod_bin_mag;
                    else
                        temp_array(j,count)=mod_bin_mag;    
                    end
                    count=count+1;
                else             
                    % this loop is for the falling edge of the triangle. Repeat the if block steps.             
                    slope=-1/(linear_upp_band_edges(j)-linear_centre_freq(j));
                    offset=(linear_upp_band_edges(j)/(linear_upp_band_edges(j)-linear_centre_freq(j)));
                    y=slope*bin_freq+offset;
                    z(i)=y;
                    mod_bin_mag=y;               
                    if j>1
                    temp_array(j,count+temp-1)=mod_bin_mag;
                    else
                    temp_array(j,count)=mod_bin_mag;    
                    end
                    count=count+1;
                end
          
            else
                
                break;
                
            end
            
        end
        
        temp=i;
        % this while loop is to bring the fft value array index back to the lower cut off of the next band as the bands are overlapping (two adjacent bands have common points).
        while linear_fft_bin_freq(temp)>=linear_centre_freq(j); %linear_Low_band_edges(j+1)
        temp=temp-1;
        end    
        temp=temp+1;
       % plot(z);              % to see the triangular filter banks
        % hold on;    
    end                         % end of computation of 40x512 temp_array matrix

end

    