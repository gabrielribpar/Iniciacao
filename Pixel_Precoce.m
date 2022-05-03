
function [t_slope_min,t_pico_min,m_min,n_min]=Pixel_Precoce(img_waveforms,mask,mask_heart)
tamanho_N=size(img_waveforms);
N=tamanho_N(3);
fs=50;
tend=(N-1)/fs;
t = 0:1/fs:tend;
m_min=1;
n_min=1;
t_pico_min=0;
t_slope=0;
t_pico=0;
t_slope_min=0;
for m=1:32
    for n=1:32
        if mask(m,n)==1 && mask_heart(m,n)==1
            y=img_waveforms(m,n,:);
            for i=1:N
                gamma_1(i)=y(1,1,i);
            end
            gamma_2=diff(gamma_1)*fs;
            gamma_3=diff(gamma_2)*fs;
            pol=polyfit(t,gamma_1,13);
            pol_1=polyfit(t(1:N-1),gamma_2,13);
            pol_2=polyfit(t(1:N-2),gamma_3,13);
            tpp=roots(pol_1);
            tms=roots(pol_2);
            for(i=1:length(tpp))
                if(tpp(i)<tend && tpp(i)>2 && imag(tpp(i))==0)
                    t_pico=tpp(i);
                end
            end
            for(i=1:length(tms))
                if(tms(i)<3 && tms(i)>1 && imag(tms(i))==0)
                    t_slope=tms(i);
                end
            end
            if t_slope_min==0
                m_min=m;
                n_min=n;
                t_slope_min=t_slope;
            elseif t_slope<t_slope_min
                m_min=m;
                n_min=n;
                t_slope_min=t_slope;
                
            end
            if t_pico_min==0
                t_pico_min=t_pico;
            elseif t_pico<t_pico_min
                t_pico_min=t_pico;
                
            end
        end
    end
end
t_pico_min=0.9*t_pico_min;
t_slope_min=0.95*t_slope_min;
end


