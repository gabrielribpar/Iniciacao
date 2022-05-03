clear all;
load('choose_OBE-2_waveforms_masks')
warning('off','all')
alfa_pulmonar=zeros(32);
beta_pulmonar=zeros(32);
k_pul=zeros(32);
t_0_pul=zeros(32);
t_pico=zeros(32);
[t_slope_min,t_pico_min,m_min,n_min]=Pixel_Precoce(img_waveforms,mask,mask_heart);
tamanho_N=size(img_waveforms);
N=tamanho_N(3);
fs=50;
tend=(N-1)/fs;
t = 0:1/fs:tend;
max_slope=zeros(32);
Pico=zeros(32);
for m=1:32
    for l=1:32
        if mask(m,l)==1
            y=-1*img_waveforms(m,l,:);
            for i=1:N
                gamma_1(i)=y(1,1,i);
            end
            if mask_heart(m,l)==0
                [gamma_pulmonar_otm,alfa_pulmonar(m,l),beta_pulmonar(m,l),k_pul(m,l),t_0_pul(m,l)]=Puramente_Pulmonar(gamma_1,t,N);
                gamma=gamma_norm(t,alfa_pulmonar(m,l),beta_pulmonar(m,l),k_pul(m,l),t_0_pul(m,l));
                Pico(m,l)=max(gamma);
                max_slope(m,l)=max(diff(gamma)*fs);
                t_pico(m,l)=beta_pulmonar(m,l)*(alfa_pulmonar(m,l)-1)+t_0_pul(m,l);
            end
        end
        
    end
end

for m=32:-1:1
    for l=32:-1:1
        if mask(m,l)==1
            y=-1*img_waveforms(m,l,:);
            for i=1:N
                gamma_1(i)=y(1,1,i);
            end
            if mask_heart(m,l)==1
                alfa_max=max([alfa_pulmonar(m+1,l) alfa_pulmonar(m+1,l+1) alfa_pulmonar(m,l+1)]);
                alfa_min=min([alfa_pulmonar(m+1,l) alfa_pulmonar(m+1,l+1) alfa_pulmonar(m,l+1)]);
                beta_max=max([beta_pulmonar(m+1,l) beta_pulmonar(m+1,l+1) beta_pulmonar(m,l+1)]);
                beta_min=min([beta_pulmonar(m+1,l) beta_pulmonar(m+1,l+1) beta_pulmonar(m,l+1)]);
                k_max=max([k_pul(m+1,l) k_pul(m+1,l+1) k_pul(m,l+1)]);
                k_min=min([k_pul(m+1,l) k_pul(m+1,l+1) k_pul(m,l+1)]);
                t_0_max=max([t_0_pul(m+1,l) t_0_pul(m+1,l+1) t_0_pul(m,l+1)]);
                t_0_min=min([t_0_pul(m+1,l) t_0_pul(m+1,l+1) t_0_pul(m,l+1)]);
                [alfa_pulmonar(m,l),beta_pulmonar(m,l),k_pul(m,l),t_0_pul(m,l),gamma_final(m,l,:),gamma_cardiaca(m,l,:),gamma_pulmonar(m,l,:)]=Pixel_Misto(gamma_1,t,N,alfa_max,alfa_min,beta_max,beta_min,k_max,k_min,t_0_max,t_0_min,t_slope_min,t_pico_min,m_min,n_min);
                gamma=gamma_norm(t,alfa_pulmonar(m,l),beta_pulmonar(m,l),k_pul(m,l),t_0_pul(m,l));
                Pico(m,l)=max(gamma);
                max_slope(m,l)=max(diff(gamma)*fs);
                t_pico(m,l)=beta_pulmonar(m,l)*(alfa_pulmonar(m,l)-1)+t_0_pul(m,l);
            end
        end
        
    end
end
hold off;
heatmap(max_slope)
colormap(hot)
grid off;
direita_slope=100*sum(sum(max_slope(:,1:16)))/sum(sum(max_slope));
esquerda_slope=100-direita_slope;


y=gamma_final(11,14,:);
y1=gamma_cardiaca(11,14,:);
y2=gamma_pulmonar(11,14,:);
for i=1:N
    gamma_final_1(i)=y(1,1,i);
    gamma_cardiaca_1(i)=y1(1,1,i);
    gamma_pulmonar_1(i)=y2(1,1,i);
end

figure;

y=-1*img_waveforms(11,14,:);
for i=1:N
       gamma_1(i)=y(1,1,i);
end
hold on;
plot(t,gamma_final_1,'r');
plot(t,gamma_1,'b-');
plot(t,gamma_cardiaca_1);
plot(t,gamma_pulmonar_1);
