
function [gamma_pulmonar_otm,alfa,beta,k,t0_pul]=Puramente_Pulmonar(gamma_1,t,N)
options = optimoptions(@lsqcurvefit,'Display','off');
fs=50;
tend=(N-1)/fs;
gamma_2=diff(gamma_1)*fs;
gamma_3=diff(gamma_2)*fs;
pol=polyfit(t,gamma_1,13);
pol_1=polyfit(t(1:N-1),gamma_2,13);
pol_2=polyfit(t(1:N-2),gamma_3,13);
tpp=roots(pol_1);
tms=roots(pol_2);
t0_pul=0;
for i=1:length(tpp)
    if(tpp(i)<tend && tpp(i)<2 && imag(tpp(i))==0 && tpp(i)>=0)
        t0_pul=tpp(i);
    end
end

for i=1:length(tpp)
    if(tpp(i)<tend && tpp(i)>t0_pul && tpp(i)>3 && imag(tpp(i))==0)
        t_pico=tpp(i);
    end
end

for i=1:length(tms)
    if(tms(i)<t_pico && tms(i)>t0_pul && imag(tms(i))==0)
        t_slope=tms(i);
    end
end
t_pico=t_pico-t0_pul;
t_slope=t_slope-t0_pul;
alfa=1+(t_pico)^2/(t_pico-t_slope)^2;
beta=((t_pico-t_slope)^2)/t_pico;
pico=exp(1-alfa)*(t_pico)^(alfa-1)/(gamma(alfa)*beta^(alfa));
k=polyval(pol,t_pico)/pico;
x1(1)=alfa;
x1(2)=beta;
x1(3)=k;
x1(4)=t0_pul;
lb = [0,0,0,0];
ub=[10,5,2,3];
F=@(x,xdata)x(3)*gampdf(xdata-x(4),x(1),x(2));
[x] = lsqcurvefit(F,x1,t,gamma_1,lb,ub,options);
gamma_pulmonar_otm=gamma_norm(t,x(1),x(2),x(3),x(4));
alfa=x(1);
beta=x(2);
t0_pul=x(4);
k=x(3);
end



