clear all
clc

% q1

tau = [0.205 0.210 0.215 0.220 0.225 0.230 0.235 0.240 0.245 0.250];
h = [0.005 0.010 0.015 0.020, 0.025 0.030 0.035 0.040 0.045 0.050];

P = tf(1,[0.2 -1]);

for i = 1:10
    for j = 1:10
        P_unc = tf(1,[tau(i) -1],'InputDelay',h(j));
        [mag, phase] = bode(P_unc - P);
        bode(P_unc-P);
        hold on;
        disp(i);
    end
end

a1 = 0.05;
a2 = 0.07;

Wa = tf([a1 0],[a2*a2 2*a2 1]);

grid on;
bode(Wa,'b--');
title("Bode Plots of Plants and Error Bound Wa");
legend('Wa');
grid on
hold off

% q2

beta = 4.6:0.05:14.55;
isRobust = zeros(size(beta));
max_mag_vec = [];

for be = 1:length(beta)
    C = tf([15 1],[beta(be) 0]);
    robust_tf = (Wa*C/(1+P*C));
    [mag2,phase2] = bode(robust_tf);
    max_mag = max(mag2);
    max_mag_vec = [max_mag_vec max_mag];
end

figure;
plot(beta,max_mag_vec);
grid on;
title('Maximum of the Robust Transfer Functions for Each Fixed \beta');
xlabel('\beta')
ylabel('Infinite Norm of TF');

% q3

P1 = tf(1,[0.25 -1],'InputDelay',0.05);
beta_best = 10;
C_best = tf([15 1],[beta_best 0]);
tf_openloop = C*P1;
nyquist(tf_openloop)
grid on
margins = allmargin(tf_openloop);

% q4

Wr = tf(1,[1 0]);
gamma = 8.4;
beta_new = zeros(1,201);
beta_index = 1;
beta_cn = 0;
w = 0:0.01:500;

for n = 4.6:0.05:14.55
   Cnew = tf([15 1],[n 0]);
   Pnew = tf(1,[0.25 -1],'inputDelay',0.05);
   Snew = 1/(1+Pnew*Cnew);
   
   [mag_r1,phase] = bode(Wr*Snew/gamma,w);
   [mag_r2,phase] = bode(Wa*Cnew*Snew,w);
   
   mag_rmax = max(squeeze(mag_r1) + squeeze(mag_r2));
      
   if(mag_rmax <= 1)
      beta_new(beta_index) = n; 
      beta_cn = beta_cn + 1;
   end
   beta_index = beta_index + 1;
end

Cnew2 = tf([15 1],[8.40 0]);
Pnew2 = tf(1,[0.25 -1],'inputDelay',0.05);
Snew2 = 1/(1+Pnew2*Cnew2);

[mag_r12,phase] = bode(Wr*Snew2/gamma,w);
[mag_r22,phase] = bode(Wa*Cnew2*Snew2,w);

mag_r2 = squeeze(mag_r12) + squeeze(mag_r22);
figure;
plot(w,mag_r2);
grid on
title('Robust Performance');
xlabel('\omega');