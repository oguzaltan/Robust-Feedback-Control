clear all;
clc;

N = tf([20 -30 10],[1 11 11 10]);
D = tf([1 -1 1],[1 1 1]);
P = minreal(N/D);

Ni = tf([1 -1.5 0.5],[1 1.5 0.5]);
No = tf([20 30 10],[1 11 11 10]);
Di = D;
Do = 1;

r = roots([1 -1 1]);
a1 = 1/evalfr(N,r(1));
a2 = 1/evalfr(N,r(2));
X = zpk([-0.3478],[-4],[-4.6]);
display("X is stable: " + isstable(X));

Y = zpk([-107.5 -0.521],[-10 -4],[1]);
display("Y is stable: " + isstable(Y));

W1 = tf([0.1 1],[1 0.1]);
W2 = tf([0.1 0.01],[1]);
W1con = tf([0.1 -1],[1 -0.1]);
W2con = tf([-0.1 0.01],[1]);

Gspec = W1*W1con + W2*W2con;
[Gspecnum, Gspecden] = tfdata(Gspec, 'v');
rnum_G = roots(Gspecnum);
rden_G = roots(Gspecden);

left_poles = real(rnum_G) < 0;
left_zeros = real(rden_G) < 0;
G = zpk(rnum_G(left_poles),rden_G(left_zeros),(1));
[~,~,kG] = zpkdata(Gspec);
G = sqrt(-kG)*G;
display("G is stable: " + isstable(G));

V = minreal(W1*W2/G);
Vcon = tf([0.01 -0.0101 0.01],[1 -4.585 10]);

Nicon = tf([1 1.5 0.5],[1 -1.5 0.5]);
Dicon = tf([1 1 1],[1 -1 1]);
W1con = tf([-0.1 1],[-1 0.1]);

Gcon = tf([1 -4.585 10],[-1 0.1]);

%Step 0
gamma_max = 100;
[mag,~] = bode(V);
gamma_min = max(mag);

while(gamma_max - gamma_min > 0.00001)
    
    %Step 1
    gamma = (gamma_min + gamma_max)/2;
    Vgamma_spectral = gamma^2 - Vcon*V;
    [Vgamma_spectral_num, Vgamma_spectral_den] = tfdata(Vgamma_spectral, 'v');
    rnum_vgamma = roots(Vgamma_spectral_num);
    rden_vgamma = roots(Vgamma_spectral_den);
    
    left_poles = real(rnum_vgamma) < 0;
    left_zeros = real(rden_vgamma) < 0;
    
    [~,~,kgamma] = zpkdata(Vgamma_spectral);
    
    Vgamma = zpk(rnum_vgamma(left_poles),rden_vgamma(left_zeros),(sqrt(kgamma)));
    display("Vgamma is stable: " + isstable(Vgamma));
    display("1/Vgamma is stable: " + isstable(1/Vgamma));
    
    %Step 2
    Rgamma = minreal((1/Vgamma)*(Nicon*Dicon*W1con/Gcon*W1-Dicon*G*No*X));
    
    %Nehari
    [Rs,Ru] = stabsep(minreal(Rgamma));
    
    [A,B,C,~] = ssdata(tf(Ru));
    
    Wc = lyap(A,-B*B');
    Wo = lyap(A',-C'*C);
    
    gamma1 = sqrt(max(eig(Wc*Wo)));
    [V_eig,D_eig] = eig(Wc*Wo);
    xmax = V_eig(:,1);
    ymax = (1/gamma1)*Wo*xmax;
    
    if (gamma1 < 1 || gamma1 == 1)
        gamma_max = gamma;
    else
        gamma_min = gamma;
    end
end

%Step 1
gamma = (gamma_min + gamma_max)/2;
Vgamma_spectral = gamma^2 - Vcon*V;
[Vgamma_spectral_num, Vgamma_spectral_den] = tfdata(Vgamma_spectral, 'v');
rnum_vgamma = roots(Vgamma_spectral_num);
rden_vgamma = roots(Vgamma_spectral_den);

left_poles = real(rnum_vgamma) < 0;
left_zeros = real(rden_vgamma) < 0;

[~,~,kgamma] = zpkdata(Vgamma_spectral);

Vgamma = zpk(rnum_vgamma(left_poles),rden_vgamma(left_zeros),(sqrt(kgamma)));
display("Vgamma is stable: " + isstable(Vgamma));
display("1/Vgamma is stable: " + isstable(1/Vgamma));

%Step 2
Rgamma = minreal((1/Vgamma)*(Nicon*Dicon*W1con/Gcon*W1-Dicon*G*No*X));

%Nehari
[Rs,Ru] = stabsep(minreal(Rgamma));

[A,B,C,~] = ssdata(tf(Ru));

Wc = lyap(A,-B*B');
Wo = lyap(A',-C'*C);

gamma1 = sqrt(max(eig(Wc*Wo)));
[V_eig,D_eig] = eig(Wc*Wo);
xmax = V_eig(:,1);
ymax = (1/gamma1)*Wo*xmax;

if (gamma1 < 1 || gamma1 == 1)
    gamma_max = gamma;
else
    gamma_min = gamma;
end

%Step 3
gamma_opt = gamma;

Q1opt = Rgamma - gamma1*(ss(A,xmax,C,0)/ss(-A',ymax,B',0));
Q1opt = zpk(minreal(Q1opt));
display("Q1opt is stable: " + isstable(Q1opt));

%Step 4
% Qopt = (Vgamma*Q1opt);
Qopt = minreal(Vgamma*Q1opt);
display("Qopt is stable: " + isstable(Qopt));

Qc = minreal(Qopt/(No*Do*G));
display("Qc is stable: " + isstable(Qc));

Copt = minreal((X+D*Qc)/(Y-N*Qc));

openlooppolezero = P*Copt;
openlooptf = minreal(P*Copt);
Sopt = minreal(1/(1+openlooptf));
display("Sopt is stable: " + isstable(Sopt));

% nyquist(openlooptf);

w = 0.001:0.01:1000;
[mag1,~] = bode(W1*Sopt,w);
[mag2,~] = bode(W2*(1-Sopt),w);
mag1 = squeeze(mag1);
mag2 = squeeze(mag2);

gamma_opt_2 = sqrt(mag1.^2 + mag2.^2);
gamma_opt_2_t = (gamma_opt_2-gamma_opt)';
figure;
plot(w,gamma_opt_2_t);
grid on;
xlabel("\omega");
ylabel("Difference");
