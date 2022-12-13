clear all;

q0min = 0.95;
q0max = 1;
q1min = 0.35;
q1max = 0.4;
r0min = 0.9;
r0max = 1;
r1min = 9;
r1max = 12;
r2min = 50;
r2max = 100;
r3min = 120;
r3max = 150;
r4min = 195;
r4max = 200;
r5min = 120;
r5max = 130;

%part a

K_interval = zeros(1,100001); %K candidates

del_max = 3.151626793; %experimental finding

K_elements = 0;
inter_index = 1;

for K = 0:0.001:100
    
    r6min = - del_max;
    r6max = + del_max;
    
    Nc = [-4 K]; 
    Dc = [1 0];

    N1 = [q0min q1min]; %Kharitonov numerator polynomials 
    N2 = [q0max q1max];
    N3 = [q0max q1min];
    N4 = [q0min q1max];

    D1 = [r0min r1min r2max r3max r4min r5min r6max]; %Kharitonov denominator polynomials 
    D2 = [r0max r1max r2min r3min r4max r5max r6min];
    D3 = [r0min r1max r2max r3min r4min r5max r6max];
    D4 = [r0max r1min r2min r3max r4max r5min r6min];
      
    XN_1 = zeros(1,8);
    XN_2 = zeros(1,8);
    XN_3 = zeros(1,8);
    XN_4 = zeros(1,8);
    
    XN_1(6:8) = conv(Nc,N1);
    XN_2(6:8) = conv(Nc,N2);
    XN_3(6:8) = conv(Nc,N3);
    XN_4(6:8) = conv(Nc,N4);

    XN_matrix = [XN_1;XN_2;XN_3;XN_4];
    
    XD_1 = conv(Dc,D1);
    XD_2 = conv(Dc,D2);
    XD_3 = conv(Dc,D3);
    XD_4 = conv(Dc,D4);
 
    X_matrix = (zeros(16,8));
    
    for i = 1:4
        X_matrix(4*i-3,:) = XN_matrix(ceil((4*i-3)/4),:) + XD_1;
        X_matrix(4*i-2,:) = XN_matrix(ceil((4*i-3)/4),:) + XD_2;
        X_matrix(4*i-1,:) = XN_matrix(ceil((4*i-3)/4),:) + XD_3;
        X_matrix(4*i,:) = XN_matrix(ceil((4*i-3)/4),:) + XD_4;
    end

    r_real_max_vec = zeros(16,1);

    %find the max real root of polynomials
    
    for i = 1:16
        r = roots(X_matrix(i,:));
        r_real_max = max(real(r));
        r_real_max_vec(i) = r_real_max;
    end
    
    %find the max K so that all 16 polynomials are stable and
    %find the interval of all Ks that the system is stable
    
    if(all(r_real_max_vec(:,1) < 0))
        K_interval(inter_index) = K; 
        K_elements = K_elements + 1;
        K_max = K;
    end    
    inter_index = inter_index + 1;
end

%part b

NF = [1 0.35];
DF = [1 10 60 140 200 121 ((-del_max / 2)-1.4) 0];
F = tf(NF,DF);

roots_K = rlocus(F,K_max);

figure
rlocus(F);
grid on
hold on
plot(roots_K,'r*');
hold off

NC_new = [-4 K_max];
DC_new = [1 0];

N_P0 = [1 0.35];
D_P0 = [1 10 60 140 200 125 (-del_max/2)];

Cc = tf(NC_new,DC_new);
Pc = tf(N_P0,D_P0);

%part c

Margins = allmargin(Cc*Pc);

Sensitivity = 1/(1 + Cc*Pc);
[mag,phase,w] = bode(Sensitivity);

disp(max(mag));

figure
plot(w,squeeze(mag));
title("Sensitivity Function");
grid on
xlabel("\omega");
ylabel("S(s)");