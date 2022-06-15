clear; clc;

%Problem 1 part a
P = tf([4,-8],[1,-2,2]);
X = tf([-2.5,3.75],[1,1]);
N = tf([4,-8],[1,3,2]);
D = tf([1,-2,2],[1,3,2]);
Y = minreal((1-N*X)/D);

%Problem 1 part b
Qc = tf([-10.8077,13.3846,-36],[1,6,9]);
C = (X+D*Qc)/(Y-N*Qc);

%Problem 2 part a
Wm = tf([1 1],1);
M = tf([1 -2 2],[1 2 2]);
N_out = 4*tf([1,2],([1,3,2]));
W = Wm*N_out*X;
W = minreal(W);
b = [evalfr(W,1+i) evalfr(W,1-i)];
a = [1+i 1-i];
[gopt,Fopt] = NevPickNew(a,b);
Qopt = minreal((-Fopt+W)/M);
deltamax = 1/gopt;
Qcopt = minreal(-(W - Fopt)/(Wm*N_out*D));
Copt = minreal((X+D*Qcopt)/(Y-N*Qcopt));

%Problem 2 part b
delta = 0.9*deltamax;
P_unc = P + P*Wm*delta;
tf_openloop = minreal(Copt*P_unc);
nyquist(tf_openloop);
pole(minreal(tf_openloop/(1+tf_openloop)))
margins = allmargin(tf_openloop);