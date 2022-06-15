clear all;
clc;

N = tf(1,[1 1]);
D = tf([1 -3 2],[1 3 2]);
P = minreal(N/D);

Nc = tf(1,[-1 1]);
Dc = tf([1 3 2],[1 -3 2]);

rhs = N*tf(1,[-1 1])+D*tf([1 3 2],[1 -3 2]);

G = tf([1 sqrt(2)],[1 1]);

Nn = minreal(N/G);
Dn = minreal(D/G);

% X(s) = (x1s+x2)/(s+1)

x1 = 1/evalfr(Nn,1);
x2 = 1/evalfr(Nn,2);

Xn = tf([5.414 -0.5858],[1 1]); % we find Xn(s)

% Yn = minreal((1 - Xn*Nn)/Dn); % we find Yn(s)
Yn = tf([1 2],[1 1]);

% Problem 2

Nnc = tf(1,[-1 1.414]);
Dnc = tf([1 3 2],[1 -3.414 2.828]);

R = Nnc*Yn-Dnc*Xn;
[Rs,Ru] = stabsep(minreal(R));

[A,B,C,~] = ssdata(tf(Ru));

Wc = lyap(A,-B*B');
Wo = lyap(A',-C'*C);

gopt = sqrt(max(eig(Wc*Wo)));
bmax = 1/gopt;
[V,D] = eig(Wc*Wo);
xmax = V(:,1);
ymax = (1/gopt)*Wo*xmax;

Qopt = R - gopt*(ss(A,xmax,C,0)/ss(-A',ymax,B',0)); % Correct implementation

Qopt = zpk(minreal(Qopt));

Qcopt = Qopt;
Copt = (Xn + Dn*Qcopt)/(Yn - Nn*Qcopt);
Copt = minreal(Copt);

[Qcnum, Qcden] = tfdata(Qcopt, 'v');
isstable(Qcopt)

r = roots(Qcden);

%% different Xn and Yn

clear all;
clc;

N = tf(1,[1 1]);
D = tf([1 -3 2],[1 3 2]);
P = minreal(N/D);

Nc = tf(1,[-1 1]);
Dc = tf([1 3 2],[1 -3 2]);

G = tf([1 sqrt(2)],[1 1]);

Nn = minreal(N/G);
Dn = minreal(D/G);

% X(s) = (x1s+x2)/(s+2)

x1 = 1/evalfr(Nn,1);
x2 = 1/evalfr(Nn,2);

Xn = tf([6.4143 0.8283],[1 2]); % we find Xn(s)

% Yn = minreal((1 - Xn*Nn)/Dn); % we find Yn(s)
Yn = tf([1],[1]); % we find Yn(s)

Nnc = tf(1,[-1 sqrt(2)])
Dnc = tf([1 3 2],[1 -(2 + sqrt(2)) sqrt(8)]);

Nnc = tf(1,[-1 sqrt(2)])
Dnc = tf([1 3 2],[1 -(2 + sqrt(2)) sqrt(8)]);

Nnc = tf(1,[-1 1.414]);
Dnc = tf([1 3 2],[1 -3.414 2.828]);

R = Nnc*Yn-Dnc*Xn;
[Rs,Ru] = stabsep(minreal(R));

[A,B,C,~] = ssdata(tf(Ru));

Wc = lyap(A,-B*B');
Wo = lyap(A',-C'*C); 

gopt = sqrt(max(eig(Wc*Wo)));
bmax = 1/gopt;
[V,D] = eig(Wc*Wo);
xmax = V(:,1);
ymax = (1/gopt)*Wo*xmax;
bmax = 1/gopt;

Qopt = R - gopt*(ss(A,xmax,C,0)/ss(-A',ymax,B',0));

Qopt = zpk(minreal(Qopt));

Qcopt = Qopt;
Copt = (Xn + Dn*Qcopt)/(Yn - Nn*Qcopt);
Copt = minreal(Copt);

[Qcnum, Qcden] = tfdata(Qcopt, 'v');
isstable(Qcopt)

r = roots(Qcden);
