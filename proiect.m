t = scope67(:,1);
u = scope67(:,2);
y1 = scope67(:,3);
y2 = scope67(:,4);

plot(t, [u, y1]); grid;
title("Datele primite:");
%%
% Pentru semnalul y1

y1_max=217;
y1_min=223;
u_max=215;
u_min=221;



Mr=(y1(y1_max)-y1(y1_min))/(u(u_max)-u(u_min));   % Mr=1.333;
tita=sqrt(2-sqrt(4-4/Mr^2))/2;                    % tita=0.4114;

k=mean(y1)/mean(u);                               % k=1.0119

i1=217;
i2=206;
Tr=t(i1)-t(i2);                                   % Tr =1.1e-04

wr = 2*pi/Tr;                                     % wr = 5.7120e+04
wn = wr/sqrt(1-2*tita^2);                         % wn=7.1056e+04

A=[0 1; -wn^2 -2*tita*wn];
B=[0; k*wn^2];
C=[1 0];
D=0;

H=ss(A,B,C,D);
ysim = lsim(A,B,C,D, u, t,[y1(1),(y1(2)-y1(1))/(t(2)-t(1))]);                

plot(t,[y1,ysim]);

xlabel('t');
ylabel( 'y(t), ysim(t)');
legend('y','ysim');


%eroarea medie patratica normalizata
eMPN=norm(y1-ysim)/norm(y1-mean(y1))
J=norm(y1-ysim)/sqrt(length(y1))

dy1=iddata(y1);
y1_simulat=iddata(ysim);
%compare(dy1, y1_simulat)


H1=tf([k*wn^2], [1 2*tita*wn wn^2]);
% faza la rezonanta



phr=rad2deg((t(215)-t(217))*wr);
%bode(H1)


%% Partea 2: Raspunsul in frecventa

%%frecvente joase 

%M1=(y1(ymax1)-y1(ymin1))/(u(umax1)-u(umin1));
%w1=pi/(t(ymin1)-t(ymax1));
%M2=(y1(ymin1)-y1(ymax2))/(u(umin1)-u(umax2));
t = scope67(:,1);
u = scope67(:,2);
y1 = scope67(:,3);
y2 = scope67(:,4);

plot(t, [u, y1]); grid;
title("Datele primite:");
%%
%mici
iu=[28, 49, 63, 75, 85, 143, 150, 226, 231, 336, 340];
iy=[33, 50, 64, 76, 87, 145, 152, 228, 233, 338, 342];

M1=(y1(33)-y1(50))/(u(28)-u(49));
M2=(y1(50)-y1(64))/(u(49)-u(63));
M3=(y1(64)-y1(76))/(u(63)-u(75));
M4=(y1(76)-y1(87))/(u(75)-u(85));


w1=pi/(t(50)-t(33));
w2=pi/(t(64)-t(50));
w3=pi/(t(76)-t(64));
w4=pi/(t(87)-t(76));

delta_t1=t(49)-t(50);   
delta_t2=t(49)-t(50);  
delta_t3=t(63)-t(64);  
delta_t4=t(75)-t(76);  

ph1=rad2deg(delta_t1*w1);
ph2=rad2deg(delta_t2*w2);
ph3=rad2deg(delta_t3*w3);
ph4=rad2deg(delta_t4*w4);

%%medii
M5=(y1(145)-y1(152))/(u(143)-u(150));
w5=pi/(t(152)-t(145));
delta_t5=t(143)-t(145);
ph5=rad2deg(delta_t5*w5);

M6=(y1(228)-y1(233))/(u(226)-u(231));
w6=pi/(t(233)-t(228));
delta_t6=t(226)-t(228);
ph6=rad2deg(delta_t6*w6);

%%mari
M7=(y1(338)-y1(342))/(u(336)-u(340));
w7=pi/(t(342)-t(338));
delta_t7=t(336)-t(338);
ph7=rad2deg(delta_t7*w7);


M8=(y1(706)-y1(709))/(u(704)-u(707));
delta_t=t(707)-t(709);
w8=pi/(t(709)-t(706));
ph8=rad2deg(delta_t*w8);


i1=990;
i2=992;
u1=988;
u2=990;

M9=(y1(i2)-y1(i1))/(u(u2)-u(u1));
w9=pi/(t(i2)-t(i1));
delta_t=t(u1)-t(i1);
ph9=rad2deg(delta_t*w9); 
%modul in jur de 1 si 0.5 dupa rezonanta, panta
%-40

H=H1;


w=logspace(3,6);
[num,den]=tfdata(H, 'v');
[M,Ph]=bode(num, den, w);

modul=[M1, M2, M3, M4, M5, Mr, M6, M7, M8, M9];
pulsatie=[w1, w2, w3, w4, w5, wr, w6, w7, w8, w9];
faza=[ph1, ph2, ph3, ph4, ph5, phr, ph6, ph7, ph8, ph9];

subplot(211);
semilogx(w, 20*log10(M), pulsatie, 20*log10(modul), "X");
grid on; ylabel('Amplitudine(dB)'); xlabel("Pulsatie(rad/s)");
title('Caracteristica de modul:'); hold on;

subplot(212);
semilogx(w, Ph); hold on;
semilogx(pulsatie, faza ,"X");
grid on; ylabel('Faza(grade)'); xlabel("Pulsatie(rad/s)")
title('Caracteristica de faza:');


%%modul in jur de 1 si 0.5 dupa rezonanta, panta -40




%% Partea 3: 2 modele pt fiecare validate prin intercorelatie, 2 prin auto, cate unul pt fiecare
%% Y1
%% ARX - validat prin intercorelatie
figure;
dt=t(2)-t(1);
dy1=iddata(y1,u,dt);
Marx=arx(dy1,[2,1,2]);

Hz_arx=tf(Marx)
compare(dy1,Marx)
figure;
resid(dy1, Marx);


%% ARMAX -trece amandoua, intercorelatia e mai buna
dt=t(2)-t(1);
dy1=iddata(y1,u,dt);
My1_armax=armax(dy1,[2,2,2,0])
Hz_armax=tf(My1_armax)
Hs=d2c(Hz_armax, 'zoh')
%%
figure;
compare(dy1,My1_armax)
figure;
resid(dy1, My1_armax);

%%
[num,den]=tfdata(Hz_armax, 'v')
[A, B, C, D]=tf2ss(num,den)
y1c=dlsim(A', C', B', D, u, [y1(1), y1(2)-A(1,1)*y1(1)-C(1,1)*u(1)])
plot(t, [y1, y1c])
eroare=norm(y1-y1c)/norm(y1-mean(y1))

%% OE-trece intercorelatia
dt=t(2)-t(1);
dy1=iddata(y1,u,dt);
M_y1_oe=oe(dy1, [2 2 0])

Hz_oe=tf(M_y1_oe)
Hs=d2c(Hz_oe, 'zoh')

figure;
compare(dy1,M_y1_oe)
figure;
resid(dy1, M_y1_oe, 5);
%%
yc=dlsim(Hz_oe, t, u);
plot(t, [y1, yc])
eroare=norm(y1-yc)/norm(y1-mean(y1))
%%
[num,den]=tfdata(Hz_oe, 'v')
[A, B, C, D]=tf2ss(num,den)
y1c=dlsim(A', C', B', D, u, [y1(1), y1(2)-A(1,1)*y1(1)-C(1,1)*u(1)])
plot(t, [y1, y1c])
eroare=norm(y1-y1c)/norm(y1-mean(y1))

%%

%% IV-trece inter
dt=t(2)-t(1);
dy1=iddata(y1,u,dt);
M_y1_iv=iv4(dy1,[2,2,1])

Hz_iv=tf(M_y1_iv)
Hs=d2c(Hz_iv, 'zoh')

figure;
compare(dy1,M_y1_iv)
figure;
resid(dy1, M_y1_iv, 5);

%%
[A, B, C, D]=tf2ss(M_y1_iv.B, M_y1_iv.A);
y1c=dlsim(A', C', B', D, u, [y1(1), y1(2)-A(1,1)*y1(1)-C(1,1)*u(1)])
plot(t, [y1, y1c]);
eroare=norm(y1-y1c)/norm(y1-mean(y1))

%% iv cu n4sid - amandoua
dt=t(2)-t(1);
dy1=iddata(y1,u,dt);
M_y1_iv=iv4(dy1,[2,1,0]);
My1_n4sid=n4sid(dy1, M_y1_iv)

figure;
resid(dy1,My1_n4sid);

figure;
compare(dy1,My1_n4sid);
%%
A=My1_n4sid.A;
B=My1_n4sid.B;
C=My1_n4sid.C;
D=My1_n4sid.D;
y1c=dlsim(A', C', B', D, u, [y1(1), y1(2)-A(1,1)*y1(1)-C(1,1)*u(1)])
plot(t, [y1, y1c]);
eroare=norm(y1-y1c)/norm(y1-mean(y1))
%% Y2
%% ARX - trece intercorelatia
figure;
dt=t(2)-t(1);
dy2=iddata(y2,u,dt);
M_y2_arx=arx(dy2,[2,2,0]);

num=M_y2_arx.B;
den=M_y2_arx.A;

Hz_arx=tf(num,den,dt)

compare(dy2,M_y2_arx)
figure;
resid(dy2, M_y2_arx, 7);


%% ARMAX - trece autocorelatia

dt=t(2)-t(1);
dy2=iddata(y2,u,dt);
M_y2_armax=armax(dy2,[2,2,2,0])

Hz_armax=tf(M_y2_armax)
figure
compare(dy2,M_y2_armax)
figure;
resid(dy2, M_y2_armax, 5);
%%
[num,den]=tfdata(Hz_armax, 'v')
%%
[A, B, C, D]=tf2ss(num,den)
y2c=dlsim(A', C', B', D, u, [y2(1), y2(2)-A(1,1)*y2(1)-C(1,1)*u(1)])

plot(t, [y2, y2c])
norm(y2-y2c)/norm(y2-mean(y2))
%% OE - niciuna
figure;
dt=t(2)-t(1);
dy2=iddata(y2,u,dt);
M_y2_oe=oe(dy2,[2,2,1])

num=M_y2_oe.B;
den=M_y2_oe.F;

Hz_oe=tf(num,den)

compare(dy2,M_y2_oe)
figure;
resid(dy2, M_y2_oe,7);

%% IV - intercorelatia

dt=t(2)-t(1);
dy2=iddata(y2,u,dt);
M_y2_iv=iv4(dy2,[2,2,0])

figure; 
compare(dy2,M_y2_iv)
figure;
resid(dy2, M_y2_iv, 7);


%%
compare(dy2,M)
figure;
resid(dy2, M, 5);
A=M.A;
B=M.B;
C=M.C;
D=M.D;
[num,den]=ss2tf(A,B,C,D);
Hz=tf(num,den,dt, "Variable", "z^-1")
%% IV cu n4sid trece intercorelatia sau chiar amandoua
figure;
dt=t(2)-t(1);
dy2=iddata(y2,u,dt);
M_y2_iv=iv4(dy2,[2,2,0])

%%
M_y2_iv_n4sid=n4sid(dy2,M_y2_iv)


Hz=tf(M_y2_iv_n4sid);
Hs=d2c(Hz, 'zoh', dt);

compare(dy2,M_y2_iv_n4sid)
figure;
resid(dy2, M_y2_iv_n4sid,'corr',7);

%%

[num,den]=tfdata(Hz, 'v')
[A, B, C, D]=tf2ss(num,den)
y2c=dlsim(A', C', B', D, u, [y2(1), y2(2)-A(1,1)*y2(1)-C(1,1)*u(1)])

plot(t, [y2, y2c])
norm(y2-y2c)/norm(y2-mean(y2))


