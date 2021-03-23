
%% Laboratory 5 - Control Engineering I
%% CONTROLLER DESIGN IN FREQUENCY DOMAIN BASED ON SECOND ORDER SYSTEM
%% GOALS

% To follow and to understand the design method steps 

% To check the resulted performance indicators

%% Problem
% For the process described by Hf(s) = 3.5 / s(0.5*s + 1) and the
% performance indicators [..] following the above described steps, design a
% P and a PI controller. Simulate the step and ramp output of the closed
% loop to highlight the performance indicators.



%% Exercise 1

% P Controller
%corner frequency = 2 rad/s
Hf = tf(3.5, [0.5 1 0])
bodemag(Hf), hold on

wf = 1/0.5 %=2 1/Tf
sigma = 0.15
zeta = abs(log(sigma)) / sqrt (pi^2 + log(sigma)^2)

A = 1/4/sqrt(2)/zeta^2


% A in decibels
AdB = 20 *log10(A)

A = tf(A, 1);
bodemag(A)

FN = -AdB + 1.84 % 1.84 from Hf(blue) at freq 2
Vr = 10^(-FN/20)

Hd = Hf*Vr
bode(Hd)


wt = 1.5 %at magnitude 0 dB de pe Hd
wn = 2 * zeta * wt;

ts = 4/zeta/wn

% from Hd at freq. 1 (Hd galben)
cvdB = 4.39; cv = 10^(4.39/20)  % >=1

deltawb = wt % 1.5 < 15

H0 = feedback(Hf*Vr, 1)

figure, step(H0)

t = 0:0.1:20;
figure, lsim(H0, t, t)

% Screenshot with values for each graphic

% Bode Diagram
 figure;
 image_bode = imread('lab5ce_ex1_bode.jpg');
 imshow(image_bode);
 
 % Step Response
 figure;
 image_step = imread('lab5ce_ex1_step.jpg');
 imshow(image_step);
 
% Linear Simulation Results
 figure;
 image_lsim = imread('lab5ce_ex1_linear.jpg');
 imshow(image_lsim);

%% Exercise 2
% Pi Controller

sigma = 0.07

Hf = tf(3.5, [0.5 1 0])
bodemag(Hf), hold on

wf = 1/0.5 %=2
zeta = abs(log(sigma)) / sqrt (pi^2 + (log(sigma)^2))

A = 1/4/sqrt(2)/zeta^2

% A in decibels
AdB = 20 *log10(A)

A = tf(A, 1);
bodemag(A)

%1.94 from Hf(blue) at freq. 2
FN = -AdB + 1.94
Vr = 10^(-FN/20)

Hd = Hf*Vr
bode(Hd)

%at frequency 1

wt = 1.05 % de pe Hd galben at magnitude 0 dB
wn = 2 * zeta * wt;

ts = 4/zeta/wn

% from Hd at freq. 1 (Hd galben)
cvdB = 0.573; cv = 10^(cvdB/20)  % >=1

deltawb = wt % 1.05 < 15

wz = 0.1 * wt;
cvstar = 5;
wp = cv/cvstar *wz;
Tz = 1/wz; Tp = 1/wp;

VrPI = cvstar / cv

 HPI = VrPI*tf([Tz 1], [Tp 1]);
 cv_new = VrPI * 3.5 * Vr;
 Hd_PI = Hf * Vr * HPI;
 bode(Hd_PI)
 H0 = feedback(Vr*HPI*Hf, 1)
 
 figure, step(H0) 
 t= 0:0.1:80;
 figure, lsim(H0, t, t)

 essv = 80 - 79.8;
 cv = 1/essv
 
 
 %Screenshot with values for each graphic

 % Bode Diagram
 figure;
  image_bode2 = imread('lab5ce_ex2_bode.jpg');
  imshow(image_bode2);
  
  % Step Response
  figure;
  image_step2 = imread('lab5ce_ex2_step.jpg');
  imshow(image_step2);
  
 % Linear Simulation Results
  figure;
  image_lsim2 = imread('lab5ce_ex2_linear.jpg');
  imshow(image_lsim2);
