clc
close all

%Created by Steven Stemmley 
%April 2016

%This takes the IV trace from the oscilloscope and reads out parameters
%based upon the trace.
%Currently only valid for a collisionless plasma.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%resistor value
R = 11e3; %Ohms

%Type of Gas
gasmass = 40; %amu


%Gas Pressure
gasp = 380; %mTorr

%100 Watts
% Current = xlsread('TEK0000.xlsx','TEK0000', 'E556:E1537')./R; %A
% Time = xlsread('TEK0000.xlsx','TEK0000', 'D556:D1537');
% Voltage = xlsread('TEK0001.xlsx','TEK0001', 'E556:E1537');

%300 Watts
% Current = xlsread('TEK0002.xlsx','TEK0002', 'E543:E1619')./R; %A
% Time = xlsread('TEK0002.xlsx','TEK0002', 'D543:D1619');
% Voltage = xlsread('TEK0003.xlsx','TEK0003', 'E543:E1619');

% %Single 100 Watt Curve
% Current = xlsread('TEK0004.xlsx','TEK0004', 'E45:E2480')./R; %A
% Time = xlsread('TEK0004.xlsx','TEK0004', 'D45:D2480');
% Voltage = xlsread('TEK0005.xlsx','TEK0005', 'E45:E2480');

% Lab Class 100 Watts
%  Current = xlsread('40war080bp70on','40war080bp70on', 'K4:K1003')./R; %A
%  Time = xlsread('40war080bp70on','40war080bp70on', 'A4:A1003');
%  Voltage = xlsread('40war080bp70on','40war080bp70on', 'C4:C1003');
 
%  Current = xlsread('40war130bp70on','40war130bp70on', 'K4:K1003')./R; %A
%  Time = xlsread('40war130bp70on','40war130bp70on', 'A4:A1003');
%  Voltage = xlsread('40war130bp70on','40war130bp70on', 'C4:C1003');

%  Current = xlsread('40war180bp70on','40war180bp70on', 'K4:K1003')./R; %A
%  Time = xlsread('40war180bp70on','40war180bp70on', 'A4:A1003');
%  Voltage = xlsread('40war180bp70on','40war180bp70on', 'C4:C1003');

%  Current = xlsread('40war230bp70on','40war230bp70on', 'K4:K1003')./R; %A
%  Time = xlsread('40war230bp70on','40war230bp70on', 'A4:A1003');
%  Voltage = xlsread('40war230bp70on','40war230bp70on', 'C4:C1003');

%  Current = xlsread('40war280bp70on','40war280bp70on', 'K4:K1003')./R; %A
%  Time = xlsread('40war280bp70on','40war280bp70on', 'A4:A1003');
%  Voltage = xlsread('40war280bp70on','40war280bp70on', 'C4:C1003');

 Current = xlsread('40war330bp70on','40war330bp70on', 'K4:K916')./R; %A
 Time = xlsread('40war330bp70on','40war330bp70on', 'A4:A916');
 Voltage = xlsread('40war330bp70on','40war330bp70on', 'C4:C916');

%  Current = xlsread('40war380bp70on','40war380bp70on', 'K4:K938')./R; %A
%  Time = xlsread('40war380bp70on','40war380bp70on', 'A4:A938');
%  Voltage = xlsread('40war380bp70on','40war380bp70on', 'C4:C938');

%  Current = xlsread('40war430bp70on','40war430bp70on', 'K4:K1003')./R; %A
%  Time = xlsread('40war430bp70on','40war430bp70on', 'A4:A1003');
%  Voltage = xlsread('40war430bp70on','40war430bp70on', 'C4:C1003');



%Sampling Frequency
Fs = length(Current)/(Time(end)-Time(1)); %Hz

%Probe Characteristics
ProbeHeight = 1.1/1000; %m
ProbeRadius = .25/2/1000; %m
ProbeArea = 2*pi*ProbeRadius*ProbeHeight+ pi*ProbeRadius^2; %m^2

e = 1.602e-19;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove High Frequency Noise
%[VF,CF] = HighFreqFilter(Fs,Voltage,Current);
VF = Voltage;
CF = Current;
%Fit the Data

% With good data use this
%  a = [1 5 5 0];
%  r = nlinfit(VF,CF,'langmuir',a);
%  l = langmuir(r,VF);
%  plot(VF,l,VF,CF)

% % %With Bad Data use this
%  r = fit(VF,CF,'poly7');
%  l = r(VF);
%  r = fit(VF,CF,'exp2');
%  l = r(VF);


%With Bad Data use this
r = fit(VF,CF,'poly9');
l = r(VF);

%Perfect Data
% VF = Voltage;% CF = Current;
% l = CF;


%plot
figure
plot(VF,CF,VF,l)
xlabel('Voltage [V]')
ylabel('Current [Amps]')
title('RF-Corrected I-V Curve ')
legend('RF-Corrected Current', 'Polynomial Fit','Location', 'southwest')

%Confirm the Fit is acceptable
prompt = 'Is this fit acceptable? (Hit Enter if Yes or Ctrl C if No):';
x = input(prompt);

clc
CF = l;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start Data Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Floating Potential

fp = mean(VF(find(min(abs(CF)) == abs(CF))));
fprintf('Floating Potential = %f Volts\n\n',fp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plasma Potential

dC = diff(CF);
dV = diff(VF);
dCdV = dC./dV; %Derivative of Current wrt Voltage

for i = 1:length(dCdV)

    if isnan(dCdV(i))==1
    dCdV(i) = 0;
    end
    if isfinite(dCdV(i)) == 0
        dCdV(i) = 0;
    end
    

end


 index1 = find(VF > fp);%find values above floating Potential
 %find the minimum value of the derivative above floating potential
 index1 = index1(1)+find(min(dCdV(index1(1):length(dCdV))) == dCdV(index1(1):length(dCdV)),1,'last');
 pp = VF(index1) ;
 fprintf('Plasma Potential = %f Volts\n\n',pp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mean Free Path
mfp = 0.061/gasp; %meters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Guessing Density

index2 = find(max(CF) == CF,1,'first');
Visat = VF(index2);
Iisat = CF(index2);
%Density
n = 1.42e15*max(CF)*sqrt(gasmass)/(ProbeArea*sqrt(pp-Visat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ion Current

mi = gasmass*1.6605402e-27 ;%ion mass, kg
Ii = (1/4).*(1.602e-19).*n.*sqrt(8*(1.602e-19)./(pi*mi)).*ProbeArea.*1.127.*sqrt(pp-VF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Electron Current

Ie = CF-Ii;

%fit Electron Current
a = [1 5 5 0];
Iefit = nlinfit(VF,Ie,'langmuir',a);
Iefitted = langmuir(Iefit,VF);
plot(VF,real(Iefitted),VF,real(Ie))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Electron Temperature

y = real(log(Ie)); % Natural Log of Electon Current.
l = diff(y)./diff(VF); %Derivative of Natural Log wrt Voltage
m = abs(diff(l)); %Getting rid of negatives
index = find(m < 0.0002); %find where the slope is relatively constant of the natural log

figure
plot(VF,real(log(Ie)),VF(index),y(index),'o')
xlabel('Probe Voltage [V]')
ylabel('Natural Log of Electron Current')
title('Natural Log of Electron Current vs. Voltage')
axis([fp inf -inf inf])



[xbounds,ybounds]=ginput(2);
xmin = find(VF <= min(xbounds),1,'last');
xmax = find(VF >= max(xbounds),1,'first');
ymin = find(y <= min(ybounds),1,'last');
ymax = find(y >= max(ybounds),1,'first');
% Slope = mean(l(ymin:ymax));
Slope = (y(ymax)-y(ymin))/(VF(xmax)-VF(xmin));

Te = 1/real(Slope);
fprintf('Electron Temperature = %f eV\n\n',Te)

close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Electron Energy Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dIdV = diff(Iefitted)./diff(VF);
dV = diff(VF);
dVF2 = (dV(1:end-1)+dV(2:end))/2;
d2IdV2 = diff(dIdV)./ dVF2;


 me = 9.109e-31; %kg mass of electron
 f = (-4/(ProbeArea*exp(2))).*sqrt(me*(pp-VF(1:length(VF)-2))./(2*e)).*d2IdV2;
 
 figure
 plot(real(VF(3:end)),real(f))
 %axis([0 inf -inf inf])
% f = f.*e;
% f= f./max(f);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find Real Density Using Laframboise Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start by Guessing
Iesat = ProbeArea*e*n*sqrt(e*Te/(2*pi*me));
Nguess = 3.73e13*Iesat/(ProbeArea*sqrt(Te));
N = Nguess;

%Normalized Potential
Xp = (pp-Visat)/Te;

for i=1:3
%Debye Length
lambD = 7430*sqrt(Te/N);

%rp/lambD
rpld = ProbeRadius/lambD;

%find jstar in lafromboise tables
if Xp < 25
[jstar]=LafrInterpolation(rpld,Xp);
else
jstar = LafrExtra(rpld,Xp);
end

N = 1.6e15*sqrt(gasmass/Te)*Iisat/(ProbeArea*jstar);

i=i+1;

end

if(isnan(N) == 1)

fprintf('Electron Density is %8.2E m^-3\n\n Lafr did not work\n\n',Nguess)

else

%Output Real Density
fprintf('Electron Density is %8.2E m^-3\n\n',N)
end

%Output Region Type

if ProbeRadius > 4*lambD && ProbeRadius<mfp
    fprintf('You are in the Collisionless Thin Sheath Regime\n\n')
elseif ProbeRadius < 4*lambD && 4*lambD<mfp
        fprintf('You are in the Collisionless Thick Sheath Regime\n\n')
elseif mfp < 4*lambD && 4*lambD < ProbeRadius
    fprintf('You are in the Collisional Thin Sheath Regime\n\n')
elseif mfp < ProbeRadius && ProbeRadius < 4*lambD
    fprintf('You are in the Collisional Thick Sheath Regime\n\n')
else
    fprintf('You are not in a Standard Regime\n\n')

end