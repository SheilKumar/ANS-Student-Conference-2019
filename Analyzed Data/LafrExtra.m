function y = LafrExtra(rpld,Xp)

xp = xlsread('LafrCurves.xlsx','Sheet1', 'A2:A418'); % eVp/KTe
jstar = xlsread('LafrCurves.xlsx','Sheet1', 'B2:B418'); % Ip/Ip0 
rp =  xlsread('LafrCurves.xlsx','Sheet1', 'C2:C418'); % r_p/lambda_D
    F = scatteredInterpolant(rp,xp,jstar);
    y = F(rpld,Xp);