function [jstar1]=LafrInterpolation(rpld,Xp)

xp = xlsread('MasterColdPlasma.xlsx','Sheet1', 'A4:A412'); % eVp/KTe
jstar = xlsread('MasterColdPlasma.xlsx','Sheet1', 'L4:S412'); % Ip/Ip0 
rp =  xlsread('MasterColdPlasma.xlsx','Sheet1', 'L3:S3'); % r_p/lambda_D

%inter = something;
 jstar1 = griddata(rp,xp,jstar,rpld,Xp);
 
 %Turn off Warning
w = warning('query','last');
id = w.identifier;
warning('off',id)
