function f=langmuir(a,x);
% This function describes any Langmuir probe characteristics using four
% free fitting parameters a(1) ... a(4) details related to this equation
% can be seen in REVIEW OF SCIENTIFIC INSTRUMENTS 79, 103501 (2008)
f=exp(a(1)*tanh((x+a(2))/a(3)))+a(4);