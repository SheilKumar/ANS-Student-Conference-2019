function [VF,CF] = HighFreqFilter(Fs,Voltage,Current)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Low Pass Filter
Wn = (2*1e4)/Fs;
[b,a] = butter(6,Wn,'low');
Current_Filter = filter(b,a,Current);

Wn = (2*1e4)/Fs;
[b,a] = butter(6,Wn,'low');
Voltage_Filter = filter(b,a,Voltage);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 index = find(max(abs(Voltage_Filter))==abs(Voltage_Filter));
 VF = Voltage_Filter(index:length(Voltage_Filter));
 CF = Current_Filter(index:length(Current_Filter));

 index1 = find(min(Voltage) >= VF,1,'last');
 VF = VF(index1:end);
 CF = CF(index1:end);

VF = Voltage_Filter;
CF = Current_Filter;
