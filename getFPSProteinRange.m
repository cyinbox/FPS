function [FPS] = getFPSProteinRange(seq,p1,p2)
% Fractional periods from integer period p1 to p2, with stepwise 0.1
%{
clear 
disp('Comparison of different protein DFT spectrum');
 
%4RUP, all alpha
seq='GSHMLSDEQMQIINSLVEAHHKTYDDSYSDFVRFRPPVREGPVTRSASRAASLHSLSDASSDSFNHSPESVDTKLNFSNLLMMYQDSGSPDSSEEDQQSRLSMLPHLADLVSYSIQKVIGFAKMIPGFRDLTAEDQIALLKSSAIEIIMLRSNQSFSLEDMSWSCGGPDFKYCINDVTKAGHTLEHLEPLVKFQVGLKKLKLHEEEHVLLMAICLLSPDRPGVQDHVRIEALQDRLCDVLQAYIRIQHPGGRLLYAKMIQKLADLRSLNEEHSKQYRSLSFQPEHSMQLTPLVLEVFGSEVS';
p1=2;
p2=6;
%}

k=10;
%Note: l will be 30,31,32,33,....60, PS(30/10,31/10,32/10,..60/10)
idx=1;
for j = p1:0.1:p2
 l=j*10;
 l=round(l); %52.000-->52
 
 profile=getProteinProfile(seq,l);
 FPS(idx)=getFPSFromProtein(profile,l,k); %Both are the same
 idx=idx+1;
end

%{
j=p1:0.1:p2;
s=size(j)

figure
fig1=plot(j,FPS)
title('Fourier spectrum at periodicities (3.1 to 4.0) of protein domain 1avyA00','FontSize',8,'FontWeight','bold');
%}
%return FPS
end

