function [ FPS] = getFPSFromSignalSpeed(signal,l,k) 
%Equation 12
%{
clear
k=7;
l=16;

seq='GSHMLSDEQMQIINSLVEAHHKTYDDSYSDFVRFRPPVREGPVTRSASRAASLHSLSDASSDSFNHSPESVDTKLNFSNLLMMYQDSGSPDSSEEDQQSRLSMLPHLADLVSYSIQKVIGFAKMIPGFRDLTAEDQIALLKSSAIEIIMLRSNQSFSLEDMSWSCGGPDFKYCINDVTKAGHTLEHLEPLVKFQVGLKKLKLHEEEHVLLMAICLLSPDRPGVQDHVRIEALQDRLCDVLQAYIRIQHPGGRLLYAKMIQKLADLRSLNEEHSKQYRSLSFQPEHSMQLTPLVLEVFGSEVS';

N = length(seq);
signal=zeros(1,N);
for i = 1:N
     signal(i) = codeAAHydrophobicity(seq(i));  
end
%}
y = getSignalProfile(signal,l);

% l is odd
isOdd=mod(l,2);
if isOdd ==1
  sa= getAutoCorr(y,0);
  sb=0;
  %sc=0;
  j= (l-1)/2 ;
  for q=1:j
    zq=getAutoCorr(y,q);
    sb=sb+zq*cos(q*2*pi*k/l);
  end

  FPS=sa+2*sb;
% l is even
else
  sa= getAutoCorr(y,0);
  sb=0;
  sc=0;
  j= l/2-1;
  for q=1:j
    zq=getAutoCorr(y,q);
    sb=sb+zq*cos(q*2*pi*k/l);
    sc=sc+y(q)*y(q+l/2);
  end
  sc=sc+y(l/2)*y(l);
 % for t=1:j+1
    %sc=sc+y(t)*y(t+l/2);
 % end
  FPS=sa+2*sb+2*cos(pi*k)*sc;
end %end of else

%FPS_check1 = getFPSFromSignalFast(signal,l,k)
%FPS_check2=getDFTPPS(signal,l/k)
end



