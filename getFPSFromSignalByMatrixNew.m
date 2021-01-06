function [ FPS] = getFPSFromSignalByMatrixNew(signal,l,k,A) 
% Computation fractional period spectrum (FPS) from congruence vector of a singal
% Inputs: 1D real number signal, two integers l and k
%         matrix A (pre-stored object for fast speed).
% Output: Power spectrum at the periodicity l/k
%v
% Note: This probram does not need the pre-build matrix A as described our paper

% Last update: 11/07/2019

y = getSignalProfile(signal,l); %Congruence derivative vector y is the signal distribution profile
sa = getAutoCorr(y,0); 
sb = 0;
sc = 0;

%use pre-populated object coefficient matrix A for  k=1

%for q=1:l
%  A(q)=cos(q*2*pi/l); 
%end

isOdd = mod(l,2);
if isOdd == 1 % when l is odd
  j = (l-1)/2 ;
  for t =1:j
    r = mod(t*k,l);
    if r>j
        %c = cos((l-r)*2*pi/l); % Or use the pre-stored matrix A object below line code
        c=A(l-r); %It is %(One-r)
    else
        %c = cos(r*2*pi/l);
         if r~=0
          c = A(r); 
         else
          c = cos(r*2*pi/l);
         end
    end
    zt = getAutoCorr(y,t);% Self shift summation
    sb = sb+zt*c;
  end
  FPS = sa+2*sb;

else % when l is even
  j = l/2-1;
  for t =1:j
     r = mod(t*k,l);
     if r>j
        % c=cos((l-r)*2*pi/l); % Or use the pre-stored matrix A object below line code
        c = A(l-r);
     else
         %c=cos(r*2*pi/l);% or use the prestored vector
          if r~=0
           c = A(r); 
          else
           c = cos(r*2*pi/l);
          end
     end
     zt = getAutoCorr(y,t);
     sb = sb+zt*c;
     sc = sc+y(t)*y(t+l/2);
  end
  sc = sc+y(l/2)*y(l);  %sc is Z(l/2,l/2)
  FPS = sa+2*sb+2*cos(pi*k)*sc;
end

end



