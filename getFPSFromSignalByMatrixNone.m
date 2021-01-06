function [ FPS] = getFPSFromSignalByMatrixNone(signal,l,k) 
% Computation fractional period spectrum (FPS) from congruence vector of a singal
% Inputs: 1D real number signal, two integers l and k
%         matrix A (pre-stored object for fast speed).
% Output: Power spectrum at the periodicity l/k
%
% Changchuan Yin
% Department of Mathematics, Statistics and Computer Science
% University of Illinois at Chicago
% Last update 03/02/2016


% Note: This contains no pre-build matrix A

y = getSignalProfile(signal,l); %Congruence derivative vector y is the signal distribution profile
sa= getAutoCorr(y,0); 
sb=0;
sc=0;

isOdd=mod(l,2);
if isOdd ==1 % when l is odd
  j= (l-1)/2 ;
  for t=1:j
    r=mod(t*k,l);
    if r>j
        c=cos((l-r)*2*pi/l); % Or use the pre-stored matrix A object below line code
        %c=A(l-r); %It is %(One-r)
    else
        c=cos(r*2*pi/l);
        % if r~=0
        %  c=A(r); 
        % else
        %  c=cos(r*2*pi/l);
        % end
    end
    zt=getAutoCorr(y,t);% Self shift summation
    sb=sb+zt*c;
  end
  FPS=sa+2*sb;

else % when l is even
  j= l/2-1;
  for t=1:j
     r=mod(t*k,l);
     if r>j
         c=cos((l-r)*2*pi/l); % Or use the pre-stored matrix A object below line code
        %c=A(l-r);
     else
         c=cos(r*2*pi/l);% or use the prestored vector
         % if r~=0
         %  c=A(r); 
         % else
         %  c=cos(r*2*pi/l);
         % end
     end
     zt=getAutoCorr(y,t);
     sb=sb+zt*c;
     sc=sc+y(t)*y(t+l/2);
  end
  sc=sc+y(l/2)*y(l);
  FPS=sa+2*sb+2*cos(pi*k)*sc;
end

end



