function [ PPS] = getPSFromSignalFraction(sigProfile,p,k) 
%
for j=1:p
    C(j)=cos(2*(j-1)*k*pi/p); 
    V(j)=sin(2*(j-1)*k*pi/p);
end

%for j=1:p
%    v2(j)=sin(2*(j-1)*k*pi/p);
%end

%m1 = v1'* v1;
%m2= v2'* v2;
%cm = m1+m2;
U=C'*C+V'* V;

%Using matrix operation, make upper diagonal be zero
for i=1:p
    for j=1:p
        if i==j
            S(i,j)=1;
        %else
        %    S(i,j)=U(i,j)+U(j,i);
        end
        
        if (i>j)
             S(i,j)=U(i,j)+U(j,i);
        end
        if (i<j)
             S(i,j)=0;
        end
    end
end

% 2. Compute power spectrum using position profile and spectral matrix
PPS=sigProfile*S*sigProfile';
end

