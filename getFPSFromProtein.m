function [ FPS] = getFPSFromProtein(profile,l,k) 
A=zeros(l,l);
%Equation (10):it is validated here.
%{
for r=1:l
    for s=1:l
         vcos=cos((r-1)*2*pi*k/l)*cos((s-1)*2*pi*k/l); 
         vsin=sin((r-1)*2*pi*k/l)*sin((s-1)*2*pi*k/l); 
         A(r,s)=vcos+vsin;
    end
end  
%}

%Equation (11):validated correctness
for r=1:l
    for s=1:l
          A(r,s)=cos((r-s)*2*pi*k/l); 
    end
end  

FPS=0;
for j=1:20
  w=profile(j,:);
  FPS=FPS+w*A*w';
end

end

