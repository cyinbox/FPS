function [ sigProfile] = congruenceDerivativeSequence(x,p)

N = length(x);
sigProfile = zeros(1,p);


for i = 1:N
   m=mod(i,p);
   
   if m ==0
       k=p; %if no reminader, =0, periodicity set to p
   else
       k=m; 
   end
   
    sigProfile(k) = sigProfile(k)+x(i);
end

end



