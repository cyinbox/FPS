function [ sigProfile] = getSignalProfile(signal,p)

N = length(signal);
sigProfile = zeros(1,p);

for i = 1:N
   m=mod(i,p);
   
   if m ==0
       k=p; %if no reminader, =0, periodicity set to p
   else
       k=m; 
   end
   
    sigProfile(k) = sigProfile(k)+signal(i);
end

end



