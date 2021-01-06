function [ sum ] = getAutoCorr(y,q)
%clear
%y=[1 2 3 4 5 6]
%q=2;

l=length(y);
sum=0;

 for t=1:l
    
   c= t+q;
   h=mod(c,l); % if q =0, then 
  
   if h ==0
       r=l;
   else
       r=h;
   end
  
   sum=sum+y(t)*y(r);
 end
 
 %sum;
 
end

