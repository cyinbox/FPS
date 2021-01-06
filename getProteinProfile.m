function [ profile] = getProteinProfile(protein,periodicity)
%  Function to get protein hydrophobicity profile according to periodicity. It is to sum
%  hydrophobicity values from peridic positions for each amino acid.
%  Inputs: protein sequence, periodicity
%  Output: periodic hydrophobicity profile
%  Changchuan Yin 01/15/2016

f = periodicity;
N = length(protein);
profile=zeros(20,f);

for i=1:N
   m=mod(i,f);
   if m==0
       p=f; 
   else
       p=m; 
   end
   
   switch upper(protein(i))
    case {'A'}
       profile(1,p)=profile(1,p)+1.8;
       
    case {'R'}
       profile(2,p)=profile(2,p)-4.5;
       
    case {'N'} 
        
       profile(3,p)=profile(3,p)-3.5;
    case {'D'}
      
       profile(4,p)=profile(4,p)-3.5;
    case {'C'}
      
       profile(5,p)=profile(5,p)+2.5;
    case {'Q'} 
       
       profile(6,p)=profile(6,p)-3.5;
    case {'E'}
       
       profile(7,p)=profile(7,p)-3.5;
    case {'G'}
       
       profile(8,p)=profile(8,p)-0.4;
    case {'H'} 
       
       profile(9,p)=profile(9,p)-3.2;
    case {'I'}
       
       profile(10,p)=profile(10,p)+4.5;
    case {'L'}
       
       profile(11,p)=profile(11,p)+3.8;
    case {'K'} 
      
       profile(12,p)=profile(12,p)-3.9;
    case {'M'} 
      
       profile(13,p)=profile(13,p)+1.9;
    case {'F'}
     
       profile(14,p)=profile(14,p)+2.8;
    case {'P'}
     
       profile(15,p)=profile(15,p)-0.9;
    case {'S'} 
      
       profile(16,p)=profile(16,p)-0.8;
    case {'T'}
     
      profile(17,p)=profile(17,p)-0.7;
    case {'W'}
      
       profile(18,p)=profile(18,p)-0.9;
   case{'Y'}
       
       profile(19,p)=profile(19,p)-1.3;
   case{'V'}
       profile(20,p)=profile(20,p)+4.2;
       
   otherwise
       y = 0.0000;%  assign pad zero for non existing amino acid (padding zeros)
  end 
  
end

end
