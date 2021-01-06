% Numerical Representations of Protein sequence using Kyte and Doolittle Hydrophobicity
% Changchuan Yin
% Last updated on 12/05/2012
%{
 @article{kyte1982simple,
  title={A simple method for displaying the hydropathic character of a protein},
  author={Kyte, Jack and Doolittle, Russell F},
  journal={Journal of molecular biology},
  volume={157},
  number={1},
  pages={105--132},
  year={1982},
  publisher={Elsevier}
}
Ile	4.5	
Val	4.2	
Leu	3.8	
Phe	2.8	
Cys	2.5	
Met	1.9	
Ala	1.8	
Gly	-0.4	
Thr	-0.7	
Ser	-0.8	
Trp	-0.9	
Tyr	-1.3	
Pro	-1.6	
His	-3.2	
Glu	-3.5	
Gln	-3.5	
Asp	-3.5	
Asn	-3.5	
Lys	-3.9	
Arg	-4.5
%}

function y = codeAAHydrophobicity(aa)
  switch upper(aa)
    case {'A'}
       y = 1.8;
    case {'R'}
       y = -4.5;
    case {'N'} 
       y = -3.5;
    case {'D'}
       y = -3.5;
    case {'C'}
       y = 2.5;
    case {'Q'} 
       y = -3.5;
    case {'E'}
       y = -3.5;
    case {'G'}
       y = -0.4;
    case {'H'} 
       y = -3.2;
    case {'I'}
       y = 4.5;
    case {'L'}
       y = 3.8;
    case {'K'} 
       y = -3.9;
    case {'M'} 
       y = 1.9;
    case {'F'}
       y = 2.8;
    case {'P'}
       y = -0.9;
    case {'S'} 
       y = -0.8;
    case {'T'}
       y = -0.7;
    case {'W'}
       y = -0.9;
   case{'Y'}
       y = -1.3;
   case{'V'}
       y = 4.2;
   otherwise
       y = 0.0000;%  Assign pad zero for non existing amino acid (padding zeros)
  end 
end