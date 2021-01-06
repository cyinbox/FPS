
% Changchuan Yin
% Department of Mathematics, Statistics and Computer Science
% University of Illinois at Chicago
%
% Citation:
% Wang, J., & Yin, C. (2020). A fast algorithm for computing the Fourier spectrum of a fractional period. Journal of Computational Biology.
% Last update on 12/22/2020

clear

%4RUP, all alpha
seq='GSHMLSDEQMQIINSLVEAHHKTYDDSYSDFVRFRPPVREGPVTRSASRAASLHSLSDASSDSFNHSPESVDTKLNFSNLLMMYQDSGSPDSSEEDQQSRLSMLPHLADLVSYSIQKVIGFAKMIPGFRDLTAEDQIALLKSSAIEIIMLRSNQSFSLEDMSWSCGGPDFKYCINDVTKAGHTLEHLEPLVKFQVGLKKLKLHEEEHVLLMAICLLSPDRPGVQDHVRIEALQDRLCDVLQAYIRIQHPGGRLLYAKMIQKLADLRSLNEEHSKQYRSLSFQPEHSMQLTPLVLEVFGSEVS';
%for a=1:2
%    seq=[seq seq];
%

%Good data for 3.6 periods
%>4GAX:A|PDBID|CHAIN|SEQUENCE
seq='GSHMASMTGGQQMGRGSMSLTEEKPIRPIANFSPSIWGDQFLIVDNQVEQGVEQIVKDLKKEVRQLLKEALDIPMKHANLLKLVDEIQRLGISYLFEQEIDHALQHIYETYGDNWSGARSSLWFRLMRKQGYFVTCDVFNNHKDESGVFKQSLKNHVEGLLELYEATSMRVPGEIILEDALVFTQSHLSIIAKDTLSINPALSTEIQRALKKPLWKRLPRIEAVQYIPFYEQQDSHNKTLIKLAKLEFNLLQSLHREELSQLSKWWKAFDVKNNAPYSRDRIVECYFWALASRFEPQYSRARIFLAKVIALVTLIDDIYDAYGTYEELKIFTEAIERWSITCLDMIPEYMKPIYKLFMDTYTEMEEILAKEGKTNIFNCGKEFVKDFVRNLMVEAQWANEGHIPTTEELDSVAVITGGANLLTTTCYLGMSDIVTKEAFEWAVSEPPLLRYKGILGRRLNDLAGHKEEQERKHVSSSVESYMKEYNVSEEYAKNLLYKQVEDLWKDINREYLITKTIPRPLLVAVINLVHFLDVLYAAKDAFTAMGEEYKNLVKSLLVYPMSI';

%http://www.pnas.org/content/95/15/8580.full.pdf following sequence is from here
%seq='ATKAVCVLKGDGPVQGTIHFEAKGDTVVVTGSITGLTEGDHGFHVHQFGDNTQGCTSAGPHFNPLSKKHGGPKDEERHVGDLGNVTADKNGVAIVDIVDPLISLSGEYSIIGRTMVVHEKPDDLGRGGNEESTKTGNAGSRLACGVIGIAK'
%multiple peaks!

%seq='SGFEFHGYARSGVIMNDSGASTKSGAYITPAGETGGAIGRLGNQADTYVEMNLEHKQTLDNGATTRFKVMVADGQTSYNDWTASTSDLNVRQAFVELGNLPTFAGPFKGSTLWAGKRFDRDNFDIHWIDSDVVFLAGTGGGIYDVKWNDGLRSNFSLYGRNFGDIDDSSNSVQNYILTMNHFAGPLQMMVSGLRAKDNDERKDSNGNLAKGDAANTGVHALLGLHNDSFYGLRDGSSKTALLYGHGLGAEVKGIGSDGALRPGADTWRIASYGTTPLSENWSVAPAMLAQRSKDRYADGDSYQWATFNLRLIQAINQNFALAYEGSYQYMDLKPEGYNDRQAVNGSFYKLTFAPTFKVGSIGDFFSRPEIRFYTSWMDWSKKLNNYASDDALGSDGFNSGGEWSFGVQMETWF'


N = length(seq)
signal=zeros(1,N);
for i = 1:N
     signal(i) = codeAAHydrophobicity(seq(i));  
end
signal


%-------------------------TEST-----------------
k=5;
l=18;

%k=10;
%l=36;

y = getSignalProfile(signal,l);
PS_36ByFPS = getPSFromSignalFraction(y,l,k) 

%or
%PS_36ByFPS=getFPSFromSignalCongruence(y,l,k) 

p=3.6;

k=5;
l=18; %Example 36/10=3.6

%use pre-populated object coefficient matrix A for  k=1
for q=1:l
  A(q)=cos(q*2*pi/l); 
end

FPSFast_P36=getFPSFromSignalFast(signal,l,k) %Good validated, =2.1864e+04
FPS_Value_Speed = getFPSFromSignalSpeed(signal,l,k) %Good, validated,
FPS_Value_MatrixNone=getFPSFromSignalByMatrixNone(signal,l,k) %Good 
FPS_Value_MatrixNewCorrected=getFPSFromSignalByMatrixNewCorrected(signal,l,k,A) %Good  with corrected 

%Need to perform comlexity test here
%{

tStart2 = tic; 
 for k=1:l
  FPSFast=getFPSFromSignalFast(signal,l,k);
 end
tElapse_Fast = toc(tStart2)

tStart3 = tic; 
 for k=1:l
  FPSSpeed=getFPSFromSignalSpeed(signal,l,k);
 end
tElapsed_Speed = toc(tStart3)
%}

tStart4 = tic; 
 for k=1:l
  FPSSuper = getFPSFromSignalByMatrixNone(signal,l,k);
 end
tElapsed_Super = toc(tStart4)

%use pre-stored matrix A: This algorithm is described in the paper
tStart5 = tic; 
 for k=1:l
  FPSSuper2 = getFPSFromSignalByMatrixNewCorrected(signal,l,k,A);
 end
tElapsed_Super2 = toc(tStart5)


%---------------------------------------------------------------------
p1=2;
p2=10;

%FPS = getFPSProteinRange(seq,p1,p2);
%j=p1:0.1:p2;
%s=size(j);


k=10;
%Note: l will be 30,31,32,33,....60, PS(30/10,31/10,32/10,..60/10)
idx=1;
for j = p1:0.1:p2
 l=j*10;
 l=round(l); %52.000-->52
 profile=getProteinProfile(seq,l);
 FPS(idx)=getFPSFromSignalFast(signal,l,k); %Both are the same
 idx=idx+1;
end


%---------------------------------------------------------------------
j=p1:0.1:p2;
figure
%fig1=plot(j,FPS)
fig1=stem(j,FPS,'filled')
%title('Fractional power spectrum of protein all alpha','FontSize',8,'FontWeight','bold');
set(fig1                        , ...
  'LineWidth'       , 1.5 );
hXLabel = xlabel('periodicity');
hYLabel = ylabel('power spectrum');
set([hXLabel, hYLabel]  , ...
    'FontName'   , 'AvantGarde', ...
    'FontSize'   , 10, ...
    'FontWeight' , 'bold');
set(gca, ...
  'Box'         , 'off'     , ...  %No rectangle cover the figure
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1,... 
  'YColor'      , [.3 .3 .3]);


%---------------------------------------------------------------------
Ps=[1.0, 1.5, 1.0, 1.5, 2.0, 1.5, 2.0, 1.5, 1.0, 2.5, 2.0, 1.5, 3.0, 1.5, 2.0, 2.5, 2.0, 2.5, 2.0, 1.5, 1.0, 4.5, 2.0, 1.5, 3.0, 2.5, 2.0, 1.5, 5.0, 1.5, 1.0, 3.5, 4.0, 1.5, 2.0, 3.5, 2.0, 2.5, 3.0, 2.5, 3.0, 2.5, 2.0, 1.5, 1.0, 8.5, 0.0, 0.5, 5.0, 2.5, 0.0, 3.5, 3.0, 0.5, 4.0, 1.5, 3.0, 2.5, 2.0, 1.5, 9.0, 1.5, -1.0, 5.5, 2.0, -0.5, 5.0, 5.5, 2.0, 2.5, 4.0, 2.5, 3.0, 1.5, 4.0, 1.5, 2.0, 3.5, 4.0, 2.5, 2.0, 4.5, 2.0, 2.5, 4.0, 3.5, 2.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 16.5, -7.0, 3.5, 5.0, 6.5, -2.0, 3.5, 3.0, -1.5, 8.0, -0.5, 1.0, 4.5, 6.0, -1.5, 5.0, 4.5, -4.0, 5.5, 9.0, -2.5, 4.0, 9.5, 3.0, -1.5, 12.0, -2.5, 0.0, 16.5, -5.0, 8.5, 5.0, 0.5, -1.0, 8.5, -2.0, 5.5, 6.0, 5.5, 13.0, -6.5, 3.0, 8.5, -4.0, 10.5, 1.0, -1.5, 6.0, 3.5, 1.0, 4.5, 6.0, 2.5, 6.0, 4.5, 0.0, 6.5, 1.0, 3.5, 3.0, 0.5, 5.0, 4.5, -1.0, 5.5, 2.0, 4.5, 3.0, 2.5, 6.0, 0.5, 4.0, 4.5, 3.0, 1.5, 5.0, 2.5, 1.0, 5.5, 3.0, 3.5, 4.0, 3.5, 5.0, 2.5, 3.0, 3.5, 3.0, 4.5, 2.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 32.5, -24.0, -0.5, 18.0, 2.5, -14.0, 15.5, 11.0, -6.5, 6.0, 23.5, -19.0, -0.5, 9.0, 5.5, 15.0, 0.5, 4.0, 1.5, 9.0, -1.5, 8.0, 2.5, -7.0, 1.5, 12.0, 3.5, -1.0, 8.5, -4.0, -5.5, 19.0, -1.5, 2.0, 2.5, -1.0, -6.5, 6.0, -0.5, 7.0, 5.5, 4.0, -2.5, 21.0, -7.5, 4.0, 0.5, 1.0, 10.5, 12.0, -7.5, 12.0, 1.5, -4.0, 2.5, -2.0, 5.5, -2.0, 11.5, 4.0, -9.5, 23.0, 10.5, -7.0, -4.5, 15.0, 17.5, -22.0, 9.5, 8.0, 8.5, 4.0, -11.5, 14.0, 2.5, 10.0, -4.5, 17.0, -1.5, 0.0, 0.5, 12.0, 11.5, -1.0, 6.5, 5.0, 13.5, -10.0, -0.5, 11.0, -0.5, 7.0, -5.5, 9.0, 16.5, -18.0, 3.5, 9.0, 3.5, 9.0, 6.5, -10.0, 5.5, 8.0, -1.5, -4.0, 13.5, 2.0, 7.5, 4.0, -1.5, 9.0, 3.5, 0.0, 2.5, 9.0, 0.5, 5.0, -1.5, 10.0, 1.5, -2.0, 11.5, 4.0, -1.5, 2.0, 9.5, -6.0, 9.5, 10.0, -1.5, 0.0, 10.5, -1.0, 5.5, -2.0, 8.5, 7.0, -3.5, 6.0, 6.5, 1.0, 5.5, 5.0, 0.5, 7.0, 1.5, 2.0, 3.5, 4.0, 4.5, 4.0, 2.5, 4.0, 1.5, 4.0, 8.5, 1.0, 3.5, 3.0, 8.5, 1.0, 1.5, 7.0, 2.5, 1.0, 7.5, 0.0, 5.5, 5.0, 3.5, 5.0, 2.5, 5.0, 5.5, 3.0, 2.5, 4.0, 5.5, 2.0, 2.5, 4.0, 5.5, 1.0, 4.5, 2.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 64.5, -57.0, 9.5, 17.0, -5.5, 6.0, 25.5, -9.0, 0.5, 6.0, 45.5, -48.0, 3.5, 31.0, -27.5, -6.0, 38.5, -12.0, -7.5, 5.0, 26.5, 0.0, -9.5, 17.0, -7.5, 6.0, -2.5, 9.0, 14.5, -24.0, 10.5, 27.0, -3.5, 0.0, 6.5, 6.0, 12.5, 8.0, -16.5, 20.0, -25.5, 4.0, 23.5, -3.0, 26.5, -10.0, -3.5, 15.0, 14.5, -52.0, 17.5, 24.0, 21.5, -4.0, 1.5, 3.0, 10.5, 2.0, 1.5, -9.0, 4.5, 6.0, 7.5, 5.0, 6.5, -10.0, 3.5, 27.0, -22.5, 16.0, 6.5, -8.0, 24.5, 2.0, -13.5, 14.0, 8.5, 12.0, -2.5, 13.0, -3.5, -4.0, 5.5, 7.0, -12.5, 15.0, -6.5, 3.0, 1.5, 12.0, 0.5, -6.0, 0.5, 7.0, 11.5, -5.0, -0.5, -6.0, 21.5, -13.0, -13.5, 44.0, -3.5, -8.0, 27.5, 5.0, -4.5, -3.0, 6.5, 15.0, -13.5, 6.0, 20.5, -5.0, -14.5, 20.0, 7.5, -6.0, 13.5, -24.0, 21.5, 21.0, -12.5, 4.0, 4.5, -9.0, 7.5, 2.0, 0.5, 10.0, 18.5, -5.0, -5.5, 34.0, -8.5, 4.0, 16.5, 4.0, 10.5, -38.0, 41.5, -8.0, 3.5, 23.0, -11.5, 12.0, 24.5, -32.0, 18.5, 13.0, -10.5, -15.0, 7.5, 22.0, -4.5, -11.0, 14.5, 21.0, 4.5, -26.0, 30.5, 18.0, -26.5, 24.0, 4.5, -23.0, 12.5, 2.0, 6.5, 3.0, -4.5, 15.0, -21.5, 22.0, 16.5, -52.0, 33.5, 21.0, -20.5, 18.0, 19.5, -27.0, 28.5, -3.0, 6.5, -2.0, 1.5, 19.0, 2.5, 16.0, 8.5, 0.0, 16.5, 8.0, -0.5, 13.0, 14.5, -18.0, 14.5, 26.0, -19.5, -9.0, 35.5, -31.0, 13.5, 7.0, 6.5, 1.0, 4.5, -16.0, 23.5, 0.0, -2.5, 8.0, 13.5, -38.0, 18.5, 21.0, -13.5, -5.0, 12.5, 11.0, -14.5, 2.0, 25.5, -10.0, 2.5, 10.0, -9.5, 16.0, 13.5, -16.0, -0.5, 13.0, 0.5, 6.0, 7.5, 0.0, -1.5, 11.0, -1.5, 9.0, -4.5, -5.0, 24.5, 10.0, -7.5, 19.0, 7.5, -12.0, 9.5, 13.0, -9.5, 15.0, -5.5, 4.0, 1.5, 17.0, 1.5, 4.0, -9.5, 14.0, 2.5, 1.0, 9.5, -7.0, 2.5, 19.0, -10.5, 3.0, 5.5, 2.0, 5.5, 16.0, 0.5, -7.0, 11.5, 10.0, -1.5, 2.0, 11.5, -1.0, 0.5, 7.0, 10.5, -2.0, 3.5, 3.0, 8.5, 7.0, -4.5, 10.0, 7.5, -5.0, 4.5, 3.0, 5.5, -1.0, 9.5, 0.0, 5.5, 2.0, 3.5, 6.0, 2.5, 4.0, 4.5, 5.0, 4.5, 4.0, 1.5, 6.0, 9.5, -1.0, 3.5, 12.0, 3.5, 5.0, 1.5, 9.0, 6.5, 1.0, 6.5, 7.0, 3.5, 5.0, 0.5, 9.0, 5.5, 3.0, 4.5, 1.0, 12.5, -2.0, 7.5, 4.0, -2.5, 10.0, 10.5, 0.0, 1.5, 6.0, 9.5, 0.0, 3.5, 8.0, 7.5, -7.0, 11.5, 2.0, 0.5, 5.0, 4.5, 6.0, 3.5, 4.0, 2.5, 4.0, 5.5, 3.0, 4.5, 3.0, 4.5, 5.0, 3.5, 1.0, 3.5, 7.0, 6.5, -2.0, 3.5, 3.0, 6.5, 2.0, 5.5, 1.0, 4.5, 2.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 128.5, -106.0, 13.5, 66.0, -23.5, -13.0, 21.5, -8.0, 27.5, -33.0, 28.5, -19.0, -1.5, 26.0, 39.5, -28.0, 4.5, 2.0, 23.5, -4.0, -12.5, 40.0, -4.5, -28.0, 13.5, 96.0, -112.5, 20.0, 50.5, -31.0, 9.5, 17.0, 52.5, -90.0, 29.5, 50.0, -32.5, 15.0, 11.5, 19.0, 40.5, -7.0, -34.5, 34.0, 6.5, -26.0, 20.5, 22.0, 2.5, 0.0, -10.5, 23.0, -27.5, 12.0, 57.5, -39.0, -28.5, 52.0, -16.5, 3.0, 36.5, -40.0, 29.5, 25.0, -31.5, 0.0, 26.5, -6.0, 0.5, 33.0, 13.5, -11.0, 6.5, -22.0, 9.5, -11.0, 29.5, 61.0, -87.5, 17.0, 38.5, -4.0, -25.5, 8.0, 31.5, -3.0, 32.5, 4.0, 17.5, 28.0, -53.5, 5.0, 92.5, -49.0, -23.5, 83.0, -23.5, -24.0, -2.5, 24.0, 0.5, 27.0, 14.5, -39.0, 33.5, 37.0, -46.5, 42.0, 30.5, -57.0, 6.5, 29.0, 21.5, -36.0, 27.5, 27.0, -33.5, 11.0, -28.5, 39.0, -13.5, 19.0, 3.5, 7.0, 10.5, -10.0, -11.5, 34.0, 28.5, -26.0, -16.5, 36.0, 5.5, -16.0, -1.5, 31.0, -38.5, 31.0, 39.5, -33.0, 21.5, 29.0, -3.5, 23.0, 20.5, -1.0, 6.5, -16.0, 15.5, 6.0, 2.5, -25.0, 22.5, -3.0, -2.5, 5.0, -8.5, 18.0, 15.5, -21.0, 3.5, -1.0, 19.5, -16.0, -7.5, -7.0, 2.5, 9.0, 10.5, 17.0, 11.5, 7.0, -3.5, -11.0, -4.5, 2.0, 16.5, -5.0, -23.5, 0.0, 26.5, -22.0, -12.5, 43.0, -1.5, -1.0, -3.5, 24.0, -30.5, 36.0, 4.5, -13.0, -1.5, 1.0, -18.5, 25.0, -2.5, 0.0, -19.5, 76.0, -0.5, -43.0, 26.5, 42.0, -15.5, -29.0, 78.5, -35.0, 29.5, 4.0, -7.5, -20.0, 36.5, -15.0, 18.5, 4.0, -11.5, 21.0, -4.5, 13.0, 25.5, -38.0, -0.5, 0.0, -18.5, -10.0, 53.5, -24.0, -44.5, 107.0, -6.5, 0.0, 2.5];

p1=2;
p2=100;

%Note: l will be 30,31,32,33,....60, PS(30/10,31/10,32/10,..60/10)
idx=1;
for j = p1:1:p2
 l=j*10;
 l=round(l); %52.000-->52
 profile=getSignalProfile(Ps,l);
 FPS(idx)=getFPSFromProtein(profile,l,k); %Both are the same
 idx=idx+1;
end
FPS
j=p1:1:p2;

%---------------------------------------------------------------------
figure
fig1=plot(j,FPS)
%fig1=stem(j,FPS,'filled')

%title('Fractional power spectrum of protein all alpha','FontSize',8,'FontWeight','bold');
set(fig1                        , ...
  'LineWidth'       , 1.5 );
hXLabel = xlabel('periodicity');
hYLabel = ylabel('power spectrum2');
set([hXLabel, hYLabel]  , ...
    'FontName'   , 'AvantGarde', ...
    'FontSize'   , 10, ...
    'FontWeight' , 'bold');
set(gca, ...
  'Box'         , 'off'     , ...  %No rectangle cover the figure
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1,... 
  'YColor'      , [.3 .3 .3]);



