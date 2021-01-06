%Theorem 5
clear
format short g

l=18
%l=36;
%18/1, 18/2,18/3,18/
%300/1
%use pre-populated object coefficient matrix A 
for q=1:l
  A(q)=cos(q*2*pi/l); 
end

% Good data for 3.6 periods
%>4GAX:A|PDBID|CHAIN|SEQUENCE
seq = 'GSHMASMTGGQQMGRGSMSLTEEKPIRPIANFSPSIWGDQFLIVDNQVEQGVEQIVKDLKKEVRQLLKEALDIPMKHANLLKLVDEIQRLGISYLFEQEIDHALQHIYETYGDNWSGARSSLWFRLMRKQGYFVTCDVFNNHKDESGVFKQSLKNHVEGLLELYEATSMRVPGEIILEDALVFTQSHLSIIAKDTLSINPALSTEIQRALKKPLWKRLPRIEAVQYIPFYEQQDSHNKTLIKLAKLEFNLLQSLHREELSQLSKWWKAFDVKNNAPYSRDRIVECYFWALASRFEPQYSRARIFLAKVIALVTLIDDIYDAYGTYEELKIFTEAIERWSITCLDMIPEYMKPIYKLFMDTYTEMEEILAKEGKTNIFNCGKEFVKDFVRNLMVEAQWANEGHIPTTEELDSVAVITGGANLLTTTCYLGMSDIVTKEAFEWAVSEPPLLRYKGILGRRLNDLAGHKEEQERKHVSSSVESYMKEYNVSEEYAKNLLYKQVEDLWKDINREYLITKTIPRPLLVAVINLVHFLDVLYAAKDAFTAMGEEYKNLVKSLLVYPMSI';

N = length(seq);
signal = zeros(1,N);
for i = 1:N
     signal(i) = codeAAHydrophobicity(seq(i));  
end

% This is the congruence derivative sequence of size l (18)
y = getSignalProfile(signal,l);%this can be used for computing integer period power spectrum
%y = congruenceDerivativeSequence(signal,l);%this can be used for computing integer period power spectrum
%y = [13.5,4.1 -31.6,-26  -9.1, -6  -3.8,  2  -0.3    -7.9 -38.5,  9,  3 -18.6 -26.4  18.9,6.1 -21.6]

jj=floor(l/2);

FPS=zeros(1,jj);
for k=1:jj
 sa= getAutoCorr(y,0); 
 sb=0;
 sc=0;
 isOdd=mod(l,2);
 if isOdd ==1 % when l is odd
  j= (l-1)/2 ;
  for t=1:j
    r=mod(t*k,l);
    if r>j
        %c=cos((l-r)*2*pi/l); % Or use the pre-stored matrix A object below line code
        c=A(l-r); %It is 'One-r'
        
    else
        %c=cos(r*2*pi/l);
         if r~=0
          c=A(r); 
         else
          c=1; % it is 1
         end
    end
    zt=getAutoCorr(y,t);% Self shift summation
    sb=sb+zt*c;
  end
  FPS(k)=sa+2*sb;
 
 else % when l is even
  j= l/2-1;
  for t=1:j
     r=mod(t*k,l);
      if r>j
        % c=cos((l-r)*2*pi/l); % Or use the pre-stored matrix A object below line code
        c=A(l-r);
     else
         %c=cos(r*2*pi/l);% or use the prestored vector as follows
          if r~=0
           c=A(r); 
          else
           c=1; %it is one
          end
     end
     zt=getAutoCorr(y,t);
     sb=sb+zt*c;
     sc=sc+y(t)*y(t+l/2);
  end
 
  sc=sc+y(l/2)*y(l);
  FPS(k)=sa+2*sb+2*cos(pi*k)*sc;
 end
 
end

[maxFPS,idx] = max(FPS);
disp(FPS)
disp(maxFPS)
disp(idx)

%names = {'18/1'; '18/2'; '18/3'; '18/4'; '18/5';'18/6'; '18/7'; '18/8'; '18/9'};

figure
k=1:1:jj;
fig1=stem(k,FPS,'filled')
%title('Fractional power spectrum of protein all alpha','FontSize',8,'FontWeight','bold');
set(fig1                        , ...
  'LineWidth'       , 1.5 );
hXLabel = xlabel('k, (l = 18)');
%hXLabel = xlabel('k, (l = 36)');
%set(gca,'xtick',[1:jj],'xticklabel',names)
hYLabel = ylabel('power spectrum');
set([hXLabel,hYLabel]  , ...
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

