function df=Eqs_NK_NR(t,f,rBp, rNK, lNK, nMB, eBp, KBp, KBpr, KBpi)
%-----Variables-----
% f(1)=nB+
% f(2)=nTA
% f(3)=nTN

% rNK0 = initial growth rate of NKs
% lNK = apoptosis rate of NKs
% nMB = carrying capacity of B-ALLs
% eBP = rate of killing of B-ALLs by the NKs
% KBp = Michaelis constant for binding of CAR to B-ALLs
% KBpi = Michaelis constant for CAR-independent binding


%-----ODEs-----
df=zeros(2,1);


% rate of NK expansion in vivo

df(1)=rBp*f(1)*(1-f(1)/nMB)-eBp*f(1)/(f(1)+KBp)*f(2)-eBp*f(1)/(f(1)+KBpi)*f(2); 
% growth term minus two Michaelis-Menten binding-killing terms
% first term: CAR-dependent binding
% second term: CAR-independent binding

df(2)=f(1)/(f(1)+KBpr)*rNK*f(2)-lNK*f(2); 
% growth term minus apoptosis term

