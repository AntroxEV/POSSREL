function [Out]=FuzzyClustering(data1D,ncl,mPar,nAlpha,ACFoption,ARoption,galpha)
% Prof Luciano Stefanini - University of  Urbino
% modified by Dr Alessandro Tombari - University of Exeter 
% Please cite: https://doi.org/10.1016/j.compgeo.2023.105967
%% Input data
% data1D: array of data
% ncl: number of clusters
%  mPar = 1.8, 2.0 exponent for fuzzy clustering
% nAlpha: alpha discretization (number)
% ACFoption = 1 ACF-based transformatation; else AAT-based transformation
% ARoption = 0 (no expansion is performed) 
% galpha >= 1.0 is the expansion factor for ACF-based transformation
% if galpha = 1.0, no expansion is performed
%% DATA - values at given location
phi=data1D(:);
n=length(phi);
phi = sort(phi,'ascend');
%% ***** Step 1: Apply Fuzzy clustering
% Run partition with specified mPar 
% Option for full AC function (1 (Average based), 2(DAAT based)) 
%Create the fuzzy partition of phi values into ncl fuzzy classes
phiU = [];
phiCenters = [];
indCenters = [];
parCenters = [];
parSupports = [];
nMFpar = [];
parPoints = [];
parMFval = [];
nMFset = [];
setPoints = [];
setMFval = [];
parNphi = [];
parMFs = []; 
 
%Data defining the final fuzzy partition
parSupports = zeros(ncl,2); %supports of fuzzy partition  
parNphi = zeros(ncl,2); % first and last phi index in fuzzy partition
parMFs = zeros(ncl,n);   %membership values of final fuzzy partition according to fcm
setMFs = zeros(ncl,n); %membeship of phi within subsets of fuzzy partition
phiSwitch = zeros(ncl,1);
%
phiClusters = zeros(n,1); %%index of prevailing membership (0/1 clusters)
listClusters = zeros(1,ncl);
cluNphi = zeros(ncl,2); % first and last phi index in in 0/1 clusters 
%

%
 optCluster = [mPar,500,1.e-6,0];
 [phiCenters,phiU] = fcm(phi,ncl,optCluster); %use modified parameters
%
% reorder clusters by increasing values of phiCenters
 [parC,indC]=sort(phiCenters);
 parCenters = parC;
 indCenters = indC;
 phiUmax = max(phiU);

%
 %% ***** Step 2a: Fuzzy Partition - deduce clusters and separation
 % Ruspini
 % Element phi(i) is first assigned to cluster k if phiU(i,k) is maximum
 % cluster k contains elements of phi from cluNphi(k,1) to cluNphi(k,2)
 % with membership value phiUmax(i)
 listClusters(1) = 1;
 phil = 1;
 for k = 1:ncl
    kcl = indCenters(k);
    cl = find(phiU(kcl,:) == phiUmax);
    lcl = length(cl);
    phiClusters(phil:phil+lcl-1) = cl(1:lcl);
    listClusters(k) = phil+lcl-1;
    phil = phil+lcl;
    phiSwitch(k) = phi(listClusters(k));
 end
 for k = 1:ncl
    if k == 1
        cluNphi(k,1) = 1;
    else
        cluNphi(k,1) = listClusters(k-1);
    end
    cluNphi(k,2) = listClusters(k);
 end

%
%% ***** Step 2b: Fuzzy Partition - compute Empirical membership functions
% the fuzzy partition assigns each phi to exactly two clusters
% primary cluster with membership phiUmax(i) and
% secondary cluster with membership 1-phiUmax(i) 
% Final partition contains ncl partially overlapping clusters between
% successive core values. 
% A fuzzy partition of interval [phiMin,phiMax] is then composed by ncl
% fuzzy intervals 
parSupports(1,1) = min(phi);
parSupports(1,2) = parCenters(2);
parSupports(ncl,1) = parCenters(ncl-1);
parSupports(ncl,2) = max(phi);
for k = 2:(ncl-1)
    parSupports(k,1) = parCenters(k-1);
    parSupports(k,2) = parCenters(k+1);
end
parNphi(1,1) = 1;
parNphi(ncl,2) = n;
for k = 1:ncl
    if k == 1
        parNphi(k,2) = find(phi <= parCenters(k+1),1,'last');
    elseif k == ncl
        parNphi(k,1) = find(phi >= parCenters(k-1),1,'first');
    else
        parNphi(k,2) = find(phi <= parCenters(k+1),1,'last');
        parNphi(k,1) = find(phi >= parCenters(k-1),1,'first');
    end
end
% Construct final parMFs
% at points in phiSwitch primary and secondary membership 
% need to be equal (=0.5)
n1 = 1;
n2 = parNphi(1,2);
% first (left) partition need to have membership 1 on left of core
for i = n1:n2
   if phi(i) <= parCenters(1)
       parMFs(1,i) = 1.0;
   end
end
n1 = parNphi(ncl,1);
n2 = n;
% last (right) partition need to have membership 1 on right of core
for i = n1:n2
   if phi(i) >= parCenters(ncl)
       parMFs(ncl,i) = 1.0;
   end
end
% intermediate partitions meet same membership at switching points
for k = 2:ncl
    n1 = find(phi >= parCenters(k-1),1,'first');
    n2 = find(phi <= parCenters(k),1,'last');
    phihalf = min(phiUmax(n1:n2));
    for i = n1:n2        
        if phi(i) <= phiSwitch(k-1)
           parMFs(k-1,i) = 0.5*(1+(phiUmax(i)-phihalf)/(1-phihalf));
           parMFs(k,i) = 1 - parMFs(k-1,i);
        elseif phi(i) <= parCenters(k)
           parMFs(k,i) = 0.5*(1+(phiUmax(i)-phihalf)/(1-phihalf));     
           parMFs(k-1,i) = 1 - parMFs(k,i);
        end
    end
end
parMFs(1,1) = 1.0;
parMFs(ncl,n)=1.0;

%
%% ***** Step 3:  Design MF
% Compute membership functions of final partition and final sets
% store all for further use
nMFpar = nAlpha; %501
nMFset = nAlpha;
parPoints = zeros(ncl,nMFpar);
parMFval = zeros(ncl,nMFpar);
setPoints = zeros(ncl,nMFset);
setMFval = zeros(ncl,nMFset);
%nAlpha = 2001;
%Alpha = linspace(0,1,nAlpha);
for k = 1:ncl
    phiL = parSupports(k,1);
    phiR = parSupports(k,2);
    phiPk = linspace(phiL,phiR,nMFpar-1);
    phiPk(nMFpar) = parCenters(k);
    phiPk = sort(phiPk,'ascend');
    % membership of fuzzy partition to be used for collocation
    n1 = parNphi(k,1);
    n2 = parNphi(k,2);
    xf = phi(n1:n2);
    valPk = MFval(parMFs(k,n1:n2),xf,phiPk);
    valPk(1) = 0.0;
    valPk(nMFpar) = 0.0;
    if k == 1
       valPk(1) = 1.0;
    end
    if k == ncl
        valPk(nMFpar) = 1.0;
    end
    valPk = valPk/max(valPk);
    parPoints(k,:) = phiPk;
    parMFval(k,:) = valPk;
        
    % membership of fuzzy sets in partition to be merged
    phiSk = linspace(phiL,phiR,nMFset-1);
    phiSk(nMFset) = parCenters(k);
    phiSk = sort(phiSk,'ascend');
    if ACFoption == 1
      % ACF-based transformation is used
      phidens = ksdensity(xf,phiSk,'Function','pdf');
      coreMF = parCenters(k);
      [tphihist,tphipos] = ecdf(phi(n1:n2));
      tphihist = tphihist/max(tphihist);
      starL = 1-ACFval(tphihist,tphipos,coreMF);
      eMF = min(tphihist/(1-starL),(1-tphihist)/starL);
      maxMF = max(eMF);
      eMF = eMF/maxMF; 
      valSk = MFval(eMF,tphipos,phiSk);
    else
      % AAT-based transformation
      phidens = ksdensity(xf,phiSk,'Function','pdf','Bandwidth',0.5,...
          'Kernel','epanechnikov','BoundaryCorrection','reflection');
      valSk = DAATrans(phidens);
      valSk = valSk/max(valSk);
      tvalSk = MFval(valPk,phiPk,phiSk)';
      valSk = valSk.*tvalSk;
    end
    %valSk(1) = 0.0;
    %valSk(nMFset) = 0.0;
    MvalSk = max(valSk);
    valSk = valSk/MvalSk;
    if galpha > 1.0
        if ARoption == +1 % reduction
          valSk = 1.0 - (1.0 - valSk).^galpha;
        elseif ARoption == -1 % amplification
          valSk = valSk.^galpha;
        end
    end
    setPoints(k,:) = phiSk;
    setMFval(k,:) =  valSk;
end    
parMFval(1,1) = 1.0;
parMFval(ncl,nMFpar) = 1.0;
setMFval(1,1) = 0.0;
setMFval(ncl,nMFset) = 0.0;

%

%% Save DATA
Out.ncl=ncl;
Out.phi = phi;
Out.phiU = phiU;
Out.phiCenters = phiCenters;
Out.indCenters = indCenters;
Out.parCenters=parCenters;
Out.parSupports=parSupports;
Out.nMFpar=nMFpar;
Out.parPoints=parPoints;
Out.parMFval=parMFval;
Out.nMFset=nMFset;
Out.setPoints=setPoints;
Out.setMFval=setMFval;
Out.parNphi=parNphi;
Out.parMFs=parMFs;

end




%%
% ................................................................
function q = DAATrans(p)
% Membership from probability 
% by (Direct) Arising Accumulation Transformation
 [sp,isp] = sort(p,'descend');
 m = length(p);
 q = zeros(1,m);
 for i = 1:m
    ii = isp(i); 
    sqqi = i*sp(i);
    for j = (i+1):m
        sqqi = sqqi + sp(j);
    end
    q(ii) = sqqi;
 end
end %end function




%%
% ................................................................
function ACFv = ACFval(acf,xacf,x)
% Compute the AC function at pointx x by linear interpolation
% Input points of ACF are (xacf,acf), xacf assumed increasing
% ACFval is computed at point x, i.e., (x,ACFv)
nacf = length(xacf);
nx = length(x);
ACFv = zeros(nx,1);
for k = 1:nx
    xk = x(k);
    if xk <= xacf(1)
        ACFv(k) = 0.0;
    elseif xk >= xacf(nacf)
        ACFv(k) = 1.0;
    else
        kp = find(xk <= xacf,1,'first');
        if xk == xacf(kp)
            ACFv(k) = acf(kp);
        else
            km = max(1,kp-1);
            fUp = acf(kp);
            fUm = acf(km);
            ACFv(k) = fUm + (xk - xacf(km))*(fUp-fUm)/(xacf(kp)-xacf(km));
        end
    end
end
end %end function

%%
% ................................................................
function MFv = MFval(MF,xf,x)
% Compute the MF function at points x by linear interpolation
% Input points of MF are (xf,MF), xf assumed increasing
% MFval is computed at points x, i.e., (x,MFv)
nmf = length(xf);
nx = length(x);
MFv = zeros(nx,1);
for k = 1:nx
    xk = x(k);
    if xk <= xf(1)
        MFv(k) = MF(1);
    elseif xk >= xf(nmf)
        MFv(k) = MF(nmf);
    else
        kp = find(xk <= xf,1,'first');
        if isempty(kp) ; kp = 1; end
        if xk == xf(kp)
            MFv(k) = MF(kp);
        else
            km = max(1,kp-1);
            fUp = MF(kp);
            fUm = MF(km);
            MFv(k) = fUm + (xk - xf(km))*(fUp-fUm)/(xf(kp)-xf(km));
        end
    end
end
end %end function

%%

%%


