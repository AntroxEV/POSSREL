function Out = cluFuzzyMembership(CuTest,AA)
%% INPUT DATA
% Prof Luciano Stefanini - University of  Urbino
% modified by Dr Alessandro Tombari - University of Exeter 
% Please cite: https://doi.org/10.1016/j.compgeo.2023.105967

% Locate phiTest value into fuzzy partition AA
% and compute the fuzzy membership function from collocation
% phiTest is the value of phi to test (assumed single element)
 mphi = CuTest(1);
 % nm = 1;

%% ***** Load partition to locate measurements
 ncl = AA.ncl;
 parCenters = AA.parCenters;
 parSupports = AA.parSupports;
 nMFpar = AA.nMFpar;
 parPoints = AA.parPoints;
 parMFval = AA.parMFval;
 nMFset = AA.nMFset;
 setPoints = AA.setPoints;
 setMFval = AA.setMFval;
% Define levels for alpha-cuts
 nAlpha = nMFpar; %1001
 Alpha = linspace(0,1,nAlpha);    
 
%% Collocation phase
 %find membership values in left and right partitions    
 kR = find(mphi <= parCenters,1,'first');
 if isempty(kR); kR = ncl; end
 kL = find(mphi >= parCenters,1,'last');
 if isempty(kL); kL = 1; end
 mkL = kL;
 mkR = kR;
 parkL = parMFval(kL,:);
 parkR = parMFval(kR,:);
 pointL = parPoints(kL,:);
 pointR = parPoints(kR,:);
 parLalpha = MFval(parkL,pointL,mphi);
 parRalpha = MFval(parkR,pointR,mphi);

% Merge left and right cluster membership functions
 nPts = nAlpha; %1001
 kL = mkL;
 kR = mkR;
 Pts = linspace(parSupports(kL,1),parSupports(kR,2),nPts);
 if kL < kR
    inp1 = find( Pts < parCenters(kL) );
    inp2 = find((Pts >=  parCenters(kL)) & (Pts <= parSupports(kL,2)));
    inp3 = find(Pts > parCenters(kR));
    v1 = MFval(setMFval(kL,:),setPoints(kL,:),Pts(inp1));
    v1 = v1 * parLalpha;
    v2L = MFval(setMFval(kL,:),setPoints(kL,:),Pts(inp2));
    v2R = MFval(setMFval(kR,:),setPoints(kR,:),Pts(inp2));
    v3 = MFval(setMFval(kR,:),setPoints(kR,:),Pts(inp3));
    v3 = v3 * parRalpha;
    v2 = v2L*parLalpha + v2R*parRalpha;
    vPts = [v1',v2',v3'];
    vPts = vPts' / max(vPts);
 else
    vPts = MFval(setMFval(kL,:),setPoints(kL,:),Pts);
    vPts = vPts/max(vPts);
 end

% compute merged MF of phi-values
 tphimf = vPts;   
 xphidens = Pts;
 tphimf(1) = 0.0;
 tphimf(nPts) = 0.0;
 Mtphimf = max(tphimf);
 tphimf = tphimf/Mtphimf;
 Mphi = max(tphimf);
 icore = find(tphimf >= Mphi,1,'first');%deve esserci
 CoreL = xphidens(icore);
 icore = find(tphimf >= Mphi,1,'last');%deve esserci
 CoreR = xphidens(icore);
 phiACF = ACFfromMF(xphidens,tphimf,CoreL,CoreR);
 SuppL = min(xphidens);
 SuppR = max(xphidens);
% [phiUm,phiUp] = AlphaCutsACF(phiACF,xphidens,nAlpha,Alpha);
% phiMF = MFAlphaCuts(phiUm,phiUp,Alpha,xphidens);
 phiMF = 2*min(phiACF,1-phiACF);
 phiMF = phiMF/max(phiMF);
 [phiUm,phiUp] = AlphaCutsMF(phiMF,xphidens,nAlpha,Alpha);
 
 nMFfull = 11; %1001
 fullPoints = linspace(SuppL,SuppR,nMFfull);
 fullSupport = [SuppL,SuppR];
 fullMFval = MFval(phiMF,xphidens,fullPoints);
 
 fullMFval = fullMFval/max(fullMFval);
 

%% Save DATA 
 Out.fullSupport=fullSupport;
 Out.fullPoints=fullPoints;
 Out.nMFfull=nMFfull;
 Out.fullMFval=fullMFval;
 Out.nAlpha = nAlpha;
 Out.Alpha = Alpha;
 Out.phiUm = phiUm;
 Out.phiUp = phiUp;
 
end



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
% ................................................................
function [fUm,fUp] = AlphaCutsMF(mfU,xmf,Nalpha,Alpha)
fUm = zeros(1,Nalpha);
fUp = zeros(1,Nalpha);
nx = length(xmf);
for ia = 1:Nalpha
    alpha = Alpha(ia);
    if alpha == 0.0
        fUm(ia) = xmf(1);
        fUp(ia) = xmf(nx);
    else
        km = find(mfU >= alpha,1,'first');
        kp = find(mfU >= alpha,1,'last');
        if isempty(km); km = 1; end
        if isempty(kp); kp = nx; end
        if km <= 1; km = 1;end
        if kp >= nx; kp = nx;end
        fUm(ia) = xmf(km);
        fUp(ia) = xmf(kp);
        if alpha == 1
            fUm(ia) = (xmf(km)+xmf(kp))/2;
            fUp(ia) = fUm(ia);
        end
    end
end
end %end function

%%


%%
% ................................................................
function ACF = ACFfromMF(x,q,c,d)
% ACF from mf (x,q) at points x
% Interval [c,d] is the core
 m = length(x);
 ACF = zeros(1,m);
 for i = 1:m
    ACF(i) = 0.5;
    if x(i)< c
        ACF(i) = 0.5*q(i);
    elseif x(i) > d
        ACF(i) = 1.0 - 0.5*q(i);
    end
 end
 ACF(1) = 0.0;
 ACF(m) = 1.0;
end %end function
