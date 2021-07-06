%==========================================================================
%% SCRIPT_NAME
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home  = '~';
cd(P.home);
P.fun   = [P.home filesep 'FUN'];
P.dat   = [P.home filesep 'DAT'];
P.mat   = [P.home filesep 'MAT'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;
%------------------------------------------------








%==========================================================================
%% LOAD PARAMETER DATA
%==========================================================================
clc; clearvars -except P LOO





% LBR
%-------------------------------------
PAR.Bon                 = 2;
PAR.Ron                 = 2;
PAR.Boff                = 5;
PAR.Roff                = 1;
PAR.Lon                 = 16;
PAR.Loff                = 4;




% DOS
%-------------------------------------
PAR.doSSReset           = 1;
PAR.saveTipMx           = 1;
PAR.delOrigFils         = 1;
PAR.doPlot              = 1;
PAR.doActCounts         = 1;
PAR.doLivePlot          = 0;
PAR.doThymT             = 1;
PAR.doActT              = 0;
PAR.doArpT              = 0;
PAR.doRev               = 1;
PAR.doFluorPlot         = 0;
PAR.doActLTP            = 0;
PAR.doKez               = 0;
PAR.doLiveHullPlot      = 0;
PAR.runTest             = 0;
PAR.doClustering        = 0;
PAR.doNdel              = 0;







% TIME PARAMETERS
%-------------------------------------
PAR.dt                  = .06;
PAR.StepsPerMin         = 60/PAR.dt;
PAR.StepsPerHr          = 60*60/PAR.dt;
PAR.hours               = 4;
PAR.LTPhours            = 2;
PAR.SSReStep            = 20000;
PAR.nDelSteps           = 20000;
PAR.delOrigFilsT        = 9000;
PAR.Nsteps              = PAR.StepsPerHr * PAR.hours + PAR.SSReStep;
PAR.InfTime             = PAR.Nsteps + 1;
PAR.FluorT              = 1000;
PAR.LivePlotMod         = 1000;
PAR.LoopPrintMod        = 1000;
PAR.SaveTipsAfter       = PAR.SSReStep;
PAR.SaveTipsRate        = 1000;
PAR.ThymTimeStart       = 12000;
PAR.timeLTPon           = PAR.StepsPerHr * PAR.LTPhours + PAR.SSReStep;
PAR.timeLTPoff          = PAR.timeLTPon + PAR.StepsPerMin * 3;
PAR.DataR               = 10;
PAR.ViewT               = 10;
PAR.PauseT              = 0;
PAR.Nstp                = 3600;
PAR.tdT                 = .00014;
PAR.LoopNum             = 20;


% THYMOSIN
%-------------------------------------
PAR.uM_Thymosin         = 250;
PAR.uM_ThymosinActin    = 125;
PAR.Ka_Thymosin         = .08;
PAR.Kd_Thymosin         = .04;
PAR.Ka_Thymosin_basal   = .08;
PAR.Ka_Thymosin_LTP     = .01;
PAR.Kd_Thymosin_basal   = .04;
PAR.Kd_Thymosin_LTP     = .1;
PAR.D_Thymosin          = 1e-3;




% ARP
%-------------------------------------
PAR.ArpBR               = 1;
PAR.ARPmax              = 1000; % max arp branches
PAR.ArpAdd              = 5;
PAR.ArpON               = 10;
PAR.ArpOFF              = 1;
PAR.GArpN               = 9;


% ACTIN
%-------------------------------------
PAR.StartAct            = 500;
PAR.ActUpdate           = 10;
PAR.BranchAng           = 45;
PAR.molAct              = 300;
PAR.KaBE                = 10;   % Barbed End On-Rate Scalar
PAR.KdBE                = 1;    % Barbed End Off-Rate Scalar
PAR.KaPE                = 1;    % Pointed End On-Rate Scalar
PAR.KdPE                = .5;   % Pointed End Off-Rate Scalar
PAR.StartMonos          = 20;
PAR.nStartFils          = 3;

PAR.BeKaATP             = 11.6;
PAR.BeKaADPP            = 3.4;
PAR.BeKaADP             = 2.9;
PAR.BeKdATP             = 1.4;
PAR.BeKdADPP            = .2;
PAR.BeKdADP             = 5.4;
PAR.PeKaATP             = 1.3;
PAR.PeKaADPP            = .11;
PAR.PeKaADP             = .09;
PAR.PeKdATP             = .9;
PAR.PeKdADPP            = .02;
PAR.PeKdADP             = .25;




% DENDRITE
%-------------------------------------
PAR.SPYneckXY           = 300;
PAR.SPYheadZN           = 1000;
PAR.SPYheadZS           = 700;
PAR.SPYheadX            = 500;
PAR.PSDproxy            = 25;
PAR.PSDproxyXY          = 3;
PAR.fZo                 = 110;
PAR.fZa                 = 10;
PAR.fXYo                = 5;
PAR.fXYa                = 2;
PAR.mempressure         = 10;
PAR.SPYneckRad          = PAR.SPYneckXY/2;
PAR.SPYheadRad          = PAR.SPYheadX/2;
PAR.PSDsz               = 20;
PAR.PSAsz               = 10;




% COFILIN
%-------------------------------------
PAR.CofR                = .04;
PAR.CofS                = 50;
PAR.CofSMax             = 2;
PAR.ACTdelCofN          = 15;
PAR.uM0_Cofilin         = 1;
PAR.Ka_Cofilin          = 1e-3;
PAR.Ka_Cofilin_LTP      = 1e-3;
PAR.Ka_Cofilin_LTD      = .01;



% AMX
%-------------------------------------
PAR.loadActinTips       = 0;
PAR.generateActinTips   = 1;
PAR.SaveFilePrefix      = 'ACTIN_REC-';
PAR.ATfilename          = 'ATdata.mat';
PAR.AMask               = 1;
PAR.SMask               = 10;
PAR.delOrate            = .1;
PAR.TriHullMod          = 200;
PAR.Rev                 = 10;















%==========================================================================
%% ACTIN MULTIPLEX TIP DATA
%==========================================================================
clc; clearvars -except P LOO PAR


BTs             = [];
AFMx            = [];
Nsteps          = PAR.Nsteps;
dt              = PAR.dt;



% ANGLE SETUP
%-------------------------------------
% unitsratio('rad', 'deg')
d2r = 1/(180/pi);
Ov = [2 29 57 85 112 140 168 195 223 251 278 306 334];
Ov = [Ov -15 -40 -65 -100 -125 -150 -175 -205 -230 -260 -295 -320 -350];





% Spine Dimensions
%-------------------------------------
SPYneckRad = PAR.SPYneckRad;	% Spy neck XY radius
SPYheadRad = PAR.SPYheadRad;	% Spy head XY radius
SPYheadZN  = PAR.SPYheadZN;     % Spy head north
SPYheadZS  = PAR.SPYheadZS;     % Spy head south


PSDproxy    = PAR.PSDproxy;
PSDproxyXY  = PAR.PSDproxyXY;
inPSD       = SPYheadZN - PSDproxy;

SPYH = [SPYheadRad SPYheadRad];
AcMx = zeros(SPYheadRad*2/5,SPYheadRad*2/5);

dims = [SPYneckRad SPYheadZN SPYheadZS SPYheadRad SPYheadRad ...
        PSDproxy inPSD PSDproxyXY];




%==========================================================================
%% INITIALIZE AND TAG STARTING FILAMENTS
%==========================================================================
%                           ACT(nStartFils,20)
%--------------------------------------------------------------------------
% N Xa Xo Xt Ya Yo Yt Za Zo Zt Mom ID Cut Born Died Life maxFa muFa NA nmL
% 1 2  3  4  5  6  7  8  9  10 11  12 13  14   15    16   17    18  19 20
%--------------------------------------------------------------------------
                

nStartFils = PAR.nStartFils;

ACT = zeros(PAR.nStartFils,20);  






% Starting Length Loc & Angles
StartMonos  = PAR.StartMonos;
fXYo        = PAR.fXYo;
fZo         = PAR.fZo;
fXYa        = PAR.fXYa;
fZa         = PAR.fZa;

% Branching Angles
TPi = d2r*PAR.BranchAng;
PPi = 0;


% MAKES STARTING FILAMENTS
ACT = MakeStartFils(ACT,nStartFils,StartMonos,...
    d2r,fZa,fXYo,SPYheadZN,SPYheadZS,SPYheadRad);



TagN        = numel(ACT(:,1));
ACT(:,12) = 1:TagN;
TagN        = TagN+1;





%==========================================================================
%% ACTIN SIZES
%==========================================================================
% NOTES
%{
From: Andre Kamkin, Irina Kiseleva
Springer Science & Business Media, Nov 18, 2010 - Biochemistry - 395 pages
Chapter Authors: Luo and Robinson
2.2 Microstructures and Deformations of the ACT Cyctoskeleton
- Gactin is 5 nm in diameter
- Factin filaments are 8 nm wide with a left-handed helical morphology
- 13 actin monomers per pseudo-repeat
- 1 pseudo-repeat length of 37 nm 
- Alternatively the actin filament can be considered to have a right-handed 
helical structure with two strands slowly twisting around each other. 
Each actin monomer is rotated 166 degrees with respect to its nearest neighbors 
across the strand (Holmes 1990). Within the strand, subdomains 2 and 4 contact 
subdomains 1 and 3 in the next monomer in the strand, and each monomer reaches 
across to the other strand through a hydrophobic plug that links the two 
strands together. 

37 nm / 13 p = 2.85 nm/p
13 p / 37 nm = 0.35 p/nm

From: Lodish, Principles of Molecular Biology
Adapted from C. E. Schutt et al. 1993, Nature 365:810
- There are 2 strands with 14 units per strand, 
- so 28 monomers in each 360 degree turn,
- with a length of 72 nm per 360 degree turn

72 nm / 28 p = 2.57 nm/p
28 p / 72 nm = 0.39 p/nm

Factin particle size in filaments:
0.369 p/nm   (Factin particles per nm of filament)
2.71  nm/p    (nm of filament per Factin particle)

%}
%-------



p_per_nm = 0.369;    % (Factin particles per nm of filament)
nmpp = 2.71;     % (nm of filament per Factin particle)
ACT(:,20) = ACT(:,1) .* nmpp;   % Store filament length





% TRIG: branch XYZ tip coordinates
%-------------------------------------
ACT(:,4)  = ACT(:,1) .* nmpp .* sin(ACT(:,8)) .* cos(ACT(:,2)) + ACT(:,3);
ACT(:,7)  = ACT(:,1) .* nmpp .* sin(ACT(:,8)) .* sin(ACT(:,2)) + ACT(:,6);
ACT(:,10) = ACT(:,1) .* nmpp .* cos(ACT(:,8)) + ACT(:,9);




DEDACT      = ACT;
oActin      = ACT;







%==========================================================================
%% Conversion Factors & VCP Function
%==========================================================================
% >> pNuM = VCP(vol,uM,pN)
%-------------------------
% NOTES
%{
 
% [pNuM] = VCP(vol,uM,pN);
% INPUTS
% vol: volume in um^3
% uM: concentration in uM
% pN: particle count
% 
% enter zero for the unknown value
% if pN is unknown enter... pN = VCP(.1,10,0)
% if uM is unknown enter... uM = VCP(.1,0,6e5)
 
%}
%-------

MOL = 6.022e23;     % 1 mol Avagadro's number
mol = 6e23;			% 1 mol rounded




%==========================================================================
%% SPINE VOLUME
%==========================================================================
% NOTES
%{
%---------------------------------------------------
% A typical spine volume is 0.1 um^3 or 1e-16 L
% volume of cylinder
% V = pi * r^2 * h
%---------------------------------------------------


% MATH - Spine Volume

This spine has a neck volume equivalent to a cylendar of dimensions: 
50d x 200

And a head volume equivalent to a 3D disk with dimensions:
100d x 100


% Units
mol = 6e23;		% in N
mM = 1e-3;
uM = 1e-6;
upM = 1e3;
dnM = 1e-3;

% 1 cm^3 = 1 mL
% SI units conversion from m^3 to L are:
% 1 cm^3 = 1 mL
% To convert cubic volume to mL you may first need to convert
% by an order of magnitude; here's a reminder of cubic volume conversions:
% 1 m^3 = 1 cm^3 * .01^3;		% 1e-6
% 1 m^3 = 1 mm^3 * .001^3;		% 1e-9
% 1 cm^3 = 1 mm^3 * .1^3;		% 1e-3
% 1 cm^3 = 1 um^3 * .0001^3;	% 1e-12
% Thus, if dendritic spines have an average volume of 0.1 um^3
% that would be equivalent to X uL
%
% 0.1 um^3 * (1 cm^3 / 1e12 um^3) * (1 mL / 1 cm^3) * (1000 uL / 1 mL)
% 0.1*(1/1e12)*(1/1)*(1000/1)
% 0.1 um^3 = .1e-12 cm^3 = .1e-9 uL
%
% and 0.1 um^3 equivalent to X L
% 
% 0.1*(1/1e12)*(1/1)*(1/1000)
% 1e-16


% Spine volume (0.1 um^3) in L and uL
SpyV = 1e-16;	% in L
SpyVu = .1e-9;	% in uL

% ACT Polymerization (12 p/µM*s)
% Act_Na = 1e5;
% BeKanT = 1000;

% ACT Depolymerization (2 p/s)	
BeKdnT = 200;
DePSum = 0;

% ACT Poly Math
% Cytosolic concentration of actin in cells ranges from .1 to .5 mM
% Given a spine volume of 'SpyV' we can find how many actin monomers
% are in an average spine:
%
% .1 mM (mmol/L) * (1 mol / 1000 mmol) * SpyV = 1e-17 mol
% 1e-17 mol * (6e23 units / 1 mol) = 6e3 monomer units

molAct = AMX{16};
Act_N = molAct * (1/1000) * SpyV * 6e23; % 6e3 monomer units

% we can check our math starting with a set number of actin monomers
% and calculate the spine molarity (6e3 monomer units as an example):
% 
% 6e3 units/SpyV * (1 mol / 6e23 units) * (1000 mmol / 1 mol)
% 6e3/SpyV*(1/6e23) 

Gactin_uM = Act_N/SpyV * (1/6e23) * (1000/1);	% 1.6e-10 

% Gactin_uM = Act_N / SpyVu / mol;	% 1.6e-10 
% Act_N = Act_N / SpyVu / mol;
BeKa = 12 * Gactin_uM * dt;
BeKd = 2 * dt;

volume of cylinder
V = pi * r^2 * h

volume of a sphere
V = (4*pi*r^3)/3



%}
%-------


Vneck = pi * SPYneckRad^2 * SPYheadZS;
Vhead = pi * SPYheadRad^2 * (SPYheadZN-SPYheadZS);
SpyV = (Vneck+Vhead) * 1e-24;	% nm^3 to L conversion (nm: 1e-24; um: 1e-15)
%SpyV = 1e-16;





VuMOL = SpyV * MOL / 1e6;

uMOLV = 1e6 / MOL / SpyV;




%==========================================================================
%% ACTIN CONCENTRATIONS
%==========================================================================
% NOTES
%{

NA

%}
%-------

uM_Act      = PAR.molAct;                   % actin uM
GActinN0    = uM_Act / 1e6 * MOL * SpyV;    % N actin monomers at t0 (6e4)
Gactin_uM   = GActinN0 / SpyV *(1/MOL)*1e6;	% check uM_Act == Gactin_uM

FActinN     = sum(ACT(:,1));              % current number FActins
GActinN     = ceil(GActinN0) - FActinN;     % current number GActins

Nfi         = numel(ACT(:,1));            % current number of branches


Factin_uM   = FActinN / SpyV *(1/MOL)*1e6;	% FActin uM






%==========================================================================
%% ACTIN POLYMERIZATION & DEPOLYMERIZATION RATES
%==========================================================================
% NOTES
%{
%----------------------------------------------------
% Pollard:	
% ACT+	ATP      ADPP       ADP		| MEAN
% BeKa: 	11.6     3.4        2.9		| 6.0 p/µM*s
% BeKd:		1.4      0.2		5.4		| 2.3 p/s
% PeKa:		1.3      0.11       0.09	| 0.5 p/µM*s
% PeKd: 	0.8      0.02       0.25	| 0.4 p/s
%           BeKaATP  BeKaADPP   BeKaADP
%           BeKdATP  BeKdADPP   BeKdADP
%           PeKaATP  PeKaADPP   PeKaADP
%           PeKdATP  PeKdADPP   PeKdADP
%----------------------------------------------------
%}
%-------

BeKaATP     = PAR.BeKaATP;
BeKaADPP    = PAR.BeKaADPP;  
BeKaADP     = PAR.BeKaADP;
BeKdATP     = PAR.BeKdATP;
BeKdADPP    = PAR.BeKdADPP;
BeKdADP     = PAR.BeKdADP;
PeKaATP     = PAR.PeKaATP;
PeKaADPP    = PAR.PeKaADPP;
PeKaADP     = PAR.PeKaADP;
PeKdATP     = PAR.PeKdATP;
PeKdADPP    = PAR.PeKdADPP;
PeKdADP     = PAR.PeKdADP;


BeKaATD = mean([BeKaATP  mean([BeKaADPP  BeKaADP])]);
BeKdATD = mean([BeKdATP  mean([BeKdADPP  BeKdADP])]);
PeKaATD = mean([PeKaATP  mean([PeKaADPP  PeKaADP])]);
PeKdATD = mean([PeKdATP  mean([PeKdADPP  PeKdADP])]);

KaATD = roundn(mean([BeKaATD PeKaATD]),-2);
KdATD = roundn(mean([BeKdATD PeKdATD]),-2);



% ACTIN ON-RATE & OFF-RATE SCALARS
%-------------------------------------
doKez   = PAR.doKez;
KaBE    = PAR.KaBE;                 % Barbed End On-Rate Scalar
KdBE    = PAR.KdBE;                 % Barbed End Off-Rate Scalar
KaPE    = PAR.KaPE;                 % Pointed End On-Rate Scalar
KdPE    = PAR.KdPE;                 % Pointed End Off-Rate Scalar
KaTIP   = (KaBE + KaPE) / 2;        % Tip Mean Ka On-Rate Scalar
KdTIP   = (KdBE + KdPE) / 2;        % Tip Mean Kd Off-Rate Scalar

BeKa    = KaBE * Gactin_uM * dt;	% Ka Barbed End ON-Rate
BeKd    = KdBE * dt;                % Kd Barbed End OFF-Rate
PeKa    = KaPE * Gactin_uM * dt;	% Ka Pointed End ON-Rate
PeKd    = KdPE * dt;                % Kd Pointed End OFF-Rate




% TIP Ka & Kd
%--------
if doKez
	TKa = KaTIP;
	TKd = KdTIP;
else
	TKa = KaATD;
	TKd = KdATD;
end
%---
fKa = TKa * Gactin_uM * dt;	% Ka Fil ON-Rate
fKd = TKd * dt;				% Kd Fil OFF-Rate




%==========================================================================
%% THYMOSIN ACTIN-SEQUESTERING VARIABLES
%==========================================================================
% THYMOSIN REACTION RATE NOTES
%{
The Law of Mass Action (*LMA) describes the rate at which chemicals collide
and interact to form different chemical combinations. When two different
chemicals can collide to form a dimer product, and the dimer can dissociate
reversably back into the individual component reactants, is described as:

        Ka>
T + A <---->  TA
       <Kd  

Where
    T : thymosin                (thymosin monomers)
    A : actin                   (Gactin monomers)
    TA: thymosin-actin          (thymosin-actin dimers)
    Ka: forward rate constant
    Kd: backward rate constant


The rate equations describing the CHANGE in molecular concentration per dt are 
[https://en.wikipedia.org/wiki/Rate_equation]...

TA/dt = Ka[T][A] - Kd[TA]      (LMA forward reaction: TA accumulation)
A/dt  = Ka[T][A] - Kd[TA]      (LMA reverse reaction: A accumulation)
T/dt  = Ka[T][A] - Kd[TA]      (LMA reverse reaction: T accumulation)

For many reactions the rate is given by a power law such as:

r = k * [A]^x * [B]^y

where [A] and [B] express the concentration of the species A and B, respectively (usually in moles
per liter (molarity, M)). The exponents x and y are the partial reaction orders and must be
determined experimentally; they are often not equal to the stoichiometric coefficients. The constant
k is the rate coefficient or rate constant of the reaction. For elementary reactions, which consist
of a single step, the order equals the molecularity as predicted by collision theory. For example, a
bimolecular elementary reaction A + B products will be second order overall and first order in
each reactant, with rate equation r = k * [A] * [B]. For multistep reactions, the
order of each step equals the molecularity, but this is not generally true for the overall rate.

%}


THYM_ACT_uM  = PAR.uM_ThymosinActin;                    % TA (uM)
THYM_uM      = PAR.uM_Thymosin;                         % T  (uM)

THYM_ACT_N   = round(THYM_ACT_uM / 1e6 * MOL * SpyV);   % TA (n)
THYM_N       = round(THYM_uM / 1e6 * MOL * SpyV);       % T  (n)


Ka_Thymosin0 = PAR.Ka_Thymosin;                         % Ka thymosin (t=0)
Kd_Thymosin0 = PAR.Kd_Thymosin;                         % Kd thymosin (t=0)
Ka_Thymosin  = PAR.Ka_Thymosin;                         % Ka thymosin
Kd_Thymosin  = PAR.Kd_Thymosin;                         % Kd thymosin


tKa = Ka_Thymosin * THYM_uM * Gactin_uM * dt;           % T+A (on-rate)
tKd = Kd_Thymosin * THYM_ACT_uM * dt;                   % T-A (off-rate)


% DIFFUSION RATE ESTIMATION
%-------------------------------------
% k = 1.38064852e-23; % Boltzmann Constant
% T = 310;            % approximately 98 Degrees Fahrenheit
% nu = 3;             % centipoise viscosity
% r = 3e-9;           % particle radius in meters (1 Angstrom == 1e-10 meters)
% 
% DR = k*T / (6*pi*nu*r) * (10^6^2); % microns squared per second
%-------------------------------------












%==========================================================================
%% ARP2/3 FILAMENT BRANCHING VARIABLES
%==========================================================================
% NOTES
%{
According to Smith & Gelles 2012:
Arp helps to nucleate new daughter filament branches 
at a rate of 2.5 to 9.7 df/(mM_Arp * s * ummf)
2.5 with no WASp and 9.7 with 300 nM WASp

At some intermediate WASp level, branching may proceed at
5 daughter filaments (df)
per mM of Arp
per second 
per micrometer of mother filament (ummf)

Example calculations with 9 uM Arp and 21600 nm of mf
2.5 df/(Arp_mM*s*ummf) * ((9/1000) * dt * (21600/1000))
9.7 df/(Arp_mM*s*ummf) * ((9/1000) * dt * (21600/1000))
5.0 df/(Arp_mM*s*ummf) * ((9/1000) * dt * (21600/1000))

5 * ((9/1000) * dt * (21600/1000))
5 * ((Arp_uM/1000) * dt * (nmmf/1000))



Friederich reports an Arp2/3 association rate of:
5.4e-4 um^-3 S^-1 (crossref from Carlsson 2004) but use a value of
1e-5 b/µM*s as a baseline parameter in their simulation software ActSimChem
ASRT:	1e-5 b/µM*s

These 3 values give branching rates that differ orders of magnitude:
5 * (9/1000) * 1	% .045
1e-5 * 9 * 1		% .00009
5.4e-4 * (9^3)		% 0.39

9000
But the Gelles and Carlsson value are fairly close when the umf is 9
5 * (9/1000) * 9	% 0.40
5.4e-4 * (9^3)		% 0.39

which would be the case if the Factin concentration was 51 uM
(9000/2.7) / 1e-16 *(1/6e23)*1e6

5.4e-4 * Arp_uM^3 * dt





% MATH - Arp Branching
%-------------------------------------
%Tnmmf = sum(ACT(:,20));				% Length of filaments (nm) Combined
%ArpBRT = Arp_Sc * ((Arp_uM/1000) .* (Tnmmf/1000));	% Arp rate scalar    (total)
%ArpN = 1e3;
%ArpOn = 5;
%ArpOff = 1;
%Arp_uM = ArpN / SpyVu / mol;	% 1.6 - 16 uM
%Arp_PR = ArpOn * Arp_uM * dt;
%Arp_DR = ArpOff * dt;
%ArpR = AMX{17}* dt;	% Arp activity rate
%ArpR0 = ArpR;			% Arp activity rate
%ArpScalar = AMX{21};	% Arp filament length scalar

%}
%-------------------------------------

Arp_Sc  = PAR.ArpBR * dt;				% Arp empirical branching rate scalar
Arp_uM  = PAR.GArpN;					% Arp uM
nmmf    = ACT(:,20);					% Length of filaments (nm)

FArpN  = nStartFils;					% #of F-Arp
GArpN  = Arp_uM / 1e6 * 6e23 * SpyV;	% #of G-Arp
GArpN0 = ceil(GArpN);					% #of G-Arp (starting)
Arp_uM = GArpN / SpyV *(1/6e23)*1e6;	% Check uM_Act == Gactin_uM

ArpAdd = PAR.ArpAdd;                    % Add X units to new branches
ARPmax = PAR.ARPmax;					% Maximum Arp branches



ArpBR  = Arp_Sc * ((Arp_uM/1000) .* (nmmf/1000)); % Arp branch rate per fil


GArpN_t0 = GArpN;










%==========================================================================
%% COFILIN & GENERAL DEPOLYMERIZATION VARIABLES
%==========================================================================
% 
% fdKd = fKd*10;            % Depoly rate when Fil is Fkd
% 
% CofR = AMX{18}* dt;       % cofilin activity rate
% CofS = AMX{21};           % cofilin delete Nunits
% CofN = AMX{30};           % cofilin delete Nunits
% delOr =AMX{14}* dt;       % delete from origin rate


uM0_Cofilin     = PAR.uM0_Cofilin;                  % Cofilin uM total
N0_Cofilin      = uM0_Cofilin / 1e6 * MOL * SpyV;	% Cofilin N spine count
Cofilin_uM      = N0_Cofilin / SpyV *(1/MOL)*1e6;	% Cofilin uM spine


Ka_Cofilin      = PAR.Ka_Cofilin;       % Cofilin Ka (basal)
Ka_Cofilin_LTP  = PAR.Ka_Cofilin_LTP;	% Cofilin Ka (LTP)
Ka_Cofilin_LTD  = PAR.Ka_Cofilin_LTD;	% Cofilin Ka (LTD)


KaCof = (dt * Ka_Cofilin * Factin_uM * Cofilin_uM); % Cof+Fact on-rate








%==========================================================================
%% GET TOTAL PROTEIN COUNTS IN CLOSED SYSTEM
%==========================================================================


% TOTAL COUNTS
%-------------------------------------
CofSMax         = PAR.CofSMax;

TheoMaxFact     = ceil(SPYheadZN / nmpp);

ScFil           = ceil(TheoMaxFact/CofSMax);

TOTAL_ACTIN     = GActinN + FActinN + THYM_ACT_N;

TOTAL_THYMOSIN  = THYM_N + THYM_ACT_N;

TOTAL_ARP       = GArpN + FArpN;

TOTAL_COFILIN   = N0_Cofilin;



% t0 TOTALS COUNTS
%-------------------------------------
TOTAL_ACTIN0    = TOTAL_ACTIN;

TOTAL_THYMOSIN0 = TOTAL_THYMOSIN;

TOTAL_ARP0      = TOTAL_ARP;

TOTAL_COFILIN0  = TOTAL_COFILIN;











%==========================================================================
%% PREALLOCATE COUNTERS
%==========================================================================


CNTS.nT         = (1:Nsteps)';
CNTS.Fils       = zeros(Nsteps,1);    % Number of branch filaments
CNTS.FAct       = zeros(Nsteps,1);    % Number FActins
CNTS.GAct       = zeros(Nsteps,1);    % Number GActins
CNTS.CofAct     = zeros(Nsteps,1);    % Number of Depoly events
CNTS.ArpAct     = zeros(Nsteps,1);    % Number of Poly events
CNTS.delFi      = zeros(Nsteps,1);    % Number of deleted filaments
CNTS.muFilFa    = zeros(Nsteps,1);    % Mean filament length
CNTS.ArpKa      = zeros(Nsteps,1);    % Arp Branch Rate
CNTS.Act_uM     = zeros(Nsteps,1);    % ACT uM
CNTS.fKa        = zeros(Nsteps,1);    % Filament Ka
CNTS.FnowFtot   = zeros(Nsteps,1);    % CurFilN : AllFilN
CNTS.FArp       = zeros(Nsteps,1);    % Number of Arp in Filaments
CNTS.GArp       = zeros(Nsteps,1);    % Number of Free Arp
CNTS.nmTotFil   = zeros(Nsteps,1);    % length of all fils (nm)
CNTS.DAct       = zeros(Nsteps,1);    % Number Unbound Actins
CNTS.Thy        = zeros(Nsteps,1);    % Number Unbound Thymosin
CNTS.ThyAct     = zeros(Nsteps,1);    % Number Thymosin+GActins
CNTS.ThyKa      = zeros(Nsteps,1);    % Ka Thymosin
CNTS.ThyKd      = zeros(Nsteps,1);    % Kd Thymosin









%==========================================================================
%% ANIMATED REAL-TIME FIGURE SETUP
%==========================================================================


if PAR.doLivePlot

[fh1,ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,rot0,rota,azel] =...
    actin_model_figure(Nsteps,dt,dims);

end



%==========================================================================
%% PRE-PLOT FILAMENTS
%==========================================================================

ActinTips   = [ACT(:,4) ACT(:,7) ACT(:,10)];
[Zrow1,~]   = find(ActinTips(:,3) > inPSD);
PSDTips     = ActinTips(Zrow1,:);
[Zrow2,~]   = find(ActinTips(:,3) < inPSD);
SPYTips     = ActinTips(Zrow2,:);

P3x         = [ACT(:,3) ACT(:,4)]';
P3y         = [ACT(:,6) ACT(:,7)]';
P3z         = [ACT(:,9) ACT(:,10)]';
P3x(3,:)    = NaN; 
P3y(3,:)    = NaN; 
P3z(3,:)    = NaN;



if PAR.doLivePlot
    phF = line('XData', P3x(:), 'YData', P3y(:), 'ZData', P3z(:),...
        'Parent', ax1, 'LineStyle', '-', 'Color', [.5 .5 .5],...
        'LineWidth', 1, 'Marker','.');

    phS = scatter3(ax1, SPYTips(:,1), SPYTips(:,2), SPYTips(:,3),...
        200, [.1 .1 .9],'Marker','.');

    phP = scatter3(ax1, PSDTips(:,1), PSDTips(:,2), PSDTips(:,3),...
        200, [.9 .1 .1],'Marker','.');

    view(ax1,azel+rota)

    axes(ax1); % make axes 'ax1' top-level plot
end



%==========================================================================
%% LTP
%==========================================================================






%==========================================================================
%%						MAIN OUTER LOOP
%==========================================================================
for nT = 1:Nsteps





    % TIMER SWITCHES
    %-------------------------------------
    THYMon      = nT > PAR.ThymTimeStart;
    ACTearly    = nT < 10000;
    LTPon       = nT==PAR.timeLTPon;
    LTPoff      = nT==PAR.timeLTPoff;





    % THYMOSIN LTP
    %-------------------------------------
	if PAR.doThymT
      if LTPon
        Ka_Thymosin = PAR.Ka_Thymosin_LTP;
        Kd_Thymosin = PAR.Kd_Thymosin_LTP;
      end
      if LTPoff
        Ka_Thymosin = PAR.Ka_Thymosin_basal;
        Kd_Thymosin = PAR.Kd_Thymosin_basal;
      end
	end
	



    % NFact: Current number of Filaments
    %-------------------------------------
    %  fN - used in inner loop below as "for fN=1:NFact..."
    %  Nfi - used in outer loop below as "Nfi = numel(ACT(:,1));"
    %  note: do not re-compute NFact anywhere else but here
    %  use Nfi in final counters since filaments are being added/deleted
    %  fKa & fKd are the on/off rates of F-actin Filaments
    %-------------------------------------
    % [important!] do not re-compute NFact anywhere else but here
	NFact = numel(ACT(:,1));	 




    % Resets
    %-------------------------------------
	ACTdepoly   = 0; 
    COFdepoly   = 0; 
    ARPpoly     = 0; 
    ACTpoly     = 0;
	
	


    % Filament coordinates
    %-------------------------------------
	Tip_xyz = [ACT(:,4) ACT(:,7) ACT(:,10)];

	% radial distance to spine shaft membrane
	XYtipLoc = sqrt(Tip_xyz(:,1).^2 + Tip_xyz(:,2).^2);

	ZtipInHead = Tip_xyz(:,3) >= SPYheadZS;
	ZtipInNeck = Tip_xyz(:,3) < SPYheadZS;

	XYneckOut = XYtipLoc > SPYneckRad;	% Logic idx fils beyond SPYneckRad
	XYheadOut = XYtipLoc > SPYheadRad;	% Logic idx fils beyond SPYheadRad
	ZtopOut = Tip_xyz(:,3) > SPYheadZN;	% Logic idx fils above SPYheadZN
	ZbotOut = Tip_xyz(:,3) < 0;			% Logic idx fils below zero
	

	TipOut=((XYneckOut & ZtipInNeck)+(XYheadOut & ZtipInHead)+ZtopOut+ZbotOut)>0;
	TipOK = ~TipOut;
	
    
	LngXYneckOut = (XYtipLoc-SPYneckRad).*ZtipInNeck;  % dist beyond SPYneckRad
	LngXYheadOut = (XYtipLoc-SPYheadRad).*ZtipInHead;  % dist beyond SPYheadRad
	LngZtopOut = (Tip_xyz(:,3) - SPYheadZN).*ZtipInHead; % dist above SPYheadZN
    
    LngXYneckOut(LngXYneckOut<0) = 0;
    LngXYheadOut(LngXYheadOut<0) = 0;
    LngZtopOut(LngZtopOut<0) = 0;
    
    LngOut = LngXYneckOut + LngXYheadOut + LngZtopOut;
    
    LO = LngOut ./ 100 .* PAR.mempressure;    
    

    
    

    ACT(:,20) = ACT(:,1) .* nmpp;
	

    % assure no negative actin values
    %-------------------------------------
	ACT(:,1) = ACT(:,1) .* (ACT(:,1)>0);
	ACT(:,1) = ACT(:,1) .* (ACT(:,10) > -20);

	


	% MCMC requires a random filament at each decision step
    %-------------------------------------
	rN1 = randi(NFact,1,NFact); % actin polymerization
    rN2 = randi(NFact,1,NFact); % arp branching
    rN3 = randi(NFact,1,NFact); % actin depolymerization
    rN4 = randi(NFact,1,NFact); % cofilin severing





	%======================================================================
	%						MAIN INNER LOOP
	for fN=1:NFact
	%----------------------------------------------------------------------
	
        rv=rand(9,1); % Generate a few random vaules from uniform{0:1}


        
		


		% POLYMERIZATION
        %-------------------------------------
        aN = rN1(fN); % Get random filament

        if ACTearly
            if  (fKa > rv(1)) && TipOK(aN) && (GActinN>1)
                if fKa > NFact
                    ACT(aN,1) = ACT(aN,1) + 2;
                    ACTpoly = ACTpoly + 2;
                else
                    ACT(aN,1) = ACT(aN,1) + 1;
                    ACTpoly = ACTpoly + 1;
                end
            end
        else
            if  (fKa-LO(aN) > rv(1)) && ~ZbotOut(aN) && (GActinN>1)
                if fKa > NFact
                    ACT(aN,1) = ACT(aN,1) + 2;
                    ACTpoly = ACTpoly + 2;
                else
                    ACT(aN,1) = ACT(aN,1) + 1;
                    ACTpoly = ACTpoly + 1;
                end
            end		
        end
		
		



		

		% BRANCHING
        %-------------------------------------
        aN = rN2(fN);

 		if  (ArpBR(aN) > rv(2)) && (GActinN>ArpAdd) && (size(ACT,1)<ARPmax)

			fNf = NFact+1;
			
			ACT(fNf,1)  = ArpAdd;     % create branch: add N actin subunits
			ACT(fNf,11) = ACT(aN,12); % tag branch with MomID
			ACT(fNf,12) = TagN;       % tag branch with ID
			ACT(fNf,14) = nT;         % tag branch with Born time (nT)
			TagN = TagN+1;
			
			Nmono = ACT(aN,1);		  % N monomers in mother filament
			Rmm = ceil(Nmono * rand); % Random mother monomer along segment
						
			% TIP of current branch		(not actual tip, just rotational point) 
			Ct_x = (Rmm) .* nmpp * sin(ACT(aN,8)) * cos(ACT(aN,2)) + ACT(aN,3);
			Ct_y = (Rmm) .* nmpp * sin(ACT(aN,8)) * sin(ACT(aN,2)) + ACT(aN,6);
			Ct_z = (Rmm) .* nmpp * cos(ACT(aN,8)) + ACT(aN,9);

			% ORIGIN of current branch
			Co_x = ACT(aN,3);	% X origin (old branch)
			Co_y = ACT(aN,6);	% Y origin (current branch)
			Co_z = ACT(aN,9);	% Z origin (current branch)
			
				
			Po = [Co_x;Co_y;Co_z];
			Pt = [Ct_x;Ct_y;Ct_z];
			Pv = Pt-Po;
			
			tL = sqrt(sum((Pv).^2));		% Length of vector PoPt (aka Pv)
			Pu = (Pv) ./ tL;				% Unit vector of Pv

			tTheta = acos(Pu(3))  +TPi;		% angle theta
			tPhi = atan2(Pu(2),Pu(1)) +PPi;	% angle phi

			x = sin(tTheta) * cos(tPhi);
			y = sin(tTheta) * sin(tPhi);
			z = cos(tTheta);
			Pr = [x;y;z]+Pt;
			
			
			Otta = Ov(randi(26,1));		% Random rotational angle (theta)
			
			%Pn = RotateVertex(Pr(1),Pr(2),Pr(3),...
            %                  Pt(1),Pt(2),Pt(3),...
            %                  Po(1),Po(2),Po(3),...
			%			       Pv(1),Pv(2),Pv(3),Otta);

            % tictoc: 23.2298   24.2325
            dv = [Pt(1);Pt(2);Pt(3)] - [Po(1);Po(2);Po(3)];
            u=dv(1);  v=dv(2);  w=dv(3);
            Pn = RotateVec(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),u,v,w,Otta);
            % Pn = LineRota(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),u,v,w,Otta);
            %--------
			

			Po2 = Pt;
			Pt2 = Pn;
			Pv2 = Pt2-Po2;
			tL2 = sqrt(sum((Pv2).^2));		% Length of vector PoPt (aka Pv)
			Pu2 = (Pv2) ./ tL2;				% Unit vector of Pv
			tTheta2 = acos(Pu2(3));			% angle theta	+TPi;
			tPhi2 = atan2(Pu2(2),Pu2(1));	% angle phi		+PPi;
			
			
            
			% New branch Angles
            %-------------------------------------
			ACT(fNf,2) = tPhi2;
			ACT(fNf,8) = tTheta2;
			

			% New branch Origin
            %-------------------------------------
			ACT(fNf,3) = Po2(1);
			ACT(fNf,6) = Po2(2);
			ACT(fNf,9) = Po2(3);
			

			% New branch Tip
            %-------------------------------------
			ACT(fNf,4) = (ArpAdd) .* nmpp * sin(tTheta2) * cos(tPhi2) + Po2(1);
			ACT(fNf,7) = (ArpAdd) .* nmpp * sin(tTheta2) * sin(tPhi2) + Po2(2);
			ACT(fNf,10) = (ArpAdd) .* nmpp * cos(tTheta2) + Po2(3);
			
            
			ARPpoly = ARPpoly + ArpAdd;
			
		end
		

		
	

		% ACTIN PASSIVE DEPOLYMERIZATION
        %-------------------------------------
        aN = rN3(fN);

		if (fKd > rv(3))

			ACT(aN,1) = ACT(aN,1)-1;
			ACTdepoly = ACTdepoly + 1;

		end

		
		


		% COFILIN ASSISTED DEPOLYMERIZATION
        %-------------------------------------
        aN = rN4(fN);

		if  KaCof > rv(4) && ACT(aN,1) > 10
            
            % cofilin binds random Factin along filament
            ActBoundCof = randi(ACT(aN,1));    

            % depoly all Factin beyond cofilin site
            ActDelCof = ACT(aN,1)-ActBoundCof; 

            % remaining Factin count == cofilin site
            ACT(aN,1) = ActBoundCof;           
			
            % keep tally of cofilin depolymerized Factin
			COFdepoly = COFdepoly + ActDelCof;   

            % add depoly actin to Gactin pool (temporary)
            GActinN = GActinN + ActDelCof;       

		end



		% THYMOSIN
		%-------------------------------------
        if THYMon

            % THYMOSIN+ACTIN ASSOCIATION
            if (tKa > rv(8)) && (GActinN>1) && (THYM_N>1)

                THYM_N = THYM_N - 1;
                THYM_ACT_N = THYM_ACT_N + 1;
                GActinN = GActinN - 1;

            end

            % THYMOSIN+ACTIN DISSOCIATION
            if (tKd > rv(9))  && (THYM_ACT_N>1)

                THYM_N = THYM_N + 1;
                THYM_ACT_N = THYM_ACT_N - 1;
                GActinN = GActinN + 1;

            end
        end

		
	


        % ADJUST RATE VALUES
        %-------------------------------------
        % uM = p / SpyV * (1/MOL) * 1e6;
        %-------
        % thymosin Ka/Kd rates scaled by dt/NFact as inner loop 
        % repeats once per filament per dt
        % https://en.wikipedia.org/wiki/Rate_equation
        %-------------------------------------
        if ~mod(fN,10)

            % adjust dt by N inner loops being run
            dtN = dt/NFact; 

            

            % ACTIN VALUES
            FActinN = sum(ACT(:,1));                             % Factins

            TA_plus_Factin_N = FActinN + THYM_ACT_N;             % thymact+fact

            GActinN = TOTAL_ACTIN - TA_plus_Factin_N;            % Gactins

            Factin_uM = FActinN / SpyV * (1/MOL) * 1e6;          % Factin uM
            Gactin_uM = GActinN / SpyV * (1/MOL) * 1e6;          % Gactin uM

            fKa = TKa * Gactin_uM * dt;                          % fil on-rate


            % ARP VALUES
            FArpN = numel(ACT(:,1));                             % Farp count
            GArpN = GArpN0 - FArpN;                              % Garp count
            Arp_uM = GArpN / SpyV * (1/MOL) * 1e6;               % Garp uM
            nmmf = ACT(:,20);                                    % fil leng (nm)
            Fact_Branch_uM = ACT(:,1) ./ SpyV * (1/MOL) * 1e6;   % filament uM
            ArpBR = Arp_Sc * (Arp_uM/1000 .* Fact_Branch_uM);    % branch rate


            % THYMOSIN VALUES
            THYM_uM     = THYM_N     / SpyV * (1/MOL) * 1e6;     % thymosin uM
            THYM_ACT_uM = THYM_ACT_N / SpyV * (1/MOL) * 1e6;     % thymact uM

            tKa = Ka_Thymosin * THYM_uM * Gactin_uM * dtN;       % uM on  per sec
            tKd = Kd_Thymosin * THYM_ACT_uM * dtN;               % uM off per sec

            TA_Non  = round(tKa / 1e6 * MOL * SpyV);
            TA_Noff = round(tKd / 1e6 * MOL * SpyV);


            % THYM_uM = THYM_N / SpyV *(1/MOL)*1e6;
            % THYM_ACT_uM = THYM_ACT_N / SpyV *(1/MOL)*1e6;
            % TA_on  = THYM_uM * Gactin_uM * TKp * dt / NFact;   % uM on  per sec
            % TA_off = THYM_ACT_uM * TKm * dt / NFact;           % uM off per sec


            % COFILIN VALUES (Cofilin+ACT Association)
            KaCof = Ka_Cofilin * Factin_uM * Cofilin_uM * dtN;

        end
	
    


	%----------------------------------------------------------------------
	end %end main inner loop
	%----------------------------------------------------------------------
	% ADJUST RATE VALUES
	%----------------------------------------------------------------------

    
    

    % ACTIN VALUES
    %-------------------------------------
	FActinN = sum(ACT(:,1));                    % Factin count

    TA_plus_Factin_N = FActinN + THYM_ACT_N;    % thymact + Factin

    GActinN = TOTAL_ACTIN - TA_plus_Factin_N;   % Gactin count

    Factin_uM = FActinN / SpyV * (1/MOL) * 1e6; % Factin uM

	Gactin_uM = GActinN / SpyV * (1/MOL) * 1e6;	% Gactin uM

	fKa = TKa * Gactin_uM * dt;                 % filament on-rate
    




    % ARP VALUES
    %-------------------------------------
    FArpN = numel(ACT(:,1));                % Farp count
    GArpN = GArpN0 - FArpN;                 % Garp count
    Arp_uM = GArpN / SpyV * (1/MOL) * 1e6;  % Garp uM
    nmmf = ACT(:,20);                       % Length of each filament (nm)


    % Arp Branching Rate depends on uM Arp and lengths of existing filaments
    Fact_Branch_uM = ACT(:,1) ./ SpyV * (1/MOL) * 1e6;   % uM of each filament
    ArpBR = Arp_Sc * (Arp_uM/1000 .* Fact_Branch_uM);
    


    
       
    % COFILIN VALUES
    %-------------------------------------
    KaCof = (dt/NFact * Ka_Cofilin * Factin_uM * Cofilin_uM);  % Cof+Fact Ka

 
    
    

    % THYMOSIN VALUES
    %-------------------------------------
    THYM_uM     = THYM_N     / SpyV * (1/MOL) * 1e6; % Thymosin uM
    THYM_ACT_uM = THYM_ACT_N / SpyV * (1/MOL) * 1e6; % Thymosin+ACT uM
    
    tKa = Ka_Thymosin * THYM_uM * Gactin_uM * dt;    % Thymosin+Gactin on-rate
    tKd = Kd_Thymosin * THYM_ACT_uM * dt;            % Thymosin+Gactin off-rate
    
    



    % THYMOSIN DIFFUSION DURING LTP
    %-------------------------------------
    if PAR.doThymT

        if nT==PAR.timeLTPon

            % Before we can replace thymosin inside spines, we need to
            % know the ratio between bound and unbound thymosin in the ES.
            % We are going to assume the T:TA ratio in the ES is equal to
            % the T:TA ratio inside spines during steady-state conditions.
            % Therefore, right before LTP initiates, we will capture the
            % percent of bound and unbound thymosin inside the spine...

            ES_pctT  = THYM_N     / (THYM_ACT_N + THYM_N);
            ES_pctTA = THYM_ACT_N / (THYM_ACT_N + THYM_N);

        end

        if (nT > PAR.timeLTPon) && (nT < PAR.timeLTPoff)

        % Next we need to compute how much bound and unbound thymosin
        % inside spines will diffuse out to the ES. The amount that
        % leaves the spine each second is a constant: PAR.D_Thymosin.
        % Here we compute the number of Thym, ThymAct, and Total
        % molecules will exit the spine at this moment...

            Dout_N_Thym     = round( THYM_N     * PAR.D_Thymosin );
            Dout_N_ThymAct  = round( THYM_ACT_N * PAR.D_Thymosin );
            Dout_N_Total    = Dout_N_Thym + Dout_N_ThymAct;


        % Next we need to compute how much bound and unbound thymosin
        % diffuses from the ES into a spine. The total amount that
        % enters the spine each second is a constant: PAR.D_Thymosin.
        % Here we compute the number of Thym, ThymAct, and Total
        % molecules will enter the spine at this moment. The total
        % amount of bount + unbound thymosin that exits the spine
        % should exactly match the amount enters the spine...

            Din_N_Thym      = round( ES_pctT  * Dout_N_Total );
            Din_N_ThymAct   = round( ES_pctTA * Dout_N_Total );


        % Finally we need to update the amount of bound and unbound
        % thymosin that is currently in the spine, along with the 
        % total amount of actin in the spine...

            THYM_N      = THYM_N     - Dout_N_Thym    + Din_N_Thym;
            THYM_ACT_N  = THYM_ACT_N - Dout_N_ThymAct + Din_N_ThymAct;
            TOTAL_ACTIN = GActinN + FActinN + THYM_ACT_N;

        end
    end


    
    


    % DELETE ORIGINAL FILAMENTS
    %-------------------------------------
    if PAR.delOrigFils && (nT == PAR.delOrigFilsT)
        ACT(1:nStartFils,:) = [];
    end



	% PROCESS STORE AND REMOVE DEL BRANCHES
	%-------------------------------------
	
	% Filament Lifetime
	FLif = ACT(:,16) + 1;
	ACT(:,16) = FLif;
	
	% Mean Length
	ACT(:,18) = (ACT(:,18).*(FLif-1) + ACT(:,1))./FLif;
	
	%Longest Ever Length
	MaxL = ACT(:,1) > ACT(:,17);
	ACT(MaxL,17) = ACT(MaxL,1); 
	
	
	% Delete Filaments With No ACT
	delFi = (ACT(:,1)<1);
	ACT(delFi,15) = nT;	% Death Timestamp
    
    DEDACT = cat(1,DEDACT,ACT(delFi,:));
	
    ACT(delFi,:) = [];
	
	Nfi = numel(ACT(:,1));	% Current #of Filaments ("NFact" used above loop)


	
	

	% COMPUTE XYZ TIP LOCATION
	%-------------------------------------
	% MATH - branch XYZ tip coordinates
	ACT(:,4)  = ACT(:,1) .* nmpp .* sin(ACT(:,8)) .* cos(ACT(:,2)) + ACT(:,3);
	ACT(:,7)  = ACT(:,1) .* nmpp .* sin(ACT(:,8)) .* sin(ACT(:,2)) + ACT(:,6);
	ACT(:,10) = ACT(:,1) .* nmpp .* cos(ACT(:,8)) + ACT(:,9);

	
	


	% SAVE TipMatrix
	%-------------------------------------
	if PAR.saveTipMx && (nT >PAR.SaveTipsAfter) && (mod(nT,PAR.SaveTipsRate) == 0)

        ActMx = TipMatrix(nT,ACT,dims,AcMx,SPYH);

        BTs{numel(BTs)+1}   = ActMx;
        AFMx{numel(AFMx)+1} = ACT;

	end

	
	
	

	% Counters
	%-------------------------------------
	if PAR.doActCounts
		
		CNTS.Fils(nT)        = Nfi;					% Num of branch fils
		CNTS.FAct(nT)        = FActinN;				% Num FActins
		CNTS.GAct(nT)        = GActinN;				% Num GActins
		CNTS.CofAct(nT)      = COFdepoly+ACTdepoly;	% Num of Depoly events
		CNTS.ArpAct(nT)      = ACTpoly+ARPpoly;		% Num of Poly events
		CNTS.delFi(nT)       = sum(delFi);
        CNTS.muFilFa(nT)     = mean(ACT(:,1));
		CNTS.ArpKa(nT)       = sum(ArpBR);
		CNTS.Act_uM(nT)      = Gactin_uM;
		CNTS.fKa(nT)         = fKa;
		CNTS.FnowFtot(nT)    = Nfi / ACT(end,12);
		CNTS.FArp(nT)        = FArpN;
		CNTS.GArp(nT)        = GArpN;
		CNTS.nmTotFil(nT)    = sum(nmmf);
        CNTS.DAct(nT)        = GActinN+FActinN;      % Number free Actins
		CNTS.Thy(nT)         = THYM_N;               % Number free Thymosin
		CNTS.ThyAct(nT)      = THYM_ACT_N;           % Number Thymosin+Gactin
        CNTS.ThyKa(nT)       = tKa;                  % Ka Thymosin
		CNTS.ThyKd(nT)       = tKd;                  % Kd Thymosin
	end
	
    
    
    
    

	% LIVE PLOTS
	%-------------------------------------
    if PAR.doLivePlot && ~mod(nT,PAR.LivePlotMod)

        AcTips  = [ACT(:,4) ACT(:,7) ACT(:,10)];
        PSDTips = AcTips(AcTips(:,3) > inPSD,:);
        SPYTips = AcTips(AcTips(:,3) < inPSD,:);

        P3x = [ACT(:,3) ACT(:,4)]';
        P3y = [ACT(:,6) ACT(:,7)]';
        P3z = [ACT(:,9) ACT(:,10)]';
        P3x(3,:)=NaN; P3y(3,:)=NaN; P3z(3,:)=NaN;

        set(phS,'XData',SPYTips(:,1),'YData',SPYTips(:,2),'ZData',SPYTips(:,3));
        set(phP,'XData',PSDTips(:,1),'YData',PSDTips(:,2),'ZData',PSDTips(:,3));
        set(phF,'XData',P3x(:),'YData',P3y(:),'ZData',P3z(:));
        view(ax1,azel+rota)

        rota = rota + rot0;
    
        if PAR.doSSReset
        if nT == PAR.SSReStep

            ax2.XLim = [PAR.SSReStep*dt/60 Nsteps*dt/60];
            ax3.XLim = [PAR.SSReStep*dt/60 Nsteps*dt/60];
            ax4.XLim = [PAR.SSReStep*dt/60 Nsteps*dt/60];
            ax5.XLim = [PAR.SSReStep*dt/60 Nsteps*dt/60];
            ax6.XLim = [PAR.SSReStep*dt/60 Nsteps*dt/60];
            ax7.XLim = [PAR.SSReStep*dt/60 Nsteps*dt/60];
            ax8.XLim = [PAR.SSReStep*dt/60 Nsteps*dt/60];
            ax9.XLim = [PAR.SSReStep*dt/60 Nsteps*dt/60];
            % ax10.XLim = [10 Nsteps*dt/60];
        end
        end

        T = round(1/dt*60);
        X = 1:numel(CNTS.FAct(1:T:nT));

        % PLOT ACTIN LEVELS
        %-------------------------------------
        line(X,CNTS.GAct(1:T:nT),'Color','b','Parent',ax3,'LineWidth',3);
        line(X,CNTS.FAct(1:T:nT),'Color','r','Parent',ax4,'LineWidth',3);

        % PLOT ARP LEVELS
        %-------------------------------------
        line(X,CNTS.GArp(1:T:nT),'Color','b','Parent',ax5,'LineWidth',3);
        line(X,CNTS.FArp(1:T:nT),'Color','r','Parent',ax6,'LineWidth',3);

        % PLOT THYMOSIN LEVELS
        %-------------------------------------
        line(X,CNTS.Thy(1:T:nT),'Color','r','Parent',ax7,'LineWidth',3);
        line(X,CNTS.ThyAct(1:T:nT),'Color','b','Parent',ax8,'LineWidth',3);

        % PLOT FILAMENT NUMBER & LENGTH
        %-------------------------------------
        line(X,CNTS.Fils(1:T:nT),'Color','r','Parent',ax9,'LineWidth',3);

        % PLOT FILAMENT LENGTH
        %-------------------------------------
        line(X,CNTS.muFilFa(1:T:nT),'Color','r','Parent',ax2,'LineWidth',3);


        delete(ax10.Children);
        YL  = ax10.YLim(2);
        s10 = sprintf(' Step:       % 3.6g  Min:  % 2.6g  ',[nT nT*dt/60]);
        s0  = sprintf(' Tot ACT:    % 9.6g  ',TOTAL_ACTIN);
        s1  = sprintf(' GActin:     % 9.6g  ',GActinN);
        s2  = sprintf(' FActin:     % 9.6g  ',FActinN);
        s3  = sprintf(' Thymosin:   % 9.6g  ',THYM_N);
        s4  = sprintf(' ThymAct:    % 9.6g  ',THYM_ACT_N);
        s5  = sprintf(' GArpN:      % 9.6g  ',GArpN);
        s6  = sprintf(' FArpN:      % 9.6g  ',FArpN);
        s7  = sprintf(' Ka ACT:     % 9.5f  ',fKa);
        s8  = sprintf(' Ka Thymo:   % 9.3f  ',tKa);
        s9  = sprintf(' Kd Thymo:   % 9.3f  ',tKd);
        text(.1,(YL*.92),s0,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        text(.1,(YL*.84),s1,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        text(.1,(YL*.76),s2,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        text(.1,(YL*.68),s3,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        text(.1,(YL*.60),s4,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        text(.1,(YL*.52),s5,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        text(.1,(YL*.44),s6,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        text(.1,(YL*.36),s7,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        text(.1,(YL*.28),s8,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        text(.1,(YL*.20),s9,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        
        text(.1,(YL*.05),s10,'FontSize',14,'FontName','FixedWidth','Parent',ax10);
        pause(.001)
    end


    if ~mod(nT,PAR.LoopPrintMod)
        spfd1 = sprintf(' nT: %6.0f ',                nT);
        spfd2 = sprintf(' GActinN: %6.0f ',      GActinN);
        spfd3 = sprintf(' FActinN: %6.0f ',      FActinN);
        spfd4 = sprintf(' NFils: %4.0f ',          NFact);
        spfd5 = sprintf(' GArpN: %4.0f ',          GArpN);
        spfd6 = sprintf(' FArpN: %4.0f ',          FArpN);
        spfd7 = sprintf(' ThymN: %6.0f ',         THYM_N);
        spfd8 = sprintf(' ThymActN: %6.0f ',  THYM_ACT_N);
        disp([spfd1 spfd2 spfd3 spfd4 spfd5 spfd6 spfd7 spfd8]);
    end


% if nT > 10000; keyboard; end;
end % end main outer loop
%==========================================================================
%% POST-LOOP CALCULATIONS
%==========================================================================



% Lifetime of remaining filaments
ACT(:,15) = nT + (nT-ACT(:,14));      % Artificial Death Time
ACT(:,16) = ACT(:,15) - ACT(:,14);  % Lifetime (nTdied - nTborn)


% DEDACT INCLUDES ALL ACTINS THAT EVER WERE
DEDACT = cat(1,DEDACT,ACT);






%==========================================================================
%% DELETE UNWANTED STEPS
%==========================================================================

if PAR.doNdel

        Ndel    = PAR.nDelSteps;    % Deletes steps 1:SSReStep
        nT0     = nT;               % Save original value of Nt
        Ns0     = Nsteps;           % Save original value of Nsteps

        nT      = nT - Ndel;        % Reduce Nt by num of deleted steps
        Nsteps	= nT;               % Reduce Nsteps by num of deleted steps

        CNTS.Fils(1:Ndel)       = [];
        CNTS.FAct(1:Ndel)       = [];
        CNTS.GAct(1:Ndel)       = [];
        CNTS.CofAct(1:Ndel)     = [];
        CNTS.ArpAct(1:Ndel)     = [];
        CNTS.delFi(1:Ndel)      = [];
        CNTS.ArpKa(1:Ndel)      = [];
        CNTS.Act_uM(1:Ndel)     = [];
        CNTS.fKa(1:Ndel)        = [];
        CNTS.FnowFtot(1:Ndel)   = [];
        CNTS.FArp(1:Ndel)       = [];
        CNTS.GArp(1:Ndel)       = [];
        CNTS.nmTotFil(1:Ndel)   = [];
        CNTS.DAct(1:Ndel)       = [];
        CNTS.Thy(1:Ndel)        = [];
        CNTS.ThyAct(1:Ndel)     = [];
end
%--------------------------------------------------





%==========================================================================
%% PREP FILAMENT TIP DATA FOR OUTPUT
%==========================================================================



ActinTips = [ACT(:,4) ACT(:,7) ACT(:,10)];


[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);


[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);


AxLP = {nT,ACT,inPSD};
Ax = {0,nT,ACT,dims,AcMx,SPYH};






%==========================================================================
%%                         ACT                 AFMx
%==========================================================================
%                       ACT = N x 20       AFMx{nT} = ACT
%--------------------------------------------------------------------------
% N Xa Xo Xt Ya Yo Yt Za Zo Zt Mom ID Cut Born Died Life maxFa muFa NA nmL
% 1 2  3  4  5  6  7  8  9  10 11  12 13  14   15    16   17    18  19 20
%--------------------------------------------------------------------------

% DISPLAY WHAT IS CONTAINED IN EACH CELL OF AFMx
% The first row in this preview table is the column number
% those values aren't actually part of the Nx20 matrix

Tnams={'N';'Xa';'Xo';'Xt';'Ya';'Yo';'Yt';'Za';'Zo';'Zt';...
       'MomID';'ID';'Fkd';'Born';'Died';'Lif';'MaxL';'MeanL';'Null';'Lgth'};

l=[1:20; ACT(50:60,:)]; 
g=1:8;

T=table(l(g,1),l(g,2),l(g,3),l(g,4),l(g,5),l(g,6),l(g,7),l(g,8),l(g,9),...
    l(g,10),l(g,11),l(g,12),l(g,13),l(g,14),l(g,15),l(g,16),l(g,17),...
    l(g,18),l(g,19),l(g,20),'VariableNames',Tnams(1:20));
disp(T);





%==========================================================================
%% GATHER AND FORMAT FILAMENT DATA
%==========================================================================

PSD_Zdim    = dims(7);
SPYtips     = {};
PSDtips     = {};
ActXYZ      = {};

for nn = 1:numel(AFMx)

    Ain = AFMx{nn};

    ActXYZ{nn}  = [Ain(:,3) Ain(:,4) Ain(:,6) Ain(:,7) Ain(:,9) Ain(:,10)];

    ActinTips   = [ActXYZ{nn}(:,2) ActXYZ{nn}(:,4) ActXYZ{nn}(:,6)];

    [Zrow,Zcol] = find(ActinTips(:,3) < PSD_Zdim);

    SPYtps      = ActinTips(Zrow,:);

    SPYtips{nn} = SPYtps;

    [Zrow,Zcol] = find(ActinTips(:,3) >= PSD_Zdim);

    PSDtps      = ActinTips(Zrow,:);

    PSDtips{nn} = PSDtps;

end







%==========================================================================
%% SAVE QUANTITATIVE DATA TO DISK
%==========================================================================





% ActinData
%-------------------------------------
ActinData.PAR       = PAR;      % Input parameters
ActinData.Nsteps    = Nsteps;   % Simulated Steps
ActinData.dt        = dt;       % Real time per step
ActinData.DEDACT    = DEDACT;   % ACT matrix from dead fils
ActinData.ACT       = ACT;      % ACT(n,20) matrix on last step
ActinData.CNTS      = CNTS;     % count data
ActinData.AFMx      = AFMx;     % ACT(n,20) matrix snapshots



% ActinMx
%-------------------------------------
ActinMx.ActXYZ  = ActXYZ;
ActinMx.SPYtips = SPYtips;
ActinMx.PSDtips = PSDtips;
ActinMx.AFMx    = AFMx;     % ACT(n,20) matrix snapshots
ActinMx.AxLP    = AxLP;     % AxLP{nT} = {nT,ACT,inPSD}
ActinMx.Ax      = Ax;       % Ax{nT} = {0,nT,ACT,dims,AcMx,SPYH}
ActinMx.BTs     = BTs;      % BTs{nT} = TipMatrix(nT,ACT,dims,AcMx,SPYH)










%==========================================================================
%% CREATE TABLES WITH VALUES FOR PLOTTING
%==========================================================================
clc; clearvars -except P LOO PAR ActinData ActinMx CNTS



step_vec    = (1:PAR.Nsteps)';          % 1:Nsteps
dt_minutes  = PAR.dt / 60;              % scales step_vec in minutes
time        = step_vec .* dt_minutes;   % scaled time vector
ts          = 20;                       % plot start time in minutes
timeSSidx   = time > ts;                % logical idx for times after 'ts'
timeSSvec   = time(time > ts) - ts;     % array of time values after 'ts'



METRICS_ALL = table();
METRICS     = table();

METRICS_ALL.time     = time          ; % t  cum(dt)
METRICS_ALL.Fils     = CNTS.Fils     ; % n  fils
METRICS_ALL.FAct     = CNTS.FAct     ; % n  Fa monos
METRICS_ALL.GAct     = CNTS.GAct     ; % n  Ga monos
METRICS_ALL.CofAct   = CNTS.CofAct   ; % n  cof cuts
METRICS_ALL.ArpAct   = CNTS.ArpAct   ; % n  arp-FaFaFa
METRICS_ALL.delFi    = CNTS.delFi    ; % n  dead fils
METRICS_ALL.FArp     = CNTS.FArp     ; % n  arp in fils
METRICS_ALL.GArp     = CNTS.GArp     ; % n  free arp
METRICS_ALL.DAct     = CNTS.DAct     ; % n  free Ga
METRICS_ALL.Thy      = CNTS.Thy      ; % n  free thymo
METRICS_ALL.ThyAct   = CNTS.ThyAct   ; % n  thym-Ga
METRICS_ALL.ArpKa    = CNTS.ArpKa    ; % n  arp branchs
METRICS_ALL.Act_uM   = CNTS.Act_uM   ; % uM Ga+Fa
METRICS_ALL.fKa      = CNTS.fKa      ; % ka fils
METRICS_ALL.FnowFtot = CNTS.FnowFtot ; % ra curr:all fils
METRICS_ALL.nmTotFil = CNTS.nmTotFil ; % nm filnet
METRICS_ALL.ThyKa    = CNTS.ThyKa    ; % ka thym
METRICS_ALL.ThyKd    = CNTS.ThyKd    ; % kd thym
METRICS_ALL.muFilFa  = CNTS.muFilFa  ; % mu Fa per fil


METRICS.time     = timeSSvec                ; % t  cum(dt)
METRICS.Fils     = CNTS.Fils     (timeSSidx); % n  fils
METRICS.FAct     = CNTS.FAct     (timeSSidx); % n  Fa monos
METRICS.GAct     = CNTS.GAct     (timeSSidx); % n  Ga monos
METRICS.CofAct   = CNTS.CofAct   (timeSSidx); % n  cof cuts
METRICS.ArpAct   = CNTS.ArpAct   (timeSSidx); % n  arp-FaFaFa
METRICS.delFi    = CNTS.delFi    (timeSSidx); % n  dead fils
METRICS.FArp     = CNTS.FArp     (timeSSidx); % n  arp in fils
METRICS.GArp     = CNTS.GArp     (timeSSidx); % n  free arp
METRICS.DAct     = CNTS.DAct     (timeSSidx); % n  free Ga
METRICS.Thy      = CNTS.Thy      (timeSSidx); % n  free thymo
METRICS.ThyAct   = CNTS.ThyAct   (timeSSidx); % n  thym-Ga
METRICS.ArpKa    = CNTS.ArpKa    (timeSSidx); % n  arp branchs
METRICS.Act_uM   = CNTS.Act_uM   (timeSSidx); % uM Ga+Fa
METRICS.fKa      = CNTS.fKa      (timeSSidx); % ka fils
METRICS.FnowFtot = CNTS.FnowFtot (timeSSidx); % ra curr:all fils
METRICS.nmTotFil = CNTS.nmTotFil (timeSSidx); % nm filnet
METRICS.ThyKa    = CNTS.ThyKa    (timeSSidx); % ka thym
METRICS.ThyKd    = CNTS.ThyKd    (timeSSidx); % kd thym
METRICS.muFilFa  = CNTS.muFilFa  (timeSSidx); % mu Fa per fil





% COMPILE DATA FROM THE ACTIN FILAMENT MATRIX
%--------------------------------------------------------------------------
% N Xa Xo Xt Ya Yo Yt Za Zo Zt Mom ID Cut Born Died Life maxFa muFa NA nmL
% 1 2  3  4  5  6  7  8  9  10 11  12 13  14   15    16   17    18  19 20
%--------------------------------------------------------------------------

CNTS.FilSurvivorNum = [];
CNTS.FilSurvivorPct = [];

if PAR.saveTipMx

    ACTcel = ActinData.AFMx';

    ACTmu = cell2mat(cellfun(@mean, ACTcel, 'UniformOutput',false ));


    ACTmx = cell2mat(ACTcel);
    ACTmx(:,21) = NaN;


    ACTsz= cell2mat(cellfun(@size, ACTcel, 'UniformOutput',false ));
    ACTsz(:,2) = [];
    aj  = cumsum(ACTsz);
    ai  = [0; aj(1:end-1)] + 1;

    ACTmx(1:ACTsz(1),21) = 0;
    for i = 1:numel(ACTsz)

        ACTmx(ai(i):aj(i),21) = i;

    end

    n = numel(ACTcel);

    MXo = ACTcel{1};
    IDo = MXo(:,12);

    StillAlive = zeros(n,1);
    StillAlive(1) = numel(IDo);

    for i = 2:n

        StillAlive(i) = sum(sum(IDo == ACTcel{i}(:,12)' ,2)>0);

    end

    FilSurvivorNum = StillAlive;
    FilSurvivorPct = StillAlive ./ StillAlive(1) .* 100;

    CNTS.FilSurvivorNum = FilSurvivorNum;
    CNTS.FilSurvivorPct = FilSurvivorPct;

end


%==========================================================================
%% SAVE MAT FILE TO DISK
%==========================================================================
clc; clearvars -except P LOO PAR ActinData ActinMx CNTS METRICS_ALL METRICS



save([P.mat P.f PAR.SaveFilePrefix getdt '.mat'],...
    'PAR','ActinData','CNTS','METRICS_ALL','METRICS','ActinMx');






%==========================================================================
%% END MEGALOOP
%==========================================================================
%end; return







%==========================================================================
%% PLOT FINAL VALUES
%==========================================================================
clc; clearvars -except P LOO PAR ActinData ActinMx CNTS METRICS_ALL METRICS
close all;





%--------------------------------------------------------------------------
set(groot,'defaultGraphplotInterpreter','none')
set(groot,'defaultTextInterpreter','none')
set(groot,'defaultAxesTickLabelInterpreter','none')
fh=figure('Units','pixels','Position',[20 40 1400 750],'Color','w');


tiledlayout(4,5);

nexttile; plot(METRICS.time, METRICS.Fils     ); ylabel('n filaments');
nexttile; plot(METRICS.time, METRICS.FAct     ); ylabel('n f-actin monomers')
nexttile; plot(METRICS.time, METRICS.GAct     ); ylabel('n g-actin monomers')
nexttile; plot(METRICS.time, METRICS.CofAct   ); ylabel('n cofilin events')
nexttile; plot(METRICS.time, METRICS.ArpAct   ); ylabel('n seed arp-actins')
nexttile; plot(METRICS.time, METRICS.delFi    ); ylabel('n deleted fils')
nexttile; plot(METRICS.time, METRICS.FArp     ); ylabel('n fil arp')
nexttile; plot(METRICS.time, METRICS.GArp     ); ylabel('n free arp')
nexttile; plot(METRICS.time, METRICS.DAct     ); ylabel('n free gactin')
nexttile; plot(METRICS.time, METRICS.Thy      ); ylabel('n free thym')
nexttile; plot(METRICS.time, METRICS.ThyAct   ); ylabel('n thym-actin dimers')
nexttile; plot(METRICS.time, METRICS.ArpKa    ); ylabel('arp branch rate')
nexttile; plot(METRICS.time, METRICS.Act_uM   ); ylabel('actin uM')
nexttile; plot(METRICS.time, METRICS.fKa      ); ylabel('fil on-rate')
nexttile; plot(METRICS.time, METRICS.FnowFtot ); ylabel('current:all fils')
nexttile; plot(METRICS.time, METRICS.nmTotFil ); ylabel('tot filnet length')
nexttile; plot(METRICS.time, METRICS.ThyKa    ); ylabel('thym Ka')
nexttile; plot(METRICS.time, METRICS.ThyKd    ); ylabel('thym Kd')
nexttile; plot(METRICS.time, METRICS.muFilFa  ); ylabel('mean fil monomers')
nexttile; plot(CNTS.FilSurvivorPct            ); ylabel('survivor pct.')





%==========================================================================
%% DETERMINE HOW MANY TIPS ARE IN THE PSD
%==========================================================================
clc; clearvars -except P LOO PAR ActinData ActinMx CNTS METRICS_ALL METRICS
close all;




PSDtips = ActinMx.PSDtips;
Ntips = cell2mat(cellfun(@size, PSDtips, 'UniformOutput',false )');
Ntips = Ntips(:,1);
Ntime = (1:numel(Ntips))';


% SST PRISM PREP
%-------------------------------------
NN = 1000;
A1PRISM             = METRICS(1000:NN:end,:);
A1PRISM.FAct1k      = A1PRISM.FAct ./ 1000;
A1PRISM.hours       = A1PRISM.time ./ 60;
A1PRISM.Ntips       = Ntips;
A1PRISM = movevars(A1PRISM,{'Fils','muFilFa','Ntips'},'After','time');



%==========================================================================
% FIGURE_DESCRIPTION
%-------------------------------------
close all;
fh1 = figure('Units','pixels','Position',[30 40 1400 400],'Color','w');
ax1=axes('Units','normalized','Position',[.04 .16 .28 .72],'Color','none');
ax2=axes('Units','normalized','Position',[.37 .16 .28 .72],'Color','none');
ax3=axes('Units','normalized','Position',[.70 .16 .28 .72],'Color','none');

axes(ax1); xlabel('time'); title('total fil count');
ARGS.XP = -.09; AXFORMAT(ax1,ARGS); hold on;

axes(ax2); xlabel('time'); title('mean fil length');
ARGS.XP = -.09; AXFORMAT(ax2,ARGS); hold on;

axes(ax3); xlabel('time'); title('receptors');
ARGS.XP = -.09; AXFORMAT(ax3,ARGS); hold on;

fh1.Renderer = 'painters'; 
%-------------------------------------




plot(ax1,A1PRISM.time, A1PRISM.Fils     , 'LineWidth', 2 );
plot(ax2,A1PRISM.time, A1PRISM.muFilFa  , 'LineWidth', 2 );
plot(ax3,A1PRISM.time, A1PRISM.Ntips    , 'LineWidth', 2 );













%==========================================================================
%% GET LOOP DATA
%==========================================================================
clc; clearvars -except P




% GET PATH TO ACTINSST... MAT FILES
%---------------------------------------------
FILES.w = what([P.mat '/4HR_TOTAL_2HR_LTP']);
%FILES.w = what([P.mat '/4HR_TOTAL_2HR_LTP_ARPLIMIT']);
%FILES.w = what(P.mat);
FILES.finfo = dir(FILES.w.path);
FILES.finames = {FILES.finfo.name};
c=~cellfun(@isempty,regexp(FILES.finames,'(ACTIN_REC(\S)+(\.mat+))'));
FILES.finames = string(FILES.finames(c)');
FILES.folder = FILES.finfo.folder;
FILES.fipaths = fullfile(FILES.folder,FILES.finames);
FILES.nfiles = size(FILES.fipaths,1);
disp(FILES.fipaths);
FPATHS = FILES.fipaths;
P.FPATHS = FPATHS;



for i = 1:10
    DATA{i} = load(P.FPATHS(i));
end

clc; clearvars -except P DATA



%% PRISM PREP
%-------------------------------------
clc; clearvars -except P DATA

ActinMx = DATA{1}.ActinMx;
METRICS = DATA{1}.METRICS;
PSDtips = ActinMx.PSDtips;
Ntips = cell2mat(cellfun(@size, PSDtips, 'UniformOutput',false )');
Ntips = Ntips(:,1);
Ntime = (1:numel(Ntips))';
NN = 1000;
TAB             = METRICS(1000:NN:end,:);
TAB.FAct1k      = TAB.FAct ./ 1000;
TAB.hours       = TAB.time ./ 60;
TAB.Ntips       = Ntips;
TAB = movevars(TAB,{'Fils','muFilFa','Ntips'},'After','time');
VNAMES = string(TAB.Properties.VariableNames);


% PRISM PREP
%-------------------------------------
%TAB = table();
for i = 2:10
ActinMx = DATA{i}.ActinMx;
METRICS = DATA{i}.METRICS;

PSDtips = ActinMx.PSDtips;
Ntips = cell2mat(cellfun(@size, PSDtips, 'UniformOutput',false )');
Ntips = Ntips(:,1);
Ntime = (1:numel(Ntips))';

NN = 1000;
TBL             = METRICS(1000:NN:end,:);
TBL.FAct1k      = TBL.FAct ./ 1000;
TBL.hours       = TBL.time ./ 60;
TBL.Ntips       = Ntips;
TBL = movevars(TBL,{'Fils','muFilFa','Ntips'},'After','time');

for j = 1:numel(VNAMES)
    TAB.(VNAMES(j))(:,i) = TBL.(VNAMES(j));
end
end




%==========================================================================
%% DETERMINE HOW MANY TIPS ARE IN THE PSD
%==========================================================================
clc; clearvars -except P DATA TAB


TBL = TAB;

TIME = (1:240)';

tLTP    = 120;
Randy   = 5;
MM      = 3;


T1 = TBL.Ntips(1:tLTP,:);

for i = 1:10
    T1(:,i) = T1(randperm(tLTP),i);
end
TBL.Ntips(1:tLTP,:) = T1;



Recep = TBL.Ntips;
Recep = circshift(Recep,[3,0]);
Recep(1:3,:) = Recep(4:6,:);
Recep = movmean(Recep,[MM,0]);
r1 = round((rand(size(Recep)) - .5) .* Randy) + 2;
Recep = round(Recep ./ 2 + r1);
TBL.Recep = Recep;


%==========================================================================
% CREATE FIGURE WINDOW
%-------------------------------------
close all;
f=figure('Units','pixels','Position',[20 40 1400 750],'Color','w');
ax1=axes('Units','pixels','Position',[100 440 550 290],'Color','none');
ax2=axes('Units','pixels','Position',[780 440 550 290],'Color','none');
ax3=axes('Units','pixels','Position',[100 75 550 290],'Color','none');
ax4=axes('Units','pixels','Position',[780 75 550 290],'Color','none');
%---AXFORMAT()
axes(ax1); xlabel('time'); title('total fil count');
ARGS.XP = -.09; AXFORMAT(ax1,ARGS); hold on;

axes(ax2); xlabel('time'); title('mean fil length');
ARGS.XP = -.09; AXFORMAT(ax2,ARGS); hold on;

axes(ax3); xlabel('time'); title('psd filaments');
ARGS.XP = -.09; AXFORMAT(ax3,ARGS); hold on;

axes(ax4); xlabel('time'); title('receptors');
ARGS.XP = -.09; AXFORMAT(ax4,ARGS); hold on;

fh1.Renderer = 'painters'; 
%-------------------------------------




plot(ax1, mean(TBL.time ,2),  mean(TBL.Fils ,2)    , 'LineWidth', 2 );
plot(ax2, mean(TBL.time ,2),  mean(TBL.muFilFa ,2) , 'LineWidth', 2 );
plot(ax3, mean(TBL.time ,2),  mean(TBL.Ntips ,2)   , 'LineWidth', 2 );
plot(ax4, mean(TBL.time ,2),  mean(TBL.Recep ,2)   , 'LineWidth', 2 );


plot(ax1, TBL.time, TBL.Fils,...
    'LineStyle','none','Marker','.','MarkerSize',8,...
    'MarkerEdgeColor',[.1 .3 .7],'MarkerFaceColor',[.1 .3 .7]);
%ax1.YLim = [300, 400];

plot(ax2, TBL.time, TBL.muFilFa,...
    'LineStyle','none','Marker','.','MarkerSize',8,...
    'MarkerEdgeColor',[.1 .3 .7],'MarkerFaceColor',[.1 .3 .7]);
%ax2.YLim = [40, 100];

plot(ax3, TBL.time, TBL.Ntips,...
    'LineStyle','none','Marker','.','MarkerSize',8,...
    'MarkerEdgeColor',[.1 .3 .7],'MarkerFaceColor',[.1 .3 .7]);
%ax3.YLim = [0, 150];


plot(ax4, TBL.time, TBL.Recep,...
    'LineStyle','none','Marker','.','MarkerSize',8,...
    'MarkerEdgeColor',[.1 .3 .7],'MarkerFaceColor',[.1 .3 .7]);
%ax4.YLim = [0, 100];









%==========================================================================
%% ADVANCED PARAMETERS (DO NOT DELETE)
%==========================================================================
%{

		Abbreviations
-------------------------------------
(+)end: Barbed end of filament
(-)end: Pointed end of filament
mono: 	monomer, free protein monomer not associated with a filament
proto: 	protomer, polymerized protein or protein associated with a filament 
Pi: 	Phosphate
ATG: 	G-actin + ATP (same as ATM)
ADG: 	G-actin + ADP (same as ADM)
ATM: 	G-actin + ATP
ADM: 	G-actin + ADP 
ATF: 	F-actin proto + ATP
APF: 	F-actin proto + ADP-Pi
ADF: 	F-actin proto + ADP
FTB: 	(+)ends terminating in ATP-actin
FPB: 	(+)ends terminating in ADP-Pi-actin
FDB: 	(+)ends terminating in ADP-actin
FTP: 	(-)ends terminating in ATP-actin
FPP: 	(-)ends terminating in ADP-Pi-actin
FDP: 	(-)ends terminating in ADP-actin
CBM: 	(+)end capping mono
CBF: 	(+)end capper proto
CPM, 	(-)end capping mono
CPF: 	(-)end capper proto
FOM: 	Formin mono
FOF: 	Formin proto at (+)end
ARM: 	Arp in mono
ARF: 	Arp proto
FRP: 	Arp proto at (-)end (complex has bound actin at growth end)
FRB: 	Arp proto at (+)end (complex has no bound actin at growth end)



SNUC:	Spontanious Nucleation
FNUC:	Formin nucleation
CBNU:	Nucleation by barbed cap
CPNU:	Nucleation by pointed cap
ASTB:	Poly on rate of ATG at (+)end
ASDB:	Poly on rate of ADG at (+)end
ASTP:	Poly on rate of ATG at (-)end
ASDP:	Poly on rate of ATG at (-)end
DITB:	DPoly off rate of ATF at (+)end
DIPB:	DPoly off rate of APF at (+)end
DIDB:	DPoly off rate of ADF at (+)end
DITP:	DPoly off rate of ATF at (-)end
DIPP:	DPoly off rate of APF at (-)end
DIDP:	DPoly off rate of ADF at (-)end
TTOP:	ATP hydrolysis (ATF -> APF)
PTOD:	P release (APF -> ADF)
DTOT:	Recharge of G-actin by ATF
ASRT:	Arp association rate
DIRP:	Arp dissociation rate
FRGM:	Spontanious fragmentation
ANNL:	Annealing
ASFB:	Association of formin to (+)end
DIFB:	Dissocaition of formin
FASB:	Formin-aided filament growth
ASCB:	Capping of (+)end
DICB:	Uncapping of (-)end
ASCP:	Capping of (+)end
DICP:	Uncapping of (-)end


		Default Rates
-------------------------------------
SNUC:	1e-8	p/(uM^2*s^1)
FNUC:	7e-5	p/(uM^3*s^1)
CBNU:	1e-5	p/(uM^3*s^1)
CPNU:	1e-5	p/(uM^3*s^1)
ASTB: 	11.5	p/µM*s
ASDB: 	3.8		p/µM*s
ASTP: 	1.3		p/µM*s
ASDP: 	0.16	p/µM*s
DITB:	0.90	p/s
DIPB:	0.90	p/s
DIDB:	1.50	p/s
DITP:	0.19	p/s
DIPP:	0.19	p/s
DIDP:	0.26	p/s
FRGM:	1.8e-8	p/s
ANNL:	1e-8	p/µM*s
TTOP:	0.30	p/s
PTOD:	2.6e-3	p/s
DTOT:	1e-2	p/s
ASRT:	1.0e-5	p/µM*s
DIRP:	1.0e-3	p/s
ASCB:	3.0		p/µM*s
DICB:	4e-4	p/s
ASCP:	1.0		p/µM*s
DICP:	1e-2	p/s
ASFB:	3.0		p/µM*s
DIFB:	1e-4	p/s
FASB:	110		p/µM*s



		Rate Calculation Notes
-------------------------------------
exponents of µ^-1 == 1/µ
exponents of µ^-2 == 1/µ/µ
exponents of µ^-2 == 1/µ/µ/µ
[note these are sequential divisions, not 1/(µ/µ)]
Thus if something has a rate of [0.1 µM^-3 s^-1]
to calculate the current rate do [p=.1	uM=uM	s=dt]

	p/(uM^-3)/(s^-1)   or   p*(uM^3)*(s^1)
-------------------------------------



		Default Rates
-------------------------------------
NUCLEATIONS
SNUC:	1e-8	p/(uM^2*s^1)
FNUC:	7e-5	p/(uM^3*s^1)
CBNU:	1e-5	p/(uM^3*s^1)
CPNU:	1e-5	p/(uM^3*s^1)

ACTIN ASSOCIATIONS
ASTB: 	11.5	p/µM*s
ASDB: 	3.8		p/µM*s
ASTP: 	1.3		p/µM*s
ASDP: 	0.16	p/µM*s

ACTIN DISSOCIATIONS
DITB:	0.90	p/s
DIPB:	0.90	p/s
DIDB:	1.50	p/s
DITP:	0.19	p/s
DIPP:	0.19	p/s
DIDP:	0.26	p/s

FRAGMENTATION/ANNEALING
FRGM:	1.8e-8	p/s
ANNL:	1e-8	p/µM*s

ATP/ADP
TTOP:	0.30	p/s
PTOD:	2.6e-3	p/s
DTOT:	1e-2	p/s

ARP2/3
ASRT:	1.0e-5	p/µM*s
DIRP:	1.0e-3	p/s

CAPPING
ASCB:	3.0		p/µM*s
DICB:	4e-4	p/s
ASCP:	1.0		p/µM*s
DICP:	1e-2	p/s

FORMIN
ASFB:	3.0		p/µM*s
DIFB:	1e-4	p/s
FASB:	110		p/µM*s








Abbreviations Full Description: 
F-actin, Filamentous actin; 
nSRF model, non-structurally-resolved filament model; 
SRF model, structurally-resolved filament model; 
Pi, Phosphate; 
ATM: Globular actin (monomeric form) with incorporated ATP; 
ADM: Globular actin (monomeric form) with incorporated ADP; 
ATF: Filamentous actin protomer (F-actin) with incorporated ATP; 
APF: Filamentous actin protomer (F-actin) with incorporated ADP-Pi; 
ADF: Filamentous actin protomer (F-actin) with incorporated ADP; 
FTB: Barbed ends of filaments, terminating by ATP-actin; 
FPB: Barbed ends of filaments, terminating by ADP-Pi-actin; 
FDB: Barbed ends of filaments, terminating by ADP-actin; 
FTP: Pointed ends of filaments, terminating by ATP-actin; 
FPP: Pointed ends of filaments, terminating by ADP-Pi-actin; 
FDP: Pointed ends of filaments, terminating by ADP-actin; 
CBM: Barbed-end capping protein (capper) in monomer (free) form; 
CBF: Barbed-end capper bound to filament; CPM, Pointed-end capper in monomer form; 
CPF: Pointed-end capper bound to filament; 
FOM: Formin in monomer (free) form; 
FOF: Formin, bound to filament barbed ends; 
ARM: Arp2/3 in monomer (free) form; 
ARF: Arp2/3 associated with filament; 
FRP: Arp2/3 associated with filament pointed end (fil terminates at ARP2/3); 
FRB: Arp2/3 associated with filament with no bound actins (fil terminates at ARP2/3); 

%}




%==========================================================================
%% MATH NOTES
%==========================================================================
%{

% Equation to linear transform range
%-------------------------------------
% [a b] to range [c d] and solve for [x]:
% F(x) = c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a))
%-------------------------------------
global sja;
global sjb;
global sjc;
global sjd;

SJKb = StartMonos;

sja=AMX{23};
sjb=SJKb*AMX{24};
sjc=AMX{25};
sjd=AMX{26};

linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));

global ska;
global skb;
global skc;
global skd;
ska=AMX{48};
skb=AMX{49};
skc=AMX{50};
skd=AMX{51};

sigsc = @(jx) (1/(1+ exp(-((skc*(1-(jx-ska)/(skb-ska)) + skd*((jx-ska)/(skb-ska)))))));




miscfun = @(xyz,tta) ([1 0 0; 0 cos(tta) sin(tta); 0 -sin(tta) cos(tta)] * xyz);


%%
close all
clear ArpCurve1
rbg1=.01;rbg2=.99;rbg3=.01;

TheoMaxFact = ceil(SPYheadZN / nmpp);
ScFil = ceil(TheoMaxFact/2);
Nac = TheoMaxFact;
CofR = .002;
CofS = 30;


figure(44)
for CofS = 10:10:100

	nn = 1;
	for Flength = 0:1:Nac
	ArpCurve1(nn) = CofR * exp((Flength-CofS)/ScFil);
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.08;
rbg2=rbg2-.08;
rbg3=rbg3+.08;
drawnow; 
hold on
pause(.2);
end






%%
close all
clear ArpCurve1
rbg1=.01;rbg2=.99;rbg3=.01;

ska=0;
skb=700;
skc=-7;
skd=7;

sigsc = @(jx) (1/(1+ exp(-((skc*(1-(jx-ska)/(skb-ska)) + skd*((jx-ska)/(skb-ska)))))));

ArpR = ArpBR(1);

figure(44)
for ArpScalr = 10:20:200

	nn = 1;
	for Flength = 0:1:skb
	ArpCurve1(nn) = ArpR * sigsc(Flength+ArpScalr);
	%ArpCurve1(nn) = ArpR * (1/(exp(1/170*(1190-14*(Flength+ArpScalr)))+1));
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.08;
rbg2=rbg2-.08;
rbg3=rbg3+.08;
drawnow; 
hold on
pause(.2);
end


%%

close all
clear ArpCurve1
rbg1=.01;rbg2=.99;rbg3=.01;

ska=0;
skb=700;
skc=-7;
skd=7;

figure(46)
for ArpScalr = 10:20:200

	nn = 1;
	for Flength = 0:1:skb
	ArpCurve1(nn) = ArpR * (1/(1+ exp(-(7*(Flength+ArpScalr)/300-7))));
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.08;
rbg2=rbg2-.08;
rbg3=rbg3+.08;
drawnow; pause(.2);
end







%%
% (ArpR * (1/(exp(1/170*(1190-14*(Flength+ArpScalr)))+1))) > rv(2) && ...



%-------------------------------------------%

clear ArpCurve1
close all;
rbg1=.01;rbg2=.99;rbg3=.01;
sja=0;
sjb=300*1;
sjc=0;
sjd=1;
linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));

figure(44)
for CofSc = -30:10:30

	nn = 1;
	for Flength = 0:1:sjb
	ArpCurve1(nn) = CofR * exp(linsc(Flength-CofSc));
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.08;
rbg2=rbg2-.08;
rbg3=rbg3+.08;
drawnow; pause(.5);
end

%-------------------------------------------%
clear ArpCurve1


rbg1=.1;rbg2=.9;rbg3=.1;
figure(43)
for sj = 0:.5:3
	
	linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));

	nn = 1;
	for FLEG = 1:1:300
	ArpCurve1(nn) = ArpR * exp(linsc(FLEG-ArpScalar));
	nn=nn+1;
	end


ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.1;
rbg2=rbg2-.1;
rbg3=rbg3+.1;

end


%-------------------------------------------%
% sjd
% as sjd increases from 0 to 4, the curve shifts from linear to exponential
% as such, a good value for sjd = 3

% sjc
% as sjc decreases from 0 to -4, the curve also becomes more exponential
% a good value for sjc = -1

% logistic sigmoid function (from -5 to 5)
% f(x) = 1 / (1 + e^-x)


clear ArpCurve1
rbg1=.1;rbg2=.9;rbg3=.1;
sja= 0;
sjb= 300;
sjc= -5;
sjd= 5;

figure(43)
for sj = 0:.5:3


linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));

	nn = 1;
	for FLEG = 1:2:300
	ArpCurve1(nn) = (ArpR*(1/(1+ exp(-linsc(FLEG)))));
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.1;
rbg2=rbg2-.1;
rbg3=rbg3+.1;
end




% 5e8 in cells; 1e4 synapses per neuron; 5e8/1e4 = 5e4

Polymerization Rate
(+)end: .012 N/µM*ms
** (12 N/mM*ms) (12 N/µM*s)
** thus at 1 µM free ATP-actin, .012 subunits will be added to the (+)end per ms
** at .1 mM free ATP-actin, 1.2 subunits will be added to the (+)end per ms

Depolymerization Rate
(+)end: 1.4 N/s
(-)end: 0.8 N/s
** dissociation is independent of free actin concentration

To find the critical concentration (Cc) for growth we set the two rate 
equations equal to each other:

12.0/µM*s = 1.4/s
12/1.4 = 1/µM
(+)Cc = .12 µM
(-)Cc = .6 µM





ACT monomer size:
- 5.5 nm x 5.5 nm x 3.5 nm

Factin filaments
- ?-helix composed of 2 strands of subunits
- 28 subunits (14 in each strand) in 1 full 360? turn 
- 180? turn: 36 nm
- 360? turn: 72 nm
- 28 subunits per 72 nm
--------
* For 1 filament to span a 1000 nm spine would require:
-- (mean spine length: 1.0 µm or 1000 nm)
-- each monomer spans 5.1 nm
-- every 5.1 nm requires 2 monomers (cuz double helix)
-- 13.89 turns (~14 turns)
-- 388.89 actin monomers (~400 monomers)
-- 195 monomers per strand
--------

The ATP-binding cleft of actin faces the (+) end of the filament
(+)end grows
(-)end shrinks

Polymerization Rate
(+)end: ~12.0 subunits/µM*s
(-)end: ~1.3 subunits/µM*s
** thus if there is 1 µM of free ATP-Gactin then 12 subunits will be added 
to the (+)end per second and 1.3 subunits will be added to the (-)end every second

Depolymerization Rate
(+)end: ~1.4 subunits/s
(-)end: ~0.8 subunits/s
** dissociation is independent of free Gactin concentration

To find the critical concentration (Cc) for growth we set the two rate 
equations equal to each other:

12.0/µM*s = 1.4/s
12/1.4 = 1/µM
(+)Cc = .12 µM
(-)Cc = .6 µM

Thus when the free actin concentration >.12 µM filaments will grow at the (+) end 
and when >.6 µM filaments will grow at the (-) end too.


----
There are 2 isoforms of actin
- ?-actin & ?-actin

?-actin is enriched in dendritic spines and builds filament stress fibers

G-actin = monomeric "globular" actin
F-actin = filamentous actin

Each actin molecule contains a Mg2+ ion complexed with ATP or ADP

----
Cytosolic concentration of actin in cells ranges from .1 to .5 mM

ACT makes up 1-5% of all cellular protein
(this suggests there is 10 mM total proteins in cells)


A typical cell contains around 5e8 actin molecules

The average number of synapses per neuron was 1e4
http://www.ncbi.nlm.nih.gov/pubmed/2778101

5e8 / 1e4 = 5e4
A typical spine contains 5e4 actin molecules


Total spines volume averaged 0.09 ?m^3 and ranged from 0.01 to 0.38 ?m^3 (n = 133). 
The mode peak value was 0.06 ?m^3. 
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518053/

0.09 ?m^3 = 9e-17 L
5e4 / 9e-17 L = 5.55e20
5.55 N/L * (1 L / 6e23 N) = 0.000925 M = .9 mM
.5 - .9 mM actin / spine

%----------------------------------------
%		STATIC FINAL FIGURES SETUP
%----------------------------------------
%Fh1 = FigSetup(33);
%}







%==========================================================================
%% TIP DENSITY ASSESSMENT
%==========================================================================
%{
clear in3db; clear xxyyzz; clear in3dn;
ActinTips = [ACT(:,4) ACT(:,7) ACT(:,10)];

res3d = 10;
Xdim = dims(4);
Ydim = dims(4);
Zdim = dims(2);

X3d = linspace(0,Xdim,res3d);
Y3d = linspace(0,Ydim,res3d);
Z3d = linspace(0,Zdim,res3d);


for zz = 1:(res3d-1)
	xxyyzz = 1;
	zL = Z3d(zz); zH = Z3d(zz+1);
	for xx = 1:(res3d-1)
		xL = X3d(xx); xH = X3d(xx+1); 
		for yy = 1:(res3d-1)
			yL = Y3d(yy); yH = Y3d(yy+1); 

		in3db{xxyyzz,zz} = in3Dbox(ActinTips,xL,xH,yL,yH,zL,zH);
		xxyyzz = xxyyzz+1;
		end
	end
end

%==================================================%
%%

XYZsz = size(in3db);
XYsz = XYZsz(1);
Xsz = sqrt(XYsz);
Ysz = Xsz;
Zsz = XYZsz(2);

in3dn = zeros(Xsz,Ysz,Zsz);
for SZz = 1:Zsz
	for XZz = 1:Zsz
		for YZz = 1:Zsz
	in3dn(XZz,YZz,SZz) = sum(in3db{XZz+YZz-1,SZz});
		end
	end
end
sum(sum(in3dn))

[xm,ym,zm] = meshgrid(1:9,1:9,1:9);
% scatter3(xm(:),ym(:),zm(:),5,in3dn(:))

figure
set(gcf,'OuterPosition',[300 200 400 700])
scatter3(ActinTips(:,1),ActinTips(:,2),ActinTips(:,3),50,'filled')
xlim([-250 250]);ylim([-250 250]);zlim([0 1100])
view([-27.5 6]);
hold on

%%
clear cmx szmx xmx ymx zmx
for xyzm = 1:9
cmx = in3dn(:,:,xyzm);
szmx = (in3dn(:,:,xyzm)+1) .^6;
xmx = xm(:,:,xyzm);
ymx = ym(:,:,xyzm);
zmx = zm(:,:,xyzm);
scatter3(xmx(:),ymx(:),zmx(:),szmx(:),cmx(:),'filled')
hold on;
end


%%
figure
set(gcf,'OuterPosition',[300 200 400 700])
scatter3(ActinTips(:,1),ActinTips(:,2),ActinTips(:,3),50,'filled')
xlim([-250 250]);ylim([-250 250]);zlim([0 1100])
view([-27.5 6]);
% axis vis3d

% in3dn = zeros(XYsz,Zsz)
% for szz = 1:Zsz
% 	for sxy = 1:XYsz
% 	in3dn(sxy,szz) = sum(in3db{sxy,szz});
% end
% end
% sum(sum(in3dn))

%}






%==========================================================================
%% TRIANGULATION AND CONVEX HULL
%==========================================================================
%{

%%
%==================================================%
Ori_xyz = [ACT(:,3) ACT(:,6) ACT(:,9)];
Tip_xyz = [ACT(:,4) ACT(:,7) ACT(:,10)];
ATPxyz = [Ori_xyz; Tip_xyz];

% radial distance to spine shaft membrane
XYtipLoc = sqrt(ATPxyz(:,1).^2 + ATPxyz(:,2).^2);

ZtipInHead = (ATPxyz(:,3) <= SPYheadZN) & (ATPxyz(:,3) >= SPYheadZS);
ZtipInNeck = ATPxyz(:,3) < SPYheadZS;

XYneckOut = XYtipLoc > SPYneckRad;	% Logical array of filaments beyond SPYneckRad
XYheadOut = XYtipLoc > SPYheadRad;  % Logical array of filaments beyond SPYheadRad
ZtopOut = ATPxyz(:,3) > SPYheadZN;	% Logical array of filaments above SPYheadZN
ZbotOut = ATPxyz(:,3) < 0;		% Logical array of filaments below zero

TipOut = ((XYneckOut & ZtipInNeck) + (XYheadOut & ZtipInHead) + ZtopOut + ZbotOut)>0;

ATPxyz(TipOut,:) = [];

ATPxyz(end:end+4,:) = 0;

szATP = size(ATPxyz,1);

ATPxyz(szATP+1,:) = [SPYneckRad/2 SPYneckRad/2 0];
ATPxyz(szATP+2,:) = [-SPYneckRad/2 SPYneckRad/2 0];
ATPxyz(szATP+3,:) = [SPYneckRad/2 -SPYneckRad/2 0];
ATPxyz(szATP+4,:) = [-SPYneckRad/2 -SPYneckRad/2 0];

%==================================================%
fig1 = figure;
set(fig1,'OuterPosition',[500 100 800 900])

AOTxyz = ATPxyz;
minAOT = -min(AOTxyz); 
AOTxyz(:,1) = AOTxyz(:,1) + minAOT(1);
AOTxyz(:,2) = AOTxyz(:,2) + minAOT(2);
AOTxyz(:,3) = AOTxyz(:,3);

AOTxyzU = unique(AOTxyz,'rows');

Uxyz = AOTxyzU;

TinH = Uxyz(:,3) >= (SPYheadZS);
HTxyz = Uxyz(TinH,:);

TinN = Uxyz(:,3) < (SPYheadZS);
NTxyz = Uxyz(TinN,:);

TinNH = (Uxyz(:,3) >= (SPYheadZS-80)) & (Uxyz(:,3) <= (SPYheadZS+50));
NHTxyz = Uxyz(TinNH,:);

DTriTip = delaunayTriangulation(HTxyz);
[FBtri,FBpoints] = freeBoundary(DTriTip);

trimesh(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

NTriTip = delaunayTriangulation(NTxyz);
[NFBtri,NFBpoints] = freeBoundary(NTriTip);

trimesh(NFBtri,NFBpoints(:,1),NFBpoints(:,2),NFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

NHTriTip = delaunayTriangulation(NHTxyz);
[NHFBtri,NHFBpoints] = freeBoundary(NHTriTip);

trimesh(NHFBtri,NHFBpoints(:,1),NHFBpoints(:,2),NHFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.3);
hold on


scatter3(AOTxyz(:,1),AOTxyz(:,2),AOTxyz(:,3),20,'fill','r')
axis vis3d
view([-48 8]);
%==================================================%
%%



















%==================================================%
fig1 = figure;
set(fig1,'OuterPosition',[500 100 800 900])

Ori_xyz = [ACT(:,3) ACT(:,6) ACT(:,9)];
Tip_xyz = [ACT(:,4) ACT(:,7) ACT(:,10)];
AOTxyz = [Ori_xyz; Tip_xyz];
minAOT = -min(AOTxyz); 
AOTxyz(:,1) = AOTxyz(:,1) + minAOT(1);
AOTxyz(:,2) = AOTxyz(:,2) + minAOT(2);
AOTxyz(:,3) = AOTxyz(:,3) + minAOT(3);

AOTxyzU = unique(AOTxyz,'rows');
DTriTip = delaunayTriangulation(AOTxyzU);
[FBtri,FBpoints] = freeBoundary(DTriTip);

trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on
scatter3(AOTxyz(:,1),AOTxyz(:,2),AOTxyz(:,3),20,'fill','r')
axis vis3d
view([-36.5 6]);
%==================================================%

%==================================================%
Ori_xyz = [ACT(:,3) ACT(:,6) ACT(:,9)];
Tip_xyz = [ACT(:,4) ACT(:,7) ACT(:,10)];
AOTxyz = [Ori_xyz; Tip_xyz];
minAOT = -min(AOTxyz); 
AOTxyz(:,1) = AOTxyz(:,1) + minAOT(1);
AOTxyz(:,2) = AOTxyz(:,2) + minAOT(2);
AOTxyz(:,3) = AOTxyz(:,3) + minAOT(3);

DTriTip = delaunayTriangulation(AOTxyz);
[DThull DThullV] = convexHull(DTriTip);

trisurf(DThull,DTriTip.Points(:,1),DTriTip.Points(:,2),DTriTip.Points(:,3),...
       'FaceColor','cyan')

trimesh(DThull,DTriTip.Points(:,1),DTriTip.Points(:,2),DTriTip.Points(:,3))
axis vis3d

%==================================================%



%%
%==================================================%
fig1 = figure;
set(fig1,'OuterPosition',[500 100 800 900])

Ori_xyz = [ACT(:,3) ACT(:,6) ACT(:,9)];
Tip_xyz = [ACT(:,4) ACT(:,7) ACT(:,10)];
AOTxyz = [Ori_xyz; Tip_xyz];
minAOT = -min(AOTxyz); 
AOTxyz(:,1) = AOTxyz(:,1) + minAOT(1);
AOTxyz(:,2) = AOTxyz(:,2) + minAOT(2);
AOTxyz(:,3) = AOTxyz(:,3) + minAOT(3);

AOTxyzU = unique(AOTxyz,'rows');

Uxyz = AOTxyzU;


TinH = Uxyz(:,3) >= (SPYheadZS);
HTxyz = Uxyz(TinH,:);

TinN = Uxyz(:,3) < (SPYheadZS-10);
NTxyz = Uxyz(TinN,:);

TinNH = (Uxyz(:,3) >= (SPYheadZS-50)) & (Uxyz(:,3) <= (SPYheadZS+30));
NHTxyz = Uxyz(TinNH,:);

DTriTip = delaunayTriangulation(HTxyz);
[FBtri,FBpoints] = freeBoundary(DTriTip);

trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

NTriTip = delaunayTriangulation(NTxyz);
[NFBtri,NFBpoints] = freeBoundary(NTriTip);

trisurf(NFBtri,NFBpoints(:,1),NFBpoints(:,2),NFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

NHTriTip = delaunayTriangulation(NHTxyz);
[NHFBtri,NHFBpoints] = freeBoundary(NHTriTip);

trisurf(NHFBtri,NHFBpoints(:,1),NHFBpoints(:,2),NHFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

scatter3(AOTxyz(:,1),AOTxyz(:,2),AOTxyz(:,3),20,'fill','r')
axis vis3d
view([-36.5 6]);
%==================================================%
%%
%}







%% EOF