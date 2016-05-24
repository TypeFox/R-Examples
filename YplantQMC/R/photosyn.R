photosynfun <- function(
PAR = 1000,
TLEAF = 20,
CS = 380,
RH = 0,
VPD = 1.5,
JMAX25 = 145.4,
VCMAX25 = 89.5,
RD0 = 0.9,
Q10F = 0.67,
RTEMP = 25,
MODELGS = 4,
EMAXLEAF = 400,
KTOT = 2,
WEIGHTEDSWP = 0,
PSILIN = 0,
PSIL = 0,
GSMIN = 0.01,
G0 = 0.03,
D0L = 5,
GAMMA = 0,
G1 = 4.3,
GK = 0.3,
TUZFUN=1,
SF = 3.2,
PSIV = -1.9,
FSOIL = 0,
IECO = 1,
EAVJ = 37259,
EDVJ = 200000,
DELSJ = 640.02,
THETA = 0.4,
HMSHAPE = 0.999,
AJQ = 0.324,
GMESO = -1.0,
EAVC = 47590,
EDVC = 0.0,
DELSC = 0.0,
TVJUP = -100,
TVJDN = -100,
DAYRESP = 1.0,
TBELOW = -100.0,
PATM=101000,
psilmethod=2,
returncols = c('PAR','TLEAF','VPD','PSILIN','PSIL','GS','ALEAF','ELEAF','RD','CI')
){

VPD <- VPD * 1000

# Outputs need to be initialized.
GS <- -999
ALEAF <- -999
RD <- -999
ALEAFHM <- -999
CI <- -999

if(any(is.na(c(VPD,TLEAF,PAR,RH))))return(rep(NA,length(returncols)))

if(MODELGS==6){
	PSILIN <- getpsil(PAR,TLEAF,CS,RH,VPD/1000,JMAX25,VCMAX25,RD0,Q10F,RTEMP,MODELGS,EMAXLEAF,
	KTOT,WEIGHTEDSWP,PSILIN,PSIL,GSMIN,G0,D0L,GAMMA,G1,GK,TUZFUN,SF,PSIV,HMSHAPE,FSOIL,GS,
	ALEAF,ALEAFHM,RD,IECO,EAVJ,EDVJ,DELSJ,THETA,AJQ,EAVC,EDVC,DELSC,TVJUP,
	TVJDN,DAYRESP,TBELOW,PATM,METHOD=psilmethod)  # Note: 
										 # METHOD=1 uses lame fortran optimizer (Slow!),
	                                     # METHOD=2 uses ZBRENT but does not work yet!
										 # METHOD=3 uses R-based 'optimize'. Works fine.
	
}

f <- .Fortran("photosyn", as.double(PAR),
                     as.double(TLEAF),
					 as.double(CS),
					 as.double(RH),
					 as.double(VPD),
					 as.double(PATM),
					 as.double(JMAX25),
					 as.integer(IECO),
					 as.double(EAVJ),
					 as.double(EDVJ),
					 as.double(DELSJ),
					 as.double(VCMAX25),
					 as.double(EAVC),
					 as.double(EDVC),
					 as.double(DELSC),
					 as.double(TVJUP),
					 as.double(TVJDN),
					 as.double(THETA),
					 as.double(AJQ),
					 as.double(RD0),
					 as.double(Q10F),
					 as.double(RTEMP),
					 as.double(DAYRESP),
					 as.double(TBELOW),
					 as.integer(MODELGS),
					 as.double(EMAXLEAF),
					 as.double(KTOT),
					 as.double(WEIGHTEDSWP),
					 as.double(PSILIN),
					 as.double(PSIL),
					 as.double(GSMIN),
					 as.double(G0),
					 as.double(D0L),
					 as.double(GAMMA),
					 as.double(G1),
					 as.double(GK),
					 as.integer(TUZFUN),
					 as.double(SF),
					 as.double(PSIV),
					 as.double(HMSHAPE),
					 as.double(FSOIL),
					 as.double(GS),
					 as.double(ALEAF),
					 as.double(ALEAFHM),
					 as.double(RD),
					 as.double(CI),
					 PACKAGE='YplantQMC')
		 
vec <- Reduce(c,f)
names(vec) <- c('PAR','TLEAF','CS','RH','VPD','PATM','JMAX25','IECO','EAVJ','EDVJ','DELSJ',
 'VCMAX25','EAVC','EDVC','DELSC','TVJUP','TVJDN','THETA','AJQ','RD0','Q10F','RTEMP','DAYRESP',
 'TBELOW','MODELGS','EMAXLEAF','KTOT','WEIGHTEDSWP','PSILIN','PSIL','GSMIN','G0','D0L','GAMMA',
 'G1','GK','TUZFUN','SF','PSIV','HMSHAPE','FSOIL','GS','ALEAF','ALEAFHM','RD','CI')

vec['ELEAF'] <- 1000*1.6*vec['GS']*vec['VPD']/vec['PATM']
 
vec['VPD'] <- vec['VPD'] / 1000

return(vec[returncols])
}

photosyn <- function(returncols=c('PAR','TLEAF','VPD','PSILIN','PSIL','GS','ALEAF','ALEAFHM','ELEAF','RD','CI'),...){
	oldform <- formals(photosynfun)
	formals(photosynfun)$returncols <- returncols
	dfr <- as.data.frame(t(mapply(photosynfun, ...)))
	formals(photosynfun) <- oldform
	return(dfr)
}

