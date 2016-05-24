getpsil <- function(
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
PSILIN = -999,
PSIL = -111,
GSMIN = 0.01,
G0 = 0.03,
D0L = 5,
GAMMA = 0,
G1 = 4.3,
GK = 0.3,
TUZFUN=1,
SF = 3.2,
PSIV = -1.9,
HMSHAPE = 0.999,
FSOIL = 0,
GS = 0,
ALEAF = 0,
ALEAFHM = 0,
RD = 0,
IECO = 1,
EAVJ = 37259,
EDVJ = 200000,
DELSJ = 640.02,
THETA = 0.4,
AJQ = 0.324,
EAVC = 47590,
EDVC = 0.0,
DELSC = 0.0,
TVJUP = -100,
TVJDN = -100,
DAYRESP = 1.0,
TBELOW = -100.0,
PATM=101000,
METHOD=1
){

VPD <- VPD * 1000
if(METHOD != 1){
	METHOD <- 1
	warning("METHOD > 1 obsolete. Using method 1.")

}
if(any(is.na(c(VPD,TLEAF,PAR,RH))))return(NA)

# if(METHOD == 1){
f <- .Fortran("psiltuzet", as.double(PAR),
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
					 as.double(PSIL), # nr 29. Can add pars below, not above!!
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
					 PACKAGE='YplantQMC')
	return(f[[29]])
# }

# if(METHOD == 2){

# f <- .Fortran("psilfind", as.double(PAR),
                     # as.double(TLEAF),
					 # as.double(CS),
					 # as.double(RH),
					 # as.double(VPD),
					 # as.double(PATM),
					 # as.double(JMAX25),
					 # as.integer(IECO),
					 # as.double(EAVJ),
					 # as.double(EDVJ),
					 # as.double(DELSJ),
					 # as.double(VCMAX25),
					 # as.double(EAVC),
					 # as.double(EDVC),
					 # as.double(DELSC),
					 # as.double(TVJUP),
					 # as.double(TVJDN),
					 # as.double(THETA),
					 # as.double(AJQ),
					 # as.double(RD0),
					 # as.double(Q10F),
					 # as.double(RTEMP),
					 # as.double(DAYRESP),
					 # as.double(TBELOW),
					 # as.integer(MODELGS),
					 # as.double(EMAXLEAF),
					 # as.double(KTOT),
					 # as.double(WEIGHTEDSWP),
					 # as.double(PSIL), # nr 29. Can add pars below, not above!!
					 # as.double(GSMIN),
					 # as.double(G0),
					 # as.double(D0L),
					 # as.double(GAMMA),
					 # as.double(G1),
					 # as.double(GK),
					 # as.integer(TUZFUN),
					 # as.double(SF),
					 # as.double(PSIV),
					 # as.double(HMSHAPE),
					 # as.double(FSOIL),
					 # as.double(GS),
					 # as.double(ALEAF),
					 # as.double(ALEAFHM),
					 # as.double(RD),
					 # PACKAGE='GasExchangeR')
	# return(f[[29]])

# }

# if(METHOD == 3){

# objfun <- function(P){
# CI <- -999
# f <- .Fortran("photosyn", as.double(PAR),
                     # as.double(TLEAF),
					 # as.double(CS),
					 # as.double(RH),
					 # as.double(VPD),
					 # as.double(PATM),
					 # as.double(JMAX25),
					 # as.integer(IECO),
					 # as.double(EAVJ),
					 # as.double(EDVJ),
					 # as.double(DELSJ),
					 # as.double(VCMAX25),
					 # as.double(EAVC),
					 # as.double(EDVC),
					 # as.double(DELSC),
					 # as.double(TVJUP),
					 # as.double(TVJDN),
					 # as.double(THETA),
					 # as.double(AJQ),
					 # as.double(RD0),
					 # as.double(Q10F),
					 # as.double(RTEMP),
					 # as.double(DAYRESP),
					 # as.double(TBELOW),
					 # as.integer(MODELGS),
					 # as.double(EMAXLEAF),
					 # as.double(KTOT),
					 # as.double(WEIGHTEDSWP),
					 # as.double(P),
					 # as.double(PSIL),    # nr 30. Can add pars below, not above!!
					 # as.double(GSMIN),
					 # as.double(G0),
					 # as.double(D0L),
					 # as.double(GAMMA),
					 # as.double(G1),
					 # as.double(GK),
					 # as.integer(TUZFUN),
					 # as.double(SF),
					 # as.double(PSIV),
					 # as.double(HMSHAPE),
					 # as.double(FSOIL),
					 # as.double(GS),
					 # as.double(ALEAF),
					 # as.double(ALEAFHM),
					 # as.double(RD),
					 # as.double(CI),
					 # PACKAGE='GasExchangeR')
	
	# o <- (P - f[[30]])^2
	# return(o)					 
# }

	# opt <- optimize(objfun, interval=c(-50,0))
    # return(opt$minimum)

# }


}