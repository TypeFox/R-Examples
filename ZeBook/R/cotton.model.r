################################################################################
# Working with dynamic models for agriculture - R script for pratical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-06-17
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################ FUNCTIONS #####################################
#' @title The Cotton model (dynamic for numbers of Cotton fruiting points).
#' @description \strong{Model description.} TO COMPLETE
#' @param TESQ : TO COMPLETE
#' @param PMAX : TO COMPLETE
#' @param AFL : TO COMPLETE
#' @param AL : TO COMPLETE
#' @param AOP : TO COMPLETE
#' @param P1 : TO COMPLETE
#' @param P2 : TO COMPLETE
#' @param P3 : TO COMPLETE
#' @param P4 : TO COMPLETE
#' @param P5 : TO COMPLETE
#' @param PF : TO COMPLETE
#' @param PSF : TO COMPLETE
#' @param TSQ : TO COMPLETE
#' @param P : TO COMPLETE
#' @param PR : TO COMPLETE
#' @param PT : TO COMPLETE
#' @param tend : TO COMPLETE
#' @return data.frame with daily state variable
#' @export
cotton.model<-function(TESQ,PMAX,AFL,AL,AOP,P1,P2,P3,P4,P5,PF,
	PSF,TSQ,P,PR,PT,tend) {
# calculate secondary parameters, which involve fixed constants
AF<-AFL+5
AM<-8
ANO<-AOP-5
TF<-TSQ+AF
AW<-5
# time is in pysiological days. Goes from 1 to tend
# initialize vectors
DN<-rep(0,tend+1)
DNM<-rep(0,tend+1)
DNSET<-rep(0,tend+1)

N<-rep(0,tend+1)
NM<-rep(0,tend+1)
NSET<-rep(0,tend+1)
NFLOWER<-rep(0,tend+1)
NSQUARE<-rep(0,tend+1)
NLARGE<-rep(0,tend+1)
NSMALL<-rep(0,tend+1)
NOPEN<-rep(0,tend+1)

# loop over physiological days
for (t in 1:tend)
{
# calculate daily changes
# Number of squares formed. Starts at t=TSQ
DN[t]<-0
if (t>TSQ & t<=TESQ & NSET[t]<PMAX )
{
if ( t <= TF) DN[t]<-P1*(t-TSQ)
if (t > TF ) DN[t]<-P2*(1-NSET[t]/PF)
}
 
# Number of squares marked for shedding. First one is marked at t=TSQ+AM
DNM[t]<-0
if (t > TSQ+AM)
{
if (t <=TF) DNM[t]<-DN[t-AM]*P3
if (t>TF & t<=(TF+AM)) DNM[t]<-DN[t-AM]*P4
if (t > (TF+AM) &  NSET[t] < PMAX) DNM[t]<- DN[t-AM]*(P4+(1-P4)*(1-NSET[t]/PF))
if (NSET[t] >= PMAX) DNM[t]<-DN[t-AM]
}

# Number of bolls not marked for shedding. First one marked at time t=TSQ+AF 
DNSET[t]<-0
if (t > TSQ+AF & NSET[t] <= PMAX)
{
if (t > AM  ) DNSET[t]<-DN[t-AF]-DNM[t-AF+AM]*P5*(1-NSET[t]/PSF)
}

# Update state variables 
N[t+1] <- N[t]+DN[t]
NM[t+1] <- NM[t]+DNM[t]
NSET[t+1] <- NSET[t]+DNSET[t]

# Calculate measureable variables
if (t>TSQ+AFL) NFLOWER[t]<-N[t-AFL]-NM[t-AFL+AM]
if (t>=TSQ) NSQUARE[t]<-N[t]-NFLOWER[t]-NM[t-AW]
if (t>TSQ+AL) NLARGE[t]<-NSET[t-AL]
if (t> TSQ+AFL) NSMALL[t]<-NFLOWER[t]-NLARGE[t]-(NFLOWER[t-(AF-AFL)-AW]-NSET[t-AW])
if (t>TSQ+AF+AOP) NOPEN[t]<-NSET[t-AOP]
}  # end of loop over physiological days

return(list(t=1:(tend+1),NSQUARE=NSQUARE,NSMALL=NSMALL,NLARGE=NLARGE,
	NOPEN=NOPEN))
}

#end of file