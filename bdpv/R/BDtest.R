BDtest <-
function(xmat, pr, conf.level=0.95)
{

if(!is.matrix(xmat)){stop("Argument xmat must be a matrix!")}
if(!ncol(xmat)==2 | !nrow(xmat)==2){stop("Test results (Argument xmat) must be a 2x2 matrix!")}
if(any((xmat-as.integer(xmat))!=0)){stop("Test results (Argument xmat) should contain integer values only!")}
if(pr<=0 | pr >=1){stop("Prevalence must be a number greater than 0 and smaller than 1!")}

# Refine input data:

INDAT<-as.data.frame(xmat)
colnames(INDAT)<-c("True positive", "True negative")
attr(INDAT, which="caption") <- "Input data set with columns representing the true property of the compounds and rows representing the result of a binary diagnostic test."

rownames(INDAT)<-c("Test positive", "Test negative")

# CI sensitivity and specificity:

x11<-xmat[1,1]
x01<-xmat[2,1]
x10<-xmat[1,2]
x00<-xmat[2,2]

SECIl <- binom.test(x=x11, n=x11+x01, alternative="greater", conf.level=conf.level)
SECIts <- binom.test(x=x11, n=x11+x01, alternative="two.sided", conf.level=conf.level)

SPCIl <- binom.test(x=x00, n=x00+x10, alternative="greater", conf.level=conf.level)
SPCIts <- binom.test(x=x00, n=x00+x10, alternative="two.sided", conf.level=conf.level)

lperc<-as.character(signif(conf.level*100,4))
tsperc<-as.character(signif((1-(1-conf.level)/2)*100,4))

SESPDAT <- data.frame(
c(SECIl$estimate, SPCIl$estimate),
c(SECIl$conf.int[1], SPCIl$conf.int[1]),
c(SECIts$conf.int[1], SPCIts$conf.int[1]), 
c(SECIts$conf.int[2], SPCIts$conf.int[2])
)
colnames(SESPDAT)<-c("Estimate",
paste("Lower ", lperc, "% limit", sep=""),
paste("Lower ", tsperc, "% limit", sep=""),
paste("Upper ", tsperc, "% limit", sep="")
)
rownames(SESPDAT)<-c("Sensitivity", "Specificity")

attr(SESPDAT, which="caption") <- "Estimates and exact confidence limits for assay sensitivity and specificity. "

# CI for predictive values:


x1 <- c(x11, x01)
x0 <- c(x10, x00)

if(any(xmat==0)){
CIPPVl<-CIlppvak(x0=x0, x1=x1, p=pr, conf.level=conf.level,
 alternative="greater")
CINPVl<-CIlnpvak(x0=x0, x1=x1, p=pr, conf.level=conf.level,
 alternative="greater")
CIPPVts<-CIlppvak(x0=x0, x1=x1, p=pr, conf.level=conf.level,
 alternative="two.sided")
CINPVts<-CIlnpvak(x0=x0, x1=x1, p=pr, conf.level=conf.level,
 alternative="two.sided")
}
else{
CIPPVl<-CIlppv(x0=x0, x1=x1, p=pr, conf.level=conf.level,
 alternative="greater")
CINPVl<-CIlnpv(x0=x0, x1=x1, p=pr, conf.level=conf.level,
 alternative="greater")
CIPPVts<-CIlppv(x0=x0, x1=x1, p=pr, conf.level=conf.level,
 alternative="two.sided")
CINPVts<-CIlnpv(x0=x0, x1=x1, p=pr, conf.level=conf.level,
 alternative="two.sided")
}

PPVNPVDAT<-data.frame(
c(CINPVl$estimate, CIPPVl$estimate),
c(CINPVl$conf.int[1], CIPPVl$conf.int[1]),
c(CINPVts$conf.int[1], CIPPVts$conf.int[1]),
c(CINPVts$conf.int[2], CIPPVts$conf.int[2]))

colnames(PPVNPVDAT)<-c( "Estimate",
paste("Lower ", lperc, "% limit", sep=""),
paste("Lower ", tsperc, "% limit", sep=""),
paste("Upper ", tsperc, "% limit", sep="")
)

rownames(PPVNPVDAT)<-c("NPV","PPV")

attr(PPVNPVDAT, which="caption") <- paste("Estimates and asymptotic confidence limits for predictive values. The prevalence is assumed to be ", 
signif(pr, 4), ".", sep="")

npvobs<-PPVNPVDAT[1,1]
ppvobs<-PPVNPVDAT[2,1]

npvl<-PPVNPVDAT[1,2]
ppvl<-PPVNPVDAT[2,2]

if(ppvobs<=pr){warning("The data indicate that, for the given prevalence, the binary diagnostic test is not usefull (observed PPV < specified prevalence).
Please check, whether the data have been entered correctly!")}
#else{if(ppvl<=pr){warning("NOTE: For the given data and prevalence, the binary diagnostic test might be not useful (lower limit of PPV < specified prevalence).
#Please check, whether the data have been entered correctly!")}}

if(npvobs<=(1-pr)){warning("The data indicate that, for the given prevalence, the binary diagnostic test is not useful (observed NPV < 1-specified prevalence).
Please check, whether the data have been entered correctly!")}
#else{if(npvl<=(1-pr)){warning("NOTE: For the given data and prevalence, the binary diagnostic test might be not useful (lower limit of NPV < 1-specified prevalence).
#Please check, whether the data have been entered correctly!")}}

RESULT<-list(INDAT=INDAT, SESPDAT=SESPDAT, PPVNPVDAT=PPVNPVDAT)
class(RESULT)<-"BDtest"
return(RESULT)

}

