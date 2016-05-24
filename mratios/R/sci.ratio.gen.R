"sci.ratio.gen" <-
function(Y, X, Num.Contrast, Den.Contrast, alternative = "two.sided", conf.level = 0.95, method="Plug") {

if(all(c("Plug", "Bonf", "MtI", "Unadj")!=method))
 {stop("argument method mis-specified")}

if(all(c("two.sided", "less", "greater")!=alternative))
 {stop("argument alternative mis-specified")}

if(nrow(Num.Contrast)!=nrow(Den.Contrast))
 {stop("Num.Contrast and Den.Contrast must have same number of rows")}

#NCrowsum<-apply(Num.Contrast, MARGIN=1, FUN=sum)
#DCrowsum<-apply(Den.Contrast, MARGIN=1, FUN=sum)
#if( !all( NCrowsum == DCrowsum ) )
# { cat("Warning: Check whether denominator and numerator contrast matrices are appropriately defined!", "\n") }

NC0<-apply(X=Num.Contrast, MARGIN=1, function(x){all(x==0)})
DC0<-apply(X=Den.Contrast, MARGIN=1, function(x){all(x==0)})
if(any(c(NC0,DC0)))
 {cat("Warning: At least one row of the numerator or denominator contrast matrices is a vector with all components equal to zero","\n")}

CMat <- Num.Contrast
DMat <- Den.Contrast

n.comp <- nrow(Num.Contrast)

Lm.Fit <- lm(Y~X-1)

Names.Coeff<-names(Lm.Fit$coef)

colnames(CMat)<-colnames(DMat)<-Names.Coeff

if(is.null(rownames(CMat)) && is.null(rownames(DMat)) )
 { compnames <- paste( "C", 1:nrow(Num.Contrast), sep="") }

    else{ if(any(rownames(CMat)!=rownames(DMat)) )
     {compnames <- paste(rownames(CMat), rownames(DMat), sep="/")}

    else{compnames <- rownames(CMat)}}

Beta.Coeff <- Lm.Fit$coef

if(any( c( ncol(Num.Contrast), ncol(Den.Contrast) )!=length(Beta.Coeff) ))
 {stop("Num.Contrast and Den.Contrast must have same number of columns as parameters are fitted in the linear model")}

Degree.f <- Lm.Fit$df 

Pooled.Var <- sum(Lm.Fit$residual^2)/Degree.f

gammaC.vec <- (CMat%*%Beta.Coeff)/DMat%*%Beta.Coeff          #  MLE of the ratios

M <- solve(t(X)%*%X)

CorrMat.plug <- matrix(as.numeric(rep(NA,n.comp*n.comp)), nrow=n.comp)
    for(i in 1:n.comp) {
        for(j in 1:n.comp) {
            CorrMat.plug[i,j] <- (gammaC.vec[i]*DMat[i,] - CMat[i,])%*%M%*%(gammaC.vec[j]*DMat[j,] - CMat[j,])/
            (sqrt((gammaC.vec[i]*DMat[i,] - CMat[i,])%*%M%*%(gammaC.vec[i]*DMat[i,] - CMat[i,]))*
            sqrt((gammaC.vec[j]*DMat[j,] - CMat[j,])%*%M%*%(gammaC.vec[j]*DMat[j,] - CMat[j,])))
        }
    }


 ##
 ####
 #### Solution of quadratic inequality
 ##
Quad.root <- function(Aj, Bj, Cj){
        Discrimi <- Bj^2 - 4*Aj*Cj   
        if ((Aj > 0)&(Discrimi >= 0)) Limit.s <- (-Bj + plus.minus*sqrt(Discrimi))/(2*Aj)
        else  Limit.s <- "NSD"
        return(Limit.s)}


switch(method,

# UNADJUSTED CI

Unadj={

if (alternative=="two.sided"){ 
    side <- 2
    plus.minus <- c(-1,1)
 
    cpUAd <- qt(1- (1-conf.level)/(side), Degree.f, lower.tail = TRUE)
 
   } # End of two-sided CI
    
if ((alternative=="less")|(alternative=="greater")){
    side <- 1
    if (alternative=="less") plus.minus <- 1
    else plus.minus <- -1
    cpUAd <- qt(1- (1-conf.level)/(side), Degree.f, lower.tail = TRUE)
   
    } # End of one-sided CI    
    



UAdCL <- matrix(as.numeric(rep(NA,side*n.comp)),nrow=n.comp)
for(j in 1:n.comp)
 {
                  AjUAd <- (DMat[j,]%*%Beta.Coeff)^2 - (cpUAd^2)*Pooled.Var*DMat[j,]%*%M%*%DMat[j,]
                  BjUAd <- -2*((CMat[j,]%*%Beta.Coeff)*(DMat[j,]%*%Beta.Coeff) - (cpUAd^2)*Pooled.Var*CMat[j,]%*%M%*%DMat[j,])
                  CjUAd <- (CMat[j,]%*%Beta.Coeff)^2 - (cpUAd^2)*Pooled.Var*CMat[j,]%*%M%*%CMat[j,]
    UAdCL[j,]  <- Quad.root(AjUAd, BjUAd,  CjUAd)
    }
sci.table <- data.frame(UAdCL)
df <- Degree.f; critp <- cpUAd
},

# BONFERRONI-ADJUSTED CI

Bonf={



if (alternative=="two.sided"){ 
    side <- 2
    plus.minus <- c(-1,1)

    cpBon <- qt(1- (1-conf.level)/(side*n.comp), Degree.f, lower.tail = TRUE)

   } # End of two-sided CI
    
if ((alternative=="less")|(alternative=="greater")){
    side <- 1
    if (alternative=="less") plus.minus <- 1
    else plus.minus <- -1

    cpBon <- qt(1- (1-conf.level)/(side*n.comp), Degree.f, lower.tail = TRUE)

    } # End of one-sided CI    
    


BonCL <- matrix(as.numeric(rep(NA,side*n.comp)), nrow=n.comp)
for(j in 1:n.comp)
 {
                  AjBon <- (DMat[j,]%*%Beta.Coeff)^2 - (cpBon^2)*Pooled.Var*DMat[j,]%*%M%*%DMat[j,]
                  BjBon <- -2*((CMat[j,]%*%Beta.Coeff)*(DMat[j,]%*%Beta.Coeff) - (cpBon^2)*Pooled.Var*CMat[j,]%*%M%*%DMat[j,])
                  CjBon <- (CMat[j,]%*%Beta.Coeff)^2 - (cpBon^2)*Pooled.Var*CMat[j,]%*%M%*%CMat[j,]
    BonCL[j,]  <- Quad.root(AjBon, BjBon,  CjBon)
    }
sci.table <- data.frame(BonCL)
df <- Degree.f; critp <- cpBon
},

# MtI-ADJUSTED CI (SIDAK, resp SLEPIAN)

MtI={

if (alternative=="two.sided"){ 
    side <- 2
    plus.minus <- c(-1,1)

    cpMtI <- qmvt(conf.level, interval=c(0,10),df=as.integer(Degree.f),corr=diag(n.comp),delta=rep(0,n.comp), tail="both", abseps=1e-05)$quantile

   } # End of two-sided CI
    
if ((alternative=="less")|(alternative=="greater")){
    side <- 1
    if (alternative=="less") plus.minus <- 1
    else plus.minus <- -1

    cpMtI <- qmvt(conf.level, interval=c(0,10),df=as.integer(Degree.f),corr=diag(n.comp),delta=rep(0,n.comp), 
tail="lower.tail", abseps=1e-05)$quantile

    } # End of one-sided CI    
    


MtICL <- matrix(as.numeric(rep(NA,side*n.comp)),nrow=n.comp)
for(j in 1:n.comp) {   
                  AjMtI <- (DMat[j,]%*%Beta.Coeff)^2 - (cpMtI^2)*Pooled.Var*DMat[j,]%*%M%*%DMat[j,]
                  BjMtI <- -2*((CMat[j,]%*%Beta.Coeff)*(DMat[j,]%*%Beta.Coeff) - (cpMtI^2)*Pooled.Var*CMat[j,]%*%M%*%DMat[j,])
                  CjMtI <- (CMat[j,]%*%Beta.Coeff)^2 - (cpMtI^2)*Pooled.Var*CMat[j,]%*%M%*%CMat[j,]
    MtICL[j,]  <- Quad.root(AjMtI, BjMtI,  CjMtI)
    }
sci.table <- data.frame(MtICL)
df <- as.integer(Degree.f); critp <- cpMtI
},

# PLUG-IN-CI

Plug={


if (alternative=="two.sided"){ 
    side <- 2
    plus.minus <- c(-1,1)
 
    Cplug <- qmvt(conf.level, interval=c(0,10),df=as.integer(Degree.f),corr=CorrMat.plug,delta=rep(0,n.comp), tail="both", abseps=1e-05)$quantile
    
   } # End of two-sided CI
    
if ((alternative=="less")|(alternative=="greater")){
    side <- 1
    if (alternative=="less") plus.minus <- 1
    else plus.minus <- -1

    Cplug <- qmvt(conf.level, interval=c(0,10),df=as.integer(Degree.f),corr=CorrMat.plug,delta=rep(0,n.comp), 
tail="lower.tail", abseps=1e-05)$quantile
   
    } # End of one-sided CI    
   

PlugCL <- matrix(as.numeric(rep(NA,side*n.comp)),nrow=n.comp) 
for(j in 1:n.comp) { 
                  AjPlug <- (DMat[j,]%*%Beta.Coeff)^2 - (Cplug^2)*Pooled.Var*DMat[j,]%*%M%*%DMat[j,]
                  BjPlug <- -2*((CMat[j,]%*%Beta.Coeff)*(DMat[j,]%*%Beta.Coeff) - (Cplug^2)*Pooled.Var*CMat[j,]%*%M%*%DMat[j,])
                  CjPlug <- (CMat[j,]%*%Beta.Coeff)^2 - (Cplug^2)*Pooled.Var*CMat[j,]%*%M%*%CMat[j,]
   PlugCL[j,] <- Quad.root(AjPlug, BjPlug,  CjPlug)

    } 
sci.table <- data.frame(PlugCL)
df <- as.integer(Degree.f); critp <- Cplug
}  
)  
# end of switch method


if (alternative=="two.sided")
{
names(sci.table) <- c("lower","upper")
}

if (alternative=="less")
{
names(sci.table) <- c("upper")
}

if (alternative=="greater")
{
names(sci.table) <- c("lower")
}

if( any(CorrMat.plug<0) && method=="MtI" && alternative!="two.sided")
 {
  warning("At least one element of the estimated correlation matrix is negative,
  therefore, according to Slepian inequality, the MtI method might yield incorrect estimates.")
 }


if (sum(sci.table=="NSD")>0){NSD <- TRUE}
 else{NSD <- FALSE}


 rownames(sci.table)<-compnames
 rownames(gammaC.vec)<-compnames


if(method=="Unadj")
{
methodname<-paste( signif(conf.level*100,2), "% confidence intervals", sep="")
}
else{
methodname<-paste("Simultaneous", signif(conf.level*100,2), "% simultaneous confidence intervals", sep="")
}


out<-list(
estimate=gammaC.vec,
CorrMat.est=CorrMat.plug,
Num.Contrast=CMat,
Den.Contrast=DMat,
conf.int=sci.table,
compnames=compnames,
NSD=NSD,
method=method,
methodname=methodname,
alternative=alternative,
conf.level=conf.level,
type="User defined",
Y=Y,
X=X,
fit=Lm.Fit,
df=df,
quantile=critp
)

class(out)<-c("sci.ratio.gen", "sci.ratio")

return(out)

}  # END of function sci.ratio.gen 




