"simtest.ratio" <-
function(formula, data, type="Dunnett", base=1, alternative = "two.sided", Margin.vec=NULL, FWER=0.05, Num.Contrast=NULL, Den.Contrast=NULL, 
 names=TRUE)
 {


alternative <- match.arg(alternative, choices=c("two.sided","less","greater"))

if(length(formula)!=3){stop("formula mis-specified")}

mf <- model.frame( formula, data)
 if(ncol(mf) != 2) {stop("Specify one response and only one class variable in the formula")}

if (is.numeric(mf[,1])==FALSE){stop("Response variable must be numeric")}

 Response <- mf[,1] 
 Treatment <- droplevels(as.factor(mf[,2]))

varnames <- levels(Treatment)

k <- length(varnames)

splitdat <- split(Response,Treatment)

ni <- as.numeric(lapply(splitdat, FUN=length))


if(any(ni<2))
 {stop("the number of observations in each group should be at least 2")}


#  check appropriateness of user-defined contrasts:


if(is.null(Num.Contrast)==FALSE || is.null(Num.Contrast)==FALSE)
 {
  if(is.null(Den.Contrast==TRUE) && is.null(Den.Contrast==FALSE))
   {stop("Num.Contrast is specified, but Den.Contrast is missing")}

  if(is.null(Den.Contrast==FALSE) && is.null(Den.Contrast==TRUE))
   {stop("Den.Contrast is specified, but Num.Contrast is missing")}

  if(is.null(Den.Contrast)==FALSE && is.null(Num.Contrast)==FALSE)
   {
    if(nrow(Den.Contrast)!=nrow(Num.Contrast))
     {stop("number of rows in Num.Contrast and Den.Contrast should be the same")}

    if(ncol(Den.Contrast)!=k || ncol(Num.Contrast)!=k)
     {stop("number of columns in Num.Contrast or Den.Contrast should be the same as number of groups")}

    # only warnings:
    # for the case that sum of single contrasts are not 1 or 0
    # and the ratio sum(num)/sum(den) is not 1 


   # NCrowsum<-apply(Num.Contrast, MARGIN=1, FUN=sum)
   # DCrowsum<-apply(Den.Contrast, MARGIN=1, FUN=sum)
   # if( !all( NCrowsum == DCrowsum ) )
   #  { cat("Warning: Check whether denominator and numerator contrast matrices are appropriately defined!", "\n") }

    NC0<-apply(X=Num.Contrast, MARGIN=1, function(x){all(x==0)})
    DC0<-apply(X=Den.Contrast, MARGIN=1, function(x){all(x==0)})
    if(any(c(NC0,DC0)))
     {cat("Warning: At least one row of the numerator or denominator contrast matrices is a vector with all components equal to zero","\n")}

     Num.C <- Num.Contrast
     Den.C <- Den.Contrast
     type <- "User defined"
     
    if(is.null(rownames(Num.C)) && is.null(rownames(Den.C)) )
     {compnames <- paste( "C", 1:nrow(Num.C), sep="")}

    else{ if(any(rownames(Num.C)!=rownames(Den.C)) )
     {compnames <- paste(rownames(Num.C), rownames(Den.C), sep="/")}
    else{compnames <- rownames(Num.C)}}
   }

 }
else
{

# Use Contrasts defined by contrMatRatio otherwise: 

type <- match.arg(type, choices=c("Dunnett", "Tukey", "Sequen", "AVE", "GrandMean", "Changepoint", "Marcus", "McDermott", "Williams", "UmbrellaWilliams") )
if(names==TRUE)
 {names(ni)<-varnames}

Cmat<-contrMatRatio(n=ni, type=type, base=base)

Num.C <- Cmat$numC
Den.C <- Cmat$denC
compnames <- Cmat$rnames
}

if(is.null(Margin.vec))
 {Margin.vec <- rep(1,nrow(Num.C))}
else
 {
  if(is.numeric(Margin.vec) && length(Margin.vec)<=nrow(Num.C))
   {Margin.vec <- cbind(Margin.vec,Num.C)[,1]}
  else{
    Margin.vec <- Margin.vec[1:nrow(Num.C)]
    warning( paste("Margin.vec has more elements than there are comparisons. Only the first ", nrow(Num.C)," elements are used!") )
    }
 }

out<-simtest.ratioI(Response=Response, Treatment=Treatment, alternative=alternative, Margin.vec=Margin.vec, FWER=FWER, Num.Contrast=Num.C, Den.Contrast=Den.C)

out$type<-type
out$compnames<-compnames

colnames(out$Num.Contrast)<-varnames
colnames(out$Den.Contrast)<-varnames

 names(out$p.value.raw)<-compnames
 names(out$p.value.adj)<-compnames
 names(out$estimate)<-compnames
 names(out$teststat)<-compnames

out$methodname<-"Tests for ratios of means assuming homogeneous variances \n"

class(out) <- "simtest.ratio"

return(out)

} # END OF sci.ratio

