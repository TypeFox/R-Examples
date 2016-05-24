"sci.ratio" <-
function(formula, data, type="Dunnett", base=1, method="Plug", Num.Contrast=NULL, Den.Contrast=NULL, 
alternative = "two.sided", conf.level = 0.95, names=TRUE)
 {

method<-match.arg(method, choices=c("Plug", "MtI", "Bonf", "Unadj"))

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
     {stop("number of rows in Num.Contrast and Den.Contrast is not the same")}

    if(ncol(Den.Contrast)!=k || ncol(Num.Contrast)!=k)
     {stop("number of columns in Num.Contrast or Den.Contrast is not the same as number of groups")}

    # only warnings:

  #  NCrowsum<-apply(Num.Contrast, MARGIN=1, FUN=sum)
  #  DCrowsum<-apply(Den.Contrast, MARGIN=1, FUN=sum)
  #  if( !all( NCrowsum == DCrowsum ) )
  #   { cat("Warning: Check whether denominator and numerator contrast matrices are appropriately defined!", "\n") }

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

out<-sci.ratioI(Response=Response, Treatment=Treatment, Num.Contrast=Num.C, Den.Contrast=Den.C, 
alternative=alternative,conf.level=conf.level, method=method)

out$type<-type
out$compnames<-compnames

# # #

colnames(out$Num.Contrast) <- colnames(out$Den.Contrast) <- varnames
rownames(out$conf.int)<-compnames
# # #

if(type=="User defined")
{
 rownames(out$estimate)<-compnames
}

if(method=="Unadj")
{
methodname<-paste( round(conf.level*100,2), "-% confidence intervals", sep="")
}
else{
methodname<-paste("Simultaneous ", round(conf.level*100,4), "% confidence intervals", sep="")
}

out$methodname<-methodname

class(out) <- "sci.ratio"

return(out)

} # END OF sci.ratio

