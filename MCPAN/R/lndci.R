lndci <-
function(x, f, type="Dunnett", cmat=NULL, alternative=c("two.sided","less","greater"), conf.level=0.95, method=c("GPQ","AN"), ...)
{
# additional arguments
aargs<-list(...)

# check arguments
if(!is.numeric(conf.level) || (length(conf.level)!=1| conf.level>=1| conf.level<0.5)){stop("conf.level should be a single numeric value smaller than 1 and larger than 0.5")}

alternative <- match.arg(alternative)
method <- match.arg(method)

if(!is.null(aargs$B) & method=="AN"){warning("Additional argument B will be ignored with method='AN'.")}
if(!is.null(aargs$dist) & method=="GPQ"){warning("Additional argument dist will be ignored with method='GPQ'.")}
wc <- which(names(aargs) %in% c("dist","B"))
if(length(aargs)>length(wc)){
if(length(wc)==0){warning(paste("Argument(s)", paste(names(aargs), collapse=", ") , "will be ignored."))}
if(length(wc)>0){warning(paste("Argument(s)", paste(names(aargs[-wc]), collapse=", ") , "will be ignored."))}
}

if(!is.null(aargs$B) & method=="GPQ"){
if(!is.numeric(aargs$B) & !is.integer(aargs$B)){stop("B must be a single integer value.")}
if(length(aargs$B)!=1){stop("B must be a single integer value.")}
}

ff <- droplevels(as.factor(f))
ni <- tapply(ff,ff, length)
ngroups<-length(ni)

if(any(ni<=3)){warning("At least one group has only 3 or less observations! Methods may be unreliable.")}

if(any(x<=0)){warning("Values of the response variable are less or equal 0!")}

    if (is.null(cmat)) {
      if(type=="Dunnett") {
        if(is.null(aargs$base)){base<-1}
        else{base<-aargs$base}
        cmat <- contrMat(n=ni, type=type, base=base)
       }
       else{cmat <- contrMat(n = ni, type = type)}
    }
    else {
        if (!is.matrix(cmat) || ncol(cmat) != ngroups)
         {stop("cmat must be a matrix with number of columns = number of groups")}
    }

EST <- lndest(x=x, f=f)
LEST <-  cmat %*% EST$estimate

# calculate SCI acc. to diff methods

switch(method,
"AN"={
if(is.null(aargs$dist)){dist<-"MVN"}else{dist<-aargs$dist}
SCI <- Waldci(cmat=cmat, estp=EST$estimate, varp=EST$varest, varcor=EST$varest, alternative = alternative, conf.level = conf.level, dist = dist) 
},
"GPQ"={
if(is.null(aargs$B)){B<-10000}else{B<-aargs$B}

lx<-log(x)
ldatlist<-split(x=lx, f=ff)

ni <- unlist(lapply(ldatlist, length))
mli <- unlist(lapply(ldatlist, mean))
varli <- unlist(lapply(ldatlist, var))

mi <- exp(mli + 0.5 * varli)
mimat<-matrix(mi, nrow=length(mi))

estimate<-cmat%*%mimat

PARMAT<-cbind(ni, mli, varli)

PQMAT<-apply(X=PARMAT, MARGIN=1, FUN=function(x){
PQln(nx=x[1], mlx=x[2], varlx=x[3], B=B)})

ePQMAT<-exp(PQMAT)

PQdiffs<-apply(X=ePQMAT, MARGIN=1, FUN=function(x){cmat %*% x})

SCI<-SCSrank(x=t(PQdiffs), alternative=alternative, conf.level=conf.level)

})

OUT<-c(list(estimate=LEST), SCI)

return(OUT)
}

