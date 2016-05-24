`poly3estf` <-
function(time, status, tmax=NULL, f, method="BP", k=NULL)
{

if(all(c("BP","BW","ADD1","ADD2")!=method))
 {stop("argument method must be one of 'BP','BW','ADD1','ADD2'")}

Ntot <- length(time)
#faclev <- levels(f)

if(length(status)!=Ntot || length(f)!=Ntot)
 {stop("time, status, f must not differ in length")}

if(is.null(tmax))
 {tmax <- max(time)}

if(is.null(k))
 {k<-3}

if(length(k)!=1)
 {k<-k[1]
  warning("Note: Only the first value of k has been used in the calculations!")
 }


timelist <- split(x=time, f=f, drop=TRUE)
statuslist <- split(x=status, f=f, drop=TRUE)

group <- names(timelist)

Ind<-length(timelist)

 Y <- numeric(length=Ind)
 n <- numeric(length=Ind)
 estp <- numeric(length=Ind)
 nadj <- numeric(length=Ind)
 varp <- numeric(length=Ind)
 varcor <- numeric(length=Ind)
 estimate <- numeric(length=Ind)
for (i in 1:length(timelist))
 {
 temp <- poly3est(time=timelist[[i]], status=statuslist[[i]], tmax=tmax, method=method, k=k)

 Y[i] <- temp$Y
 n[i] <- temp$n
 estp[i] <- temp$estp
 nadj[i] <- temp$nadj
 varp[i] <- temp$varp
 varcor[i] <- temp$varcor
 estimate[i] <- temp$estimate
 }

# H0pool <- poly3est(time=time, status=status, tmax=tmax)

# additional information

# earliest death
    mintime<-unlist(lapply(timelist,min))
# median time of death 
    medtime<-unlist(lapply(timelist, median))
# proportion death before tmax
    ntimebe<-lapply(timelist,function(x){length(which(x<tmax))})
    proptimebe<-unlist(ntimebe)/n

names(Y) <- names(n) <- names(estp) <- names(varp) <- names(varcor) <- group

out<-list(
Y=Y,
n=n,
estimate=estimate,
estp=estp,
nadj=nadj,
varp=varp,
varcor=varcor,
k=k,
names=group, 
method=method,
mintime=mintime,
medtime=medtime,
proptimebe=proptimebe
)

class(out)<-"poly3est"
return(out)

}

# Summary function

summary.poly3est<-function(object, ...)
{

args<-list(...)

if(is.null(args$digits) || !is.integer(args$digits))
 {digits<-4}
 else
  {digits<-args$digits}

    if (object$method == "BP") {
        methvar <- "Bailer-Portier"
    }
    if (object$method == "BW") {
        methvar <- "Bieler-Williams"
    }
    if (object$method == "ADD1") {
        methvar <- "Add-1"
    }
    if (object$method == "ADD2") {
        methvar <- "Add-2"
    }


 sample.estimate <- rbind(object$Y, object$n,
  object$Y/object$n,
 round(object$nadj, digits=digits),
 round(object$estimate, digits=digits),
 object$mintime,
 object$medtime,
 round(object$proptimebe,digits=digits) )

    rownames(sample.estimate) <- c("Number of tumours", "Number of animals",
 "Tumour proportion", "Adjusted number of animals",
 "Adjusted proportion",
 "First death",
 "Median time of death",
 "Proportion dead before end")

    cat("Raw and poly-", object$k, "-adjusted sample estimates:", 
        "\n", sep = "")

    print(sample.estimate)
    cat(" ", "\n")
    invisible(object)

}

