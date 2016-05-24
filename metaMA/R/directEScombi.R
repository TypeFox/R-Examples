`directEScombi` <-
function(ES,varES,BHth=0.05,useREM=TRUE)
{
listres=vector("list",2)
f.Q.NA=function (dadj, varadj) 
{
    w <- 1/varadj
    tmp1 <- w * dadj
    mu <- rowSums(tmp1,na.rm=TRUE)/rowSums(w,na.rm=TRUE)
    Q <- rowSums(w * (dadj - mu)^2,na.rm=TRUE)
}
tau2.NA <- function(Q, num.studies, my.weights) {
        vwts <- rowSums(my.weights,na.rm=TRUE)
        tmp2 <- rowSums(my.weights^2,na.rm=TRUE)
        tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
        return(tau2)
    }
num.studies <- dim(ES)[2]
Qvals <- f.Q.NA(ES, varES)
if (useREM) { varES <- varES + tau2.NA(Qvals, num.studies, my.weights = 1/varES)}
    wt <- 1/varES
    MUvals <- rowSums(ES * wt,na.rm=TRUE)/rowSums(wt,na.rm=TRUE)
    MUsES <- sqrt(1/rowSums(wt,na.rm=TRUE))
    zSco <- MUvals/MUsES
rpvalESc=2*(1-pnorm(abs(zSco)))
res=which(p.adjust(rpvalESc,method="BH")<=BHth)
listres[[1]]=res
listres[[2]]=zSco 
names(listres)=c("DEindices","TestStatistic")
listres
}

