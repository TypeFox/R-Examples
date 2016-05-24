"rsextrobj" <-
function(RSobj,start=1,end=8192)
{
    obj.name <- deparse(substitute(RSobj))
    if (!(class(RSobj)=="RSmpl" || class(RSobj)=="RSmplext")){
         err.text<-paste(obj.name,"not a sample object - see help(\"rsextrobj\")",sep ="",collapse="")
         stop(err.text)
    }

    n_tot  <- RSobj$n_tot
    if (end>n_tot) end<-n_tot
    n      <- RSobj$n
    k      <- RSobj$k
    nwords <- c(trunc((k+31)/32))

    objnew <- RSobj
    l_one_mat <- n*nwords
    b <- (start-1)*l_one_mat+1
    e <- end*l_one_mat
    objnew$outvec <- RSobj$outvec[b:e]
    objnew$n_tot <- end-start+1
    if (start==1) {
         objnew$n_eff <- objnew$n_tot - 1
    } else {
         objnew$n_eff <- objnew$n_tot
    }
    class(objnew)="RSmplext"

    RET<-objnew
    RET
}

