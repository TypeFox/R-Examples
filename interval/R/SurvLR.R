`SurvLR` <-
function(x){
    ## change Surv object to data.frame with L and R columns
    ## type=right, left or interval are allowed, type=counting is not
    type<-attr(x,"type")
    if (type=="right"){
        L<-R<-x[,1]
        R[x[,2]==0]<-Inf
    } else if (type=="counting") {
        stop("Surv object type='counting' not supported")
    } else if (type=="left"){
        L<-R<-x[,1]
        L[x[,2]==0]<-0
    } else if (type=="interval"){
        L<-R<-x[,1]
        R[x[,3]==0]<-Inf
        L[x[,3]==2]<-0
        R[x[,3]==3]<-x[x[,3]==3,2]
    } else { stop(paste("Surv obj type='",type,"' unrecognized",sep=""))
    }
    ## add Lin and Rin to work with right censored data
    Lin<-rep(FALSE,length(L))
    Rin<-rep(TRUE,length(L))
    Lin[L==R]<-TRUE
    Rin[R==Inf]<-FALSE
    out<-data.frame(L=L,R=R,Lin=Lin,Rin=Rin)
    return(out)
}

