`methodRuleIC1` <-
function(x,group,exact, Nbound=c(20)){
    n<-length(x)

    if (is.null(exact)){
        if (n<=Nbound) method<-"exact.network"
        else method<-"pclt"
    } else { 
        if (exact==FALSE){ method<-"pclt"
        } else if (n<=Nbound){  method<-"exact.network"
        } else { method<-"exact.mc" }
    }

    if (length(unique(group))!=2 && method=="exact.network") method<-"exact.mc"
 
    return(method)
}
