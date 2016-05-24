`methodRule1` <-
function(W,Z,exact, parms){

    below.bound<-function(W,Z,maxCalc=parms[1]){
        M<-length(W[W!=0])
        N<-length(W)
        n1<-sum(Z)
        k<- N-M
        num.calcs<- sum( choose(M,max(0,n1-k):min(n1,M)) )
        if (num.calcs<=maxCalc) BB<-TRUE
        else BB<-FALSE
        return(BB)
    }

    if (is.null(exact)){
        if (below.bound(W,Z)) method<-"exact"
        else method<-"approx"
    } else { 
        if (exact==FALSE) method<-"approx"
        else {
            if (below.bound(W,Z)) method<-"exact"
            else method<-"exactMC"
        }
    }
    return(method)
}

