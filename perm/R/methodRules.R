`methodRuleTS1` <-
function(x,group,exact, Nbound=c(1000,200,100,50,16)){
    ug<-sort(unique(group))
    if (length(ug)!=2) stop("method rule function works for two groups only")

    ## use a different Nbound depending on the size of the minimum of the sample size of the two groups
    below.bound<-function(x,Nbound){
        n<-length(x)
        kmin<-min(length(unique(x[group==ug[1]])),length(unique(x[group==ug[2]])))
        last.Nbound<-length(Nbound)
        if (kmin-1<=last.Nbound){
            if (n<=Nbound[kmin-1]) BB<-TRUE
            else BB<-FALSE
        } else{
            if (n<=min(Nbound)) BB<-TRUE
            else BB<-FALSE
        }  
        return(BB)
    }

    if (is.null(exact)){
        if (below.bound(x,Nbound)) method<-"exact.network"
        else method<-"pclt"
    } else { 
        if (exact==FALSE) method<-"pclt"
        else {
            if (below.bound(x,Nbound)) method<-"exact.network"
            else method<-"exact.mc"
        }
    }
    return(method)
}

`methodRuleTREND1` <-
function(x,y,exact, Nbound=c(20)){

    below.bound<-function(x,nbound=Nbound){
        n<-length(x)
        if (n<=nbound) BB<-TRUE
        else BB<-FALSE
        return(BB)
    }

    if (is.null(exact)){
        if (below.bound(x,Nbound)) method<-"exact.mc"
        else method<-"pclt"
    } else { 
        if (exact==FALSE){ method<-"pclt"
        } else  method<-"exact.mc"
    }
    return(method)
}

`methodRuleKS1` <-
function(x,group,exact, Nbound=c(5)){

    below.bound<-function(x,g,nbound=Nbound){
        xg<-split(x,g)
        N<-sapply(xg,length)
        
        if (min(N)<=nbound) BB<-TRUE
        else BB<-FALSE
        return(BB)
    }

    if (is.null(exact)){
        if (below.bound(x,group,Nbound)) method<-"exact.mc"
        else method<-"pclt"
    } else { 
        if (exact==FALSE){ method<-"pclt"
        } else  method<-"exact.mc"
    }
    return(method)
}


