`cm.compact` <-
function(n0,n1,M){
    K<-n1+n0-M
    ## instead of using h=number of zeros in vaccine group, we use
    ## g=number of non-zeros in vaccine group, so 
    ## g=n1-h
    ## Since h goes from max(0,n1-M) to min(n1,K)
    ## g=n1-h goes from n1-min(n1,K)  = max(0,n1-K)
    ##               to n1-max(0,n1-M)= min(n1,M)
    lowest<- max(0,n1-K)
    highest<-min(n1,M)
    gvals<-lowest:highest
    ng<-length(gvals)
    for (g in 1:ng){
         if (g==1){ 
            cmout<-chooseMatrix(M,gvals[g])
            #weight<-rep(choose(K,n1-gvals[g]),dim(cmout)[1]) 
            weight<-rep(dhyper(gvals[g],M,K,n1)/choose(M,gvals[g]),dim(cmout)[1]) 
            #gout<-rep(gvals[g],dim(cmout)[1])
        }
        else{ 
            cm.temp<-chooseMatrix(M,gvals[g])
            #weight<-c(weight,
            #    rep(choose(K,n1-gvals[g]),dim(cm.temp)[1]) )
            weight<-c(weight,
                rep(dhyper(gvals[g],M,K,n1)/choose(M,gvals[g]),dim(cm.temp)[1]) )
            #gout<-c(gout, rep(gvals[g],dim(cm.temp)[1]))
            cmout<-rbind(cmout,cm.temp) 
        }
    }
    out<-list(weight=weight,cm=cmout)
    #out<-list(g=gout,weight=weight,cm=cmout)
    out
}

