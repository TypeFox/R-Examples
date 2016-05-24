Wmatgen <-
function(mu,B,Link,Dist){
        if (Dist=="Normal")     Vmat<-rep(1,length(mu))
        if (Dist=="Poisson")    Vmat<-mu
        if (Dist=="Binomial")   Vmat<-(B-mu)*(mu/B)         # In binomial models mu=p*B therefore the transformation is used g(mu/B)=eta #
        if (Dist=="Gamma")      Vmat<-mu^2  
        if (Dist!="Binomial")   B<-1                        # This makes sure offset is not used here if distribution is different than binomial #
                                                            # Include B everywhere and set it to one for different then binomial distribution 3
        #if (Link=="Inverse")    Wvec<-(1/Vmat)
        #if (Link=="Log")        Wvec<-(1/Vmat)*(mu^2)
        #if (Link=="Identity")   Wvec<-(1/Vmat)*rep(1,length(mu))
        #if (Link=="Logit")      Wvec<-(1/Vmat)*
        Wmat<-Vmat
        Wmat
    }
