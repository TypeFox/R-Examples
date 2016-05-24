dWdmugen <-
function(mu,B,Link,Dist){
        mu1<-mu/B
        if (Dist=="Normal") {
            Vmat<-rep(1,length(mu))
            dVmatdmu<-rep(0,length(mu))
        }
        if (Dist=="Poisson") {
            Vmat<-mu
            dVmatdmu<-rep(1,length(mu))
        }
        if (Dist=="Binomial") {
            Vmat<-(B-mu)*(mu/B)
            dVmatdmu<-1-2*(mu/B)
        }
        if (Dist=="Gamma") {
            Vmat<-mu^2
            dVmatdmu<-2*mu
        }
        if (Dist!="Binomial") B<-1
        
        if (Link=="Inverse")    {
            detadmu <- 1/(mu^2)
            d2etadmu2 <- -2/(mu^3)
        }    
        if (Link=="Log")        {
            detadmu<-1/mu
            d2etadmu2<--1/(mu^2)
            
        }
        if (Link=="Identity")   {
            detadmu<-rep(1,length(mu))
            d2etadmu2<-rep(0,length(mu))
            
        }    
        if (Link=="Logit")      {
            detadmu<-1/(mu*(1-mu1))
            d2etadmu2<--(1-2*mu1)/((mu*(1-mu1))^2)
        }
        
        dWdmu<--(1/Vmat^2)*dVmatdmu*((1/detadmu)^2)+2*(1/Vmat)*(1/detadmu)*(-1/detadmu^2)*d2etadmu2       
        dWdmu
    }
