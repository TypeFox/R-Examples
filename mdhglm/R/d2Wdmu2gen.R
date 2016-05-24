d2Wdmu2gen <-
function(mu,B,Link,Dist){
        mu1<-mu/B
        if (Dist=="Normal") {
            Vmat<-rep(1,length(mu))
            dVmatdmu<-rep(0,length(mu))
            d2Vmatdmu2<-rep(0,length(mu))
        }
        if (Dist=="Poisson") {
            Vmat<-mu
            dVmatdmu<-rep(1,length(mu))
            d2Vmatdmu2<-rep(0,length(mu))
        }
        if (Dist=="Binomial") {
            Vmat<-(B-mu)*(mu/B)
            dVmatdmu<-1-2*(mu/B)
            d2Vmatdmu2<--2*(1/B)
        }
        if (Dist=="Gamma") {
            Vmat<-mu^2
            dVmatdmu<-2*mu
            d2Vmatdmu2<-2
        }
        if (Dist!="Binomial") B<-1
        
        if (Link=="Inverse")    {
            detadmu <- 1/(mu^2)
            d2etadmu2 <- -2/(mu^3)
            d3etadmu3 <- 6/(mu^4)
        }    
        if (Link=="Log")        {
            detadmu<-1/mu
            d2etadmu2<--1/(mu^2)
            d3etadmu3<-2/(mu^3)
            
        }
        if (Link=="Identity")   {
            detadmu<-rep(1,length(mu))
            d2etadmu2<-rep(0,length(mu))
            d3etadmu3<-rep(0,length(mu))
            
        }    
        if (Link=="Logit")      {
            detadmu<-1/(mu*(1-mu1))
            d2etadmu2<--(1-2*mu1)/((mu*(1-mu1))^2)
            d3etadmu3<-((2/B)*((mu*(1-mu1))^2)+2*(1-2*mu1)*mu*(1-2*mu1)*(1-mu1))/(mu*(1-mu1))^4
        }
        
        # Add d2Vmatdmu2 and d3etadmu3 to all the functions #
        d2Wdmu2<-2*(1/Vmat^3)*(dVmatdmu^2)*((1/detadmu)^2)-(1/Vmat^2)*d2Vmatdmu2*((1/detadmu)^2)+2*(1/Vmat^2)*dVmatdmu*((1/detadmu)^3)*(d2etadmu2)-
                    2*(1/Vmat^2)*(dVmatdmu)*(1/detadmu)*(-1/detadmu^2)*d2etadmu2-2*(1/Vmat)*(1/detadmu^2)*d2etadmu2*(-1/detadmu^2)*d2etadmu2-
                    4*(1/Vmat)*(1/detadmu)*(-1/detadmu^3)*(d2etadmu2^2)+2*(1/Vmat)*(1/detadmu)*(-1/detadmu^2)*d3etadmu3 
        return(d2Wdmu2)     
    }
