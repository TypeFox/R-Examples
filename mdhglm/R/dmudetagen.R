dmudetagen <-
function(mu,B,Link,Dist){
        if (Link=="Inverse")    dmudeta<-mu^2
        if (Link=="Log")        dmudeta<-mu
        if (Link=="Identity")   dmudeta<-rep(1,length(mu))
        if (Link=="Logit")      dmudeta<-(B-mu)*(mu/B)
        dmudeta
    }
