theta2z <-
function(theta,family){
    if (family==1){
        if (abs(theta)>0.99) theta<-sign(theta)*0.99
        out<-0.5*log((1+theta)/(1-theta))
    }
    if (family==3){
        if (theta<=0.0001) out<-log(0.0001) else
        out<-log(theta)
    }
    if (family==4){
            if (theta<1) theta<-1
        out<-log(theta-0.99999)
    }
    if (family==5){
        out<-theta
        if (abs(out)<0.0001) out<-0.0001*sign(out)
    }
    return(out)
}
