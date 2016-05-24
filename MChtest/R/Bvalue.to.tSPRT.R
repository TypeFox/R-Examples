"Bvalue.to.tSPRT" <-
function(parms,p0,TOL=10^-8){
    if (is.null(names(parms))){
        m<-parms[1];alpha<-parms[2];e0<-parms[3];e1<-parms[4]
    }
    else{   
        m<-parms["Nmax"]
        alpha<-parms["alpha"]
        e0<-parms["e0"]
        e1<-parms["e1"]
    }   
 
    p1<-p1.given.p0(p0,alpha,TOL)
    OR<- (p1*(1-p0))/(p0*(1-p1))
    C1<- qnorm(1-e0)*sqrt(m*alpha*(1-alpha))
    A<- OR^C1
    C2<-qnorm(e1)*sqrt(m*alpha*(1-alpha))
    B<- OR^C2
    C1<- log(A)/log(OR)
    C2<- log(B)/log(OR)

    a0<-(1-B)/(A-B)
    b0<- (B*(1-A))/(B-A)
     out<-c(p0,p1,A,B,alpha,C1,C2,a0,b0)
    names(out)<-c("p0","p1","A","B","C0","C1","C2","alpha0","beta0")
    out
}

