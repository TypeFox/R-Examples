"tSPRT.to.Bvalue" <-
function(parms){
   if (is.null(names(parms))){
       p0<-parms[1]
       p1<-parms[2]
       A<-parms[3]
       B<-parms[4]
       m<-parms[5]
    }
    else{
        p0<-parms["p0"]
        p1<-parms["p1"]
        m<-parms["Nmax"]
        if (is.na(parms["A"])){
            alpha0<-parms["alpha0"]
            beta0<-parms["beta0"]
            A<-(1-beta0)/alpha0
            B<-beta0/(1-alpha0)
        }
        else{
            A<-parms["A"]
            B<-parms["B"]
        }
    }
    #print(paste("A=",A,"1/A=",1/A,"B=",B))
    log.OR<- log((p1*(1-p0))/(p0*(1-p1)))
    ALPHA<- log( (1-p0)/(1-p1) )/log.OR
    if (log.OR==0){ stop("p0 cannot equal p1") } 
    else if (log.OR>0){
        C1<- log(A)/log.OR
        C2<- log(B)/log.OR
    }
    else if (log.OR<0){
        C1<- log(B)/log.OR
        C2<- log(A)/log.OR
    }
    e0<- 1- pnorm(C1/sqrt(m*ALPHA*(1-ALPHA)))
    e1<-pnorm(  C2/sqrt( m*ALPHA*(1-ALPHA) ) )
    out<-c(m,ALPHA,e0,e1)
    names(out)<-c("Nmax","alpha","e0","e1")
    out
}

