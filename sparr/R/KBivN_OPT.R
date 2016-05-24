KBivN_OPT <- function(x,WIN=NULL,mu=c(0,0),Sigma=matrix(c(1,0,0,1),2,2),sq=F,forSimul=F,y=NULL){
    #if(is.na(x[3])) return(0)
    
    if(!forSimul){
        u1 <- (x[1]-mu[1])/Sigma[1,1]
        u2 <- (x[2]-mu[2])/Sigma[2,2]
    } else {
        WIN=NULL
        u1 <- (x-mu[1])/Sigma[1,1]
        u2 <- (y-mu[2])/Sigma[2,2]
    }
            
    if(sq){
        if(is.null(WIN)){
            result <- ((1/(2*pi*Sigma[1,1]*Sigma[2,2]))*exp(-.5*(u1^2 + u2^2)))^2
        } else if (!inside.owin(x[1],x[2],WIN)){
            result <- 0
        } else {
            result <- ((1/(2*pi*Sigma[1,1]*Sigma[2,2]))*exp(-.5*(u1^2 + u2^2)))^2
        }
    } else {
        if(is.null(WIN)){
            result <- (1/(2*pi*Sigma[1,1]*Sigma[2,2]))*exp(-.5*(u1^2 + u2^2))
        } else if (!inside.owin(x[1],x[2],WIN)){
            result <- 0
        } else {
            result <- (1/(2*pi*Sigma[1,1]*Sigma[2,2]))*exp(-.5*(u1^2 + u2^2))
        }
    }
    return(list(val=result,u1=u1,u2=u2))
}