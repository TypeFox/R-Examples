fixed01 <-
function(y, n, xmu.1,p.xmu,xsum.1,p.xsum, x0.1,p.x0, x1.1,p.x1, prior1,
         prec.int,prec.DN,lambda.L1,lambda.L2,lambda.ARD,link,n.chain,inits)
{ 
  dataIn <- vector("list",14)  
  dataIn.name <- c("y","xmu.1","p.xmu","xsum.1","p.xsum", "x0.1","p.x0",
                   "x1.1","p.x1","n","zero","prior1","hyper","link")
  names(dataIn)<- dataIn.name  
  dataIn[[1]] <- y
  dataIn[[2]] <- as.matrix(xmu.1)
  dataIn[[3]] <- p.xmu
  dataIn[[4]] <- as.matrix(xsum.1)
  dataIn[[5]] <- p.xsum      
  dataIn[[6]] <- as.matrix(x0.1)
  dataIn[[7]] <- p.x0
  dataIn[[8]] <- as.matrix(x1.1)
  dataIn[[9]] <- p.x1       
  dataIn[[10]]<- n
  dataIn[[11]]<- rep(0,n) 
  dataIn[[12]]<- prior1
  dataIn[[13]]<- as.matrix(cbind(prec.int,prec.DN,lambda.L1,lambda.L2,lambda.ARD))
  dataIn[[14]]<- link
  
  init <- function( ){
    list("tmp1" = rnorm(1,0,0.1),
         "tmp2" = rnorm(1,0,0.1),
         "tmp3" = rnorm(1,0,0.1),
         "tmp4" = rnorm(1,0,0.1), 
         "b.tmp" =  matrix(rnorm((p.xmu-1)*4,0,0.1),ncol=4),
         "d.tmp" =  matrix(rnorm((p.xsum-1)*4,0,0.1),ncol=4),
         "b0.tmp" = matrix(rnorm((p.x0-1)*4,0,0.1),ncol=4),
         "b1.tmp" = matrix(rnorm((p.x1-1)*4,0,0.1),ncol=4),
         "sigmab.L1" =  runif((p.xmu-1),0,2), 
         "sigmad.L1" =  runif((p.xsum-1),0,2), 
         "sigmab1.L1" = runif((p.x1-1),0,2), 
         "sigmab0.L1" = runif((p.x0-1),0,2), 
         "taub.ARD" =  runif((p.xmu-1),0,2), 
         "taud.ARD" =  runif((p.xsum-1),0,2), 
         "taub1.ARD" = runif((p.x1-1),0,2), 
         "taub0.ARD" = runif((p.x0-1),0,2), 
         "taub.L2" =  runif(1,0,2), 
         "taud.L2" =  runif(1,0,2),
         "taub0.L2" = runif(1,0,2),
         "taub1.L2" = runif(1,0,2))}    

  inits.internal <- list(init( ));
  if(n.chain >= 2) {
    for(j in 2:n.chain) inits.internal <- c(inits.internal,list(init( ))) }   
  
  if(!is.null(inits)){
    
    for(i in 1:n.chain){
    
      if(!is.null(inits[[i]]$b)) {
        inits.internal[[i]][[1]] <- inits[[i]]$b[1]
        if(p.xmu>=2) inits.internal[[i]][[5]] <- matrix(rep(inits[[i]]$b[2:p.xmu],4), 
                                           ncol=4, byrow=FALSE)}
      if(!is.null(inits[[i]]$d)) {
        inits.internal[[i]][[2]] <- inits[[i]]$d[1]
        if(p.xsum>=2) inits.internal[[i]][[6]] <- matrix(rep(inits[[i]]$d[2:p.xsum],4), 
                                           ncol=4, byrow=FALSE)}
      if(!is.null(inits[[i]]$b0)) {
        inits.internal[[i]][[3]] <- inits[[i]]$b0[1]
        if(p.x0>=2) inits.internal[[i]][[7]] <- matrix(rep(inits[[i]]$b0[2:p.x0],4), 
                                           ncol=4, byrow=FALSE)}
      if(!is.null(inits[[i]]$b1)) {
        inits.internal[[i]][[4]] <- inits[[i]]$b1[1]
        if(p.x1>=2) inits.internal[[i]][[8]] <- matrix(rep(inits[[i]]$b1[2:p.x1],4), 
                                           ncol=4, byrow=FALSE)}
    }
  }
  op<- system.file("bugs", "fixed01.bug",package="zoib") 
  model<- jags.model(op,data=dataIn,n.adapt=0,inits=inits.internal,n.chains=n.chain)  
  return(model)
}
