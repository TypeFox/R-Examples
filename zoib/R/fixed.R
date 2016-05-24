fixed <-
function(y, n, xmu.1,p.xmu,xsum.1,p.xsum, prior1, prec.int,prec.DN, 
         lambda.L1, lambda.L2, lambda.ARD, link,n.chain, inits)
{
  dataIn <- vector("list",10)  
  dataIn.name <- c("y","xmu.1","p.xmu","xsum.1","p.xsum",
                   "n","zero","prior1","hyper","link")
  names(dataIn)<- dataIn.name  
  dataIn[[1]] <- y
  dataIn[[2]] <- as.matrix(xmu.1)
  dataIn[[3]] <- p.xmu
  dataIn[[4]] <- as.matrix(xsum.1)
  dataIn[[5]] <- p.xsum      
  dataIn[[6]] <- n
  dataIn[[7]] <- rep(0,n)   
  dataIn[[8]] <- prior1
  dataIn[[9]] <- as.matrix(cbind(prec.int,prec.DN,lambda.L1,lambda.L2,lambda.ARD))
  dataIn[[10]] <- link
  
  init <- function( ){
    list("tmp1" = rnorm(1,0,0.1),
         "tmp2" = rnorm(1,0,0.1),
         "b.tmp" = matrix(rnorm((p.xmu-1)*4,0,0.1),ncol=4),
         "d.tmp" = matrix(rnorm((p.xsum-1)*4,0,0.1),ncol=4),
         "sigmab.L1" = runif((p.xmu-1),0,1), 
         "sigmad.L1" = runif((p.xsum-1),0,1), 
         "taub.ARD" = runif((p.xmu-1),0,1), 
         "taud.ARD" = runif((p.xsum-1),0,1), 
         "taub.L2" = runif(1,0,1), 
         "taud.L2" = runif(1,0,1))}  
  
  # 1b, 2d
  inits.internal <- list(init( ));
  if(n.chain >= 2) {
    for(j in 2:n.chain) inits.internal <- c(inits.internal,list(init( ))) }   
  
  if(!is.null(inits)){
    for(i in 1:n.chain){
      
      if(!is.null(inits[[i]]$b)) {
        inits.internal[[i]][[1]] <- inits[[i]]$b[1]
        if(p.xmu>=2) inits.internal[[i]][[3]] <- matrix(rep(inits[[i]]$b[2:p.xmu],4), 
                                           ncol=4, byrow=FALSE)}
      if(!is.null(inits[[i]]$d)) {
        inits.internal[[i]][[2]] <- inits[[i]]$d[1]
        if(p.xsum>=2) inits.internal[[i]][[4]] <- matrix(rep(inits[[i]]$d[2:p.xsum],4), 
                                           ncol=4, byrow=FALSE)}
    }}
  op<- system.file("bugs","fixed.bug",package="zoib") 
  model <- jags.model(op,data=dataIn, n.adapt=0, inits=inits.internal, n.chains=n.chain)
  return(model)
}
