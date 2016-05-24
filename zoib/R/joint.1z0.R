joint.1z0 <-
function(y, n, q, xmu.1, p.xmu, xsum.1, p.xsum, x0.1, p.x0, inflate0,  
         rid, EUID, nEU, prior1, prior2, prior.beta, prior.Sigma, 
         prec.int, prec.DN, lambda.L1, lambda.L2,lambda.ARD,
         scale.unif, scale.halft, link, n.chain, inits) 
{ 
  dataIn <- vector("list",19)
  dataIn.name <- c("n","y","q","xmu.1","p.xmu","xsum.1","p.xsum","x0.1","p.x0",
                   "inflate0","zero","link","hyper","prior1","prior2","rid",
                   "EUID","nEU","hyper2")
  names(dataIn)<- dataIn.name 
  dataIn[[1]] <- n      
  dataIn[[2]] <- as.matrix(y)
  dataIn[[3]] <- q
  dataIn[[4]] <- as.matrix(xmu.1)
  dataIn[[5]] <- p.xmu
  dataIn[[6]] <- as.matrix(xsum.1) 
  dataIn[[7]] <- p.xsum      
  dataIn[[8]] <- as.matrix(x0.1)
  dataIn[[9]] <- p.x0  
  dataIn[[10]]<- inflate0
  dataIn[[11]]<- matrix(0,n,q)  
  dataIn[[12]]<- link 
  dataIn[[13]]<- abind(prec.int,prec.DN,lambda.L1,lambda.L2,lambda.ARD,along=3)    
  if(grepl("unif", prior.Sigma)) dataIn[[19]] <- scale.unif
  if(grepl("halfcauchy",prior.Sigma)) dataIn[[19]] <- scale.halft                    
  dataIn[[14]] <- prior1
  dataIn[[15]] <- prior2    
  dataIn[[16]] <- rid 
  dataIn[[17]] <- EUID
  dataIn[[18]] <- nEU
  
  init <- function( ){
    list("tmp1" = rnorm(q,0,0.1),
         "tmp2" = rnorm(q,0,0.1),
         "tmp3" = rnorm(q,0,0.1),
         
         "b.tmp" = array(rnorm((p.xmu-1)*4*q,0,0.1), c((p.xmu-1),q,4)),
         "d.tmp" = array(rnorm((p.xsum-1)*4*q,0,0.1),c((p.xsum-1),q,4)),
         "b0.tmp"= array(rnorm((p.x0-1)*4*q,0,0.1),  c((p.x0-1),q,4)),
         
         "sigmab.L1" =  matrix(runif((p.xmu-1)*q,0,2),(p.xmu-1),q), 
         "sigmad.L1" =  matrix(runif((p.xsum-1)*q,0,2),(p.xmu-1),q),  
         "sigmab0.L1" = matrix(runif((p.x0-1)*q,0,2),(p.xmu-1),q),  
         
         "taub.ARD" =  matrix(runif((p.xmu-1)*q,0,2), (p.xmu-1),q), 
         "taud.ARD" =  matrix(runif((p.xsum-1)*q,0,2),(p.xmu-1),q),  
         "taub0.ARD" = matrix(runif((p.x0-1)*q,0 ,2),(p.xmu-1),q),  
         
         "taub.L2" =  runif(q,0,2), 
         "taud.L2" =  runif(q,0,2),
         "taub0.L2" = runif(q,0,2),
         
         "sigma1" = runif(1,0.25,1),
         "scale2" = runif(1,0.25,1))}    
  
  # 1b, 2d, 3b0, 4d1, 
  # 5 sigma1 scale2
  inits.internal <- list(init( ));
  if(n.chain >= 2) {
    for(j in 2:n.chain) inits.internal <- c(inits.internal,list(init( ))) }   
  
  if(!is.null(inits)){
    
  for(i in 1:n.chain){
    
    if(!is.null(inits[[i]]$b)) {
      inits.internal[[i]][[1]] <- inits[[i]]$b[1,]
      if(p.xmu>=2) inits.internal[[i]][[4]] <- array(rep(inits[[i]]$b[2:p.xmu,],4), c((p.xmu-1),q,4))}
    if(!is.null(inits[[i]]$d)) {
      inits.internal[[i]][[2]] <- inits[[i]]$d[1,]
      if(p.xsum>=2) inits.internal[[i]][[5]] <- array(rep(inits[[i]]$d[2:p.xsum,],4), c((p.xsum-1),q,4))}
    if(!is.null(inits[[i]]$b0)) {
      inits.internal[[i]][[3]] <- inits[[i]]$b0[1,]
      if(p.x0>=2) inits.internal[[i]][[6]] <- array(rep(inits[[i]]$b0[2:p.x0,],4), c((p.x0-1),q,4))}
    if(!is.null(inits[[i]]$sigma)) {
      inits.internal[[i]][[16]]<- inits[[i]]$sigma
      inits.internal[[i]][[17]]<- inits[[i]]$sigma}
  }}
  op<- system.file("bugs", "joint_1z0.bug", package="zoib") 
  model <- jags.model(op,n.adapt=0, data=dataIn, inits=inits.internal, n.chains=n.chain)   
  return(model)
}
