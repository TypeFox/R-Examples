sep.2z0 <-
function(y, n, xmu.1, p.xmu, xsum.1, p.xsum, x0.1, p.x0, 
                    zdummy, qz,nz0, m, rid, EUID, nEU,
                    prior1, prior2, prior.beta, prior.Sigma, 
                    prec.int, prec.DN, lambda.L1, lambda.L2, lambda.ARD,
                    scale.unif, scale.halft, link, n.chain, inits) 
{
  dataIn <- vector("list",22)
  names(dataIn) <- c("n","y","xmu.1","p.xmu","xsum.1","p.xsum","x0.1",
                     "p.x0","z","nz0","qz","m","cumm", "zero","link",
                     "hyper","prior1","prior2","rid","EUID","nEU","hyper2")
  dataIn[[1]] <- n      
  dataIn[[2]] <- y
  dataIn[[3]] <- as.matrix(xmu.1)
  dataIn[[4]] <- p.xmu
  dataIn[[5]] <- as.matrix(xsum.1)
  dataIn[[6]] <- p.xsum      
  dataIn[[7]] <- as.matrix(x0.1)
  dataIn[[8]] <- p.x0
  dataIn[[9]] <- zdummy
  dataIn[[10]]<- nz0
  dataIn[[11]]<- qz
  dataIn[[12]]<- m
  dataIn[[13]]<- c(0,cumsum(m[-nz0]))   
  dataIn[[14]]<- rep(0,n)
  dataIn[[15]]<- link
  dataIn[[16]] <- as.matrix(cbind(prec.int,prec.DN,lambda.L1,lambda.L2,lambda.ARD))         
  dataIn[[17]] <- prior1
  dataIn[[18]] <- prior2    
  dataIn[[19]] <- rid 
  dataIn[[20]] <- EUID
  dataIn[[21]] <- nEU
  if(grepl("unif",prior.Sigma))  dataIn[[22]] <- scale.unif
  if(grepl("halfcauchy",prior.Sigma)) dataIn[[22]] <- scale.halft        
  
  init <- function( ){
    
    rho1 <- runif(1,-0.5,0.5)
    rho2 <- runif(1,-0.5,0.5) 
    rho3 <- runif(1, rho1*rho2 - sqrt((1-rho1^2)*(1-rho2^2)), 
                  rho1*rho2 + sqrt((1-rho1^2)*(1-rho2^2)))
    return(
    list("tmp1" = rnorm(1,0,0.1),
         "tmp2" = rnorm(1,0,0.1),
         "tmp3" = rnorm(1,0,0.1),
         
         "b.tmp"  = matrix(rnorm((p.xmu-1)*4,0,0.1),ncol=4),
         "d.tmp"  = matrix(rnorm((p.xsum-1)*4,0,0.1),ncol=4),
         "b0.tmp" = matrix(rnorm((p.x0-1)*4,0,0.1), ncol=4),
         
         "sigmab.L1" = runif((p.xmu-1),0,2), 
         "sigmad.L1" = runif((p.xsum-1),0,2), 
         "sigmab0.L1" = runif((p.x0-1),0,2), 
         
         "taub.ARD" = runif((p.xmu-1),0,2), 
         "taud.ARD" = runif((p.xsum-1),0,2), 
         "taub0.ARD" = runif((p.x0-1),0,2), 
         
         "taub.L2" = runif(1,0,2), 
         "taud.L2" = runif(1,0,2),
         "taub0.L2" = runif(1,0,2),
         
         "sigma.VC1" = runif(nz0,0.25,2),
         "t" = runif(nz0,0.25,1),         
         "scale1" = runif(qz,0.25,2),
         "scale2" = runif(qz,0.25,2),
         
         "rho1" = rho1,
         "rho2" = rho2,
         "rho3" = rho3))}  
  
  # 1b, 2d, 3b0, 4d1, 
  # 5 SigmaVC (sigma.VC1 or t),SigmaUN (scale1 or scale2),
  # 6 rho1,2,3
  inits.internal <- list(init( ));
  if(n.chain >= 2) {
    for(j in 2:n.chain) inits.internal <- c(inits.internal,list(init( ))) } 
  
  if(!is.null(inits)){
    for(i in 1:n.chain){
    
    if(!is.null(inits[[i]]$b)) {
      inits.internal[[i]][[1]] <- inits[[i]]$b[1]
      if(p.xmu>=2) inits.internal[[i]][[4]] <- matrix(rep(inits[[i]]$b[2:p.xmu,1],4), 
                                         ncol=4, byrow=FALSE)}
    if(!is.null(inits[[i]]$d)) {
      inits.internal[[i]][[2]] <- inits[[i]]$d[1]
      if(p.xsum>=2) inits.internal[[i]][[5]] <- matrix(rep(inits[[i]]$d[2:p.xsum,1],4), 
                                         ncol=4, byrow=FALSE)}
    if(!is.null(inits[[i]]$b0)) {
      inits.internal[[i]][[3]] <- inits[[i]]$b0[1]
      if(p.x0>=2) inits.internal[[i]][[6]] <- matrix(rep(inits[[i]]$b0[2:p.x0,1],4), 
                                         ncol=4, byrow=FALSE)}

    if(!is.null(inits[[i]]$sigma)) {
      inits.internal[[i]][[16]]<- inits[[i]]$sigma
      inits.internal[[i]][[17]]<- inits[[i]]$sigma
      inits.internal[[i]][[18]]<- runif(qz,0.25,2)
      inits.internal[[i]][[19]]<- runif(qz,0.25,2)
    }
    
    # check PD of the initial R matrix
    if(!is.null(inits[[i]]$R)) {
      notuse <-FALSE
      Rele <- inits[[i]]$R
      size <- (sqrt(1+8*length(Rele))-1)/2 # (# of random effects)
      R <- diag(size)
      R[upper.tri(R, diag=TRUE)] <- Rele 
      R <- R + t(R) - diag(diag(R))
      pd <- all(eigen(R)$values>0)
      if(!pd) {
        notuse <- TRUE
        warning('the specified initial correlation matrix is not positive definite')
        warning('Internal initial value are used')
        break}
      else{
        if(size==2) inits.internal[[i]][[20]] <-inits[[i]]$R[2]
        if(size==3){
          inits.internal[[i]][[20]] <-inits[[i]]$R[2]; 
          inits.internal[[i]][[21]] <-inits[[i]]$R[4]; 
          inits.internal[[i]][[22]] <-inits[[i]]$R[5]}
      }
      lower <- inits.internal[[i]][[20]]*inits.internal[[i]][[21]]-
        sqrt((1-inits.internal[[i]][[20]]^2)*(1-inits.internal[[i]][[21]]^2))
      upper <- inits.internal[[i]][[20]]*inits.internal[[i]][[21]]+
        sqrt((1-inits.internal[[i]][[20]]^2)*(1-inits.internal[[i]][[21]]^2))
      if(inits.internal[[i]][[22]]<lower | inits.internal[[i]][[22]]>upper)
        inits.internal[[i]][[22]] <- runif(1, lower, upper)
    }
  }}
  
  op<- system.file("bugs", "sep_2z0.bug",package="zoib") 
  model<- jags.model(op,data=dataIn,n.adapt=0,inits=inits.internal, n.chains=n.chain)   
  return(model)
}
