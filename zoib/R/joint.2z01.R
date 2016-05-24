joint.2z01 <-
function(y, n, q, xmu.1, p.xmu, xsum.1, p.xsum, x1.1, p.x1, x0.1, p.x0,
                       inflate0, inflate1, zdummy, qz,nz0, m, rid, EUID, nEU,
                       prior1, prior2, prior.beta, prior.Sigma, 
                       prec.int, prec.DN, lambda.L1, lambda.L2,lambda.ARD,
                       scale.unif, scale.halft, link, n.chain, inits) 
{
  dataIn <- vector("list",27)
  names(dataIn) <- c("y","n","q","xmu.1","p.xmu","xsum.1","p.xsum","x0.1",
                     "p.x0","x1.1","p.x1","inflate0","inflate1","z","nz0",
                     "qz","m","cumm", "zero","link","hyper","prior1",
                     "prior2", "rid","EUID","nEU","hyper2")
  dataIn[[1]] <- y
  dataIn[[2]] <- n      
  dataIn[[3]] <- q
  dataIn[[4]] <- as.matrix(xmu.1)
  dataIn[[5]] <- p.xmu
  dataIn[[6]] <- as.matrix(xsum.1)
  dataIn[[7]] <- p.xsum      
  dataIn[[8]] <- as.matrix(x0.1)
  dataIn[[9]] <- p.x0
  dataIn[[10]]<- as.matrix(x1.1)
  dataIn[[11]]<- p.x1
  dataIn[[12]]<- inflate0
  dataIn[[13]]<- inflate1  
  dataIn[[14]]<- zdummy
  dataIn[[15]]<- nz0
  dataIn[[16]]<- qz
  dataIn[[17]]<- m
  dataIn[[18]]<- c(0,cumsum(m[-nz0]))   
  dataIn[[19]]<- rep(0,n)
  dataIn[[20]]<- link
  dataIn[[21]]<- abind(prec.int,prec.DN,lambda.L1,lambda.L2,lambda.ARD,along=3)    
  dataIn[[22]] <- prior1
  dataIn[[23]] <- prior2 
  dataIn[[24]] <- rid
  dataIn[[25]] <- EUID
  dataIn[[26]] <- nEU
  if(grepl("unif", prior.Sigma)) dataIn[[27]] <- scale.unif
  if(grepl("halfcauchy",prior.Sigma)) dataIn[[27]] <- scale.halft                
  
  init <- function( ){
    rho1 <- runif(1,-0.5,0.5)
    rho2 <- runif(1,-0.5,0.5) 
    # to ensure R is Positive definite
    rho3 <- runif(1, rho1*rho2 - sqrt((1-rho1^2)*(1-rho2^2)), 
                  rho1*rho2 + sqrt((1-rho1^2)*(1-rho2^2)))
    
    list("tmp1" = rnorm(q,0,0.1),
         "tmp2" = rnorm(q,0,0.1),
         "tmp3" = rnorm(q,0,0.1),
         "tmp4" = rnorm(q,0,0.1),
         
         "b.tmp" = matrix(rnorm((p.xmu-1)*4,0,0.1),ncol=4),
         "d.tmp" = matrix(rnorm((p.xsum-1)*4,0,0.1),ncol=4),
         "b0.tmp" = matrix(rnorm((p.x0-1)*4,0,0.1),ncol=4),
         "b1.tmp" = matrix(rnorm((p.x1-1)*4,0,0.1),ncol=4),
         
         "sigmab.L1" = runif((p.xmu-1),0,2), 
         "sigmad.L1" = runif((p.xsum-1),0,2), 
         "sigmab0.L1" = runif((p.x0-1),0,2), 
         "sigmab1.L1" = runif((p.x1-1),0,2), 
         
         "taub.ARD" = runif((p.xmu-1),0,2), 
         "taud.ARD" = runif((p.xsum-1),0,2), 
         "taub0.ARD" = runif((p.x0-1),0,2), 
         "taub1.ARD" = runif((p.x1-1),0,2), 
         
         "taub.L2" = runif(1,0,2), 
         "taud.L2" = runif(1,0,2),
         "taub0.L2" = runif(1,0,2),
         "taub1.L2" = runif(1,0,2),
         
         "sigma.VC1" = runif(nz0,0.25,2),
         "t" = runif(nz0,0.25,1),     
         "scale1" = runif(qz,0.25,2),
         "scale2" = runif(qz,0.25,2),
         
         "rho1" = rho1,
         "rho2" = rho2,
         "rho3" = rho3 )} 
  
  # 1b, 2d, 3b0, 4d1, 
  # 5 SigmaVC (sigma.VC1 or t),SigmaUN (scale1 or scale2),
  # 6 rho1,2,3
  inits.internal <- list(init( ));
  if(n.chain >= 2) {
    for(j in 2:n.chain) inits.internal <- c(inits.internal,list(init( ))) }   
  
  if(!is.null(inits)){
    
  for(i in 1:n.chain){
    
    if(!is.null(inits[[i]]$b)) {
      inits.internal[[i]][[1]] <- inits[[i]]$b[1,]
      if(p.xmu>=2) inits.internal[[i]][[5]] <- array(rep(inits[[i]]$b[2:p.xmu,],4), c((p.xmu-1),q,4))}
    if(!is.null(inits[[i]]$d)) {
      inits.internal[[i]][[2]] <- inits[[i]]$d[1,]
      if(p.xsum>=2) inits.internal[[i]][[6]] <- array(rep(inits[[i]]$d[2:p.xsum,],4), c((p.xsum-1),q,4))}
    if(!is.null(inits[[i]]$b0)) {
      inits.internal[[i]][[3]] <- inits[[i]]$b0[1,]
      if(p.x0>=2) inits.internal[[i]][[7]] <- array(rep(inits[[i]]$b0[2:p.x0,],4), c((p.x0-1),q,4))}
    if(!is.null(inits[[i]]$b1)) {
      inits.internal[[i]][[4]] <- inits[[i]]$b1[1,]
      if(p.x1>=2) inits.internal[[i]][[8]] <- array(rep(inits[[i]]$b1[2:p.x1,],4), c((p.x1-1),q,4))}
    
    if(!is.null(inits[[i]]$sigma)) {
      inits.internal[[i]][[21]]<- inits[[i]]$sigma
      inits.internal[[i]][[22]]<- inits[[i]]$sigma
      inits.internal[[i]][[23]]<- runif(qz,0.25,2)
      inits.internal[[i]][[24]]<- runif(qz,0.25,2)
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
        if(size==2) inits.internal[[i]][[25]] <-inits[[i]]$R[2]
        if(size==3){
          inits.internal[[i]][[25]] <-inits[[i]]$R[2]; 
          inits.internal[[i]][[26]] <-inits[[i]]$R[4]; 
          inits.internal[[i]][[27]] <-inits[[i]]$R[5]}
      }
    }
    lower <- inits.internal[[i]][[25]]*inits.internal[[i]][[26]]-
      sqrt((1-inits.internal[[i]][[25]]^2)*(1-inits.internal[[i]][[26]]^2))
    upper <- inits.internal[[i]][[25]]*inits.internal[[i]][[26]]+
      sqrt((1-inits.internal[[i]][[25]]^2)*(1-inits.internal[[i]][[26]]^2))
    if(inits.internal[[i]][[27]]<lower | inits.internal[[i]][[27]]>upper)
      inits.internal[[i]][[27]] <- runif(1, lower, upper)
  }}
  
  op<- system.file("bugs", "joint_2z01.bug", package="zoib") 
  model <- jags.model(op,data=dataIn, n.adapt=0, inits=inits.internal, n.chains=n.chain)   
  return(model)
}
