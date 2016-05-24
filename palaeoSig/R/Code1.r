# change name of a and g to alpha and gamma
# alpha und gamma k?nnen f?r alle Arten unterschiedlich sein
# statt spacing k?nnte man auch ein Argument fun einf?hren, f?r regular eine Funktion schreiben und sonst runif oder rnorm

"make.species" <- function(nspp=30, Amax,srange, fun,xpar, alpha=4,gamma=4)     
{
  #if(missing(xpar)) xpar=c(-50,150)
  #maximum Ao
  if(is.numeric(Amax)) Ao=Amax
  else Ao=Amax(n=nspp)
  r <- srange
  #else gamma=gamma
  #mean position m
  if(missing(fun)) m <- seq(xpar[1],xpar[2],length.out = nspp)
  else m <- fun(nspp,xpar[1],xpar[2])            
  spp=data.frame(Ao=Ao, m=m, r=r, alpha=alpha, gamma=gamma)
  list(spp = spp)
}

species <- function(nspp=30, Amax,fun, xpar, srange, alpha = 4,gamma= 4, ndim, sdistr,  ocor, odistr){

if(length(xpar) == 2) xpar = matrix(ncol = 2, rep(xpar,ndim), byrow = TRUE)
if(length(alpha) == 1) alpha <- matrix(ncol = ndim,rep(alpha,ndim*nspp))
if(length(gamma) == 1) gamma <- matrix(ncol = ndim,rep(gamma,ndim*nspp))
if(length(srange) == 1) srange <- matrix(ncol = ndim,rep(srange,ndim*nspp))
if(length(Amax) == 1 & is.numeric(Amax)) srange <- matrix(ncol = 1,rep(Amax,nspp))

if(missing(fun)) { 
spp <- lapply(1:ndim, function(x) {                                      
  spp <- 'make.species'(nspp=nspp,Amax = Amax, srange = srange[,x], xpar=xpar[x,],alpha =  alpha[,x],gamma = gamma[,x])
  spp<-spp$spp[sample(1:nspp),]
  })} 
else{spp <- lapply(1:ndim, function(x) {                                      
  spp <- 'make.species'(nspp = nspp ,Amax = Amax,fun = fun, srange = srange[,x], xpar=xpar[x,],  alpha = alpha[,x], gamma = gamma[,x])
  }) 
}  
  
  if(!missing(ocor)) {                
            if(odistr == 'Gaussian') {
            opt <- mvrnorm(n=nspp, mu=rep(0,ndim),Sigma = ocor, empirical = FALSE)  # producing gaussian correlated variable (standard normal distribution)
            opt <- sapply(1:ndim, function(x) opt[,x]* xpar[x,2] + xpar[x,1])   #    adjusting mean and variance to desired values
            }
            
            if(odistr == 'uniform') {   
            opt <- mvrnorm(n=nspp, mu=rowMeans(xpar),Sigma = ocor, empirical = FALSE)  # producing gaussian correlated variable
            opt <- pnorm(scale(opt))   # inverse quantile function produces uniform correlated data   ~U[0,1]
            aa <- t(t(opt) * apply(xpar,1,diff))
            opt <- apply(aa,2, function(y) y + (colMeans(t(xpar)) - colMeans(aa)))
             }
             
     spp <- lapply(1:ndim, function(x) {
     spp[[x]]$m <- opt[,x]
     spp[[x]]
     }) 
  }     
  if(!missing(sdistr)){ 
    spp <- lapply(1:ndim, function(x) {
    spp[[x]]$m <- sdistr[,x]
    spp[[x]]
    })
  }   
  spp
}

"run.spp.present" <-function (spp){spp[,apply(spp,2,sum)>0]}
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
cor.mat.fun <- function(ndim, cors) {   
cor.mat <- diag(ndim)
if(!missing(cors)) {
for(i in 1:length(cors)) {
    cor.mat[cors[[i]][1],cors[[i]][2]] <-  cors[[i]][3]
    cor.mat[cors[[i]][2],cors[[i]][1]] <-  cors[[i]][3]
}
}          
cor.mat
}



make.env<-function(n,elen,emean, edistr, ecor, ndim){
  if(missing(ecor)) { 
  ecor <- diag(ndim)
  }
  if(edistr=="uniform"){
     env.norm <- mvrnorm(n=n, mu=rep(0,ndim),Sigma=ecor, empirical = FALSE) 
     env <- pnorm(scale(env.norm))
     env <- sapply(1:ndim, function(x) {
     aa <- env[,x] * elen[x]
     aa + (emean[x] - mean(aa))  
     })
     }
  if(edistr=="Gaussian"){
  env <- mvrnorm(n=n, mu=rep(0,ndim),Sigma=ecor, empirical = FALSE)
  env <- apply(env,2,function(x) (x-min(x)) / (max(x) - min(x))) 
  env <- sapply(1:ndim, function(x) {
  aa <- env[,x] * elen[x]
  aa + (emean[x] - mean(aa))  
  })  
  } 
  env<-as.data.frame(env)
  #names(env)<-c("env1", "env2", "env3")
  env
}


#-------------------------------------------------------------------------------
#make abundances

make.abundances <- function (env,param)   
{
  env<-as.matrix(env)
  #error checking
  stopifnot(ncol(env)==nrow(param))
  #if(!is.matrix(param)&!is.data.frame(param))
  param <- as.data.frame(param)
  with(param,{
  #calculate b and d
  d<-numeric(ncol(env))  
  b<-numeric(ncol(env))
  for(k in 1:ncol(env)){
  	b[k]=alpha[k]/(alpha[k]+gamma[k])
  	d[k]<-b[k]^alpha[k]*(1-b[k])^gamma[k]
  } 
  d<-prod(d)
  #find which env values are within taxon's range
  in.range<-matrix(numeric(0), nrow=nrow(env), ncol=ncol(env))
  for(k in 1:ncol(env)){
  	c1<-m[k]-r[k]*b[k]<=env[,k]
  	c2<-m[k]+r[k]*(1-b[k])>=env[,k]
  	in.range[,k]<-c1&c2
  }
  in.range<-apply(in.range,1, function(x) all(as.logical(x)))

  #calculate abundances
  A=numeric(nrow(env))
  A[!in.range]=0
  xx=env[in.range,,drop=F]
  if(nrow(xx)>0){
     h<-matrix(numeric(0),nrow(xx), ncol(xx))
     for(k in 1:ncol(xx)){
     	h[,k]=((xx[,k]-m[k])/r[k]+b[k])^alpha[k]*(1-((xx[,k]-m[k])/r[k]+b[k]))^gamma[k]
     }
     A[in.range]=Ao[1]/d*apply(h,1,prod)
  }
  A})
}

#-------------------------------------------------------------------------------
#testen von make abundances
'add.noise.bs' <- function (spp, cnt) 
{
  spp<-spp/rowSums(spp)   # percent data
  mat<-apply(spp,1,function(sam){
  	sam<- sample(1:ncol(spp),round(cnt), replace=TRUE, prob= sam)#,sample(1:ncol(spp),round(cnt/2), replace=TRUE))          # add more noise poisson 
  	table(c(sam,1:ncol(spp))) # include 1:n to ensure no taxa skipped by table as empty
  })
  t(mat-1)
}


abundances <- function(env, spp,nc) {
  ndim <- ncol(env)
  spp <- matrix(ncol = 5*ndim,unlist(spp))
  colnames(spp) <- rep(c('Ao','m','r','alpha','gamma'),ndim)
  abun <- sapply(1:nrow(spp), function(x)  {
      param <- matrix(ncol = 5, spp[x,], byrow = TRUE)
      colnames(param) <- c('Ao','m','r','alpha','gamma')
      make.abundances(env = env, param = param)
      
  })
  if(!missing(nc)) abun <- add.noise.bs(abun,nc)
  colnames(abun) <- 1:ncol(abun)
  abun  <- run.spp.present(abun)  
  abun <- abun/rowSums(abun)
  list(spp =as.data.frame(abun), env = env) 
}


#

make.set <- function(ndim, n, elen, emean, edistr, ecor, cnt, spec, env, ...) {
  if(missing(spec)) spec <- species (ndim,...)
  if(missing(env)) env <- make.env(n,elen,emean,edistr,ecor,ndim)
  dat <- abundances(env,spec,cnt)
  list(spp = dat$spp, env = dat$env, spec = spec)
}  


