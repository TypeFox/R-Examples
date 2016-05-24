
process.output <- function(x,DIC,params.omit,verbose=TRUE) {

if(verbose){cat('Calculating statistics.......','\n')}  
  
#Get parameter names
params <- colnames(x[[1]])
#Get number of chains
m <- length(x)

#Collapse mcmc.lists into matrix
mat = do.call(rbind,x)

#Get # of iterations / chain
n <- dim(mat)[1] / m

#Get parameter dimensions
dim <- get.dim(params)

#Create new parameter name vectors to handle non-scalar params
expand <- sapply(strsplit(params, "\\["), "[", 1)
params.simple <- unique(sapply(strsplit(params, "\\["), "[", 1))

#Functions for statistics
qs <- function(x,y){as.numeric(quantile(x,y))}
#Overlap 0 function
ov <- function(x){findInterval(0,sort(c(qs(x,0.025),qs(x,0.975))))==1}
#f function (proportion of posterior with same sign as mean)
gf <- function(x){if(mean(x)>=0){mean(x>=0)}else{mean(x<0)}}
#n.eff function
calcneff <- function(x,n,m){
  xp <- matrix(x,nrow=n,ncol=m)
  xdot <- apply(xp,2,mean)
  s2 <- apply(xp,2,var)
  W <- mean(s2)
  
  #Non-degenerate case
  if ((W > 1.e-8) && (m > 1)) {
    B <- n*var(xdot)
    sig2hat <- ((n-1)*W + B)/n      
    n.eff <- round(m*n*min(sig2hat/B,1),0)
    #Degenerate case
  } else {
    n.eff <- 1
  }
  n.eff
}

#Gelman diag function
gd <- function(i,hold){
  r <- try(gelman.diag(hold[,i], autoburnin=FALSE)$psrf[1], silent=TRUE)
  if(inherits(r, "try-error") || !is.finite(r)) {
    r <- NA
  }
  return(r)
}

#Make blank lists
sims.list <- means <- rhat <- n.eff <- se <- as.list(rep(NA,length(params.simple)))
q2.5 <- q25 <- q50 <- q75 <- q97.5 <- overlap0 <- f <- as.list(rep(NA,length(params.simple)))
names(sims.list) <- names(means) <- names(rhat) <- names(n.eff) <- params.simple
names(se) <- names(q2.5) <- names(q25) <- names(q50) <- names(q75) <- names(q97.5) <- params.simple
names(overlap0) <- names(f) <- params.simple

#This function modifies objects in global environment (output is discarded)
#Calculates statistics for each parameter
calc.stats <- function(i){
  
  #If parameter is not a scalar (e.g. vector/array)
  if(!is.na(dim[i][1])){
    
    #Get all samples
    sims.list[[i]] <<- mat[,expand==i]
    
    #If more than 1 chain, calculate rhat 
    #Done separately for each element of non-scalar parameter to avoid errors
    if(m > 1 && (!i%in%params.omit)){
      hold <- x[,expand==i]
      rhat.vals <- sapply(1:dim(hold[[1]])[2],gd,hold=hold)
      names(rhat.vals) <- colnames(hold[[1]])
      rhat[[i]] <<- populate(rhat.vals,dim[[i]])
    } else if (m == 1){
      hold <- x[,expand==i]
      rhat[[i]] <<- array(NA,dim=dim[[i]])
    }
    
    #Calculate other statistics
    ld <- length(dim(sims.list[[i]]))
    means[[i]] <<- populate(colMeans(sims.list[[i]]),dim[[i]])
    if(!i%in%params.omit){
    se[[i]] <<- populate(apply(sims.list[[i]],c(2:ld),sd),dim=dim[[i]])
    q2.5[[i]] <<- populate(apply(sims.list[[i]],c(2:ld),qs,0.025),dim=dim[[i]])
    q25[[i]] <<- populate(apply(sims.list[[i]],c(2:ld),qs,0.25),dim=dim[[i]])
    q50[[i]] <<- populate(apply(sims.list[[i]],c(2:ld),qs,0.5),dim=dim[[i]])
    q75[[i]] <<- populate(apply(sims.list[[i]],c(2:ld),qs,0.75),dim=dim[[i]])
    q97.5[[i]] <<- populate(apply(sims.list[[i]],c(2:ld),qs,0.975),dim=dim[[i]])
    overlap0[[i]] <<- populate(apply(sims.list[[i]],c(2:ld),ov),dim=dim[[i]])
    f[[i]] <<- populate(apply(sims.list[[i]],c(2:ld),gf),dim=dim[[i]])
    n.eff[[i]] <<- populate(apply(sims.list[[i]],c(2:ld),calcneff,n,m),dim=dim[[i]])   
    }
    
    sims.list[[i]] <<- populate(sims.list[[i]],dim=dim[[i]],simslist=T,samples=dim(mat)[1])
  
  #If parameter is a scalar
  } else {
    
    if(m > 1 && (!i%in%params.omit)){rhat[[i]] <<- gelman.diag(x[,i],autoburnin=FALSE)$psrf[1]}
    
    sims.list[[i]] <<- mat[,i]
    
    means[[i]] <<- mean(sims.list[[i]])
    if(!i%in%params.omit){
    se[[i]] <<- sd(sims.list[[i]])
    q2.5[[i]] <<- qs(sims.list[[i]],0.025)
    q25[[i]] <<- qs(sims.list[[i]],0.25)
    q50[[i]] <<- qs(sims.list[[i]],0.5)
    q75[[i]] <<- qs(sims.list[[i]],0.75)
    q97.5[[i]] <<- qs(sims.list[[i]],0.975)
    overlap0[[i]] <<- ov(sims.list[[i]])
    f[[i]] <<- gf(sims.list[[i]])
    n.eff[[i]] <<- calcneff(sims.list[[i]],n,m)}
  }
  
}

#Actually run function(nullout not used for anything)
nullout <- sapply(params.simple,calc.stats)

#Warn user if at least one Rhat value was NA
if(NA%in%unlist(rhat)&&verbose){
  options(warn=1)
  warning('At least one Rhat value could not be calculated.')
  options(warn=0,error=NULL)
}

#Do DIC/pD calculations if requested by user
if(DIC){
  dev <- matrix(data=mat[,'deviance'],ncol=m,nrow=n)   
  pd <- numeric(m)
  dic <- numeric(m)    
  for (i in 1:m){
    pd[i] <- var(dev[,i])/2
    dic[i] <- mean(dev[,i]) + pd[i]
  }    
  pd <- mean(pd)
  dic <- mean(dic)
  
  #Return this list if DIC/pD requested
  if(verbose){cat('\nDone.','\n')}
  return(list(sims.list=sims.list,mean=means,sd=se,q2.5=q2.5,q25=q25,q50=q50,q75=q75,q97.5=q97.5,overlap0=overlap0,
              f=f,Rhat=rhat,n.eff=n.eff,pD=pd,DIC=dic))
} else {
  #Otherwise return list without pD/DIC
  if(verbose){cat('\nDone.','\n')}
  return(list(sims.list=sims.list,mean=means,sd=se,q2.5=q2.5,q25=q25,q50=q50,q75=q75,q97.5=q97.5,overlap0=overlap0,
              f=f,Rhat=rhat,n.eff=n.eff))
}

}


