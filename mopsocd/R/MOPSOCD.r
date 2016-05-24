mopsocd <- function(fn,gn,varcnt,fncnt,lowerbound,upperbound,opt,popsize,maxgen,archivesize,verbosity,pMut,w,c1,c2) UseMethod("mopsocd")
mopsocd <- function(fn,gn=NULL,varcnt,fncnt,lowerbound,upperbound,opt,popsize=100,maxgen=100,archivesize=250,verbosity=0,pMut=0.5,w=0.4,c1=1.0,c2=1.0)
{

  varcnt <- as.integer(varcnt)
  fncnt <- as.integer(fncnt)
  lowerbound <- as.integer(lowerbound)
  upperbound <- as.integer(upperbound)
  opt <- as.integer(opt)
  popsize <- as.integer(popsize)
  maxgen <- as.integer(maxgen)
  archivesize <- as.integer(archivesize)
  verbosity <- as.integer(verbosity)
  pMut <- as.numeric(pMut)
  w <- as.numeric(w)
  c1 <- as.numeric(c1)
  c2 <- as.numeric(c2)
  # pareto <- mopsocdbase(fn,gn,varcnt,fncnt,lowerbound,upperbound,opt,popsize,maxgen,archivesize,verbosity,pMut,w,c1,c2)  

  # pareto$call <- match.call()
  
  # class(pareto) <- "pareto"
  # pareto
# }

# mopsocdbase <- function(fn,gn,varcnt,fncnt,lowerbound,upperbound,opt,popsize,maxgen,archivesize,verbosity,pMut,w,c1,c2)
# {
####################
# Global Variables #
####################
minvalue <- vector(length=varcnt)
maxvalue <- vector(length=varcnt)
archiveVar <- matrix(,nrow=archivesize,ncol=varcnt) 
archiveFit <- matrix(,nrow=archivesize,ncol=fncnt)
popVar <- matrix(,nrow=popsize,ncol=varcnt)
popFit <- matrix(,nrow=popsize,ncol=fncnt)
pbestsVar <- matrix(,nrow=popsize, ncol=varcnt)
pbestsFit <- matrix(,nrow=popsize, ncol=fncnt)
velocity <- matrix(,nrow=popsize, ncol=varcnt) 
crowdDist <- vector(length=archivesize)
minvaluetemp <- matrix(,nrow=popsize,ncol=varcnt)
maxvaluetemp <- matrix(,nrow=popsize,ncol=varcnt)
ndomCtr <- 0

maxvalue <- upperbound
minvalue <- lowerbound

#########################
# Initialize Population #
#########################

initialize_pop <- function() {

  for (i in 1:varcnt) {
    popVar[1:popsize,i] <<- runif(popsize,min=minvalue[i],max=maxvalue[i])
  }

}

##################################
# Initialize Particle Velocities #
##################################

initialize_vel <- function() { 
  
  velocity[1:popsize,1:varcnt] <<- 0.0  
  
}

##########################################################
# Evaluate: Compute Value of each Particle in Population #
##########################################################

evaluate <- function() {

      popFit[1:popsize,] <<- t(apply(popVar,1,fn))

}

######################################################
# Store Values and Fitness Values for Best Particles #
######################################################

store_pbests <- function() {  
  pbestsVar[1:popsize,] <<- popVar[1:popsize,]
  pbestsFit[1:popsize,] <<- popFit[1:popsize,]
}

##########################################
# Insert ndominated Particles in Archive #
##########################################

insert_nondom <- function() {

  if(!is.null(gn)) {
    z <- rowSums(t(matrix(!apply(popVar,1,gn),ncol=popsize)*1)) 
    minz <- min(z)
    ix <- which(z == minz)
    lenix <- length(ix)
  
    ixfpop <- NULL
    for (i in ix) {   
      if (opt == 1) {    
        x <- (t(matrix(rep(popFit[i,],lenix),nrow=fncnt))-popFit[ix,] < 0.0)*1
      } else {
        x <- (t(matrix(rep(popFit[i,],lenix),nrow=fncnt))-popFit[ix,] > 0.0)*1
      }         
      if((max(rowSums(x)) == fncnt) == FALSE){
        ixfpop <- c(ixfpop,i)
      }
    }  
  } else {
    ixfpop =seq(popsize)
  }
  lenfpop <- length(ixfpop)
  
  if(lenfpop != 0){
    if (ndomCtr + lenfpop <= archivesize) {
      archiveVar[(ndomCtr+1):(ndomCtr+lenfpop),] <<- popVar[ixfpop,]
      archiveFit[(ndomCtr+1):(ndomCtr+lenfpop),] <<- popFit[ixfpop,]
      ndomCtr <<- ndomCtr + lenfpop
    } else {
      cat("FATAL ERROR:\nArchive Size Value ",archivesize," too small\n")
      cat("Recommended Archive Size: 2.5*population size\n")
      cat("Press Enter to exit\n")
      scan(quiet=TRUE)
      cat("Goodbye !\n")
      stop()    
    }
  } else {
    cat("No solution (Unreasonable constraints)\n")
    cat("Press Enter to exit\n")
    scan(quiet=TRUE)
    cat("Goodbye !\n")
    stop()
  }
  
}

#############################
# Compute Crowding Distance #
#############################

crowding <- function() {
    
  crowdDist[1:ndomCtr] <<- 0.0
  
  for (i in seq(fncnt)) {
    
    sorted_archive <- sort(archiveFit[1:ndomCtr,i],index.return=TRUE)
    valsa <- sorted_archive$x
    ixsa <- sorted_archive$ix  
    crowdDist[ixsa[2:(ndomCtr-1)]] <<- crowdDist[ixsa[2:(ndomCtr-1)]] + valsa[3:ndomCtr] - valsa[1:(ndomCtr-2)]
    maxcd <- which.max(crowdDist[1:ndomCtr])
    crowdDist[ixsa[1]] <<- crowdDist[ixsa[1]] + crowdDist[maxcd]
    crowdDist[ixsa[ndomCtr]] <<- crowdDist[ixsa[ndomCtr]] + crowdDist[maxcd]   
    mincd <- which.min(crowdDist[1:ndomCtr])
    crowdDist[1:ndomCtr] <<- crowdDist[1:ndomCtr]/(crowdDist[maxcd]-crowdDist[mincd])
  }  
  
  sorted_cd <- sort(crowdDist[1:ndomCtr],decreasing=TRUE,index.return=TRUE)
  archiveVar[1:length(sorted_cd$ix),] <<- archiveVar[sorted_cd$ix,]
  archiveFit[1:length(sorted_cd$ix),] <<- archiveFit[sorted_cd$ix,]
 
}


######################################
# Update Personal Bests of Particles #
######################################

update_pbests <- function() {
  
  if(opt == 1){
    x <- rowSums((popFit-pbestsFit >  0.0)*1)
  } else {
    x <- rowSums((popFit-pbestsFit <  0.0)*1)
  }
  
  ix1 <- which(x != 0 & x != fncnt)
  ix2 <- sample(ix1,floor(length(ix1)/2),replace=FALSE)
  ix <- union(which(x == fncnt),ix2)

  pbestsVar[ix,] <<- popVar[ix,]  
  pbestsFit[ix,] <<- popFit[ix,]
  
}

############################
# Update Archive Particles #
############################

update_archive <- function() {
  if(!is.null(gn)) {
    z <- rowSums(t(matrix(!apply(popVar,1,gn),ncol=popsize)*1)) 
    minz <- min(z)
    ix <- which(z == minz)
    lenix <- length(ix)

    ixfpop <- NULL
    for (i in ix) {   
      if (opt == 1) {    
       x <- (t(matrix(rep(popFit[i,],lenix),nrow=fncnt))-popFit[ix,] < 0.0)*1
      } else {
        x <- (t(matrix(rep(popFit[i,],lenix),nrow=fncnt))-popFit[ix,] > 0.0)*1
      }         
      if((max(rowSums(x)) == fncnt) == FALSE){
       ixfpop <- c(ixfpop,i)
      }
    } 
  } else {
    ixfpop =seq(popsize)
  }
  lenfpop <- length(ixfpop)
  
  ixinspop <- NULL
  for (i in ixfpop) { 
    
    if (opt == 1) {    
      x <- (t(matrix(rep(popFit[i,],ndomCtr),nrow=fncnt))-archiveFit[1:ndomCtr,] < 0.0)*1
    } else {
      x <- (t(matrix(rep(popFit[i,],ndomCtr),nrow=fncnt))-archiveFit[1:ndomCtr,] > 0.0)*1
    }         
    if((max(rowSums(x)) == fncnt) == FALSE) {
      ixinspop <- c(ixinspop,i)
    } 
  }  
  leninspop <- length(ixinspop)
    
  if(leninspop != 0){  
    if (ndomCtr + leninspop <= archivesize) {
      archiveVar[(ndomCtr+1):(ndomCtr+leninspop),] <<- popVar[ixinspop,]
      archiveFit[(ndomCtr+1):(ndomCtr+leninspop),] <<- popFit[ixinspop,]
      ndomCtr <<- ndomCtr + leninspop
    } else {
      crowding()
      bottom <- floor(ndomCtr*0.9)      
      cnt <- min((archivesize-bottom),leninspop) 
      ixinssel <- sample(ixinspop,cnt,replace=F)
      ixreplace <- sample((archivesize-cnt):archivesize,cnt,replace=F)
      archiveVar[(archivesize-cnt+1):archivesize,] <<- popVar[ixinssel,]
      archiveFit[(archivesize-cnt+1):archivesize,] <<- popFit[ixinssel,]
      ndomCtr <<- archivesize  
    }      
  } else {
     return(0) 
  }
  
  ixarc <- seq(ndomCtr)
  ixndomarc <- NULL
  rmove <- NULL
  dup <- NULL
  for (i in ixarc) {   
    if (i %in% dup) {
      if(verbosity >= 3) cat("Duplicate Particle Found ",i,"\n")
      next
    }  
    if (opt == 1) {    
      x <- (t(matrix(rep(archiveFit[i,],ndomCtr),nrow=fncnt))-archiveFit[ixarc,] < 0.0)*1
    } else {
      x <- (t(matrix(rep(archiveFit[i,],ndomCtr),nrow=fncnt))-archiveFit[ixarc,] > 0.0)*1
    }   
    rmove <- c(setdiff(which(rowSums(x) == 0),i),rmove)
    wich1 <- setdiff(which(rowSums(x)==0),i) 
    if(length(wich1) != 0) {
     if(!(i %in% wich1)){
        for (j in wich1) {
          if(all(archiveFit[j,] == archiveFit[i,])){
            dup <- c(j,dup)
          } 
        }
     }
    }

  }  
  
  rmove <- union(rmove,dup)
  lenrem <- length(rmove)
  if(lenrem != 0){
    cor = setdiff(seq(ndomCtr),rmove)
    ndomCtr <<- length(cor)  
    archiveVar[1:ndomCtr,] <<- archiveVar[cor,]
    archiveFit[1:ndomCtr,] <<- archiveFit[cor,]
  }

}


#############################################
# Compute Velocity of Particles in Population
#############################################

compute_velocity <- function(){
  
	top <- max(floor((ndomCtr*0.1)),1)
  gBest <- sample(1:top,1,replace=TRUE)
  gBestmat <- t(matrix(rep(archiveVar[gBest,],popsize),nrow=varcnt))
	randum1 <- matrix(runif(popsize*varcnt),popsize,varcnt) 
	randum2 <- matrix(runif(popsize*varcnt),popsize,varcnt)
	velocity[1:popsize,] <<- w*velocity[1:popsize,] + c1*randum1*(pbestsVar-popVar) + c2*randum2*(gBestmat-popVar) 
	popVar[1:popsize,] <<- popVar[1:popsize,] + velocity[1:popsize,]
  
}

##############################################
# Maintain Particles within Feasibility Region
##############################################

maintain_particles <- function(){

	for (j in seq(varcnt)) {    
    wich1 <- which(popVar[,j] < minvalue[j])
    wich2 <- which(popVar[,j] > maxvalue[j])
    popVar[wich1,j] <<- minvalue[j]
    popVar[wich2,j] <<- maxvalue[j]
    wich12 <- union(wich1,wich2)
    velocity[wich12,j] <<- -1.0*velocity[wich12,j]
  }
  
}

######################################
# Mutate Function (adapted from MOPSO)
######################################

mutate <- function() { 
  
  if(t >= maxgen*pMut)
    return(0)

  rangev <- vector(length=varcnt)
  flip <- ifelse(runif(1) < (1-(t/(maxgen*pMut)))^1.5,1,0) 
  if(flip) {
    dim <- sample(1:varcnt,1,replace=TRUE)
    rangev[1:varcnt] <- (maxvalue[1:varcnt] - minvalue[1:varcnt])*0.5*(1-(t/(maxgen*pMut)))^1.5
    wich1 <- which((popVar[,dim] - rangev[dim]) < minvalue[dim])
    wich2 <- which((popVar[,dim] + rangev[dim]) > maxvalue[dim])
    if(length(wich1) != 0) {
      minvaluetemp[wich1,dim] <<- minvalue[dim]
      minvaluetemp[-wich1,dim] <<- popVar[-wich1,dim] - rangev[dim]
    } else {
      minvaluetemp[,dim] <<- popVar[,dim] - rangev[dim]      
    }
    if(length(wich2) != 0) {
      maxvaluetemp[wich2,dim] <<- maxvalue[dim]
      maxvaluetemp[-wich2,dim] <<- popVar[-wich2,dim] + rangev[dim]
    } else {
      maxvaluetemp[,dim] <<- popVar[,dim] + rangev[dim]
    } 
    popVar[,dim] <<- runif(popsize,min=minvaluetemp[,dim],max=maxvaluetemp[,dim]) 
  } 
  
}

##############
# Main Program
##############
if(verbosity >= 1) cat("Welcome to MOPSO-CD\n")

t1 <- Sys.time()
t <- 0 

initialize_pop()
initialize_vel()
evaluate()
store_pbests()
insert_nondom()

while (t < maxgen) {
  t <- t + 1
  
  if(verbosity >= 1)
  cat("Generation ",t,"\n") 

  if (ndomCtr>2*fncnt)
    crowding()
  
  compute_velocity()
  maintain_particles()
  mutate()
  evaluate()
  update_archive()
  update_pbests()

  if(verbosity >= 3){
    cat("\nVariables\n")
    print(archiveVar[1:ndomCtr,])
    cat("Objective Function Values:\n")
    print(archiveFit[1:ndomCtr,])  
    cat("\n")
  }
  
}

if(verbosity >= 2){
  cat("###############################\n")
  cat("         FINAL VALUES\n")
  cat("###############################\n")
  cat("Variables:\n")
  print(archiveVar[1:ndomCtr,])
  cat("\nObjective Function Values:\n")
  print(archiveFit[1:ndomCtr,])  
}

if(verbosity >= 1){
  cat("\nNondominated Solutions: ",ndomCtr,"\n")
  cat("Number of Generations: ",t,"\n")
  cat("Computation Time: ",Sys.time() - t1,"\n")
}

if(!is.null(gn)) {
  pareto <- list("paramvalues"=archiveVar[1:ndomCtr,],"objfnvalues"=archiveFit[1:ndomCtr,],"numsols"=ndomCtr,
     "fn"=fn,"gn"=gn,"varcnt"=varcnt,"fncnt"=fncnt,"lowerbound"=lowerbound,"upperbound"=upperbound,"opt"=opt,"popsize"=popsize,"maxgen"=maxgen,"archivesize"=archivesize,"pMut"=pMut,"w"=w,"c1"=c1,"c2"=c2)
} else {
  pareto <- list("paramvalues"=archiveVar[1:ndomCtr,],"objfnvalues"=archiveFit[1:ndomCtr,],"numsols"=ndomCtr,
       "fn"=fn,"gn"=NULL,"varcnt"=varcnt,"fncnt"=fncnt,"lowerbound"=lowerbound,"upperbound"=upperbound,"opt"=opt,"popsize"=popsize,"maxgen"=maxgen,"archivesize"=archivesize,"pMut"=pMut,"w"=w,"c1"=c1,"c2"=c2)  
}

  warnings()
  
  pareto$call <- match.call()
  class(pareto) <- "pareto"
  pareto
 
}

