################################################################################
##### function to convert cost distances into dispersal probabilities
################################################################################
#'Function to calculate dispersal distances based on cost distances
#'
#'@param x Cost (Euclidean) distance matrix
#'@param d0 dispersal distance
#'@param p of all individuals in a population
#'@return a dispersal probability matrix
#'@description converts cost distances to probabilities: to reach a certain patch[ d0 average distance of p of all individuals, for example d0=100, p =0.5 -> 50\% procent of all migrating individuals go up to 100 m.

 p2p <- function(x, d0, p)
 {
 return (exp(((x/d0)*log(p))))
 }

################################################################################
######landscape creation functions
################################################################################


#'Function to add lines to landscape
#'
#'@param r raster that represents the landscape
#'@param x1 from coordinates
#'@param x2 to coordinates
#'@param val resistance value 1 is no resistance
#'@param plot switch if landscape should be plotted
#'@return returns modified raster layer
#'@description adds a line from x1 to x2 to a raster layer

addline <- function (r,x1,x2,val, plot=FALSE)
{
  ll <- Line(rbind(x1,x2))
  S1 = Lines(list(ll), ID="a")
  lin <- SpatialLines(list(S1))
  r1 <-rasterize(lin, r, field=val, fun='last', background=NA,  mask=FALSE, update=T, updateValue='all', filename="", na.rm=TRUE)
  if (plot==TRUE) plot(r1)
  r1
}

#'Function to add a polygon to a landscape
#'
#'@param r raster that represents the landscape
#'@param pol coordinates of the polygon
#'@param val resistance value, One equals no resistance
#'@param plot switch if landscape should be plotted
#'@return returns modified raster layer
#'@description adds a polygon to a raster layer

addpoly <- function(r, pol, val, plot=T)
{
  addpoly <- SpatialPolygons(list(Polygons(list(Polygon(rbind(pol, pol[1,]))), 1)))
  r1 <- rasterize(addpoly, r, field=val, update=TRUE, updateValue='all')
  if (plot==T) plot(r1)
  r1
}

################################################################################
##### functions  of popdynamics
################################################################################
#mutation
#'Function to execute mutation on a pop data.frame
#'
#'@param x a pop object
#'@param n.allels number of alleles 
#'@param mutrate the mutation rate
#'@param mtype the mutation process (currently on n.all is implemented. (=Choose any allel from 1:n.allels))
#'@param n.cov (number of covariates, default  )
#'@return a dispersal probability matrix
#'@description mutation subprocess on single populations

mutation <- function(x, n.allels, mutrate, mtype="n.all", n.cov)
{

n.indp <- nrow(x)
n.loci <-  ( ncol(x)-n.cov)  / 2
p.mut = mutrate
noffset = n.cov


if (!is.null(n.indp) )    #some populations might have gone extinct...
{
mut.ind <- which(matrix(as.logical( rbinom(n.indp*n.loci,1,p.mut)), ncol=n.loci, nrow=n.indp)==TRUE,arr.ind=T )

if (nrow(mut.ind)>0) #are there any mutations in this population?
{
for (i in 1:nrow(mut.ind))
  {
  
  if (mtype=="n.all")
    {
    mf <- rbinom(1,1,0.5)
    x[mut.ind[i,1], noffset+mut.ind[i,2]*2-mf] <- sample(1:n.allels,1)
    }
  }
}
}
x
}

####################################################################################
#reprodution
#'Function to execute reproduction on a pop data.frame
#'
#'@param x a pop object
#'@param type type of density dependence K.limit only yet 
#'@param K Kapacity of a subpoplation
#'@param n.off number of offspring per female
#'@param n.cov (number of covariates, default 3 )
#'@return an updated pop object
#'@description reproduction subprocess on single populations

reproduction <- function(x, type="K.limit",K=n.ind, n.off, n.cov)
{
noffset <- n.cov
n.loci <-  ( ncol(x)-n.cov)  / 2
n.ind <- nrow(x)

index.male <- x$sex=="male"
index.female <- x$sex=="female"
n.female <- sum(index.female)
n.male <-  sum(index.male)
#offsprings  <- data.frame(cbind(pop=NA, sex=NA, age=NA,matrix(NA,ncol=n.loci*2,nrow=1)) )


if (n.male>0 & n.female>0)
{
loci.cols <- (noffset+1):(n.loci*2+noffset)

genes.females <- (x[index.female,loci.cols])

#for (i in 1:(n.off-1)) genes.females <- rbind(genes.females,)

genes.females <- genes.females[rep(1:nrow(genes.females),n.off),]

pick.femall <-sample(c(T,F),n.female*n.off*n.loci,replace=T)
all.fem <-  matrix(as.vector(rbind(pick.femall,!pick.femall)),nrow=n.female*n.off, ncol=n.loci*2, byrow=T)


genes.off.fem <- matrix(t(genes.females)[t(all.fem)], ncol=n.loci, nrow=n.female*n.off, byrow=T)



breeding.males <- sample(which(index.male),n.female*n.off, replace=T)
genes.males <- x[breeding.males,loci.cols]
pick.malall <-sample(c(T,F),n.female*n.off*n.loci,replace=T)

all.mal <-  matrix(as.vector(rbind(pick.malall,!pick.malall)),nrow=n.female*n.off, ncol=n.loci*2, byrow=T)
genes.off.mal <- matrix(t(genes.males)[t(all.mal)], ncol=n.loci, nrow=n.female*n.off, byrow=T)


popn <- x$pop[1]
off.sex <- factor( ifelse(runif(n.female*n.off)<0.5,"male","female"))
off.age <- 1

all.seq<- NA

for (i in seq(1,(2*n.loci-1),2))
{ 

if (i==1) all.seq[i] <- 1 else all.seq[i] <- all.seq[i-2]+1
all.seq[i+1] <- all.seq[i]+n.loci
}

off.loci <- cbind(genes.off.mal, genes.off.fem)[,all.seq]

colnames(off.loci) <-paste("locus", rep(seq(1,n.loci,1),each=2),rep(c("A","B"),n.loci),sep="")

offsprings <- data.frame(pop=popn, sex=off.sex , age=off.age, off.loci)


#combine adults and offsprings
x <- rbind(x, offsprings)
#x<- (rbindlist(list(x,offsprings)))

#cut off at K

npop <- nrow(x)
x<- x[sample(1:npop,K, replace=F),]

} else x<- NULL    #return NULL if only males or only females....

return(x)
#return(data.frame(x))

}





#survival   (cut of at K)
#surv <- function(x, s=t.surv) { 
#
#x<- x[runif(dim(x)[1])< s,]
#}
#lapply(pops, surv, s=0.1)

################################################################################
#emigration (for all populations at once xp=pops)
#'Function to execute emigration on a pops object
#'
#'@param xp all pops combined in a list
#'@param perc.mig percentage if migrating individuals 
#'@param emi.m emigration probability (normally based on cost dispersal distance)
#'@param emi.table a fixed number of migrating individuals can be specified (overrides emi.m)
#'@return a list, first entry are updated pops, second entry the number of disperserin a matrix
#'@description emigration process on all population in one go

emigration <- function(xp, perc.mig, emi.m, emi.table=NULL)
{
n.pops = length(xp)
migs <- matrix(0,nrow=n.pops, ncol=n.pops)
migrants<- NA

if (is.null(emi.table))
{
for (i in 1:n.pops)
{
pop.size <-nrow(xp[[i]])

if (!is.null(pop.size))
{
#how many migrants
migrants[i] <- sum(ifelse(runif(pop.size)<perc.mig,1,0))
#where do they go to (depends on emis)
fromto <- table(sample(1:n.pops,migrants[i],replace=T, prob=emi.m[,i]))

to <- as.numeric(names(fromto))
for (ii in 1:length(to)) migs[i,to[ii ]] <- migs[i,to[ii ]] + fromto[ii]
#migs[i, unique(to)] <- migs[i,unique(to)] + table(to)

}
}
} else migs <- emi.table #end emi.type


#deterministic emigration (needs to be passed by emi.table


#let them move (at the moment an individual can migrate more than once, need to be fixed)

for (from in 1:n.pops)
{

for (to in 1:n.pops)
  {
  #update population size during migration events!!!
  psize <- nrow(xp[[from]])
  m <- migs[from, to]
  if (m>0)
    {
    ind.from <- sample(1:psize,m,replace=F)
    xp[[to]] <- rbind(xp[[to]], xp[[from]][ind.from,])
    xp[[from]] <- xp[[from]][-ind.from,]
    xp[[to]][,"pop"] <- to
    }
  }
}
#output number of migrants to a population
#diag(migs) <- 0
#cat(paste("Migrants: ",colSums(migs),"\n"))
results <- list(xp, migs)
return(results)
}

################################################################################
#combines all pops into a single genind object.
#'Function converts pops to a genind object
#'
#'@param x pops object (a list of pop)
#'@param locs named coordinates of locations
#'@param n.cov number of covariates (defaults to 3)
#'@return a spatial genind object
#'@description converts pops into genind (to calculate Fst etc.)

pops2genind <- function(x, locs=NULL, n.cov=3)
{         
n.pops <- length(x)
n.loci <- (ncol(x[[1]])-n.cov )/2 #only diploid currently
noffset <- n.cov

##### convert pops to data.frame
#combine <-do.call(rbind.data.frame,x)
#combine <- do.call(rbind.data.frame,x)
combine <- as.data.frame(rbindlist(x))

#convert to genind
allele <-combine[,-(1:noffset)]
pop.size <- table(combine[,"pop"] )
res <- data.frame(matrix(NA, nrow=nrow(allele), ncol=n.loci))#rep(NA,sum(pop.size)))
for (i in seq(1, 2*n.loci, 2))
            res[, ceiling(i/2)] <- paste(allele[, i], allele[, i + 1], sep = "/")
pops.genind <- df2genind(res,pop=combine$pop, sep="/", ind.names= rownames(res))
#add locations
if (!is.null(locs))
  {
  x <- rep(locs[,1], pop.size)
  y <- rep(locs[,2], pop.size)
  pops.genind@other$xy <- cbind(x=x,y=y)
  popNames(pops.genind) <- row.names(locs)
  }
return(pops.genind)
}



################################################################################
### sim functions
################################################################################
#'Run forward step generations a popgen simulation
#'
#'@param simpops pops object (a list of pop)
#'@param steps the number of steps (generations)
#'@param cost.mat a cost matrix (e.g. calculated via costDistance)
#'@param n.offspring number of offsprings per female
#'@param n.ind number of individuals
#'@param mig.rate migration rate
#'@param disp.max maximal dispersal distance of disp.rate individuals
#'@param disp.rate percentage of individuals achieving disp.max
#'@param n.allels number of maximal alleles at a loci
#'@param mut.rate mutation rate
#'@param n.cov number of covariates (defaults to 3)
#'@param rec switch if emigration matrix should be recorded, either "none" or "emi"
#'@param emi.table a emigration matrix, if provide a fixed number of migration events will take place otherwise based on disp.max, mig.rate and disp.rate, events occuring are (in that order, emigration, reproduction, mutation)
#'@return an updated simpops object after steps steps
#'@description performs a time-forward, agent-based and spatiallly explicit genetic population simulation 
run.popgensim <- function(simpops, steps, cost.mat, n.offspring , n.ind,  mig.rate, disp.max, disp.rate, n.allels, mut.rate, n.cov=3, rec = "none", emi.table=NULL)
{
pops <- simpops
distances <- cost.mat
disdis <- p2p(distances,d0=disp.max, p=disp.rate)
stdemi.table= disdis/rep(colSums(disdis), each=ncol(disdis))
emi.tableint=NULL
#cl<- makePSOCKcluster(6)

for (year in 1:steps)
{

if (!is.null(emi.table))
  {
  emi.table <- stdemi.table * mig.rate * n.ind +  emi.table
  emi.tableint <- floor(emi.table)
  emi.table <- emi.table - emi.tableint
  }


#emigration
emi<- emigration(pops, mig.rate, emi.m=disdis, emi.table=emi.tableint)
pops <- emi[[1]]
#emi[[2]] table of migrants....

#reproduction
 pops <- lapply(pops, reproduction,K=n.ind, n.off=n.offspring, n.cov=n.cov)
 # pops <- parLapply(cl, pops, reproduction, K=n.ind, n.off=n.offspring, n.cov=n.cov)
#mutation
pops <- lapply(pops, mutation, n.allels, mutrate=mut.rate, n.cov=n.cov)
}  #end of every year....   time loop

#stopCluster(cl)

out <- pops #list(pops=pops)
if (rec=="emi")
{
out <- list()
out[[1]] <- pops
out[[2]] <- emi[[2]]
}
return (out)
} #end of simpop function

################################################################################
#initialise populations
################################################################################
#'Initialise a pops object fo a  a popgen simulation
#'
#'@param n.pops number of subpopulations
#'@param n.ind number of individuals per subpopulation
#'@param sex.ratio sex ratio of males and females
#'@param n.loci number of loci
#'@param n.allels number of maximal alleles per loci
#'@param locs coordinates of the subpopulations
#'@param n.cov number of covariates (defaults to 3)
#'@return a simpops object (to be used in run.popgensim )
#'@description initialises a time-forward, agent-based and spatiallly explicit genetic population simulation 
init.popgensim <- function(n.pops, n.ind, sex.ratio, n.loci, n.allels, locs, n.cov)
{
noffset <- n.cov #number of covariates before the loci
#initialise populations
empty.loci <- matrix(NA,ncol=n.loci*2,nrow=n.ind)
colnames(empty.loci) <-paste("locus", rep(seq(1,n.loci,1),each=2),rep(c("A","B"),n.loci),sep="")
pops <- list( )
for (i in 1:n.pops)
{
pops[[i]]<- data.frame(cbind(pop=NA, sex=NA, age=NA,empty.loci))
#population number
pops[[i]][,"pop"] <- i
#sex
pops[[i]][,"sex"] <- "male"
pops[[i]][1:(round(n.ind*sex.ratio,0)),"sex"] <- "female"

pops[[i]][,"sex"] <- factor(pops[[i]][,"sex"])

#uniform allel distribution
pops[[i]][,(noffset+1):(n.loci*2+noffset) ]  <- sample(1:n.allels,n.loci*2*n.ind, replace=T)
}

names(pops) <- row.names(locs)   #or from genind object....
return (pops)
}



################################################################################
### faster pairwise fst
################################################################################
#'Calculates pairwise fsts using a genind object (very fast)
#'
#'@param gsp a genind object
#'@return a pairwise fst matrix (same as hierfstat pairwise.fst)
#'@description for simulation the original pairwise.fst was much to slow. The fast version works only without NAs and diploid data (so do not use it on non-simulated data)

pairwise.fstb <- function(gsp)
{
n.pops <- length(popNames(gsp))
allPairs <- combn(1:n.pops, 2)
gen.mat2<- matrix(0, nrow=n.pops, ncol=n.pops)
allfreq <- function (x) 1-sum((x/2)^2)

for (i in 1:ncol(allPairs)){

np1 <- allPairs[1,i]
np2 <- allPairs[2,i]
p1 <- gsp@tab[gsp@pop==popNames(gsp)[np1],]
p2 <- gsp@tab[gsp@pop==popNames(gsp)[np2],]
p12 <- rbind(p1,p2)

lf <- gsp@loc.fac

n1 <- nrow(p1)
n2 <- nrow(p2)
n12 <- nrow(p12)

Hs1 <- mean( tapply(colSums(p1)/n1, lf, allfreq))
Hs2 <- mean( tapply(colSums(p2)/n2, lf, allfreq))
Hs12 <- mean( tapply(colSums(p12)/n12, lf, allfreq))
fst <- (Hs12-weighted.mean(c(Hs1, Hs2), c(n1,n2))) / Hs12
gen.mat2[np1,np2] <- fst
gen.mat2[np2,np1] <- fst

}
la <- popNames(gsp)
colnames(gen.mat2) <- rownames(gen.mat2) <- la
return(gen.mat2)
}
