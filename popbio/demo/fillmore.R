## fillmore.R

## copied from demo graphics
#if(dev.cur() <= 1) get(getOption("device"))()
#opar <- par(ask = interactive() &&
#            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))


caption<-function (x) 
{
    cat("Press <Enter> to continue...")
    readline()
    invisible()
}


options(digits=4)

# ---------------------------------------------------------------- #
# Data sets

data(aq.trans)
head2(aq.trans)
years<-unique(aq.trans$year)

# ---------------------------------------------------------------- #
# Fillmore stage vectors 

sv<-table(aq.trans$stage, aq.trans$year)

## totals
addmargins(sv)

## mean stage vector
round(apply(sv, 1, mean),0)


stage.vector.plot(sv[-1,], prop=FALSE, col=rainbow(4),
     ylab = "Total number of plants",
     main = "Fillmore Canyon stage vectors")

### plot proportions

op<-par(mar=c(5,4,2,7), xpd=TRUE)
x<-barplot(prop.table(sv[-1,],2), las=1,
xlab="Year", ylab="Proportion in stage class", main="Fillmore Canyon",
col=gray(1:4/4), ylim=c(0,1), xaxt='n', space=.5)
box()
yrs<-substr(colnames(sv),3,4)
axis(1,x, yrs)
## draw outside plot boundaries

legend(13,0.7, rev(rownames(sv)[-1]), fill=rev(gray(1:5/5)))
par(op)


# ---------------------------------------------------------------- #
# DENISTY plots at Fillmore Canyon (using 1 meter square plots)

a<-ftable( aq.trans$plot, aq.trans$stage, aq.trans$year)
dim(a)
## PLOT four stage classes: recruits, small and large veg, flowering plants
op<-par(mfrow=c(2,2), mar=c(4, 4, 2 , 2) , oma=c(1,1,2,0))
fstages<-c("Recruits", "Small vegetative", "Large vegetative", "Flowering plants")
for (i in 2:5)
{
   x<-a[seq(i,50,5),]   
   plot(years, x[1,], type='n',  ylim=c(0,max(x)), ylab=fstages[i-1], xlab="Year")
   for (j in 1:10)
   {
       lines(years, x[j,], col=rainbow(10)[j])
   }
}
title(expression(paste("Plant density in 10 demography plots (", m^-2, ")")), outer=TRUE)
par(op)


# ---------------------------------------------------------------- #
# BUILD a projection matrix for Aquilegia

## 1.  get transitions for a single projection interval (1996-1997)
x<-subset(aq.trans, year==1996)

## 2. count number of recruits in 1997 
recruits<- sv["recruit", "1997"]

## 3.  Use aq.matrix() to build matrix 
aq.96<- aq.matrix(x, recruits)
aq.96
caption("1996-1997 projection matrix and vectors from aq.matrix")


# ---------------------------------------------------------------- #
#  Population projections

p<-pop.projection(aq.96$A, aq.96$n)
p
## exclude seeds
stage.vector.plot(  p$stage.vectors[-1,1:10], col=rainbow(4), prop=FALSE,
   xlab="Years after 1996", main="Stage vector projections", log='y')


# ---------------------------------------------------------------- #
#  Eigenvalues and vectors

eig.96<- eigen.analysis(aq.96$A, zero=FALSE)
eig.96


## some Plot methods
barplot(eig.96$stable.stage, col="blue", ylim=c(0,1), ylab="Proportion", xlab="Stage class", main="Stable stage distribution from 96-97 matrix")
box()

ymax<-max(eig.96$repro.value)*1.25
barplot(eig.96$repro.value, col="green", ylim=c(0,  ymax ), xpd=FALSE, ylab="Reproductive value", xlab="Stage Class", main="Reproductive values")
box()


matplot2(eig.96$elasticities, type='b', ylab="Elasticity", legend="topleft", ltitle="Fate", main="Elasticity matrix")


image2(eig.96$sensitivities, mar=c(1,3,5,1))
title("Sensitivity matrix")


# ---------------------------------------------------------------- #
#  Bootstrap confidence intervals

# get transitions with individual fertilities
y<-aq.matrix(x, recruits, summary=FALSE)

# 200 bootstrap samples for speed (5000 is default) -
#  estimated seed bank survival and recruitment does not vary here, but could

aq.96.boot<-boot.transitions(y, 200, add=c(1,1, aq.96$seed.survival, 2,1, aq.96$recruitment.rate) )

# calculate confidence intervals using quantile() or perc.ci and other methods in boot package

ci<- quantile(aq.96.boot$lambda, c(0.025,0.975) )
aq.96$lambda
ci

# plot histogram  ... plotCI in gplots and other packages
hist(aq.96.boot$lambda, col="green", xlab="Lambda", main=paste('Bootstrap estimates of population\ngrowth rate from 1996-1997'))
abline(v=ci, lty=3)


# ---------------------------------------------------------------- #
# CALCULATE population growth rates using constant 10,000 seed bank and 12.6% seed survival


years<-1996:2002
n<-length(years)


fill1.all<- vector("list", n)
names(fill1.all)<-years

fill1<-matrix(numeric(length(years)*5), nrow=n, dimnames=list(years, c("lambda", "seed.new", "seed.b0", "seed.b1", "rec.rate") )) 
for (s.year in years)
{
   i<-s.year- (years[1]-1)   ## index for list starting at 1
   x<-aq.matrix(subset(aq.trans, year==s.year), sv['recruit', as.character(s.year+1)] )
   fill1[i,] <- c(x$lambda, x$seeds.from.plants, x$seed.bank, x$n1[1], x$recruitment.rate)
   fill1.all[[i]]<-x
 }
fill1.all[["1997"]]
fill1

plot(rownames(fill1), fill1[,1], type='b', ylim=c(0,1.4), xlab="Year", ylab="Growth rate" , main=paste("Population growth rate using constant \nseed bank size and seed survival rate"), col="red", pch=16, las=1)
 abline(h=1, lty=3)



# ---------------------------------------------------------------- #
#  Projection matrix elements....


## list only projection matrices - will be used later for stochastic growth
fill.A<-lapply(fill1.all, FUN="[[", "A")

## add column names using matrix element notation aij where i=row and j=column))
fill.vr<-matrix(unlist(fill.A), nrow= n , byrow=TRUE,
           dimnames=list(years, paste("a", 1:5 , rep(1:5,each=5), sep="")) )

## mean matrix and variation
fill.vr<-rbind(fill.vr, mean=apply(fill.vr, 2, mean), var=apply(fill.vr, 2, var) )

# TABLE using years as columns and ignoring empty elements
t( fill.vr[,c(1:2, 8,10,13:15,18:25)] )

caption("Annual, mean and variance of projection matrix elements (1996-2002)")


# ---------------------------------------------------------------- #
# Mean matrix
A.mean<-matrix(fill.vr["mean",],  nrow=5, dimnames=dimnames(fill.A[[1]]))

A.mean

eigen.analysis(A.mean)
caption("Eigen analysis of mean matrix")


# ---------------------------------------------------------------- #
#  AGE specific rates

# split A into T and F (could change aq.matrix output, but fertility always in a15 and a25)

splitA(A.mean, r=1:2)


generation.time(A.mean, r=1:2)
net.reproductive.rate(A.mean, r=1:2)
fundamental.matrix(A.mean, r=1:2)$N

caption("Age-specific traits from mean matrix")

## Check generation time for all matrices.  If F matrix is empty, generation time is Inf
sapply(fill.A, generation.time, r=1:2)

caption("Generation time for each projection matrix")




# ---------------------------------------------------------------- #
# VARY seed survival rates before matrix construction
#  --  seed survival part of three matrix elements (1 seed survival,2 fertilities in prebreeding census)


fill.seed.surv <- vector("list", n) 
names(fill.seed.surv)<- years

for (s.year in years)
{
   z<-s.year - (years[1]-1)
   fill.seed.surv[[z]]<-matrix(numeric(11*3), nrow=11,
            dimnames=list(1:11, c("survival", "lambda", "seeds")) ) 
   for (i in seq(0,100,10))
   {
      x<-aq.matrix(subset(aq.trans, year==s.year),
                   sv['recruit',as.character(s.year+1)], seed.survival=i/100)
      fill.seed.surv[[z]][(i+10)/10,]<- c(i/100, x$lambda, x$n1[1])
   }
}

fill.seed.surv[["1996"]]

plot(fill.seed.surv[["1996"]][,1], fill.seed.surv[["1996"]][,2], type="n",
     xlim=c(-0.08,1), ylim=c(0,1.5), xlab="Seed survival rate", ylab="Growth rate " ,
     main=paste("Changes in seed survival rate\n on population growth rate"))

abline(h=1, lty=2)
for (i in 1:n)
{
  lines(fill.seed.surv[[i]][,1], fill.seed.surv[[i]][,2], lwd=2, col=rainbow(11)[i], lty=i)
}

#LABELS?
y.abbr<-substr(as.character(years),3,4)

y<-numeric(n)
for (i in 1:n){y[i]<-fill.seed.surv[[i]][1,2]}

text(-.08, y, y.abbr, pos=4)


##-------------------------------------------------------------##
#  VARY seed bank size from 1000 to 20000  in 1997 only (for  demo) 

##  no recruits survive in 99, 01,02,03 so seed bank size does not alter the growth rate in these years


   fill.seed.97<-matrix(numeric(25*2), nrow=25 ) 
   for (i in seq(1000,25000,1000))
   {
      x<-aq.matrix(subset(aq.trans,  year==1997),
                   sv['recruit','1998'], seed.bank.size= i  )
      fill.seed.97[i/1000,]<- c(i, x$lambda )
   }


plot(fill.seed.97[,1], fill.seed.97[,2], 
     xlab="Seed bank size", ylab="Growth rate " , type="b", pch=16, col="green",
     main=paste("Changes in seed bank size\non population growth rate in 1997"))


# ---------------------------------------------------------------- #
#  STOCHASTIC growth rates  -this will take a while

eigen.analysis(A.mean)$lambda1


stoch.growth.rate(fill.A)


# ---------------------------------------------------------------- #
#  Stochastic population sizes

sv.mean<-c(seed=10000, round(apply(sv, 1, mean),0)[-1])


##  set 2000 above-ground plants as max density in 15 m2 plots...
aq.sp1<-stoch.projection(fill.A[1:3], sv.mean, nreps=1000,tmax=25, nmax=2000, sumweight=c(0,1,1,1,1) )



hist(apply(aq.sp1[,-1], 1, sum), col="green", xlim=c(0,2000), xlab="Final population size (excluding seeds)", main=paste('Stochastic growth projections\nafter 25 years using 96-98 matrices'), breaks=seq(0,2000, 100))

## starting size
abline(v=sum(sv.mean[-1]), lty=2)

# ---------------------------------------------------------------- #
# QUASI-extinction

aq.ex.seed<-stoch.quasi.ext(fill.A, n0=sv.mean, Nx=1, nreps=500, sumweight=c(0,1,1,1,1))

matplot(aq.ex.seed, xlab="Years", ylab="Quasi-extinction probability", type="l",
     main=paste("Time to reach a quasi-extinction threshold
of 1 above-ground individual"))


#par(opar)







