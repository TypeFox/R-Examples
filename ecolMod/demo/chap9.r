##############################
## Soetaert and Herman      ##
## ecological modelling     ##
## Figures from chapter 9   ##
## Discrete time models     ##
##############################

opar <- par()
par(ask=TRUE)
par(mfrow=c(1,1))
figNr <- 1
subtitle <- function()
{
 mtext(side=1,outer=TRUE,"Soetaert and Herman - chapter 9  ",cex=0.7,adj=1,line=-1.5)
 mtext(side=1,outer=TRUE,paste("  Fig. 9.",figNr,sep=""),cex=0.7,adj=0,line=-1.5)
 figNr <<- figNr +1
}

###############################################################################
####======================================================================#####
####                                Theory                                #####
####======================================================================#####
###############################################################################


########################################
# Figure 9.1. Discrete versus continuous time
########################################


par(mfrow=c(2,2),mar=c(2,3,3,1))
x  <- seq(0.1,12,length.out=100)
y  <- 0.1+(x^2/(x^2+3^2))*0.75
plot(x,y,xlim=c(0,15),ylim=c(0,1),xlab="",ylab="",axes=F,type="l",lwd=2,
main="continuous time")
Arrows(0,0.1,13,0.1,arr.length=0.25)
text(14,0.1,"t",cex=2)
Arrows(0,0.1,0,0.85,arr.length=0.25)
text(0,0.95,"C",cex=1.5)
textrect(mid=c(7.5,0.95),radx=4.5,rady=0.08,lab=expression(over(dC,dt)==source-sinks))
writelabel("A")
box(col="grey")

x  <- seq(0,10,2)
y  <- 0.1+x*0.07
plot(x,y,xlim=c(0,15),ylim=c(0,1),xlab="",ylab="",axes=F,type="S",lwd=2,
main="discrete time")
segments(x,0.1,x,y,lwd=2)
Arrows(0,0.1,13,0.1,arr.length=0.25)
text(14,0.1,"t",cex=2)
Arrows(0,0.1,0,0.85,arr.length=0.25)
text(0,0.95,"N",cex=1.5)
textrect(mid=c(7.5,0.95),radx=4.5,rady=0.08,lab=expression(N^{t+Delta~t}==N^t*p^t))

for (i in seq(2,8,2))
plotellipse(rx=0.8,ry=0.05,mid=c(i,0.02),dr=0.001,
            type="l",lwd=1,col=NULL,angle=0,
            from=-pi,to=0,lcol="black",arrow=TRUE,arr.adj=0,
            arr.type="triangle",arr.length=0.14)
box(col="grey")
writelabel("B")
subtitle()


##############################################################################
# Fig. 9.2 Host-parasitoid model
##############################################################################
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
openplotmat()
selfarrow (pos=c(0.8,0.7),curve=c(0.3,0.2),arr.pos=0.25)
#selfarrow (pos=c(0.8,0.7),curve=c(0.3,0.2),arr.pos=0.4)
selfarrow (pos=c(0.8,0.7),curve=c(0.3,0.2),arr.pos=0.6)
ss <- selfarrow (pos=c(0.8,0.7),curve=c(0.3,0.2),arr.pos=0.88)
text(ss[1]+0.03,ss[2]-0.02,"f",font=2)
ss <- selfarrow (pos=c(0.2,0.3),path="R",code=2,curve=c(0.3,0.2),arr.pos=0.12)
text(ss[1]+0.03,ss[2]+0.02,"1-f",font=2)
selfarrow (pos=c(0.2,0.3),path="R",code=2,curve=c(0.3,0.2),arr.pos=0.35)
selfarrow (pos=c(0.2,0.3),path="R",code=2,curve=c(0.3,0.2),arr.pos=0.75)
textellipse(c(0.8,0.7),radx=0.1,rady=0.1,lab="Host")
textellipse(c(0.2,0.3),radx=0.1,rady=0.1,lab="Parasitoid")
textrect(c(0.8,0.3),radx=0.12,rady=0.05,lab="Infected pupae")
textrect(c(0.2,0.7),radx=0.12,rady=0.05,lab="Larvae")
textrect(c(0.5,0.5),radx=0.12,rady=0.05,lab="Pupae")
#textrect(c(0.5,0.9),radx=0.12,rady=0.05,lab="Eggs")
mtext(side=3,"Host-Parasitoid",cex=1.25,line=-1)
subtitle()
##############################################################################
# Fig. 9.3 Age-structure population model
##############################################################################

par(mfrow=c(2,1),mar=c(0,0,0,0))

# Fecundity and Survival for each generation
Numgenerations   <- 6

# Original Population matrix M
DiffMat  <- matrix(data=0,nrow=Numgenerations,ncol=Numgenerations)   # declare it
AA <- as.data.frame(DiffMat)
AA[[1,4]]<- "f[3]"
AA[[1,5]]<- "f[4]"
AA[[1,6]]<- "f[5]"

AA[[2,1]]<- "s[list(0,1)]"
AA[[3,2]]<- "s[list(1,2)]"
AA[[4,3]]<- "s[list(2,3)]"
AA[[5,4]]<- "s[list(3,4)]"
AA[[6,5]]<- "s[list(4,5)]"

name  <- c("Age0","Age1","Age2","Age3","Age4","Age5")

PP <- plotmat(A=AA,pos=6,curve=0.6,name=name,lwd=2,
              box.size=0.05,arr.type="triangle",my=-0.1,dtext= 0.65,cex.txt=0)
for (i in 1:nrow(PP$arr))
  text(as.double(PP$arr[i,"TextX"]),as.double(PP$arr[i,"TextY"]),
  parse(text=as.character(PP$arr[i,"Value"])))
text(0.05,0.95,"A",cex=1.5)

# reduced population matrix
Numgenerations   <- Numgenerations-1
DiffMat       <- DiffMat[-2,-2]
AA <- as.data.frame(DiffMat)
AA[[1,3]]<- "f[3]*s[list(0,1)]"
AA[[1,4]]<- "f[4]*s[list(0,1)]"
AA[[1,5]]<- "f[5]*s[list(0,1)]"

AA[[2,1]]<- "s[list(0,2)]"
AA[[3,2]]<- "s[list(2,3)]"
AA[[4,3]]<- "s[list(3,4)]"
AA[[5,4]]<- "s[list(4,5)]"

name  <- c("Age0","Age2","Age3","Age4","Age5")

pos <- PP$comp[-1,]
PP <- plotmat(AA,pos=pos,curve=0.6,name=name,lwd=2,
              box.size=0.05,arr.type="triangle",my=-0.1,dtext= 0.65,cex.txt=0)
for (i in 1:nrow(PP$arr))
  text(as.double(PP$arr[i,"TextX"]),as.double(PP$arr[i,"TextY"]),
  parse(text=as.character(PP$arr[i,"Value"])))

text(0.05,0.95,"B",cex=1.5)
subtitle()

#############################################################################
# Figure 9.4 -  logistic model
#############################################################################

par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

# less detail than example in book

ricker  <- function(N,r) N*exp(r*(1-N)) 
rseq    <- seq(1.5,4,0.005) # sequence of r-values

plot(0,0,xlim=range(rseq),ylim=c(0,5),type="n",
     xlab="r",ylab="Nt",main="discrete logistic model")

for ( r in rseq)
 {
  N  <- runif(1)
  for (i in 1:100) N <- ricker(N,r)   # spinup 
  for (i in 1:100){N <- ricker(N,r)
               points(r,N,pch=".",cex=1.5)}
}

subtitle()

#############################################################################
# Figure 9.5 -  Parasitoid
#############################################################################

# less detail than in book 
par (mfrow=c(2,2),mar=c(5.1,4.1,4.1,2.1))

rH <- 2.82   # rate of increase
tS <- 100    # searching time
tH <- 1      # handling time
A  <- tS/tH  # attack rate
ks <- 30     # 1/tH*a

Parasite <- function(P_H,ks)
{
 P<-P_H[1] ;  H <- P_H[2]
 f  <- A*P/(ks+H)
 return(c(H*(1-exp(-f)),       
          H * exp(rH*(1-H)-f)))
}

out <- matrix(nrow=50,ncol=2)

plottraject<-function(ks)
{
P_H <- c(0.5,0.5)
for (i in 1:100) P_H <-Parasite(P_H,ks)
for (i in 1:50) {P_H <-Parasite(P_H,ks); out[i,]<-P_H }

plot (out[,1],type="l",ylim=range(out),lwd=2,xlab="t",ylab="Population",
      main=paste("ks=",ks))
lines(out[,2],lty=2)
}

#plottraject(35)


plottraject(25)
writelabel("A")
plottraject(20)
writelabel("B")
legend("topright",c("Parasitoid","Host"),lty=c(1,2),lwd=c(2,1))

ksSeq <- seq(15,35,0.1) # sequence of a-values
plot(0,0,xlim=range(ksSeq),ylim=c(0.,2),xlab="ks",ylab="Nt",main="Bifurcation diagram")

for ( ks in ksSeq)
{
P_H <- c(0.5,0.5)
for (i in 1:100)  P_H <-Parasite(P_H,ks)   # spinup steps
for (i in 1:200)  {P_H <-Parasite(P_H,ks); points(ks,P_H[2],pch=".",cex=1.5)}
}

writelabel("C")

# domain of attraction
ks   <- 23.09 
dz   <- 0.005 # 0.0025 
xlim <- c(0.001,0.5)
ylim <- c(0.001,0.5)

Initial <- expand.grid(P = seq(xlim[1],xlim[2],dz),
                       H = seq(ylim[1],ylim[2],dz))
plot(0,0,xlim=xlim,ylim=ylim,ylab="Parasitoid initial",xlab="Host initial",
     type="n",main="Domain of attraction")

PP   <- vector(length=100)

for ( ii in 1:nrow(Initial))
{
ini <- Initial[ii,]
P_H <- unlist(ini)
for (i in 1:100) P_H<-Parasite (P_H,ks)
for (i in 1:20) {P_H <-Parasite(P_H,ks); PP[i] <- P_H[1]}

Freq <- length(unique(trunc(PP*10)))
ifelse (Freq == 4,col<-"black",col<-"white")
rect(ini$P-dz/2,ini$H-dz/2,ini$P+dz/2,ini$H+dz/2,col=col,border=col)
}

writelabel("D")
subtitle()       

##############################################################################
# Fig. 9.6. Schematic diagram of Teasel stage model
##############################################################################
par(mfrow=c(1,1))
par(mar=c(2,2,2,2))
# Define the stages
Stagenames  <- c("dormant seeds 1yr","dormant seeds 2yr","small rosettes<2.5cm",
                 "medium rosettes 2.5-18.9cm","large rosettes>19cm ","flowering plants")
names       <- c("DS 1yr","DS 2yr","R small",
                      "R medium","R large","F")

NumStages   <- length(Stagenames)

# Population matrix 
A        <- matrix(nrow=NumStages,ncol=NumStages,byrow=TRUE,data = c(   # declare and fill it
                   0,      0,      0,      0,      0,      322.38,
                   0.966,  0,      0,      0,      0,      0     ,
                   0.013,  0.01,   0.125,  0,      0,      3.448 ,
                   0.007,  0,      0.125,  0.238,  0,      30.170,
                   0.008,  0,      0.038,  0.245,  0.167,  0.862 ,
                   0,      0,      0,      0.023,  0.75,   0      )  )

curves <- A
curves[]<-0
curves[3,1]<- -0.35             
curves[4,6]<-curves[6,4]<-curves[5,6]<-curves[6,5]<-0.08
curves[3,6]<-  0.35
curves[1,6]<- -0.35


plotmat(A,pos=c(3,2,1),curve=curves,name=names,lwd=1,box.lwd=2,my=-0.1,cex.txt=0.8, 
        box.cex=0.8,box.size=0.08,box.type="circle",box.prop=1, relsize=0.95,
        shadow.size = 0.01,self.cex=0.6,
        self.shiftx=c(0,0,0.125,-0.12,0.125,0),self.shifty=0)
legend("bottomright",c("DS = dormant seeds","R  = rosettes","F  = flowers"),cex=0.8)
legend("bottomleft",ncol=2,legend=c("small","medium","large","<2.5cm","2.5-18.9cm",">19cm"),cex=0.8)
mtext(side=3,"Dipsacus sylvestris",cex=1.25,line=0)
subtitle()


###############################################################################
####======================================================================#####
####                            R case studies                            #####
####======================================================================#####
###############################################################################


##############################################################################
# Fig. 9.7. Teasel stage model application
##############################################################################

#------------------------#
# Stable stage + lambda  #
#------------------------#
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

numsteps      <- 50
Population    <- rep(0,times=NumStages)
Population[1] <- 100
out           <- matrix(nrow=numsteps+1,ncol=NumStages,data=0) 
out[1,]       <- Population

for (i in 1:numsteps)
{
 Population <- A %*%Population    
 out[i+1,]  <- t(Population)  
}

#------------------------#
# PLOTTING  
#------------------------#

p   <- out/rowSums(out)
matplot(p,type="l",lty=1:6,main="Teasel stage distribution",col="black")
legend("topright", legend=Stagenames,lty=1:6)


###############################################################################
####======================================================================#####
####                               R projects                             #####
####======================================================================#####
###############################################################################


##############################################################################
# Fig. 9.8. population of long-lived individuals with 4 age classes
##############################################################################

par(mfrow=c(2,1))
par(mar=c(5.1,4.1,4.1,2.1))
# Fecundity and Survival for each generation
AgeClasses       <- 4
ReprodRate       <- c(2,3,1,0)    
Survival         <- c(0.8,0.5,0.25,0.0)          # survival from i to i+1

# Population matrix M
Leslie       <- matrix(data=0,nrow=AgeClasses,ncol=AgeClasses) # declare it
Leslie[1,]   <- ReprodRate * Survival[]                     # first row :fecundity
for (i in 1:(AgeClasses-1))  Leslie[i+1,i] <- Survival[i]              

Leslie                   # print the matrix to screen      
curves <- Leslie
curves[] <- 0.85
plotmat(Leslie,curve=curves,pos=4,arr.type="triangle",box.size=.075,self.shiftx=-0.075,mx=0.025,my=-0.01)
subtitle()

##------------------------------------------------------------------------
# Fig. 9.9. U.S. population 
##------------------------------------------------------------------------

par(mfrow=c(1,1))
par(mar=c(2,2,2,2))
# Fecundity and Survival for each generation
NumClass    <- 10
Fecundity   <- c(0,      0.00102,0.08515,0.30574,0.40002,
                 0.28061,0.1526 ,0.0642 ,0.01483,0.00089)
Survival    <- c(0.9967 ,0.99837,0.9978 ,0.99672,0.99607,
                 0.99472,0.99240,0.98867,0.98274,NA)            # survival from i to i+1    

cbind(Fecundity,Survival)

# Population matrix M
DiffMatrix       <- matrix(data=0,nrow=NumClass,ncol=NumClass)     # declare it
DiffMatrix[1,]   <- Fecundity                                      # first row: fecundity
for (i in 1:(NumClass-1))  DiffMatrix[i+1,i] <- Survival[i]              

DiffMatrix                                                         # print the matrix to screen  
names <- c("0-5yr","5-10yr","10-15yr","15-20yr","20-25yr","25-30yr","30-35yr","35-40yr","40-45yr","45-50yr")
# first generation will be positioned in middle; other generations on a circle
pos <- coordinates(NULL,N=NumClass-1)
pos <- rbind(c(0.5,0.5),pos)
curves <- DiffMatrix
curves[]   <- -0.4
curves[1, ] <- 0
curves[2,1] <- -0.125
curves[1,2] <- -0.125
plotmat(DiffMatrix,pos=pos,name=names,curve=curves, 
        box.size=0.07,arr.type="triangle",cex.txt=0.8,box.col=grey(0.95),box.prop =1)

mtext(side=3,"US population life cycle, 1966",cex=1.2)


subtitle()

par(ask=opar$ask)
par(mar=opar$mar)
par(oma=opar$oma)
par(mfrow=opar$mfrow)

