
#########################################################################
####
#### R-script used to produce figures, tables or numbers for the article:
#### Kones, J., Soetaert, K., van Oevelen D. and J.O. Owino, 2009.
#### Are network indices robust indicators of food web functioning?
#### a Monte Carlo approach
#### Ecological Modelling. In press.
####
#### Script written by Karline Soetaert
####
#########################################################################

# the packages to be loaded...
require(LIM)         # contains the food web specifications
require(NetIndices)  # network indices

# First the script in Montecarlo.r has to be run.
# This creates the MCA sample and estimates the indices.

# alternatively: make this directory the working directory, e.g.
# setwd("c:/karline/netindices/")

# then run the script by executing the following command:
# source("Montecarlo.r")



#==============================================================================
# Summary statistics: mean, standard deviations
#==============================================================================
# Function that estimates the mean, standard deviation, minimum and maximum
# of the indices from the Monte Carlo samples
summstats <- function (Indices)
 {
  DD<- data.frame(pars=Indices[1,],
                  mean=apply(Indices,2,mean),
                  sd  =apply(Indices,2,sd),
                  min =apply(Indices,2,min),
                  max =apply(Indices,2,max))
  DD$CoV                     <- DD$sd/DD$mean
  # cannot calculate CoV for too small values of mean (roundoff)
  DD$CoV[is.infinite(DD$CoV)]<-NA
  DD$CoV[DD$mean<1e-8]<-NA
  DD
 }

summTakapoto        <- summstats(Takapoto.X)
summRigaSummer      <- summstats(RigaSummer.X)
summRigaAutumn      <- summstats(RigaAutumn.X)
summRigaSpring      <- summstats(RigaSpring.X)


summIndTakapoto     <- summstats(IndTakapoto  $Sub[,1:17])
summIndRigaSummer   <- summstats(IndRigaSummer$Sub[,1:17])
summIndRigaAutumn   <- summstats(IndRigaAutumn$Sub[,1:17])
summIndRigaSpring   <- summstats(IndRigaSpring$Sub[,1:17])

summTrophTakapoto   <- summstats(IndTakapoto  $Full[,c("TL4","TL5","OI4","OI5")])
summTrophRigaSummer <- summstats(IndRigaSummer$Full[,c("TL4","TL5","OI4","OI5")])
summTrophRigaAutumn <- summstats(IndRigaAutumn$Full[,c("TL4","TL5","OI4","OI5")])
summTrophRigaSpring <- summstats(IndRigaSpring$Full[,c("TL4","TL5","OI4","OI5")])

Flows<-rbind(Taka=summary(summTakapoto$CoV),
      RigaSum=summary(summRigaSummer$CoV),
      RigaAut=summary(summRigaAutumn$CoV),
      RigaSpr=summary(summRigaSpring$CoV))
Indi<-rbind(Taka=summary(summIndTakapoto$CoV),
      RigaSum=summary(summIndRigaSummer$CoV),
      RigaAut=summary(summIndRigaAutumn$CoV),
      RigaSpr=summary(summIndRigaSpring$CoV))
Trophic<-rbind(Taka=summary(summTrophTakapoto$CoV),
      RigaSum=summary(summTrophRigaSummer$CoV),
      RigaAut=summary(summTrophRigaAutumn$CoV),
      RigaSpr=summary(summTrophRigaSpring$CoV))
      
#-------------------------------------------------------------------------------
# Table 3.
#-------------------------------------------------------------------------------
data.frame(Fmean=Flows[,"Mean"],Fmax=Flows[,"Max."],IndMean=Indi[,"Mean"],IndMax=Indi[,"Max."],
TrophMean=Trophic[,"Mean"],TrophMax=Trophic[,"Max."])

#-------------------------------------------------------------------------------
# Table 4.
#-------------------------------------------------------------------------------

isel <- c("T..","TST","Ltot","LD","C","Tijbar","TSTbar","Cbar","AMI","HR","DR",
"RU","Hc","CE","Ascendency","Capacity","Overhead","ACratio","TSTC","TSTS","FCI",
        "APL","HP","BC","ID")
        
table4<-
cbind(Taka=summstats(IndTakapoto  $Full[,isel])[,1:3],
RigaSumm  =summstats(IndRigaSummer$Full[,isel])[,1:3],
RigaAut   =summstats(IndRigaAutumn$Full[,isel])[,1:3],
RigaSpr   =summstats(IndRigaSpring$Full[,isel])[,1:3])

#ordering in table changed:
ii <- c(1:8,19:22,9:18,23:25)

(FF <- format(round(table4[ii,]*100)/100,scientific=FALSE,digits=2))


# parts of table 4 - to make it more readable....
data.frame(name=rownames(FF),FF$RigaSpr.mean, "+-",FF$RigaSpr.sd ,FF$RigaSpr.pars)
data.frame(name=rownames(FF),FF$RigaSumm.mean,"+-",FF$RigaSumm.sd,FF$RigaSumm.pars)
data.frame(name=rownames(FF),FF$RigaAut.mean, "+-",FF$RigaAut.sd ,FF$RigaAut.pars)
data.frame(name=rownames(FF),FF$Taka.mean,    "+-",FF$Taka.sd    ,FF$Taka.pars)

# fraction of network indices < MCA
Smaller<-with (as.data.frame(table4),
    sum(Taka.pars<Taka.mean)      +sum(RigaSumm.pars<RigaSumm.mean)+
    sum(RigaAut.pars<RigaAut.mean)+sum(RigaSpr.pars<RigaSpr.mean))
Larger<-with (as.data.frame(table4),
    sum(Taka.pars>Taka.mean)      +sum(RigaSumm.pars>RigaSumm.mean)+
    sum(RigaAut.pars>RigaAut.mean)+sum(RigaSpr.pars>RigaSpr.mean))


print(c("percent of MN smaller than MCA = ",Smaller/(nrow(table4)*4)*100))   # percentage MN smaller
print(c("percent of MN larger than MCA = ",Larger/(nrow(table4)*4)*100))     # larger

(1-(Smaller+Larger)/(nrow(table4)*4))*100    # and equal...

#-------------------------------------------------------------------------------
# probability of the minimum norm solution with respect to the median...
#-------------------------------------------------------------------------------
probMN <- function(indices)
{
  nc  <- ncol(indices)
  nr  <- nrow(indices)
  prob<-NULL

  for (i in 1:nc)
  {
    DD <- density(indices[,i])
    cy <- cumsum(DD$y)
    cy <- cy/max(cy)
    cx <- DD$x

    # MN is the first row of the matrix - find to which index of the density
    # values it corresponds (i.e. ii for which MN > cx[ii] and MN <cx[ii+1])
    ii <- which(indices[1,i]>=cx[1:(nr-1)] & indices[1,i]<cx[2:nr])

    # is even smaller or larger than the smallest/largest value
    #( note: this does not happen with 3000 samples, but may happen with fewer)
    if (length(ii)==0) ii <- 1
    
    # The probability associated with that index is added to the results
    prob<-c(prob,cy[ii])
  }
  prob[prob>0.5]<-1-prob[prob>0.5]
  names(prob)<-colnames(indices)
  return(prob)
}

pMN <-
cbind(Taka=probMN(IndTakapoto  $Full[,isel]),
  RigaSumm=probMN(IndRigaSummer$Full[,isel]),
  RigaAut =probMN(IndRigaAutumn$Full[,isel]),
  RigaSpri=probMN(IndRigaSpring$Full[,isel])) [ii,]

# significance level for all indices...
format(round(pMN*1000)/1000,digits=3)

# and corresponding symbols
DF <- as.data.frame(pMN)
i001<-which (DF<0.001,arr.ind=TRUE)
DF[DF<0.001] <- 4
i005<-which (DF<0.005,arr.ind=TRUE)
DF[DF<0.005] <- 3
i01<-which (DF<0.01,arr.ind=TRUE)
DF[DF<0.01] <- 2
i05<-which (DF<0.05,arr.ind=TRUE)
DF[DF<0.05] <- 1
DF[DF<1] <- 0

DF

pMNtaka <-probMN(IndTakapoto  $Full[,isel])
sum(pMNtaka<=0.05)/length(pMNtaka)

# fraction of MN network indices significantly under- or overestimated
print(c("fraction significantly different from MCA = ",sum(pMN<=0.05)/length(pMN)))

# fraction of MN network indices significantly underestimated
Small<-with (as.data.frame(table4),
    cbind(Taka.pars<Taka.mean,RigaSumm.pars<RigaSumm.mean,
    RigaAut.pars<RigaAut.mean,RigaSpr.pars<RigaSpr.mean)
)
#significantly smaller
print(c("significantly smaller at p=0.05", sum(pMN<=0.05 & Small)/length(pMN)))
print(c("significantly smaller at p=0.001", sum(pMN<=0.001 & Small)/length(pMN)))


#-------------------------------------------------------------------------------
# Table 3.
#-------------------------------------------------------------------------------

summCov <- function(web,ind,trophic)

{ c(flowmean=mean(web$CoV,na.rm=TRUE)   ,flowmax=max(web$CoV,na.rm=TRUE),
    indmean =mean(ind$CoV,na.rm=TRUE)   ,indmax =max(ind$CoV,na.rm=TRUE),
    tromean =mean(trophic$CoV,na.rm=TRUE),tromax=max(trophic$CoV,na.rm=TRUE))
}

table3<-
rbind(Takapoto  = summCov(summTakapoto  ,summIndTakapoto  ,summTrophTakapoto),
      RigaSummer= summCov(summRigaSummer,summIndRigaSummer,summTrophRigaSummer),
      RigaAutumn= summCov(summRigaAutumn,summIndRigaAutumn,summTrophRigaAutumn),
      RigaSpring= summCov(summRigaSpring,summIndRigaSpring,summTrophRigaSpring))

format(round(table3*100)/100,digits=2)

# selection for comparison of various food webs
compselTak <- c("ACratio","CE","FCI","BC","ID","TL3","TL5","OI3","OI5")
compsel    <- c("ACratio","CE","FCI","BC","ID","TL4","TL5","OI4","OI5")

##########################################
# BOX-WHISKER PLOTS                      #
##########################################
# ps: zoo is mesozoo in Takapoto
nn <- c("AC=A/DC","CE=Hc/Hmax","FCI=TSTC/TST","B/C","I/D","TLNano","TLZoo","OINano","OIZoo")
indices <- NULL
for (i in 1:length(nn))
{
  it <- compselTak[i]
  ir <- compsel   [i]

  indices  <- rbind(indices,data.frame(web="Taka"   ,index=nn[i],val=IndTakapoto  $Full[,it]))
  indices  <- rbind(indices,data.frame(web="RigaSum",index=nn[i],val=IndRigaSummer$Full[,ir]))
  indices  <- rbind(indices,data.frame(web="RigaAut",index=nn[i],val=IndRigaAutumn$Full[,ir]))
  indices  <- rbind(indices,data.frame(web="RigaSpr",index=nn[i],val=IndRigaSpring$Full[,ir]))
}
require(shape)
windows(width=8.5,height=11)   ## A4-size; Figure 1
par(mfrow=c(3,2))
for (i in 1:5)
{
 boxplot(val~web,data=indices,subset=index==nn[i],main=nn[i],
 names=c("Takapoto","GRSummer","GRAutumn","GRSpring"))
 writelabel(nr=i)
 }

windows()                      ## Figure 2

par(mfrow=c(2,2))
for (i in 1:4)
{
 boxplot(val~web,data=indices,subset=index==nn[5+i],main=nn[5+i],
 names=c("Takapoto","GRSummer","GRAutumn","GRSpring"))
 writelabel(nr=i)
 }


