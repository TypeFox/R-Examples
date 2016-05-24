### R code from vignette source 'comclim.Rnw'

###################################################
### code chunk number 1: comclim.Rnw:17-18
###################################################
library(comclim)


###################################################
### code chunk number 2: comclim.Rnw:23-47
###################################################
num_climateaxes = 3
num_regionalpool = 100
num_occurrences = 50

climateniches <- NULL
for (i in 1:num_regionalpool)
{
  randdata = NULL
	for (j in 1:num_climateaxes)
	{
		meanpos = runif(num_climateaxes,min=2,max=4)
		tcol = rnorm(num_occurrences, mean=meanpos[j]
            + runif(n=1,min=-2,max=2), sd=runif(1, 0.2,0.4))
		randdata <- cbind(randdata, tcol)
	}
	
	randdata <- as.data.frame(randdata)
	
	names(randdata) <- paste("ClimateAxis", 1:num_climateaxes, sep='')
	randdata$taxon = paste("Species", i, collapse='')
	
	climateniches <- rbind(climateniches, randdata)
}
climateniches$taxon <- factor(climateniches$taxon)


###################################################
### code chunk number 3: comclim.Rnw:52-54
###################################################
climateniches[,1:num_climateaxes] <- 
    scale(climateniches[,1:num_climateaxes], center=TRUE, scale=TRUE)


###################################################
### code chunk number 4: comclim.Rnw:59-60
###################################################
print(str(climateniches))


###################################################
### code chunk number 5: comclim.Rnw:66-81
###################################################
num_community = 5

nichedist <- do.call("rbind",by(
  climateniches[,1:num_climateaxes], 
  climateniches$taxon, function(x) { 
      cm <- colMeans(x)
      cm <- cm- rep(1, num_climateaxes); 
      return(data.frame(pos=sqrt(sum(cm^2))))
    }
  ))

# select for species on the lower edge of the climate space
whichsp <- order(nichedist,decreasing=FALSE)[1:num_community]
localcommunity <- row.names(nichedist)[whichsp]
print(localcommunity)


###################################################
### code chunk number 6: comclim.Rnw:87-89
###################################################
regionalpool <- as.character(levels(climateniches$taxon))
print(regionalpool)


###################################################
### code chunk number 7: comclim.Rnw:94-97
###################################################
observedclimate <- rep(-1, num_climateaxes)
names(observedclimate) <- paste("ClimateAxis", 1:num_climateaxes, sep='')
print(observedclimate)


###################################################
### code chunk number 8: comclim.Rnw:103-110
###################################################
cci <- inputcommunitydata(
  	localcommunity = localcommunity, 
		regionalpool = regionalpool, 
		climateniches = climateniches, 
		observedclimate = observedclimate)

summary(cci)


###################################################
### code chunk number 9: figure1
###################################################
plot(cci,cex.community=0.75)


###################################################
### code chunk number 10: comclim.Rnw:126-129
###################################################
result_community <- communityclimate(cci,
  climateaxes=c("ClimateAxis1","ClimateAxis2","ClimateAxis3"),
  numreplicates=100, verbose=F)


###################################################
### code chunk number 11: comclim.Rnw:135-136
###################################################
summary(result_community)


###################################################
### code chunk number 12: figure2
###################################################
plot(result_community)


###################################################
### code chunk number 13: figure3
###################################################
plot(result_community, deviations=TRUE)


