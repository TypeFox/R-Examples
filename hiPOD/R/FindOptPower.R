FindOptPower <-
function(cost, sample.size, MAF, OR, error, costPerExp=18915, costPerPool=970, costPerX=300, lower.P=20, upper.P=400, lower.N.p=2, upper.N.p=200, lower.Xmean=4, upper.Xmean=1280, sig.level=0.05, Number.Grids=100)
{


if(lower.P < 20 | upper.P > 400 | upper.N.p>200 | upper.Xmean > 1280) warning("Power is most accurately caculated with 20<P<400, 2<N.p<200 and 4<Xmean<1280.\n  Designs outside of these constraints may not have accurate predicted power.")

if(Number.Grids >= 200) warning("Too many number of grids! It might take too long! Best to have Number.Grids between 100 and 180.")

costPar <- c(costPerExp, costPerPool, costPerX)
constraint.set <- c(lower.P, upper.P, lower.N.p, upper.N.p, lower.Xmean, upper.Xmean)
scenario.set <- c(MAF, OR, error)
names(constraint.set) <- c("lower.P", "upper.P", "lower.N.p", "upper.N.p", "lower.Xmean", "upper.Xmean")
names(scenario.set) <- c("MAF", "OR", "error")
	
	
newdata.P <- round(seq(constraint.set["lower.P"],constraint.set["upper.P"],length.out=Number.Grids))
newdata.N.p <- round(seq(constraint.set["lower.N.p"],constraint.set["upper.N.p"],length.out=Number.Grids))

newdata <- data.frame(MAF=as.numeric(scenario.set["MAF"]), OR=as.numeric(scenario.set["OR"]), error=as.numeric(scenario.set["error"]), P=rep(newdata.P, each=length(newdata.N.p)), N.p=rep(newdata.N.p, length(newdata.P)))

Xmean.Best <- XmeanGivenCost(costPar, cost, newdata$P, constraint.set["lower.Xmean"], constraint.set["upper.Xmean"])

newdata$Xmean <- Xmean.Best[[2]]

temp <- t(apply(newdata, 1, function(x) diffVariantError(x[6], x[5], x[3], x[4]/2)))
newdata$T <- temp[,1]
newdata$prob.detect <- apply(newdata, 1, function(x) 1- (1-probDetEqual3(x[1], x[6], x[7], x[5], x[3]))^(x[4]/2) )


newdata$is.valid.design <- newdata$prob.detect>=0.9
newdata$upper.sample.good <- newdata$P*newdata$N.p<=sample.size
newdata$Xmean.good <- Xmean.Best[[1]]


newdata$pred.stat <- apply(newdata, 1, function(x) predictPoolNCP(x[1], x[2], x[3], x[4], x[5], x[6]))
newdata$pred.power <- with(newdata, pnorm(qnorm(1-sig.level), mean=pred.stat, lower.tail=FALSE))

opt.design.results <- list(cost, sample.size, constraint.set, scenario.set, newdata)

opt.design.results



}

