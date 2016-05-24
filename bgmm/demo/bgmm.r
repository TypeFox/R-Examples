#
# demo for real dataset
#
data(genotypes)
par(mfrow=c(2,3),pty="s")
modelSupervised = supervised(knowns=genotypes$knowns, class=genotypes$labels); plot(modelSupervised)
modelSemiSupervised = semisupervised(X=genotypes$X, knowns=genotypes$knowns, class = genotypes$labels); plot(modelSemiSupervised)
modelBelief = belief(X=genotypes$X, knowns=genotypes$knowns, B=genotypes$B); plot(modelBelief)
modelSoft = soft(X=genotypes$X, knowns=genotypes$knowns, P=genotypes$B); plot(modelSoft)
modelUnSupervised = unsupervised(X=genotypes$X, k=3); plot(modelUnSupervised)
par(mfrow=c(1,1))

# predictions
preds = predict(modelSoft, X = genotypes$X)
str(preds)

#
# simulated dataset 2d
#
set.seed(1313)
simulated = simulateData(d=2, k=3, n=300, m=60, cov="0", within="E", n.labels=2)

model = belief(X = simulated$X, knowns = simulated$knowns, B=simulated$B)
plot(model)

#
# list of models
# choosing structure
#
models1 = mModelList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, kList=3, mean=c("D","E"), between=c("D","E"), within=c("D","E"), cov=c("D","0"), funct=belief)

plot(models1)
par(mfrow=c(1,1))
plotGIC(models1, penalty="BIC")

#
# best model from the list
#
str(bestModel <- chooseOptimal(models1, penalty="BIC"))
                                                   
#
# choosing number of components
#
models2 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, kList=2:7, mean="D", between="D", within="E", cov="0")
plot(models2)
par(mfrow=c(1,1))
plotGIC(models2, penalty="BIC")

#
# choosing number of components and structure
#
models3 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, kList=2:7, mean="D")
plotGIC(models3, penalty="BIC", plot.it=FALSE)

models4 = chooseModels(models3, kList=2:5, struct=c("DDDD","DDED","DDE0"))
plot(models4)
par(mfrow=c(1,1))
plotGIC(models4, penalty="BIC")


#
# simulated dataset 1d
#
set.seed(1313)
simulated = simulateData(d=1, k=3, n=300, m=60, n.labels=2)

model = belief(X = simulated$X, knowns = simulated$knowns, B=simulated$B)
plot(model)

#
# list of models
# choosing structure
#
models1 = mModelList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, kList=3, mean=c("D","E"), between=c("D","E"), within=c("D","E"), cov=c("D","0"), funct=belief)

plot(models1)
par(mfrow=c(1,1))
plotGIC(models1, penalty="BIC")

#
# choosing number of components
#
models2 = beliefList(X=simulated$X, knowns=simulated$knowns, B=simulated$B, kList=2:7, mean="D", between="D", within="E", cov="0")
plot(models2)
par(mfrow=c(1,1))
plotGIC(models2, penalty="BIC")

