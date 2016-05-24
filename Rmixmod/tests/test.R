library(Rmixmod)

### Models

# Gaussian models
all <- mixmodGaussianModel()
all

# Multinomial models
all <- mixmodMultinomialModel()
all

only.free.proportions <- mixmodMultinomialModel(equal.proportions = FALSE)
list.models <- mixmodMultinomialModel(listModels = c("Binary_p_E","Binary_p_Ekjh", "Binary_pk_Ekj", "Binary_pk_Ej"))

var.independent <- mixmodMultinomialModel(variable.independency = TRUE)
var.comp.independent <- mixmodMultinomialModel(variable.independency = TRUE, component.independency = TRUE)

# Strategy
mixmodStrategy()
strategy1 <- mixmodStrategy(algo = "CEM", initMethod = "random", nbTry = 10, epsilonInInit = 0.000001)
strategy2 <- mixmodStrategy(algo = c("SEM", "EM"), nbIterationInAlgo = c(200, 100), epsilonInAlgo = c(NA, 0.0001))

###  Unsupervised classification

# Continuous variables: Geyser dataset
data("geyser")
xem.geyser <- mixmodCluster(data = geyser, nbCluster = 2:8,criterion = c("BIC", "ICL", "NEC"), models = mixmodGaussianModel(),strategy = strategy2)
xem.geyser
summary(xem.geyser)
plot(xem.geyser)
icl <- sortByCriterion(xem.geyser, "ICL")
icl["bestResult"]
#xem.geyser["results"]
#icl["results"]

# Categorical variables: Birds of different subspecies
data("birds")
xem.birds <- mixmodCluster(birds, 2)
barplot(xem.birds)

### Supervised classification

# First step: Learning
data("finance")
ratios2002 <- finance[finance["Year"] == 2002, 3:6]
health2002 <- finance[finance["Year"] == 2002, 2]
ratios2003 <- finance[finance["Year"] == 2003, 3:6]
health2003 <- finance[finance["Year"] == 2003, 2]
learn <- mixmodLearn(ratios2002, health2002)
learn["bestResult"]
plot(learn)

# Second step: Prediction
prediction <- mixmodPredict(data = ratios2003, classificationRule = learn["bestResult"])
summary(prediction)
mean(as.integer(health2003) == prediction["partition"])

# semi-supervised Learning
birds.partition<-c(rep(1,34),rep(2,35))
birds.partition[10]<-0
birds.partition[20]<-0
birds.partition[30]<-0
birds.partition[40]<-0
birds.partition[50]<-0
birds.partition[60]<-0
mixmodCluster(birds, 2, knownLabels=birds.partition)