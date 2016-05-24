cocktailEnsemble <-
function(x,y){

options(warn=-1)

#logistic regression
if (is.factor(y)==TRUE) yglm <- ifelse(as.integer(y)==1,0,1) else yglm <- y
lr <- glm(y ~ ., data=data.frame(x,y=y), family=binomial("logit"))
lrstep <- step(lr, direction="both", trace = FALSE)

#random forest
rF <- randomForest(x,as.factor(y),  ntree=1000, importance=FALSE, na.action=na.omit )

#ada boost
ab <- ada(x,as.factor(y),iter=100)

#kernelFactory
kF <- kernelFactory(x,as.factor(y))

options(warn=0)

result <- list(lrstep,rF,ab,kF)
class(result) <- "cocktailEnsemble"
result
}
