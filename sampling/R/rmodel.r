rmodel<-function(formula,weights,X)
{
cl <- match.call() 
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "weights"), names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
y <- model.response(mf, "numeric")
w <- as.vector(model.weights(mf))
x <- model.matrix(mt, mf, contrasts)
prob<-glm(y~x,family="binomial",weights=w)$fitted.values
result<-cbind.data.frame(X,prob)
names(result)<-c(names(X),"prob_resp")
result
}
