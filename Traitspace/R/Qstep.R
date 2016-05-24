Qstep <-
function(resp, pred){

# find the best predictors
qm.traitspace <- function(X, Y, log = TRUE){
X.q <- NULL
for(i in 1:ncol(X)){
X.2 <- poly(X[,i],2,raw =TRUE)
colnames(X.2) <- c(colnames(X)[i],paste(colnames(X)[i], "^2"))
X.q <- cbind(X.q, X.2)
}

# check
dat <- data.frame(na.omit(cbind(Y,X.q)))

# linear regression

if(log){
response.lm = lm(log(Y)~., data = dat)
}else{
response.lm = lm(Y~., data = dat)
}
result <- list(model.q = response.lm, pred = X.q)
return(result)
}


plotQstep <- function(resp.name, pred, model, par){
y <- model$model[[1]]
op <- par(ask=TRUE)
for(i in 1:length(pred)){
y.hat.new.PI <- par$fit[order(pred[,i]),]
pred.s <- sort(pred[,i])
x <- model$model[[i+1]]
plot(x, y, ylim=range(y.hat.new.PI),xlab=names(pred)[i],ylab=resp.name)
points(pred.s, y.hat.new.PI[,1], col="red", lwd=2)
lines(pred.s, y.hat.new.PI[,2], col="red", lty=2)
lines(pred.s, y.hat.new.PI[,3], col="red", lty=2)
}
par(op)
}

par <- list()
model <- list()
summary.lm <- list()
# regression model
for(i in 1:ncol(resp)){
qm.result <- qm.traitspace(X = pred, Y = resp[,i])
resp.lm <- qm.result$model.q
pred.temp <- data.frame(qm.result$pred)
par[i] <- list(predict.lm(resp.lm, newdata = data.frame(pred.temp), se.fit=TRUE,interval="prediction",level=0.95))
model <- c(model, resp.lm$coefficients)
plotQstep(names(resp)[i], pred.temp, resp.lm, par[[i]])                      
summary.lm <- c(summary.lm, summary(resp.lm)$ r.squared)
}
result <-list(summary.lm = summary.lm , par = par, model = model)
return(result)
}
