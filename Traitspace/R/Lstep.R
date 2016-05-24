Lstep <-
function(resp, pred){

# find the best predictors
lm.traitspace <- function(X, Y, log = TRUE){

dat <- na.omit(cbind(Y,X))
# linear regression
if(log){
response.lm = lm(log(Y)~., data = dat)
}else{
response.lm = lm(Y~., data = dat)
}
return(response.lm)
}

plotLstep <- function(resp.name, pred, model, par){
y <- model$model[[1]]
op <- par(ask=TRUE)
for(i in 1:length(pred)){
y.hat.new.PI <- par$fit[order(pred[,i]),]
pred.s <- sort(pred[,i])
x <- model$model[[i+1]]
plot(x, y, ylim=range(y.hat.new.PI),xlab=names(pred)[i],ylab=resp.name)
points(pred.s, y.hat.new.PI[,1], col="red", lwd=2)
abline(model$coefficients[1], model$coefficients[i+1], col="darkgreen")
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
resp.lm <- lm.traitspace(X = pred, Y = resp[,i]) 
par[i] <- list(predict.lm(resp.lm, newdata = data.frame(pred), se.fit=TRUE,interval="prediction",level=0.95))
model <- c(model, resp.lm$coefficients)
plotLstep(names(resp)[i], pred, resp.lm, par[[i]])          
summary.lm <- c(summary.lm, summary(resp.lm)$ r.squared)
}
result <-list(summary.lm = summary.lm , par = par, model = model)
return(result)
}
