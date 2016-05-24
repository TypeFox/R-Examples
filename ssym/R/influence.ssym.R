influence.ssym <-
function(model, ...){
mu <- list(case.weights=model$cw, response.p=model$pr)
theta <- list(case.weights=model$cw.theta, response.p=model$pr.theta)

colnames(mu$case.weights) <- c("local.influence","total.local.influence")
colnames(theta$case.weights) <- c("local.influence","total.local.influence")
colnames(mu$response.p) <- c("local.influence","total.local.influence")
colnames(theta$response.p) <- c("local.influence","total.local.influence")

par(mfrow=c(2,2))
if(model$censored==FALSE){
plot(model$cw.theta[,1], type="h", main="Case-weight perturbation", xlab="Index", ylab="Local Influence")
plot(model$cw.theta[,2], type="h", main="Case-weight perturbation", xlab="Index", ylab="Total Local Influence")
plot(model$pr.theta[,1], type="h", main="Response perturbation", xlab="Index", ylab="Local Influence")
plot(model$pr.theta[,2], type="h", main="Response perturbation", xlab="Index", ylab="Total Local Influence")}
else{
linea <- model$cw.theta[,1]
linea[model$event==1] <- 0
plot(linea, type="h", main="Case-weight perturbation", xlab="Index", ylab="Local Influence", ylim=c(0,max(model$cw.theta[,1])))
par(new=TRUE)
linea <- model$cw.theta[,1]
linea[model$event==0] <- 0
plot(linea, type="h", lty=3, main="Case-weight perturbation", xlab="Index", ylab="Local Influence", ylim=c(0,max(model$cw.theta[,1])))

linea <- model$cw.theta[,2]
linea[model$event==1] <- 0
plot(linea, type="h", main="Case-weight perturbation", xlab="Index", ylab="Total Local Influence", ylim=c(0,max(model$cw.theta[,2])))
par(new=TRUE)
linea <- model$cw.theta[,2]
linea[model$event==0] <- 0
plot(linea, type="h", lty=3, main="Case-weight perturbation", xlab="Index", ylab="Total Local Influence", ylim=c(0,max(model$cw.theta[,2])))

linea <- model$pr.theta[,1]
linea[model$event==1] <- 0
plot(linea, type="h", main="Response perturbation", xlab="Index", ylab="Local Influence", ylim=c(0,max(model$pr.theta[,1])))
par(new=TRUE)
linea <- model$pr.theta[,1]
linea[model$event==0] <- 0
plot(linea, type="h", lty=3, main="Response perturbation", xlab="Index", ylab="Local Influence", ylim=c(0,max(model$pr.theta[,1])))

linea <- model$pr.theta[,2]
linea[model$event==1] <- 0
plot(linea, type="h", main="Response perturbation", xlab="Index", ylab="Total Local Influence", ylim=c(0,max(model$pr.theta[,2])))
par(new=TRUE)
linea <- model$pr.theta[,2]
linea[model$event==0] <- 0
plot(linea, type="h", lty=3, main="Response perturbation", xlab="Index", ylab="Total Local Influence", ylim=c(0,max(model$pr.theta[,2])))
}
list(mu=mu, theta=theta)
}
