step_2A <-
function(par_2_3, N, mod,level_2,level_3, model, linear){

plot.step_2A <- function(level_2, level_3, level_2_sample, mod, model, par_2_3, N,linear){

op <- par(ask=TRUE)
for(i in 1:length(level_2)){
for(j in 1:length(level_3)){
y <- exp(level_2[,i]) #log data back transformed for plotting
x <- level_3[,j]
xu<- unique(level_3[,j])
if(!mod){
#pred.s <- x[rep(seq_len(length(x)), each=N)] #Xin's
pred.s <- xu[rep(seq_len(length(xu)), each=N)]
y.hat.new.PI <- exp(par_2_3[[i]]$fit[order(x),])
yup = max(exp(level_2[,i]),exp(level_2_sample[,i]),y.hat.new.PI)
ylo = min(exp(level_2[,i]),exp(level_2_sample[,i]),y.hat.new.PI)
plot(x, y, ylim=c(ylo,yup),col="black",lwd=2, xlab=names(level_3)[j],ylab=names(level_2)[i])
if(linear) abline(model$coefficients[1], model$coefficients[j+1], col="darkgreen")
lines(sort(x), y.hat.new.PI[,2], col="red", lty=2)
lines(sort(x), y.hat.new.PI[,3], col="red", lty=2)
points(pred.s, exp(level_2_sample[,i]), col="red", lwd=2)
}else{
#boxplot(y~x)
yup = max(exp(level_2[,i]),exp(level_2_sample[,i]))
ylo = min(exp(level_2[,i]),exp(level_2_sample[,i]))
plot(x, y, ylim=c(ylo,yup),col="black",lwd=2,xlab=names(level_3)[j],ylab=names(level_2)[i])
x <- unique(x)
pred.s <- x[rep(seq_len(length(x)), each=N)]
points(pred.s, exp(level_2_sample[,i]), col="red", lwd=2)

}
}
}
par(op)
}


if (mod){
# simulate based on the Mclust
level_2_sample <- matrix(0, nrow(unique(level_3))*N, ncol(level_2))
P_level_2_level_3.temp <- matrix(0,nrow(unique(level_3))*N,ncol(level_2))

# N simulation for each site
      for (n in 1:nrow(unique(level_3))){ 
level_2_sample[((n-1)*N+1):(n*N),] <- sim(par_2_3[[n]]$variance$modelName,par_2_3[[n]],N)[,2:(ncol(level_2)+1)]
    P_level_2_level_3.temp[((n-1)*N+1):(n*N),1] <- dens(par_2_3[[n]]$variance$modelName, data = level_2_sample[((n-1)*N+1):(n*N),], parameters=par_2_3[[n]])
}

result <- list(sample = level_2_sample,P_level_2_level_3 = P_level_2_level_3.temp[,1])

}else{
# simulate based on the Linear model
level_2_sample <- matrix(0, nrow(unique(level_3))*N, ncol(level_2))
P_level_2_level_3.temp <- matrix(0,nrow(unique(level_3))*N,ncol(level_2))
P_level_2_level_3 <- rep(0,nrow(level_2_sample))

# N simulation for each site

for (i in 1:ncol(level_2)){
      for (n in 1:nrow(unique(level_3))){  
         level_2_sample[((n-1)*N+1):(n*N),i] <- rnorm(N,unique(par_2_3[[i]]$fit[,1])[n],par_2_3[[i]]$residual.scale)# simulating trait values
# computing P_level_2_level_3
    P_level_2_level_3.temp[((n-1)*N+1):(n*N),i] <- dnorm(level_2_sample[((n-1)*N+1):(n*N),i],unique(par_2_3[[i]]$fit[,1])[n],par_2_3[[i]]$residual.scale)
  }
}

for (i in 1:nrow(level_2_sample)) {
  P_level_2_level_3[i] <- exp(sum(log(P_level_2_level_3.temp[i,])))
}
result <- list(sample = level_2_sample,P_level_2_level_3 = P_level_2_level_3)

}

plot.step_2A(log(level_2), level_3, level_2_sample, mod, model, par_2_3, N, linear)

return(result)
}
