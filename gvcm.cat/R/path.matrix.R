path.matrix <-
function(
x, 
y, 
weights,
family, 
control, 
acoefs,
lambda, 
phis, # vector length L
weight, # vector length L 
which.a, # vector length L
start,
offset,
coefficients, 
path, 
oml, 
...
)

{   
# new.lambdas <- round(c(lambda-lambda/10*(7:9), (lambda+lambda/10*(7:9))[which((lambda+lambda/10*(7:9))<control$lambda.upper)]), 2) 
nl1 <- c((lambda-lambda/10*(7:9))[which((lambda-lambda/10*(7:9))>control$lambda.lower)], lambda, (lambda+lambda/10*(7:9))[which((lambda+lambda/10*(7:9))<control$lambda.upper)])
nl2 <- if(control$tuning.criterion=="deviance") {
       seq(from=control$lambda.lower, to=control$lambda.upper, length.out=control$steps)[-c(1)]
       } else {
       seq(from=control$lambda.lower, to=control$lambda.upper, length.out=control$steps)[-c(1,control$steps)]
       }
new.lambdas <- round(c(nl1, nl2),2) # 2

new.path <- matrix(nrow=length(oml),ncol=length(new.lambdas))
 
for (i in 1:(length(new.lambdas))){

  opt <- gvcmcatfit(x, y, weights, family, control, acoefs, lambda=new.lambdas[i], 
                     phis, weight, which.a, start, offset)
  
  new.path[,i] <- opt$coefficients
  
  }
 
# add oml, coefficients
new.path <- cbind(oml, coefficients, new.path)
colnames(new.path) <- as.character(c(0, lambda, new.lambdas))
path <- if (!any(is.na(path))) cbind(new.path,path)  else new.path

# sortieren nach lambdas
path <- path[,order(as.numeric(colnames(path)))]

return(path) 
}

