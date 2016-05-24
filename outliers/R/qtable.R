"qtable" <-
function(p,probs,quants) {

quants <- quants[order(probs)];
probs <- sort(probs);
res <- vector();

for (n in 1:length(p)) {
pp <- p[n];

if (pp <= probs[1]) {
q0 <- quants[c(1,2)];
p0 <- probs[c(1,2)];
fit <- lm(q0 ~ p0)
}
else
if (pp >= probs[length(probs)]) { 
q0 <- quants[c(length(quants)-1,length(quants))];
p0 <- probs[c(length(probs)-1,length(probs))];
fit <- lm(q0 ~ p0)
}
else
{
x0 <- which(abs(pp - probs) == min(abs(pp - probs)));
x1 <- which(abs(pp - probs) == sort(abs(pp - probs))[2]);
x <- min(c(x0,x1));
if (x == 1) x <- 2;
if (x > length(probs)-2) x <- length(probs)-2;
i <- c(x-1,x,x+1,x+2)
q0 <- quants[i];
p0 <- probs[i];
fit <- lm(q0 ~ poly(p0,3));
}
res <- c(res,predict(fit,newdata=list(p0=pp)))
}
return(res)
}

