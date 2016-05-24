bootstrap <-
function(
x,
y,
indices,
family,
tuning,
weights,
offset,
start,
control,
method,
...
)
{
# definitions
B <- control$bootstrap
n <- nrow(x)
st <- matrix(nrow=ncol(x),ncol=B)

for (i in 1:B) {

bootstrap.sample <- sample(1:n, n, replace=TRUE)
x.i <- x[bootstrap.sample,]
y.i <- y[bootstrap.sample]

if (method %in% c("lqa"))
    model.i <- pest(x.i, y.i, indices, family, tuning, weights, offset, start, control, plot=FALSE)
    
if (method %in% c("AIC", "BIC"))
    model.i <- abc(x.i, y.i, indices, family, method, weights, control, plot=FALSE)

st[,i] <- model.i$coefficients
}

# var(coefs) <- 1/(B-1) sum (ceof.i- mean(coef))^2
mean.coefs <- rowSums(st)/B
var.coefs <- rowSums((st-mean.coefs)^2)/(B-1)
strd <- round(matrix(sqrt(var.coefs),ncol=1), digits=control$digits)
rownames(strd) <- rownames(model.i$coefficients)
colnames(strd) <- c("strd. error")

mins <- apply(st,1,min)
maxs <- apply(st,1,max)
minmax <- round(matrix(c(mins, maxs), ncol=2, byrow=FALSE), digits=control$digits)
rownames(minmax) <- rownames(model.i$coefficients)
colnames(minmax) <- c("min", "max")

output <- list(strd.error = strd, minmax = minmax)
return(output)
}