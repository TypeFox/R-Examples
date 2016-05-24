dummy <- function(f, name = "f"){

n <- length(f)
t <- table(f)
nlevels <- length(t)
levels <- as.numeric(names(t))
dum <- matrix(NA, ncol = nlevels - 1, nrow = n)
for (l in 2:nlevels){dum[, l - 1] <- as.numeric(f == levels[l])}

dimnames(dum) <- list(NULL, paste(name, ".", levels[-1], sep = ""))
return(dum)
}
