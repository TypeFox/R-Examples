check.mc.matrix <-
function (mc.matrix, tol = .Machine$double.eps)
{
erg <- logical(length = length(mc.matrix))
for (i in 1: length(mc.matrix)) {
if (all(dim(mc.matrix[[i]]) == c(2, 2))) {
erg[i] <- (mc.matrix[[i]][1, 1] + mc.matrix[[i]][2, 2] > 1)
} else {
ev <- eigen(mc.matrix[[i]])
evalue <- ev[["values"]]
evectors <- ev[["vectors"]]
d <- diag(log(evalue))
mc <- zapsmall(evectors %*% d %*% solve(evectors))
diag(mc) <- 1
erg[i] <- all(mc > -tol)
}
}
return(erg)
}

