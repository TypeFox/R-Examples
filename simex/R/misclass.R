misclass <-
function (data.org, mc.matrix, k = 1)
{
if (!is.list(mc.matrix))
stop("mc.matrix must be a list", call. = FALSE)
if (!is.data.frame(data.org))
stop("data.org must be a dataframe", call. = FALSE)
if (!all(names(mc.matrix) %in% colnames(data.org)))
stop("Names of mc.matrix and colnames of data.org do not match",
call. = FALSE)
if (k < 0)
stop("k must be positive")
data.mc <- data.org
factors <- lapply(data.org, levels)
ev <- lapply(mc.matrix, eigen)
data.names <- colnames(data.org)
for (j in data.names) {
evalue <- ev[[c(j,"values")]]
evectors <- ev[[c(j,"vectors")]]
d <- diag(evalue)
mc <- zapsmall(evectors %*% d^k %*% solve(evectors))
dimnames(mc) <- dimnames(mc.matrix[[j]])
for (i in factors[[j]]) {
data.mc[[j]][data.org[[j]] == i] <- sample(x = factors[[j]],
size = length(data.org[[j]][data.org[[j]] == i]),
prob = mc[, i], replace = TRUE)
}
}
return(data.mc)
}

