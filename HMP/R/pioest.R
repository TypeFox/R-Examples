pioest <-
function(group.data){
if(missing(group.data))
stop("data.groups missing.")

par.groups <- lapply(group.data, function(x){
p <- DM.MoM(x)
p$reads <- rowSums(x)
return(p) 
})

K <- length(par.groups[[1]]$pi)
n.groups <- length(par.groups)
group.parameter <- lapply(par.groups, function(x){c(x$pi, x$theta, x$reads)})
index <- as.matrix(seq(1:n.groups))

Xscg <- apply(index, 1, function(x) {
pi <- group.parameter[[x]][1:K]
theta <- group.parameter[[x]][K+1]
P <- length(group.parameter[[x]])
nreads.data <- group.parameter[[x]][(K+2):P]
N_1Cj <- (theta * (sum(nreads.data^2)-sum(nreads.data)) + sum(nreads.data))/(sum(nreads.data))^2
Out1 <- c(pi, N_1Cj)
return(Out1)
})

Xscg <- Xscg[rowSums(Xscg)!=0,]
K <- nrow(Xscg)-1
pi0 <- (Xscg[1:K,]%*%as.matrix(1/Xscg[K+1,]))/sum(1/Xscg[K+1,])

return(pi0)
}
