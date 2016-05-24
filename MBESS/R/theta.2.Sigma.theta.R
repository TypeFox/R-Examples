`theta.2.Sigma.theta` <-
function (model, theta, latent.vars){    
if(is.null(names(theta))) stop("The elements in 'theta' must have names")
if(sum(is.na(names(theta))) !=0) stop ("Some of the elements in 'theta' do not have names")
if(length(unique(names(theta))) != length(theta)) stop("More than one element in 'theta' has the same name")

t<- length(theta)
r<- length(latent.vars)
T <- rep(0, t)
R <- rep(0, r)

for(i in 1:t){
    T[i]<- sum(na.omit(unique(model[,2])==names(theta)[i]))
    }
if(sum(T) != t) stop ("Some elements in 'theta' have not been included in 'model'")

    parse.path <- function(path) {                                           
        path.1 <- gsub('-', '', gsub(' ','', path))
        direction <- if (regexpr('<>', path.1) > 0) 2 
            else if (regexpr('<', path.1) > 0) -1
            else if (regexpr('>', path.1) > 0) 1
            else stop(paste('ill-formed path:', path))
        path.1 <- strsplit(path.1, '[<>]')[[1]]
        list(first=path.1[1], second=path.1[length(path.1)], direction=direction)
        }
    if ((!is.matrix(model)) | ncol(model) != 3) stop ("'model' must be a 3-column matrix in RAM manner")
    startvalues <- as.numeric(model[,3])
    par.names <- model[,2]
    n.paths <- length(par.names)
    heads <- from <- to <- rep(0, n.paths)
    for (p in 1:n.paths){
        path <- parse.path(model[p,1])
        heads[p] <- abs(path$direction)
        to[p] <- path$second
        from[p] <- path$first
        if (path$direction == -1) {
            to[p] <- path$first
            from[p] <- path$second
            }
        }
    
all.vars <- unique(c(to, from))

for(i in 1:r){
    R[i] <- sum(all.vars== latent.vars[i])
    }
if(sum(R)!=r) stop ("Some elements in 'latent.vars' have not been included in the 'model'")
    
obs.vars <- setdiff(all.vars, latent.vars)        
vars<- c(obs.vars, latent.vars)
pars <- names(theta)

    ram <- matrix(0, p, 5)    
    ram[,1] <- heads
    ram[,2] <- apply(outer(vars, to, '=='), 2, which)
    ram[,3] <- apply(outer(vars, from, '=='), 2, which)   
    par.nos <- apply(outer(pars, par.names, '=='), 2, which)
    if (length(par.nos) > 0)
        ram[,4] <- unlist(lapply(par.nos, function(x) if (length(x) == 0) 0 else x))
    ram[,5]<- startvalues
    colnames(ram) <- c('heads', 'to', 'from', 'parameter', 'start')

    n <- length(obs.vars)
    m <- length(all.vars)
    t <- length(pars)

result<- list()
result$ram<- ram
result$t <- t
result$m <- m
result$n <- n
result$all.vars <- vars
result$obs.vars <- obs.vars
result$latent.vars <- latent.vars
result$pars <- pars

one.head <- ram[,1] == 1
ram.A <- ram[one.head, c(2:5), drop=FALSE]
ram.P <- ram[!one.head, c(2:5), drop=FALSE]

P<- A<- matrix(0, m, m)

L<- dim(ram.A)[1]
for(l in 1:L){
    if(ram.A[l, "parameter"]==0) A[ram.A[l, "to"], ram.A[l, "from"]] <- ram.A[l, "start"]
    else A[ram.A[l, "to"], ram.A[l, "from"]] <- theta[ram.A[l, "parameter"]]
    }

L<- dim(ram.P)[1]
for(l in 1:L){
    if(ram.P[l, "parameter"]==0) P[ram.P[l, "to"], ram.P[l, "from"]] <- ram.P[l, "start"]
    else P[ram.P[l, "to"], ram.P[l, "from"]] <- theta[ram.P[l, "parameter"]]
    }
P <- P + t(P) - diag(diag(P))

J<- matrix(0, n, m)
J[cbind(1:n, 1:n)]<-1
IA.inv<- solve(diag(m)-A)
Sigma.theta <- J %*% IA.inv %*% P %*% t(IA.inv) %*% t(J)

rownames(Sigma.theta)<- colnames(Sigma.theta) <- obs.vars 

result$P <- P
result$A <- A
result$Sigma.theta <- Sigma.theta

return(result)
} # end of theta.2.Sigma <- function()

