f.post.poo.diff <- function(coeff, covar){
##
## FOR THE POO EFFECTS, SUBTRACT MATERNALLY INHERITED FROM PATERNALLY INH.
## TO TEST FOR THE ACTUAL POO EFFECT
##
## MERK: BLIR DETTE RETT HVIS F.EKS. reference = "population"?
##
.names <- rownames(coeff[[1]])
#
## FIND NAMES/POSITIONS OF RELEVANT PARAMETERS. NOTE: \\< INSISTS ON START OF WORD, SO THAT, FOR INSTANCE, "cm1" ISN'T PICKED UP BY "m"
.tmp <- f.coefnames(.names)
.cm <- .tmp$child.poo.m
.cf <- .tmp$child.poo.f
#
##
.D0 <- diag(length(coeff[[1]]))
dimnames(.D0) <- list(.names, .names)
#
.Dm <- .D0[.cm, , drop = F]
.Df <- .D0[.cf, , drop = F]
.Dm[, .cf] <- -.Df[, .cf]
# rownames(.Dm) <- paste(.cm, .cf, sep = "_")
rownames(.Dm) <- paste("cm_", .cf, sep = "") # avoid number within
#
.D <- rbind(.D0, .Dm)
#
## RECOMPUTE COEFFICIENT AND COVARIANCE MATRICES
.coef <- lapply(coeff, function(x) .D %*% x)
.cov <- lapply(covar, function(x) .D %*% x %*% t(.D))
#
## 
.vis <- F
if(.vis){
	f.vis(.D, vis = .vis)
	f.vis(coeff[[1]], vis = .vis)
	f.vis(.D %*% coeff[[1]], vis = .vis)
	f.vis(.coef, vis = .vis)
	f.vis(.cov, vis = .vis)
}#
#
##
return(list(coeff = .coef, covar = .cov))
}
