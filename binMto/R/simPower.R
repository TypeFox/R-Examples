"simPower" <-
function(H0diff, pH1, n, n.sim=1000, conf.level=0.95, alternative="two.sided", method="Add4", adj="Dunnett")
{

# include some warnings:

alternative <- match.arg(arg=alternative, choices=c("two.sided", "less", "greater"))
method <- match.arg(arg=method, choices=c("Add4", "Add2", "NHS", "Wald"))
adj <- match.arg(arg=adj, choices=c("Dunnett", "Dunnettappr", "Bonf", "Unadj"))




if(is.matrix(H0diff)==FALSE)
 {H0diff<-matrix(H0diff, ncol=length(H0diff))}

if(is.matrix(pH1)==FALSE)
 {pH1<-matrix(pH1, ncol=length(pH1))}

if(is.matrix(n)==FALSE)
 {n<-matrix(n, ncol=length(n))}


if( (ncol(H0diff)+1) != ncol(pH1) || (ncol(H0diff)+1) != ncol(n) )
 {stop("the length of vectors respective ncol of matrices specified in H0diff, pH1, n must be the same")}



n.set<-max(nrow(H0diff), nrow(pH1), nrow(n))

if(nrow(H0diff)<n.set) {H0diff <- matrix(rep(H0diff[1,], times=n.set), byrow=TRUE, nrow=n.set )}
if(nrow(pH1)<n.set) {pH1 <- matrix(rep(pH1[1,], times=n.set), byrow=TRUE, nrow=n.set )}
if(nrow(n)<n.set) {n <- matrix(rep(n[1,], times=n.set), byrow=TRUE, nrow=n.set )}

mat.set<-cbind(H0diff, pH1, n)

coverage <- numeric(length=n.set)
anypp <- numeric(length=n.set)

for(e in 1:n.set)
 {
 temp <- simPowerI(H0diff=H0diff[e,], pH1=pH1[e,], n=n[e,], n.sim=n.sim , conf.level=conf.level,  
alternative=alternative, method=method, adj=adj)
 coverage[e] <- temp$coverage
 anypp[e] <-temp$any_pair_power
 }


mat.set<-cbind(H0diff, pH1, n, anypp, coverage)

k<-ncol(H0diff)



colnames(mat.set)<-c(paste("H0:p",1:k,"-p0", sep=""),
 paste("H1:p",0:k, sep=""),
 paste("n", 0:k, sep=""),
 "power", "coverage")


return(mat.set)

}

