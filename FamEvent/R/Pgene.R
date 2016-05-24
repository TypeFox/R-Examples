Pgene <-
function(g, pg, a.freq=0.0014){
#P(g| pg:parents' genotype), a.freq: allele freq.
qAA <- a.freq^2
qAa <- 2*a.freq*(1-a.freq)
qaa <- (1-a.freq)^2

if(length(g)==1) g <- rep(g,length(pg))
if(length(pg) == 1) pg <- rep(pg, length(g))
re <- 0
re[g==1] <- cbind(qAA, 1, 0.5, 0, 0.5, 0.25, 0,0,0,0)[pg[g==1]+1]
re[g==2] <- cbind(qAa, 0, 0.5, 1, 0.5, 0.5, 0.5, 1, 0.5, 0)[pg[g==2]+1]
re[g==3] <- cbind(qaa, 0,0,0,0, 0.25,0.5, 0, 0.5,1)[pg[g==3]+1]
return(re)

}
