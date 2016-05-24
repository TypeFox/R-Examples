INCAtest <- function(d, pert, d_test, np=1000, alpha=0.05, P=1){
#################### Application of INCA test #################################
# Input:
# d: distance matrix between individuals (nxn) 
# pert: integer vector indicating the group each individual belongs to.
# d_test: vector of length n with distancies from the specific individual to the #         individuals of different groups.
# np: bootstrap sample size for null distribution of W
# alpha: fixed level of the test
# P: the procedure is repeted 10*P times
#
# output: value of INCA statistic (W), pvalue and values of proj. (U1,..., Uk)
################################################################################

d <- as.matrix(d)
n<-dim(d)[1]
pert <- as.integer(pert)
k<-max(pert)
# populations must be named with numbers from 1 to k
if (length(tabulate(as.factor(pert))) != k)
  stop("Partitions must be named by factors or with numbers from 1 to k.")
# 0 can not be a partitions name
if (any(pert==0))
  stop("pert contains 0 named individuals.Partitions must be named by factors or with numbers from 1 to k.")

# Calculate  W and projections U
vg <- vgeo(d, pert);
delta <- deltas_simple(d,vg,pert)
phi <- proxi_simple(d_test,vg,pert)
ama <- estW_simple(vg,delta, phi)
W0 <-  ama$Wvalue
U <- ama$Uvalue

###### Calculate the null distribution. It is repeated 10*P times
frec <- tabulate(pert)  # vector of frecuencies of individuals in each population


# Acumulate frecuencies
frecacum <- matrix(0,k+1,1);
for (pob in 2:(k+1)){
     frecacum[pob] <- frecacum[pob-1]+frec[pob-1];
}

P <- 10*P
percentage  <- 0
for (r in 1:P){
	TW <- distrW(d, pert, np, frec, frecacum)  
	partial.p.value <- 0
	for (b in 1:np){
		if (TW[b] > W0){partial.p.value <- partial.p.value +1}
	}
	partial.p.value <- partial.p.value/np
	if (partial.p.value < alpha){percentage <- percentage + 1}
}

percentage  <- percentage /P*100

out <- list(StatisticW0= W0, ProjectionsU=U, Percentage_under_alpha=percentage, alpha=alpha )
class(out) <- "incat"
return(out)
}
