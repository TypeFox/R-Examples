library(hmmm)
 
# see Section 4.2, pg. 350,
# "Multinomial-Poisson homogeneous models for contingency tables",
# Lang, J.B.
# The Annals of Statistics, (2004)
# 
# Table 2 - 1999 statistics journals citation pattern counts (n_3=225)
# Citing x Cited statistics journals: JASA, BMCS, ANNS

y <- matrix(c(104,24,65,76,146,30,50,9,166),9,1)

# population matrix: 3 strata with 3 observations each
Zmat <- kronecker(diag(3),matrix(1,3,1))

# the 3rd stratum sample size is fixed
ZFmat <- kronecker(diag(3),matrix(1,3,1))[,3]

########################################################################
# Let (i,j) be a cross-citation, where i is the citing journal and j is 
# the cited journal. Let m_ij be the expected counts of cross-citations.
# The Gini concentrations of citations for each of the journals are: 
# G_i = sum_j=1_3 (m_ij/m_i+)^2  for i=1,2,3.

 
 Gini<-function(m) {
 A<-matrix(m,3,3,byrow=TRUE)
 GNum<-rowSums(A^2)
 GDen<-rowSums(A)^2
 G<-GNum/GDen
 c(G[1],G[3])-c(G[2],G[1])
 }

####################################################################

# Example of MPH model subject to equality constraints: 

## 
## --> h = c(G1,G3)-c(G2,G1) = 0
## HYPOTHESIS: G1 = G2 = G3 
## 

mod_eq <- mphineq.fit(y,Z=Zmat,ZF=ZFmat,h.fct=Gini)

print(mod_eq)



# Example of MPH model subject to inequality constraints: 

##
## --> d = c(G1,G3)-c(G2,G1) >= 0
## HYPOTHESIS: G1 > G2, G3 > G1 
##

mod_ineq <- mphineq.fit(y,Z=Zmat,ZF=ZFmat,d.fct=Gini)


# Reference model (sat_mod):
# --> model (MPH model subject to inequality constraints) without inequalities 
# ==============
# NB sat_mod --> in this case saturated model is properly "the" saturated model
# (Gsq=0), model with no constraints  

mod_sat <-mphineq.fit(y,Z=Zmat,ZF=ZFmat)


# HYPOTHESES TESTED:
# NB: testA --> H0=(mod_eq) vs H1=(mod_ineq model)
#     testB --> H0=(mod_ineq model) vs H1=(sat_mod model)

hmmm.chibar(nullfit=mod_eq,disfit=mod_ineq,satfit=mod_sat)

