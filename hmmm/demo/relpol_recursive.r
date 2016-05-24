library(hmmm)

data(relpol)
  
y<-getnames(relpol,st=12,sep=";")

names<-c("Religion","Politics")


# variable 1: Religion
# variable 2: Politics

#the lower the variable number is the faster the variable sub-script changes in the vectorized table

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# MODELS with EQUALITY CONSTRAINTS
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#=====================================================
# MODEL of STOCHASTIC INDEPENDENCE (Recursive approach)
#=====================================================

# univariate and bivariate marginals
marginals<-marg.list(c("r-m","m-r","r-r"),mflag="m")

# NB: in this case the univariate marginals can be of any type
# as they are not constrained. 

R1<-matrix(c(-1,-1,1,
             -1,1,0),2,3,byrow=T)
# logits of recursive (or nested) type for variable 1 RELIGION:
# log[p(N)/p(P,C)]; log[p(C)/p(P)];
# Note the comparison between religious and non-religious citizens

R2<-matrix(c(-1,-1,-1, 0, 1, 1, 1,
             -1,-1,-1, 1, 0, 0, 0,
              1,-1,-1, 0, 0, 0, 0,
              0,-1, 1, 0, 0, 0, 0,
              0, 0, 0, 0,-1,-1, 1,
              0, 0, 0, 0,-1, 1, 0),6,7,byrow=T)           
# logits of recursive (or nested) type for variable 2 POLITICS:
# log[p(SC,C,EC)/p(EL,L,SL)]; log[p(M)/p(EL,L,SL)]; log[p(EL)/p(L,SL)];
# log[p(SL)/p(L)]; log[p(EC)/p(SC,C)]; log[p(C)/p(SC)];
# Note the aggregation of categories following an ideological similarity
# criterion: 'Liberals' and 'Conservatives'; moreover note that moderate
# orientations are compared to extreme attitudes

rec<-recursive(R1,R2)

# =====
# if the approach is, for example, Baseline-Recursive instead of Recursive:
#
# marginals<-marg.list(c("b-m","m-r","b-r"),mflag="m")
#
# R2<-matrix(c(-1,-1,-1, 0, 1, 1, 1,
#              -1,-1,-1, 1, 0, 0, 0,
#               1,-1,-1, 0, 0, 0, 0,
#               0,-1, 1, 0, 0, 0, 0,
#               0, 0, 0, 0,-1,-1, 1,
#               0, 0, 0, 0,-1, 1, 0),6,7,byrow=T)           
# logits of recursive (or nested) type for variable 2 POLITICS:
# log[p(SC,C,EC)/p(EL,L,SL)]; log[p(M)/p(EL,L,SL)]; log[p(EL)/p(L,SL)];
# log[p(SL)/p(L)]; log[p(EC)/p(SC,C)]; log[p(C)/p(SC)];
# 
# rec<-recursive(0,R2)
# =====


# Hypothesis of stochastic independence: all recursive-log odds ratios are null 
# NB: sel=c(9:20) --> positions of the zero-constrained interactions

# definition of the model
model<-hmmm.model(marg=marginals,lev=c(3,7),sel=c(9:20),cocacontr=rec,names=names)

print(model)

# estimation of the model
est_model<-hmmm.mlfit(y,model)

#print(est_model2,printflag=T)

summary(est_model)



#==========================================================
# MODEL 1: Religion _||_ Politics (Slightly, Normal, Extreme)
# when only Catholics and Protestants of the same political 
# orientation (Conservatives or Liberals) are considered
#                  - Recursive approach -
#==========================================================


# Hypothesis: the following four recursive-log odds ratios are null 
# log{[p(P,(L,SL))*p(C,EL)]/[p(P,EL)*p(C,(L,SL))]}; 
# log{[p(P,L)*p(C,SL)]/[p(P,SL)*p(C,L)]}; 
# log{[p(P,(SC,C))*p(C,EC)]/[p(P,EC)*p(C,(SC,C))]}; 
# log{[p(P,(SC)*p(C,C)]/[p(P,C)*p(C,SC)]}; 

sel<-c(14,16,18,20) # positions of the zero-constrained interactions

# definition of the model
model1<-hmmm.model(marg=marginals,lev=c(3,7),sel=sel,cocacontr=rec,names=names)

# estimation of the model
est_model1<-hmmm.mlfit(y,model1)

print(est_model1,printflag=T)

#summary(est_model1)



#==========================================================
# MODEL 2: Religion (religious vs non religious citizens) 
# _||_ Politics (Slightly, Normal, Extreme) when only 
# Liberals are considered
#                  - Recursive approach -
#==========================================================


# Hypothesis: the following two recursive-log odds ratios are null 
# log{[p((P,C),(L,SL))*p(N,EL)]/[p((P,C),EL)*p(N,(L,SL))]}; 
# log{[p((P,C),L)*p(N,SL)]/[p((P,C),SL)*p(N,L)]}; 

# sel=c(13,15) --> positions of the zero-constrained interactions

# definition of the model
model2<-hmmm.model(marg=marginals,lev=c(3,7),sel=c(13,15),cocacontr=rec,names=names)

# estimation of the model
est_model2<-hmmm.mlfit(y,model2)

#print(est_model2,printflag=T)

summary(est_model2)


#==========================================================
# MODEL 3: Religion (religious vs non religious citizens) 
# _||_ Politics (Slightly, Normal, Extreme) when only 
# Conservatives are considered
#                  - Recursive approach -
#==========================================================


# Hypothesis: the following two recursive-log odds ratios are null 
# log{[p((P,C),(SC,C))*p(N,EC)]/[p((P,C),EC)*p(N,(SC,C))]}; 
# log{[p((P,C),SC)*p(N,C)]/[p((P,C),C)*p(N,SC)]}; 

sel<-c(17,19) # positions of the zero-constrained interactions

# definition of the model
model3<-hmmm.model(marg=marginals,lev=c(3,7),sel=sel,cocacontr=rec,names=names)

# estimation of the model
est_model3<-hmmm.mlfit(y,model3)

print(est_model3,printflag=T)

#summary(est_model3)


#==========================================================
# MODEL 4: intersection of the hypotheses of Models 1, 2, 3
#                  - Recursive approach -
#==========================================================


# Hypothesis: the following eight recursive-log odds ratios are null 
# log{[p(P,(L,SL))*p(C,EL)]/[p(P,EL)*p(C,(L,SL))]}; 
# log{[p(P,L)*p(C,SL)]/[p(P,SL)*p(C,L)]}; 
# log{[p(P,(SC,C))*p(C,EC)]/[p(P,EC)*p(C,(SC,C))]}; 
# log{[p(P,(SC)*p(C,C)]/[p(P,C)*p(C,SC)]}; 
# log{[p((P,C),(L,SL))*p(N,EL)]/[p((P,C),EL)*p(N,(L,SL))]}; 
# log{[p((P,C),L)*p(N,SL)]/[p((P,C),SL)*p(N,L)]}; 
# log{[p((P,C),(SC,C))*p(N,EC)]/[p((P,C),EC)*p(N,(SC,C))]}; 
# log{[p((P,C),SC)*p(N,C)]/[p((P,C),C)*p(N,SC)]}; 

sel<-c(13,14,15,16,17,18,19,20) # positions of the zero-constrained interactions

# definition of the model
model4<-hmmm.model(marg=marginals,lev=c(3,7),sel=sel,cocacontr=rec,names=names)

# estimation of the model
est_model4<-hmmm.mlfit(y,model4)

#print(est_model4,printflag=T)

summary(est_model4)


