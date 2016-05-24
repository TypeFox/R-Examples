library(hmmm)

data(kentucky)

y<-getnames(kentucky,st=4)

names<-c("Injury","Restraint","Year")

# variable 1: Injury status 
# variable 2: Restraint used 
# variable 3: Year 

# the lower the variable number is the faster the variable sub-script changes in the vectorized table

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# HMM MODELS with INEQUALITY CONTRAINTS on second and third order interactions
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#================================================
# Hypothesis H_1: 
# the distribution of Injury status is stochastically greater (Simple stochastic
# order) for non-user of Restraint than for user in every Year and the strenght 
# of this positive association increases over years
#================================================

marglist<-marg.list(c("m-m-l","m-l-l","g-l-l"),mflag="m")

# interactions involved in the inequality constraints
inelist<-list(marg=c(1,2,3),int=list(c(1,2),c(1,2,3)),types=c("g","l","l"))


# definition of the model with inequalities 

model<-hmmm.model(marg=marglist,lev=c(5,2,5),dismarg=list(inelist),names=names)

# definition of the model with inequalities turned into equalities, 

sel<-c(14:17,34:49) # positions of the zero-constrained interactions
modeli<-hmmm.model(marg=marglist,lev=c(5,2,5),dismarg=list(inelist),sel=sel,names=names)

# estimation of the models

# Saturated model
modnoineq<-hmmm.mlfit(y,model,noineq=TRUE)

# model with non-negative global-local log-odds ratios of Injury and Restraint
# and non-negative third-order interactions (Injury and Restraint over Years)
modineq<-hmmm.mlfit(y,model,noineq=FALSE)

# model with inequalities turned into equalities,
# no EQUALITY contraints on the univariate marginal logits
modin<-hmmm.mlfit(y,modeli)


# -------------------
#  hypotheses tested
# -------------------

# NB: testA --> H0=(modin model) vs H1=(modineq model)
#     testB --> H0=(modineq model) vs H1=(modnoineq model)
#               that is: [H_1 simple] vs [no ineq. model]

P<-hmmm.chibar(nullfit=modin,disfit=modineq,satfit=modnoineq)

print(P)


#====================================================================================
# Hypothesis H_2: 
# Hypothesis H_1 tested simultaneously with the hypothesis that the distribution of 
# Injury status becomes stochastically greater (Simple stochastic order) over years
# both for users and non-user of Restraint 
#====================================================================================

marglist<-marg.list(c("m-m-l","m-l-l","g-l-l"),mflag="m")

# interactions involved in the inequality constraints
inelist<-list(marg=c(1,2,3),int=list(c(1,2),c(1,3),c(1,2,3)),types=c("g","l","l"))

# definition of the model with inequalities, 

model<-hmmm.model(marg=marglist,lev=c(5,2,5),dismarg=list(inelist),names=names)

# definition of the model with inequalities turned into equalities, 

sel<-c(14:17,34:49) # positions of the zero-constrained interactions
modeli<-hmmm.model(marg=marglist,lev=c(5,2,5),dismarg=list(inelist),sel=sel,names=names)

# estimation of the models

# Saturated model
modnoineq<-hmmm.mlfit(y,model,noineq=TRUE)

# model with non-negative global-local log-odds ratios of Injury and Restraint;
# non-negative third-order interactions (Injury and Restraint over Years);
# non-negative global-local log-odds ratios of Injury and Year (both in user and non-user population)
modineq<-hmmm.mlfit(y,model,noineq=FALSE)

# model with inequalities turned into equalities,
# no EQUALITY contraints on the univariate marginal logits
modin<-hmmm.mlfit(y,modeli)


# -------------------
#  hypotheses tested
# -------------------

# NB: testA --> H0=(modin model) vs H1=(modineq model)
#     testB --> H0=(modineq model) vs H1=(modnoineq model)
#               that is: [H_1 & H_2 simple] vs [no ineq. model]

p<-hmmm.chibar(nullfit=modin,disfit=modineq,satfit=modnoineq)

print(p)


