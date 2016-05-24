library(hmmm)


#=========================================================
# MODELS with INEQUALITY CONTRAINTS on the log-odds ratios
#   - Double Monotone Dependence Hypotheses -
#=========================================================


# two way contingency table --> Customer Service x Training,
# 377 customers, users of a machine tool.
# (4 categories each, levels of satisfaction regarding 'Customer Service'
# and 'Training' given by the company that sells the machine tool:
# U=unsatisfied, S=satisfied, RS=really satisfied, ES=extremely satisfied).

y <- c(
 9, 10, 22, 25,
 3,  4, 10,  4,
12, 15, 81, 37,
11,  6, 31, 97)

# variable 1: Customer Service
# variable 2: Training

#the lower the variable number is the faster the variable sub-script changes in the vectorized table


##################################################
# Double Simple monotone dependence:
# global-local & local-global log-odds ratios >= 0
##################################################


# NB: no constraints on the marginal logits 

# univariate and bivariate marginals
marginals<-marg.list(c("l-m","m-l","l-l"),mflag="m")

# marginal list involved in the INEQUALITIES
marginal12dis<-list(marg=c(1,2),int=list(c(1,2)),types=c("g","l"))
marginal12bis<-list(marg=c(1,2),int=list(c(1,2)),types=c("l","g"))
dism<-list(marginal12dis,marginal12bis)

# NB: all the non-negative constraints on the log-odds ratios
# are non-redundant constraints

# definition of the model
model<-hmmm.model(marg=marginals,dismarg=dism,lev=c(4,4))
                                                         
# estimation of the models

# SATURATED model
asat<-hmmm.mlfit(y,model)

# print(asat)

# model with INEQUALITIES on the g-l log-odds ratios and l-g log-odds ratios,
# no EQUALITY constraints on the univariate marginal logits: 
# "Double Simple monotone dependence model"
a<- hmmm.mlfit(y,model,noineq=FALSE)

# print(a)

# model with INEQUALITIES turned into EQUALITIES, no EQUALITY constraints
# on the univariate marginal logits: "Stochastic independence model";
# sel=c(7:15) --> positions of the zero-constrained interactions 
models0<-hmmm.model(marg=marginals,lev=c(4,4),sel=c(7:15))

anull<- hmmm.mlfit(y,models0)

# print(anull)

# HYPOTHESES TESTED:
# NB: testA --> H0=(anull model) vs H1=(a model)
#     testB --> H0=(a model) vs H1=(asat model)


P<-hmmm.chibar(nullfit=anull,disfit=a,satfit=asat)

print(P)






