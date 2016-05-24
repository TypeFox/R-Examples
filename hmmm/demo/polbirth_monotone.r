library(hmmm)

data(polbirth)
       
y<-getnames(polbirth,st=12,sep=";")

names<-c("Politics","Birth")

# variable 1: Politics
# variable 2: Birthcontrol

#the lower the variable number is the faster the variable sub-script changes in the vectorized table

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# MODELS with EQUALITY and INEQUALITY CONSTRAINTS
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# ----------------------------------------------------------------------
# Hypothesis of Uniform monotone dependence of Politics on Birthcontrol:
#          -- non-negative continuation-local log-odds ratios --
# ----------------------------------------------------------------------

# univariate and bivariate marginals
marginals<-marg.list(c("c-m","m-l","c-l"),mflag="m")

# NB: in this case the univariate marginals can be of any type
# as they are not constrained. 

# all the counts come from the same and only stratum
Z<-ZF<-matrix(1,28,1)
            
# interactions involved in the inequality constraints 
dism<-list(marg=c(1,2),int=list(c(1,2)),types=c("c","l"))

# definition of the model with inequalities, no EQUALITY constraints on the univariate marginal logits
model<-hmmm.model(marg=marginals,dismarg=list(dism),lev=c(7,4),names=names)

# estimation of the models

# SATURATED model
msat<-hmmm.mlfit(y,model,noineq=TRUE)

# model with INEQUALITIES on the continuation-local log-odds ratios: "Uniform monotone dependence model"
# --> model with non-negative continuation-local log-odds ratios
mu <- hmmm.mlfit(y,model,noineq=FALSE)

# model with INEQUALITIES turned into EQUALITIES, no EQUALITY constraints
# on the univariate marginal logits: "Stochastic independence model"
# --> model with null continuation-local log-odds ratios
# sel=c(10:27) --> positions of the zero-constrained interactions

model0<-hmmm.model(marg=marginals,lev=c(7,4),sel=c(10:27))
mnull <- hmmm.mlfit(y,model0)

# HYPOTHESES TESTED:
# NB: testA --> H0=(mnull model) vs H1=(mu model)
#     testB --> H0=(mu model) vs H1=(msat model)

PP<-hmmm.chibar(nullfit=mnull,disfit=mu,satfit=msat)

summary(PP)


# ---------------------------------------------------------------------
# Hypothesis of Simple monotone dependence of Politics on Birthcontrol:
#           -- non-negative global-local log-odds ratios --
# ---------------------------------------------------------------------

# univariate and bivariate marginals
marginals<-marg.list(c("g-m","m-l","g-l"),mflag="m")

# NB: in this case the univariate marginals can be of any type
# as they are not constrained. 

# all the counts come from the same and only stratum
Z<-ZF<-matrix(1,28,1)

# interactions involved in the inequality constraints 
dism<-list(marg=c(1,2),int=list(c(1,2)),types=c("g","l"))

# definition of the model with inequalities, no EQUALITY constraints on the univariate marginal logits
model<-hmmm.model(marg=marginals,dismarg=list(dism),lev=c(7,4),names=names)

# estimation of the models

# SATURATED model
msat<-hmmm.mlfit(y,model,noineq=TRUE)

# model with INEQUALITIES on the global-local log-odds ratios: "Simple monotone dependence model"
# --> model with non-negative global-local log-odds ratios
ms <- hmmm.mlfit(y,model,noineq=FALSE)

# model with INEQUALITIES turned into EQUALITIES, no EQUALITY constraints
# on the univariate marginal logits: "Stochastic independence model"
# --> model with null global-local log-odds ratios
# sel=c(10:27) --> positions of the zero-constrained interactions

model0<-hmmm.model(marg=marginals,lev=c(7,4),sel=c(10:27))

mnull <- hmmm.mlfit(y,model0)

# HYPOTHESES TESTED:
# NB: testA --> H0=(mnull model) vs H1=(ms model)
#     testB --> H0=(ms model) vs H1=(msat model)

P<-hmmm.chibar(nullfit=mnull,disfit=ms,satfit=msat)

summary(P)

