library(hmmm)

data(accident)

y<-getnames(accident,st=9,sep=";")


# RESPONSE VARIABLES:
# variable 1: Type
# variable 2: Time

# COVARIATES:
# variable 3: Age
# variable 4: Hour

#the lower the variable number is the faster the variable sub-script changes in the vectorized table

#================================================================
# Modelling the effect of covariates on local-global interactions
#================================================================

# univariate and bivariate marginals FOR THE RESPONSE VARIABLES
marginals<-marg.list(c("l-m","m-g","l-g"),mflag="m")
names<-c("Type","Time")

# --------------------------------------------------------
# unconstrained marginals and Age, Hour additive effect on
# the log-odds ratios between Type and Time

al<-list(
Type=~Type*Age*Hour,
Time=~Time*Age*Hour,
Type.Time=~Type.Time*(Age+Hour)
)
                                                              
# estimation of the models                                                 
                                                              
model1<-hmmm.model.X(marg=marginals,lev=c(3,4),names=names,Formula=al,strata=c(3,2),fnames=c("Age","Hour"))
mod1<-hmmm.mlfit(y,model1,y.eps=0.00001)
                                                  
# H0=(mod1) vs H1=(unrestricted)

summary(mod1)

# --------------------------------------------------------
# Age, Hour additive effect on the marginal logits
# of the response variables (Type, Time) and on the log-odds ratios between Type and Time

al<-list(
Type=~Type*(Age+Hour),
Time=~Time*(Age+Hour),
Type.Time=~Type.Time*(Age+Hour)
)

model2<-hmmm.model.X(marg=marginals,lev=c(3,4),names=names,Formula=al,strata=c(3,2),fnames=c("Age","Hour"))
mod2<-hmmm.mlfit(y,model2,y.eps=0.00001)

# H0=(mod2) vs H1=(unrestricted)

print(mod2)

# --------------------------------------------------------
# Age, Hour additive effect on the marginal logits,
# stochastic independence between Type and Time
# in each sub-table identified by the levels of the covariates Age and Hour

alind<-list(
Type=~1+Type*(Age+Hour),
Time=~1+Time*(Age+Hour),
Type.Time="zero"
)

modelind<-hmmm.model.X(marg=marginals,lev=c(3,4),names=names,Formula=alind,strata=c(3,2),fnames=c("Age","Hour"))
modind<-hmmm.mlfit(y,modelind,y.eps=0.00001)

# H0=(modind) vs H1=(unrestricted)

print(modind,printflag=TRUE)

