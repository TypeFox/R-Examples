############################################################
##                                                        ##
## Example 1: Adjacent categories logit model for CEMS    ##
##                                                        ##
############################################################

############################################################
# Reproduce results from Table 3 of Dittrich et al. (2001)
############################################################

# Get the CEMS data and generate design matrix
example(wide2long, package="ordBTL", echo=FALSE)
des1 <- design(CEMSlong, var1="object1", var2="object2", 
              use.vars="Y", reference="Stockholm")

# Fit the adjacent categories model, which corresponds to 
# the log-linear BTL model (see Agresti, 1992)
mod1 <- ordBTL(Y~., data=des1, family="acat", 
               family.control=list(reverse=TRUE))

# We get the same results from Table 3 of Dittrich et al (2001).
# Since Stockholm is the reference university, its estimate 
# is set to zero (due to identifiability)
getRank(mod1)

############################################################
# Reproduce results from Table 6 of Dittrich et al. (2001)
############################################################

# Generate design matrix and specify model formula
des2 <- design(CEMSlong, var1="object1", var2="object2", 
               use.vars="ALL", reference="Stockholm")
form2 <- Y~GAMMA.London + GAMMA.Paris + GAMMA.Milano + 
  GAMMA.StGallen + GAMMA.Barcelona + WOR +
  SEX + WOR:GAMMA.Paris + WOR:GAMMA.Milano +
  WOR:GAMMA.Barcelona + DEG:GAMMA.StGallen +
  STUD:GAMMA.Paris + STUD:GAMMA.StGallen +
  ENG:GAMMA.StGallen + FRA:GAMMA.London + 
  FRA:GAMMA.Paris + SPA:GAMMA.Barcelona +
  ITA:GAMMA.London + ITA:GAMMA.Milano +
  SEX:GAMMA.Milano

# Fit the adjacent categories model with symmetric 
# constraint for covariable WOR and SEX 
mod2 <- ordBTL(form2, data=des2, family="acat", 
               family.control=list(reverse=TRUE),
               restrict=c("WOR", "SEX"))

# We get the same results from Table 6 of Dittrich et al. (2001)
getRank(mod2)

# Notice that the change in sign for (Intercept), WOR and SEX
# is because we use here a different "coding".

############################################################
##                                                        ##
## Example 2: Fitting models from Agresti (1992)          ##
##                                                        ##
############################################################

# Data from Table 1 of Agresti (1992)
data(ribbon)

# design matrix 
des3 <- design(ribbon, var1="obj1", var2="obj2", use.vars="ALL")

# Note that Agresti (1992) used the constraint that the object 
# parameters sum up to 1. To get the same results, we use the model
form3 <- cbind(V1,V2,V3,V4,V5,V6,V7)~I(GAMMA.1-GAMMA.5)+
  I(GAMMA.2-GAMMA.5)+I(GAMMA.3-GAMMA.5)+I(GAMMA.4-GAMMA.5)

# Fit the adjacent categories logit model
ac <- ordBTL(form3, data=des3, family="acat", 
             family.control=list(reverse=TRUE))

# Fit the cumulative logit model
clm.logit <- ordBTL(form3, data=des3)

# Fit the cumulative probit model
clm.probit <- ordBTL(form3, data=des3,
                     family.control=list(link="probit"))

# Parameter estimates
coefs <- t(rbind("Adjacent categories logit"=coefficients(ac), 
                 "Cumulative probit"=coefficients(clm.probit),
                 "Cumulative logit"=coefficients(clm.logit)))
coefs <- rbind(coefs, "GAMMA.5"=0-colSums(coefs[4:7,]))
coefs


############################################################
##                                                        ##
## Example 3: Fitting models for Bundesliga 2005/2006     ##
##                                                        ##
############################################################
# real ranking can be obtained from:
# http://fussballdaten.sport.de/bundesliga/2006

# load data
example(design, package="ordBTL", echo=FALSE)

# Model without home advantage
des.nohome <- design(buli0506, var1="Heim", var2="Gast", 
                     use.vars="Y3", home.advantage="no", 
                     reference="GAMMA.MSV.Duisburg")
mod.nohome <- ordBTL(Y3~., data=des.nohome)
# team 'abilities' (should be approximately the ranking of the final standings)
getRank(mod.nohome, prefix="GAMMA", reference="GAMMA.MSV.Duisburg")

# Model with home advantage
des.onehome <- design(buli0506, var1="Heim", var2="Gast",
                      use.vars="Y3", home.advantage="yes", 
                      reference="GAMMA.MSV.Duisburg")
mod.onehome <- ordBTL(Y3~., data=des.onehome)
# team 'abilities'
getRank(mod.onehome, prefix="GAMMA", reference="GAMMA.MSV.Duisburg")
# home advantage
getRank(mod.onehome, prefix="ALPHA")

# Model with team-specific home advantage
des.teamhome <- design(buli0506, var1="Heim", var2="Gast",
                      use.vars="Y3", home.advantage="specific", 
                      reference="GAMMA.MSV.Duisburg")
mod.teamhome <- ordBTL(Y3~., data=des.teamhome)
# team 'abilities' (should be approximately the ranking for the away table)
getRank(mod.teamhome, prefix="GAMMA", reference="GAMMA.MSV.Duisburg")
# team-specific home advantages
getRank(mod.teamhome, prefix="ALPHA")
