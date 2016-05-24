# Get the CEMS data and generate design matrix
example(wide2long, package="ordBTL", echo=FALSE)

des2 <- design(CEMSlong[-which(is.na(CEMSlong$Y)),], 
               var1="object1", var2="object2", 
               use.vars="ALL", reference="Stockholm")

# Formula for full model considering all subject-object interactions
form2 <- Y ~ 
  (GAMMA.London+GAMMA.Paris+GAMMA.Milano+GAMMA.StGallen+GAMMA.Barcelona)+
  (GAMMA.London+GAMMA.Paris+GAMMA.Milano+GAMMA.StGallen+GAMMA.Barcelona):
  (WOR+SEX+DEG+STUD+ENG+FRA+SPA+ITA)

## Not run: 
# Exemplatory boosting call with mstop=5
#BoostDev <- BTLboost(form2, data=des2, groupVars=c("WOR","DEG","SEX","STUD"), 
#                     selection="DEVIANCE", mstop=5)
## End(Not run)
