# Get the CEMS data and generate design matrix
example(wide2long, package="ordBTL", echo=FALSE)
des1 <- design(CEMSlong, var1="object1", var2="object2", 
               use.vars="Y", reference="Stockholm")

# Fit the adjacent categories model, which corresponds to 
# the log-linear BTL model (see Agresti, 1992)
mod1 <- ordBTL(Y~., data=des1, family="acat", 
               family.control=list(reverse=TRUE))

# Extract all parameter estimates
getRank(mod1)

# Extract all parameter estimates and add parameter for
# reference object (which is set to zero)
getRank(mod1, reference="GAMMA.Stockholm")

# Extract only parameter estimates that include the 
# string "Intercept"
getRank(mod1, prefix="Intercept")

# Extract only parameter estimates that include the 
# string "GAMMA" (which will be the object parameters)
getRank(mod1, prefix="GAMMA")
