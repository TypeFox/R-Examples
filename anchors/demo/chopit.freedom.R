#######################################################################
##
## Author:    Jonathan Wand <wand(at)stanford.edu>
## Created:   2007-02-01
##
## Modified :  2008-05-09 : JW
## - anchors 3.0 syntax
##
#######################################################################
cat("Repl Wand et al (2007) Table 4 (chopit, linear)\n")
data(freedom)

cat("\nSpecify list of model components\n")
## Step 2. List of names of columns from dataset which will be used in analysis
fo <- list(self =  self ~ sex + age + educ + factor(country)  ,
           vign = cbind(vign1,vign2,vign3,vign4,vign5,vign6) ~ 1  ,
           tau  =       ~ sex + age + educ  + factor(country) )

cat("\nDefault invocation of chopit\n")
## Step 3. Invoke chopit() function
out0  <- chopit( fo, freedom, options=anchors.options(verbose=TRUE))
summary(out0)
