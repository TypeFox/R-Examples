#######################################################################
## R Code to estimate CHOPIT model of political efficacy
## Example: default model
## 
## Author:    Jonathan Wand <wand(at)stanford.edu>
## Created:   2002-08-01
##
## Modified :  2008-05-09 : JW
## - anchors 3.0 syntax
##
#######################################################################

cat("\n\n annchor(...,method='chopit') Demo\n\n")

## Step 1. Get library and data
# library(anchors)
data(mexchn)

cat("\nSpecify list of model components\n")
## Step 2. List of names of columns from dataset which will be used in analysis
fo <- list(self =  xsayself ~ china + age + male + educyrs  ,
           vign = cbind(xsay1,xsay2,xsay3,xsay4,xsay5) ~ 1  ,
           tau  =           ~ china + age + male + educyrs  )

cat("\nDefault invocation of chopit\n")
## Step 3. Invoke chopit() function
out0  <- chopit( fo, mexchn)

cat("\nSummary of default chopit\n")
summary(out0)

cat("\nTake a look at some additional information\n")
cat("Gradients\n")
print(out0$gr)
cat("Computation time\n")
print(out0$time)


