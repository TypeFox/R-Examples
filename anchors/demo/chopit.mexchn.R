#######################################################################
## R Code to estimate CHOPIT model of political efficacy
##
## Example: estimation with non-linear taus
## 
## Author:    Jonathan Wand <wand(at)stanford.edu>
## Created:   2002-08-01
##
## Modified
## - 2008-05-09 : JW
##   anchors 3.0 syntax
##
#######################################################################
cat("Replication of Table 2 in King et al 2004\n")

## Step 1. Get data
data(mexchn)

## Step 2. List of names of columns from dataset which will be used in analysis
fo <- list(self =  xsayself ~ china + age + male + educyrs,
           vign = cbind(xsay1,xsay2,xsay3,xsay4,xsay5) ~ 1         ,
           tau  =           ~ china + age + male + educyrs )

## Step 3. Run model
out1  <- chopit( fo, mexchn,
                 options=anchors.options(linear=FALSE,vign.var="homo"))
summary(out1)
