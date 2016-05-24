
### Regression summaries - multiple formats:


book4<-XLwriteOpen("attenu.xls") 

quakenames=c("Magnitude (Richter), per unit","Distance (log km), per x10")

# Ground acceleration as a function of magnitude and distance, all log scale.
quakemod1=summary(lm(log10(accel)~mag+log10(dist),data=attenu))


## Model-scale summaries; we don't care for the intercept.

XLregresSummary(book4,"ModelScale",varnames=quakenames,
                betas=quakemod1$coef[-1,1],SE=quakemod1$coef[-1,2],
                ,title="Log-Ground Acceleration Effects",
                confac=qt(0.975,179),pfun=function(x) 2*pt(-abs(x),df=179))

## Same thing, but using matrix input; no need to provide SE and names. 
## It is arguably still nicer to provide your own names - but could be a reproducibility risk. 
## Also, increasing the p-value resolution by changing 'pround'.

XLregresSummary(book4,"ModelScale",betas=quakemod1$coef[-1,],
                pround=6,title="Log-Ground Acceleration Effects",
                confac=qt(0.975,179),pfun=function(x) 2*pt(-abs(x),df=179),row1=8)

## Effects are arguably more meaningful as percent change. 
## So... still same model, but different summaries. 
## Also, note the combination of matrix input with names over-written via 'varnames':

XLregresSummary(book4,"PercentChange",varnames=quakenames,
                betas=quakemod1$coef[-1,],
                roundig=1,pround=6,title="Relative Ground Acceleration Effects",
                transfun=function(x) 100*(10^x-1),
                effname="Percent Change",confac=qt(0.975,179),pfun=function(x) 2*pt(-abs(x),df=179))

cat("Look for",paste(getwd(),"attenu.xls",sep='/'),"to see the results!\n")

### lm() does not take account of station or event level grouping.
### So we use a mixed model, losing 16 data points w/no station data:
### Run this on your own... and ask the authors of "lme4" about p-values at your own risk :)

# library(lme4)
# quakemod2=lmer(log10(accel)~mag+log10(dist)+(1|event)+(1|station),data=attenu)
# 
# XLregresSummary(book4,"MixedModel",varnames=quakenames,betas=fixef(quakemod2)[-1],
# SE=sqrt(diag(vcov(quakemod2)))[-1],
# roundig=1,pround=6,
# title="Relative Ground Acceleration Effects",
# transfun=function(x) 100*(10^x-1),effname="Percent Change",
# confac=qt(0.975,160),pfun=function(x) 2*pt(-abs(x),df=160))

