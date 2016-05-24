cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++                   Begin testRtlbDisplay                 +++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

source("../R/rtlbDisplay.R")


######################################################################
###                           Dispatching                          ###
###                               S4                               ###
######################################################################

setGeneric(name="r2lBiv",def=function(y,x,tabTitle,textBefore="",textAfter="",graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide"){standardGeneric("r2lBiv")})

for(i in 1:5){
    for(j in 1:5){
        cat(paste("f",i,j,"<- function(y,x,textBefore='',textAfter='',graphDir='graphBiv',
graphName='V',type='png',out='latex',displayStyle='wide')cat('f",i,"-",j,"; Y=',class(y),' X=',class(x),'Display=',displayStyle,'\\n')\n",sep=""))
    }
}

setMethod(f="r2lBiv",signature=c("logical","logical"),def=f11)#r2lBivContinuousLogical)
setMethod(f="r2lBiv",signature=c("logical","factor"),def=f12)#r2lBivContinuousFactor)
setMethod(f="r2lBiv",signature=c("logical","ordered"),def=f13)#r2lBivContinuousOrdered)
setMethod(f="r2lBiv",signature=c("logical","discrete"),def=f14)#r2lBivContinuousDiscrete)
setMethod(f="r2lBiv",signature=c("logical","continuous"),def=f15)#r2lBivContinuousContinuous)

setMethod(f="r2lBiv",signature=c("factor","logical"),def=f21)#r2lBivContinuousLogical)
setMethod(f="r2lBiv",signature=c("factor","factor"),def=f22)#r2lBivContinuousFactor)
setMethod(f="r2lBiv",signature=c("factor","ordered"),def=f23)#r2lBivContinuousOrdered)
setMethod(f="r2lBiv",signature=c("factor","discrete"),def=f24)#r2lBivContinuousDiscrete)
setMethod(f="r2lBiv",signature=c("factor","continuous"),def=f25)#r2lBivContinuousContinuous)

setMethod(f="r2lBiv",signature=c("ordered","logical"),def=f31)#r2lBivContinuousLogical)
setMethod(f="r2lBiv",signature=c("ordered","factor"),def=f32)#r2lBivContinuousFactor)
setMethod(f="r2lBiv",signature=c("ordered","ordered"),def=f33)#r2lBivContinuousOrdered)
setMethod(f="r2lBiv",signature=c("ordered","discrete"),def=f34)#r2lBivContinuousDiscrete)
setMethod(f="r2lBiv",signature=c("ordered","continuous"),def=f35)#r2lBivContinuousContinuous)

setMethod(f="r2lBiv",signature=c("discrete","logical"),def=f41)#r2lBivContinuousLogical)
setMethod(f="r2lBiv",signature=c("discrete","factor"),def=f42)#r2lBivContinuousFactor)
setMethod(f="r2lBiv",signature=c("discrete","ordered"),def=f43)#r2lBivContinuousOrdered)
setMethod(f="r2lBiv",signature=c("discrete","discrete"),def=f44)#r2lBivContinuousDiscrete)
setMethod(f="r2lBiv",signature=c("discrete","continuous"),def=f45)#r2lBivContinuousContinuous)

setMethod(f="r2lBiv",signature=c("continuous","logical"),def=f51)#r2lBivContinuousLogical)
setMethod(f="r2lBiv",signature=c("continuous","factor"),def=f52)#r2lBivContinuousFactor)
setMethod(f="r2lBiv",signature=c("continuous","ordered"),def=f53)#r2lBivContinuousOrdered)
setMethod(f="r2lBiv",signature=c("continuous","discrete"),def=f54)#r2lBivContinuousDiscrete)
setMethod(f="r2lBiv",signature=c("continuous","continuous"),def=f55)#r2lBivContinuousContinuous)




####################
#   Latex output   #
####################

ds <- list(lim=4,wide=c("f1","o2"),long=c("f2","y2"))
ds2 <- list(wide=c("f1","o2"),long=c("f2","y2"))

### Logical
r2lBiv(f1,f1)
r2lBiv(f1,f2)
r2lBiv(f1,f3)

r2lBiv(f1,o1)
r2lBiv(f1,o2)
r2lBiv(f1,o3)

r2lBiv(f1,d1)
r2lBiv(f1,d2)
r2lBiv(f1,d3)

r2lBiv(f1,c1)
r2lBiv(f1,c2)
r2lBiv(f1,c3)

### factor
r2lBiv(f2,f1)
r2lBiv(f2,f2)
r2lBiv(f2,f3)

r2lBiv(f2,o1)
r2lBiv(f2,o2)
r2lBiv(f2,o3)

r2lBiv(f2,d1)
r2lBiv(f2,d2)
r2lBiv(f2,d3)

r2lBiv(f2,c1)
r2lBiv(f2,c2)
r2lBiv(f2,c3)




cat("---------------------------------------------------------------
---                    End testRtlbDisplay                  ---
---------------------------------------------------------------\n")

