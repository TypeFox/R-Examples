######################################################################
###                           Continuous                           ###
######################################################################

#for(i in c("Logical","Factor","Ordered","Discrete","Continuous")){
#    for(j in c("Logical","Factor","Ordered","Discrete","Continuous")){
#        cat('
#r2lBiv',i,j,' <- function(y,x,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
#    if(displayStyle=="long"){
#        r2lBiv',i,j,'Wide(y=y,x=x,graphDir=graphDir,graphName=graphName,type=type,out=out)
#    }else{
#        r2lBiv',i,j,'Long(y=y,x=x,graphDir=graphDir,graphName=graphName,type=type,out=out)
#    }
#}
#
#',sep="")
#    }
#}

#r2lBivLogicalLogical <- function(y,x,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide")








#r2lBivFactorLogical <- function(y,x,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide")


r2lBivFactorFactor <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivFactorFactorWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivFactorFactorLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


r2lBivFactorOrdered <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivFactorOrderedWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivFactorOrderedLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


r2lBivFactorDiscrete <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivFactorDiscreteWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivFactorDiscreteLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


r2lBivFactorContinuous <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivFactorContinuousWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivFactorContinuousLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


#r2lBivOrderedLogical <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide")

r2lBivOrderedFactor <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivOrderedFactorWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivOrderedFactorLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


r2lBivOrderedOrdered <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivOrderedOrderedWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivOrderedOrderedLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


r2lBivOrderedDiscrete <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivOrderedDiscreteWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivOrderedDiscreteLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


r2lBivOrderedContinuous <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivOrderedContinuousWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivOrderedContinuousLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


#r2lBivDiscreteLogical <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide")

r2lBivDiscreteFactor <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivDiscreteFactorWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivDiscreteFactorLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


r2lBivDiscreteOrdered <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivDiscreteOrderedWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivDiscreteOrderedLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


r2lBivDiscreteDiscrete <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivDiscreteDiscreteWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivDiscreteDiscreteLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


r2lBivDiscreteContinuous <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivDiscreteContinuousWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivDiscreteContinuousLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


#r2lBivContinuousLogical <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {

r2lBivContinuousFactor <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivContinuousFactorWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivContinuousFactorLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}

r2lBivContinuousOrdered <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivContinuousOrderedWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivContinuousOrderedLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}

r2lBivContinuousDiscrete <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivContinuousDiscreteWide(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivContinuousDiscreteLong(y=y,x=x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}

#r2lBivContinuousContinuous <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {


