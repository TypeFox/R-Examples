cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++                    Begin testDiscrete                   +++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

source("../R/discrete.R")


####################
#   Latex output   #
####################

# Type Postscript

sink(file="bivDisLogical.tex")
r2lBivDiscreteLogical(d2,f1,tabTitle="title",graphDir="graphBiv",graphName="V1",out="latex")
sink()

sink(file="bivDisFactor.tex")
r2lBivDiscreteFactorWide(d2,f1,tabTitle="title",graphDir="graphBiv",graphName="V3a",out="latex")
r2lBivDiscreteFactorMixed(d2,f2,tabTitle="title",graphDir="graphBiv",graphName="V3b",out="latex")
r2lBivDiscreteFactorLong(d2,f3,tabTitle="title",graphDir="graphBiv",graphName="V3c",out="latex")
sink()

#r2lBivDiscreteFactor(d2,f2,tabTitle="title",graphDir="graphBiv",graphName="V3",out="latex",displayStyle=3)
#r2lBivDiscreteFactor(d2,f3,tabTitle="title",graphDir="graphBiv",graphName="V3",out="latex",displayStyle=10)


sink(file="bivDisOrdered.tex")
r2lBivDiscreteOrderedWide(d2,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",out="latex")
r2lBivDiscreteOrderedMixed(d3,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",out="latex")
r2lBivDiscreteOrderedLong(d3,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",out="latex")
sink()

#r2lBivDiscreteOrdered(y2,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",out="latex")
#r2lBivDiscreteOrdered(y2,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",out="latex")
#r2lBivDiscreteOrdered(y2,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",out="latex")


sink(file="bivDisDiscrete.tex")
r2lBivDiscreteDiscreteWide(d2,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",out="latex")
r2lBivDiscreteDiscreteMixed(d3,d2,tabTitle="title",graphDir="graphBiv",graphName="V8",out="latex")
r2lBivDiscreteDiscreteLong(d3,d2,tabTitle="title",graphDir="graphBiv",graphName="V9",out="latex")
sink()

#r2lBivDiscreteDiscrete(d2,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",out="latex")
#r2lBivDiscreteDiscrete(y2,d2,tabTitle="title",graphDir="graphBiv",graphName="V8",out="latex")
#r2lBivDiscreteDiscrete(y3,d3,tabTitle="title",graphDir="graphBiv",graphName="V9",out="latex")

sink(file="bivDisContinuous.tex")
r2lBivDiscreteContinuousWide(d2,c1,tabTitle="title",graphDir="graphBiv",graphName="V10",out="latex")
r2lBivDiscreteContinuousMixed(d2,c2,tabTitle="title",graphDir="graphBiv",graphName="V10",out="latex")
r2lBivDiscreteContinuousLong(d2,c3,tabTitle="title",graphDir="graphBiv",graphName="V10",out="latex")
sink()


###################
#   Html output   #
###################

# Type PNG
sink(file="bivDisLogical.html")
r2lBivDiscreteLogical(d2,f1,tabTitle="title",graphDir="graphBiv",graphName="V1",type="png",out="html")
sink()

sink(file="bivDisFactor.html")
r2lBivDiscreteFactorWide(d2,f2,tabTitle="title",graphDir="graphBiv",graphName="V2",type="png",out="html")
r2lBivDiscreteFactorLong(d2,f3,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html")
sink()


#r2lBivDiscreteFactor(d2,f1,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html")
#r2lBivDiscreteFactor(d2,f2,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html")
#r2lBivDiscreteFactor(d2,f3,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html")

#r2lBivDiscreteFactor(d2,f2,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html",displayStyle=3)
#r2lBivDiscreteFactor(d2,f3,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html",displayStyle=10)

sink(file="bivDisOrdered.html")
r2lBivDiscreteOrderedWide(d3,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",type="png",out="html")
r2lBivDiscreteOrderedMixed(d3,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",type="png",out="html")
r2lBivDiscreteOrderedLong(d3,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",type="png",out="html")
sink()

#r2lBivDiscreteOrdered(y2,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",type="png",out="html")
#r2lBivDiscreteOrdered(y2,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",type="png",out="html")
#r2lBivDiscreteOrdered(y2,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",type="png",out="html")

sink(file="bivDisDiscrete.html")
r2lBivDiscreteDiscreteWide(d2,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",type="png",out="html")
r2lBivDiscreteDiscreteMixed(d2,d3,tabTitle="title",graphDir="graphBiv",graphName="V8",type="png",out="html")
r2lBivDiscreteDiscreteLong(d3,d2,tabTitle="title",graphDir="graphBiv",graphName="V9",type="png",out="html")
sink()

#r2lBivDiscreteDiscrete(d2,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",type="png",out="html")
#r2lBivDiscreteDiscrete(y2,d2,tabTitle="title",graphDir="graphBiv",graphName="V8",type="png",out="html")
#r2lBivDiscreteDiscrete(y3,d3,tabTitle="title",graphDir="graphBiv",graphName="V9",type="png",out="html")

sink(file="bivDisContinuous.html")
r2lBivDiscreteContinuousWide(d2,c1,tabTitle="title",graphDir="graphBiv",graphName="V10",type="png",out="html")
r2lBivDiscreteContinuousMixed(d2,c2,tabTitle="title",graphDir="graphBiv",graphName="V10",type="png",out="html")
r2lBivDiscreteContinuousLong(d2,c3,tabTitle="title",graphDir="graphBiv",graphName="V10",type="png",out="html")
sink()


cat("---------------------------------------------------------------
---                       End testDiscrete                  ---
---------------------------------------------------------------\n")
