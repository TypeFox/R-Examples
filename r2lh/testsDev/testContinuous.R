
cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++                   Begin testContinuous                  +++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

source("../R/continuous.R")



####################
#   Latex output   #
####################

# Type Postscript

sink(file="bivContLogical.tex")
r2lBivContinuousLogical(c1,f1,tabTitle="title",graphDir="graphBiv",graphName="V1",out="latex")
sink()

sink(file="bivContFactor.tex")
r2lBivContinuousFactorWide(c1,f1,tabTitle="title",graphDir="graphBiv",graphName="V3a",out="latex")
#r2lBivContinuousFactorMixed(c1,f2,tabTitle="title",graphDir="graphBiv",graphName="V3b",out="latex")
r2lBivContinuousFactorLong(c1,f3,tabTitle="title",graphDir="graphBiv",graphName="V3c",out="latex")
sink()

#r2lBivContinuousFactor(c1,f2,tabTitle="title",graphDir="graphBiv",graphName="V3",out="latex",displayStyle=3)
#r2lBivContinuousFactor(c1,f3,tabTitle="title",graphDir="graphBiv",graphName="V3",out="latex",displayStyle=10)


sink(file="bivContOrdered.tex")
r2lBivContinuousOrderedWide(c2,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",out="latex")
#r2lBivContinuousOrderedMixed(c2,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",out="latex")
r2lBivContinuousOrderedLong(c2,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",out="latex")
sink()

#r2lBivContinuousOrdered(c2,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",out="latex")
#r2lBivContinuousOrdered(c2,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",out="latex")
#r2lBivContinuousOrdered(c2,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",out="latex")


sink(file="bivContDiscrete.tex")
r2lBivContinuousDiscreteWide(c1,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",out="latex")
#r2lBivContinuousDiscreteMixed(c2,d2,tabTitle="title",graphDir="graphBiv",graphName="V8",out="latex")
r2lBivContinuousDiscreteLong(c3,d3,tabTitle="title",graphDir="graphBiv",graphName="V9",out="latex")
sink()

#r2lBivContinuousDiscrete(c1,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",out="latex")
#r2lBivContinuousDiscrete(c2,d2,tabTitle="title",graphDir="graphBiv",graphName="V8",out="latex")
#r2lBivContinuousDiscrete(c3,d3,tabTitle="title",graphDir="graphBiv",graphName="V9",out="latex")

sink(file="bivContContinuous.tex")
r2lBivContinuousContinuous(c1,c2,tabTitle="title",graphDir="graphBiv",graphName="V10",out="latex")
sink()


###################
#   Html output   #
###################

# Type PNG
sink(file="bivContLogical.html")
r2lBivContinuousLogical(c1,f1,tabTitle="title",graphDir="graphBiv",graphName="V1",type="png",out="html")
sink()

sink(file="bivContFactor.html")
r2lBivContinuousFactorWide(c1,f2,tabTitle="title",graphDir="graphBiv",graphName="V2",type="png",out="html")
r2lBivContinuousFactorLong(c1,f3,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html")
sink()


#r2lBivContinuousFactor(c1,f1,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html")
#r2lBivContinuousFactor(c1,f2,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html")
#r2lBivContinuousFactor(c1,f3,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html")

#r2lBivContinuousFactor(c1,f2,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html",displayStyle=3)
#r2lBivContinuousFactor(c1,f3,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html",displayStyle=10)

sink(file="bivContOrdered.html")
r2lBivContinuousOrderedWide(c2,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",type="png",out="html")
#r2lBivContinuousOrderedMixed(c2,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",type="png",out="html")
r2lBivContinuousOrderedLong(c2,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",type="png",out="html")
sink()

#r2lBivContinuousOrdered(c2,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",type="png",out="html")
#r2lBivContinuousOrdered(c2,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",type="png",out="html")
#r2lBivContinuousOrdered(c2,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",type="png",out="html")

sink(file="bivContDiscrete.html")
r2lBivContinuousDiscreteWide(c1,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",type="png",out="html")
#r2lBivContinuousDiscreteMixed(c2,d2,tabTitle="title",graphDir="graphBiv",graphName="V8",type="png",out="html")
r2lBivContinuousDiscreteLong(c3,d3,tabTitle="title",graphDir="graphBiv",graphName="V9",type="png",out="html")
sink()

#r2lBivContinuousDiscrete(c1,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",type="png",out="html")
#r2lBivContinuousDiscrete(c2,d2,tabTitle="title",graphDir="graphBiv",graphName="V8",type="png",out="html")
#r2lBivContinuousDiscrete(c3,d3,tabTitle="title",graphDir="graphBiv",graphName="V9",type="png",out="html")

sink(file="bivContContinuous.html")
r2lBivContinuousContinuous(c1,c2,tabTitle="title",graphDir="graphBiv",graphName="V10",type="png",out="html")
sink()


cat("---------------------------------------------------------------
---                      End testContinous                  ---
---------------------------------------------------------------\n")
