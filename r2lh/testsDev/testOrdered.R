
cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++                    Begin testOrdered                    +++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

source("../R/ordered.R")



####################
#   Latex output   #
####################

# Type Postscript

sink(file="bivOrdLogical.tex")
r2lBivOrderedLogical(o1,f1,tabTitle="title",graphDir="graphBiv",graphName="V1",out="latex")
sink()

sink(file="bivOrdFactor.tex")
r2lBivOrderedFactorWide(o2,f1,tabTitle="title",graphDir="graphBiv",graphName="V3a",out="latex")
r2lBivOrderedFactorMixed(o2,f2,tabTitle="title",graphDir="graphBiv",graphName="V3b",out="latex")
r2lBivOrderedFactorLong(o3,f3,tabTitle="title",graphDir="graphBiv",graphName="V3c",out="latex")
sink()

sink(file="bivOrdOrdered.tex")
r2lBivOrderedOrderedWide(o2,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",out="latex")
r2lBivOrderedOrderedMixed(o3,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",out="latex")
r2lBivOrderedOrderedLong(o2,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",out="latex")
sink()

sink(file="bivOrdDiscrete.tex")
r2lBivOrderedDiscreteWide(o2,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",out="latex")
r2lBivOrderedDiscreteMixed(o3,d2,tabTitle="title",graphDir="graphBiv",graphName="V8",out="latex")
r2lBivOrderedDiscreteLong(o3,d3,tabTitle="title",graphDir="graphBiv",graphName="V9",out="latex")
sink()

sink(file="bivOrdContinuous.tex")
r2lBivOrderedContinuousWide(o2,c2,tabTitle="title",graphDir="graphBiv",graphName="V10",out="latex")
r2lBivOrderedContinuousMixed(o3,c2,tabTitle="title",graphDir="graphBiv",graphName="V11",out="latex")
r2lBivOrderedContinuousLong(o3,c2,tabTitle="title",graphDir="graphBiv",graphName="V12",out="latex")
sink()


###################
#   Html output   #
###################

# Type PNG
sink(file="bivContLogical.html")
r2lBivOrderedLogical(o1,f1,tabTitle="title",graphDir="graphBiv",graphName="V1",type="png",out="html")
sink()


#r2lBivOrderedFactorWide(o1,f2,tabTitle="title",graphDir="graphBiv",graphName="V2",type="png",out="html")
#r2lBivOrderedFactorLong(o1,f3,tabTitle="title",graphDir="graphBiv",graphName="V3",type="png",out="html")

sink(file="bivContFactor.html")
r2lBivOrderedFactorWide(o2,f1,tabTitle="title",graphDir="graphBiv",graphName="V3a",type="png",out="html")
r2lBivOrderedFactorMixed(o3,f2,tabTitle="title",graphDir="graphBiv",graphName="V3b",type="png",out="html")
r2lBivOrderedFactorLong(o2,f3,tabTitle="title",graphDir="graphBiv",graphName="V3c",type="png",out="html")
sink()

sink(file="bivContOrdered.html")
r2lBivOrderedOrderedWide(o2,o1,tabTitle="title",graphDir="graphBiv",graphName="V4",type="png",out="html")
r2lBivOrderedOrderedMixed(o2,o2,tabTitle="title",graphDir="graphBiv",graphName="V5",type="png",out="html")
r2lBivOrderedOrderedLong(o3,o3,tabTitle="title",graphDir="graphBiv",graphName="V6",type="png",out="html")
sink()

sink(file="bivContDiscrete.html")
r2lBivOrderedDiscreteWide(o2,d1,tabTitle="title",graphDir="graphBiv",graphName="V7",type="png",out="html")
r2lBivOrderedDiscreteMixed(o2,d2,tabTitle="title",graphDir="graphBiv",graphName="V8",type="png",out="html")
r2lBivOrderedDiscreteLong(o3,d3,tabTitle="title",graphDir="graphBiv",graphName="V9",type="png",out="html")
sink()

sink(file="bivOrdContinuous.html")
r2lBivOrderedContinuousWide(o2,c2,tabTitle="title",graphDir="graphBiv",graphName="V10",type="png",out="html")
r2lBivOrderedContinuousMixed(o3,c2,tabTitle="title",graphDir="graphBiv",graphName="V11",type="png",out="html")
r2lBivOrderedContinuousLong(o3,c2,tabTitle="title",graphDir="graphBiv",graphName="V12",type="png",out="html")
sink()


cat("---------------------------------------------------------------
---                     End testOrdered                     ---
---------------------------------------------------------------\n")
