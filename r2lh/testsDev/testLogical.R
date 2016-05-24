cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++                    Begin testLogical                    +++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")


source("../R/logical.R")



cat("####################
#   Latex output   #
####################\n")

# Type Postscript

sink(file="bivLogicalLogical.tex")
r2lBivLogicalLogical(concours,sexe,tabTitle="concours $\\sim$ sexe",graphDir="graphBiv",graphName="Va1",out="latex")
r2lBivLogicalLogical(malade,traitement,tabTitle="malade $\\sim$ traitement",graphDir="graphBiv",graphName="Va2",out="latex")
sink()

sink(file="bivLogicalFactor.tex")
cat("\\paragraph{Wide}")
r2lBivLogicalFactorWide(concours,bacType,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vb1",out="latex")
cat("\\paragraph{Mixed}")
r2lBivLogicalFactorMixed(concours,bacType,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vb2",out="latex")
cat("\\paragraph{Mixed}")
r2lBivLogicalFactorMixed(concours,bacType2,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vb3",out="latex")
cat("\\paragraph{Long}")
r2lBivLogicalFactorLong(concours,bacType2,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vb4",out="latex")
sink()

sink(file="bivLogicalOrdered.tex")
cat("\\paragraph{Wide}")
r2lBivLogicalOrderedWide(concours,bacMension,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc1",out="latex")
cat("\\paragraph{Mixed}")
r2lBivLogicalOrderedMixed(concours,bacMension,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc2",out="latex")
cat("\\paragraph{Mixed}")
r2lBivLogicalOrderedMixed(concours,bacMension2,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc3",out="latex")
cat("\\paragraph{Long}")
r2lBivLogicalOrderedLong(concours,bacMension2,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc4",out="latex")
sink()

sink(file="bivLogicalDiscrete.tex")
cat("\\paragraph{Wide}")
r2lBivLogicalDiscreteWide(concours,nbRedoublement,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc1",out="latex")
cat("\\paragraph{Mixed}")
r2lBivLogicalDiscreteMixed(concours,nbRedoublement,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc2",out="latex")
cat("\\paragraph{Mixed}")
r2lBivLogicalDiscreteMixed(concours,nbRedoublement,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc3",out="latex")
cat("\\paragraph{Long}")
r2lBivLogicalDiscreteLong(concours,nbRedoublement,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc4",out="latex")
sink()

sink(file="bivLogicalContinuous.tex")
cat("\\paragraph{Wide}")
r2lBivLogicalContinuous(concours,taille,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc1",out="latex")
cat("\\paragraph{Wide}")
r2lBivLogicalContinuous(concours,poids,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc1",out="latex")
cat("\\paragraph{Long}")
r2lBivLogicalContinuousLong(concours,taille,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc2",out="latex")
cat("\\paragraph{Long}")
r2lBivLogicalContinuousLong(concours,poids,tabTitle="var1 $\\sim$ var2",graphDir="graphBiv",graphName="Vc2",out="latex")
sink()



cat("###################
#   Html output   #
###################\n")


# Type PNG
sink(file="bivLogicalLogical.html")
r2lBivLogicalLogical(concours,sexe,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Va1",type="png",out="html")
r2lBivLogicalLogical(malade,traitement,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Va2",type="png",out="html")
sink()

sink(file="bivLogicalFactor.html")
cat("<H3>Wide</H3>")
r2lBivLogicalFactorWide(concours,bacType,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vb1",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivLogicalFactorMixed(concours,bacType,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vb2",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivLogicalFactorMixed(concours,bacType2,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vb3",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivLogicalFactorLong(concours,bacType2,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vb4",type="png",out="html")
sink()

sink(file="bivLogicalOrdered.html")
cat("<H3>Wide</H3>")
r2lBivLogicalOrderedWide(concours,bacMension,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc1",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivLogicalOrderedMixed(concours,bacMension,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc2",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivLogicalOrderedMixed(concours,bacMension2,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc3",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivLogicalOrderedLong(concours,bacMension2,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc4",type="png",out="html")
sink()

sink(file="bivLogicalDiscrete.html")
cat("<H3>Wide</H3>")
r2lBivLogicalDiscreteWide(concours,nbRedoublement,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc1",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivLogicalDiscreteMixed(concours,nbRedoublement,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc2",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivLogicalDiscreteMixed(concours,nbRedoublement,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc3",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivLogicalDiscreteLong(concours,nbRedoublement,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc4",type="png",out="html")
sink()

sink(file="bivLogicalContinuous.html")
cat("<H3>Wide</H3>")
r2lBivLogicalContinuous(concours,taille,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc1",type="png",out="html")
cat("<H3>Wide</H3>")
r2lBivLogicalContinuous(concours,poids,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc1",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivLogicalContinuousLong(concours,taille,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc2",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivLogicalContinuousLong(concours,poids,tabTitle="var1 ~ var2",graphDir="graphBiv",graphName="Vc2",type="png",out="html")
sink()




cat("---------------------------------------------------------------
---                     End testLogical                     ---
---------------------------------------------------------------\n")
