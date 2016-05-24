cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++                     Begin testFactor                    +++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

source("../R/factor.R")



####################
#   Latex output   #
####################

# Type Postscript

sink(file="bivFactorLogical.tex")
r2lBivFactorLogical(bacType,sexe,tabTitle="title",graphDir="graphBiv",graphName="Va1",out="latex")
r2lBivFactorLogical(bacType2,sexe,tabTitle="title",graphDir="graphBiv",graphName="Va2",out="latex")
sink()

sink(file="bivFactorFactor.tex")
cat("\\paragraph{Wide}")
r2lBivFactorFactorWide(bacType,cheveux,tabTitle="title",graphDir="graphBiv",graphName="Vb1",out="latex")
cat("\\paragraph{Mixed}")
r2lBivFactorFactorMixed(bacType,bacType2,tabTitle="title",graphDir="graphBiv",graphName="Vb2",out="latex")
cat("\\paragraph{Mixed}")
r2lBivFactorFactorMixed(bacType2,cheveux,tabTitle="title",graphDir="graphBiv",graphName="Vb3",out="latex")
cat("\\paragraph{Long}")
r2lBivFactorFactorLong(bacType2,bacType,tabTitle="title",graphDir="graphBiv",graphName="Vb4",out="latex")
sink()

sink(file="bivFactorOrdered.tex")
cat("\\paragraph{Wide}")
r2lBivFactorOrderedWide(bacType,bacMension,tabTitle="title",graphDir="graphBiv",graphName="Vc1",out="latex")
cat("\\paragraph{Mixed}")
r2lBivFactorOrderedMixed(bacType,bacMension2,tabTitle="title",graphDir="graphBiv",graphName="Vc2",out="latex")
cat("\\paragraph{Mixed}")
r2lBivFactorOrderedMixed(bacType2,bacMension,tabTitle="title",graphDir="graphBiv",graphName="Vc3",out="latex")
cat("\\paragraph{Long}")
r2lBivFactorOrderedLong(bacType2,bacMension2,tabTitle="title",graphDir="graphBiv",graphName="Vc4",out="latex")
sink()

sink(file="bivFactorDiscrete.tex")
cat("\\paragraph{Wide}")
r2lBivFactorDiscreteWide(bacType,nbRedoublement,tabTitle="title",graphDir="graphBiv",graphName="Vc1",out="latex")
cat("\\paragraph{Mixed}")
r2lBivFactorDiscreteMixed(bacType,nbRedoublement2,tabTitle="title",graphDir="graphBiv",graphName="Vc2",out="latex")
cat("\\paragraph{Mixed}")
r2lBivFactorDiscreteMixed(bacType2,nbRedoublement,tabTitle="title",graphDir="graphBiv",graphName="Vc3",out="latex")
cat("\\paragraph{Long}")
r2lBivFactorDiscreteLong(bacType2,nbRedoublement2,tabTitle="title",graphDir="graphBiv",graphName="Vc4",out="latex")
sink()


sink(file="bivFactorContinuous.tex")
cat("\\paragraph{Wide}")
r2lBivFactorContinuousWide(bacType,taille,tabTitle="title",graphDir="graphBiv",graphName="Vc1",out="latex")
cat("\\paragraph{Wide}")
r2lBivFactorContinuousWide(bacType,poids,tabTitle="title",graphDir="graphBiv",graphName="Vc1",out="latex")
cat("\\paragraph{Long}")
r2lBivFactorContinuousLong(bacType2,taille,tabTitle="title",graphDir="graphBiv",graphName="Vc2",out="latex")
cat("\\paragraph{Long}")
r2lBivFactorContinuousLong(bacType2,poids,tabTitle="title",graphDir="graphBiv",graphName="Vc2",out="latex")
sink()



###################
#   Html output   #
###################


# Type PNG
sink(file="bivFactorLogical.html")
r2lBivFactorLogical(bacType,sexe,tabTitle="title",graphDir="graphBiv",graphName="Va1",type="png",out="html")
#r2lBivFactorLogical(bacType2,traitement,tabTitle="title",graphDir="graphBiv",graphName="Va2",type="png",out="html")
sink()

sink(file="bivFactorFactor.html")
cat("<H3>Wide</H3>")
r2lBivFactorFactorWide(bacType,bacType,tabTitle="title",graphDir="graphBiv",graphName="Vb1",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivFactorFactorMixed(bacType2,bacType,tabTitle="title",graphDir="graphBiv",graphName="Vb2",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivFactorFactorMixed(bacType,bacType2,tabTitle="title",graphDir="graphBiv",graphName="Vb3",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivFactorFactorLong(bacType2,bacType2,tabTitle="title",graphDir="graphBiv",graphName="Vb4",type="png",out="html")
sink()

sink(file="bivFactorOrdered.html")
cat("<H3>Wide</H3>")
r2lBivFactorOrderedWide(bacType,bacMension,tabTitle="title",graphDir="graphBiv",graphName="Vc1",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivFactorOrderedMixed(bacType,bacMension,tabTitle="title",graphDir="graphBiv",graphName="Vc2",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivFactorOrderedMixed(bacType2,bacMension2,tabTitle="title",graphDir="graphBiv",graphName="Vc3",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivFactorOrderedLong(bacType2,bacMension2,tabTitle="title",graphDir="graphBiv",graphName="Vc4",type="png",out="html")
sink()

sink(file="bivFactorDiscrete.html")
cat("<H3>Wide</H3>")
r2lBivFactorDiscreteWide(bacType,nbRedoublement,tabTitle="title",graphDir="graphBiv",graphName="Vc1",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivFactorDiscreteMixed(bacType,nbRedoublement,tabTitle="title",graphDir="graphBiv",graphName="Vc2",type="png",out="html")
cat("<H3>Mixed</H3>")
r2lBivFactorDiscreteMixed(bacType2,nbRedoublement,tabTitle="title",graphDir="graphBiv",graphName="Vc3",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivFactorDiscreteLong(bacType2,nbRedoublement,tabTitle="title",graphDir="graphBiv",graphName="Vc4",type="png",out="html")
sink()

sink(file="bivFactorContinuous.html")
cat("<H3>Wide</H3>")
r2lBivFactorContinuousWide(bacType,taille,tabTitle="title",graphDir="graphBiv",graphName="Vc1",type="png",out="html")
cat("<H3>Wide</H3>")
r2lBivFactorContinuousWide(bacType,poids,tabTitle="title",graphDir="graphBiv",graphName="Vc1",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivFactorContinuousLong(bacType2,taille,tabTitle="title",graphDir="graphBiv",graphName="Vc2",type="png",out="html")
cat("<H3>Long</H3>")
r2lBivFactorContinuousLong(bacType2,poids,tabTitle="title",graphDir="graphBiv",graphName="Vc2",type="png",out="html")
sink()




cat("---------------------------------------------------------------
---                      End testFactor                     ---
---------------------------------------------------------------\n")
