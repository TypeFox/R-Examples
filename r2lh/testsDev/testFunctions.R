cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++                   Begin testFunctions                   +++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

source("../R/functions.R")

cat("###############################################################
###      Utility functions to write out pieces of text      ###
###                    Both univ and biv                    ###
###############################################################\n")

cat("####################
#   Latex output   #
####################\n")

cleanProg(r2lEol)
r2lEol("my line")

cleanProg(r2lComment)
r2lComment("To comment",out="latex")

cleanProg(r2lStartTable)
r2lStartTable(attrs="|c|cc|",out="latex")

cleanProg(r2lEndTable)
r2lEndTable(out="latex")

cleanProg(r2lBuildRow)
r2lBuildRow(c("case1","case2"),span=c(1,1),out="latex")
r2lBuildRow(c("case1","case2"),span=c(1,2),out="latex")
r2lBuildRow(c("case1","case2"),span=c(2,3),hline=FALSE,out="latex")

cleanProg(r2lBuildColumnTitle)
r2lBuildColumnTitle("MyTitle",span=1,out="latex")
r2lBuildColumnTitle(c("MyTitle","Boxplot","Barplot"),span=c(1,2,1),out="latex")

cleanProg(r2lPercentRow)
r2lPercentRow(1/3,out="latex")

cleanProg(r2lIncludeGraphics)
r2lIncludeGraphics("toto.eps",width=300,out="latex")

cleanProg(r2lBold)
r2lBold("ert",out="latex")

cleanProg(r2lUnderline)
r2lUnderline("ert",out="latex")

cleanProg(r2lItal)
r2lItal("aze",out="latex")

cleanProg(r2lStartBlock)
r2lStartBlock("zer",out="latex")

cleanProg(r2lEndBlock)
r2lEndBlock(out="latex")

cleanProg(r2lEndStruct)
r2lEndStruct(out="latex")


cat("###################
#   Html output   #
###################\n")

cleanProg(r2lEol)
r2lEol("my line")

cleanProg(r2lComment)
r2lComment("To comment",out="html")

cleanProg(r2lStartTable)
r2lStartTable(attrs="|c|cc|",out="html")

cleanProg(r2lEndTable)
r2lEndTable(out="html")

cleanProg(r2lBuildRow)
r2lBuildRow(c("case1","case2"),span=c(1,1),out="html")
r2lBuildRow(c("case1","case2"),span=c(1,2),out="html")
r2lBuildRow(c("case1","case2"),span=c(2,3),hline=FALSE,out="html")

cleanProg(r2lPercentRow)
r2lPercentRow(1/3,out="html")

cleanProg(r2lBuildColumnTitle)
r2lBuildColumnTitle("MyTitle",span=1,out="html")
r2lBuildColumnTitle(c("MyTitle","Boxplot","Barplot"),span=c(1,2,1),out="html")

cleanProg(r2lIncludeGraphics)
r2lIncludeGraphics("toto.eps",width=200,out="html")

cleanProg(r2lBold)
r2lBold("ert",out="html")

cleanProg(r2lUnderline)
r2lUnderline("ert",out="html")

cleanProg(r2lItal)
r2lItal("aze",out="html")

cleanProg(r2lStartBlock)
r2lStartBlock("zer",out="html")

cleanProg(r2lEndBlock)
r2lEndBlock(out="html")

cleanProg(r2lEndStruct)
r2lEndStruct(out="html")



cat("###############################################################
###                    Functions for Univ                   ###
###############################################################\n")

cat("####################
#   Latex output   #
####################\n")

cleanProg(r2lUnivBeginStruct)
cat(r2lUnivBeginStruct(f1,"titre",nbColum=3,tabSpec="|c|cc|",out="latex"))

cleanProg(r2lUnivFrequency)
cat(r2lUnivFrequency(f1,out="latex"))
cat(r2lUnivFrequency(f2,out="latex"))

cleanProg(r2lUnivSummary)
cat(r2lUnivSummary(n1,out="latex"))


cat("###################
#   Html output   #
###################\n")

cleanProg(r2lUnivBeginStruct)
cat(r2lUnivBeginStruct(f1,tabTitle="titi",nbColum=3,tabSpec="|c|cc|",out="html"))

cleanProg(r2lUnivFrequency)
cat(r2lUnivFrequency(f1,out="html"))
cat(r2lUnivFrequency(f2,out="html"))

cleanProg(r2lUnivSummary)
cat(r2lUnivSummary(n1,out="html"))



cat("###############################################################
###                     Functions for Biv                   ###
###############################################################\n")

cat("####################
#   Latex output   #
####################\n")

cleanProg(r2lBivBeginStruct)
cat(r2lBivBeginStruct(n1,f1,tabTitle="TITI $\\sim$ toto",nbColumn=3,tabSpec="ccc"))

cleanProg(r2lBivBannerStr)
cat(r2lBivBannerStr(n1))

cleanProg(r2lBivSummaryArray)
cat(r2lBivSummaryArray(n1,f1))
cat(r2lBivSummaryArray(n1,f2))

cleanProg(r2lBivSummary)
cat(r2lBivSummary(n1,f1))
cat(r2lBivSummary(n1,f2))

cleanProg(r2lBivQuartilesArray)
r2lBivQuartilesArray(o1,f1)
r2lBivQuartilesArray(o1,f2)

cleanProg(r2lBivQuartilesTable)
cat(r2lBivQuartilesTable(o1,f1))
cat(r2lBivQuartilesTable(o1,f2))

cleanProg(r2lBivContingencyTable)
cat(r2lBivContingencyTable(f1,f2,out="latex"))



cat("###################
#   Html output   #
###################\n")

cleanProg(r2lBivBeginStruct)
cat(r2lBivBeginStruct(n1,f1,tabTitle="TITI $\\sim$ toto",nbColumn=3,tabSpec="|c|cc|",out="html"))

cleanProg(r2lBivSummary)
cat(r2lBivSummary(n1,f1,out="html"))
cat(r2lBivSummary(n1,f2,out="html"))

cleanProg(r2lBivContingencyTable)
cat(r2lBivContingencyTable(f1,f2,out="html"))

cleanProg(r2lBivQuartilesTable)
cat(r2lBivQuartilesTable(o1,f1,out="html"))
cat(r2lBivQuartilesTable(o1,f2,out="html"))


cat("###############################################################
###                           tests                         ###
###############################################################\n")

cleanProg(r2lSignificance)
cat(r2lSignificance(0.2))
cat(r2lSignificance(0.005))

cleanProg(r2lPValueStr)
cat(r2lPValueStr(0.2))
cat(r2lPValueStr(0.005))
cat(r2lPValueStr(0.001))
cat(r2lPValueStr(0.0001))
cat(r2lPValueStr(0.000012345))
cat(r2lPValueStr(0.00000012345))
cat(r2lPValueStr(0.000000000000000012345))

cat("####################
#   Latex output   #
####################\n")

cleanProg(r2lBivTestCorPearson)
cat(r2lBivTestCorPearson(n1,n2,out="latex"))
cat(r2lBivTestCorPearson(n2,n3,out="latex"))

cleanProg(r2lBivTestCorSpearman)
cat(r2lBivTestCorSpearman(n1,n2,out="latex"))
cat(r2lBivTestCorSpearman(n2,n3,out="latex"))

cleanProg(r2lBivTestStudent)
cat(r2lBivTestStudent(n1,f1,out="latex"))

cleanProg(r2lBivTestWilcoxon)
cat(r2lBivTestWilcoxon(n1,f1,out="latex"))

cleanProg(r2lBivTestAnova)
cat(r2lBivTestAnova(n1,f2,out="latex"))

cleanProg(r2lBivTestKruskalWallis)
cat(r2lBivTestKruskalWallis(n1,f2,out="latex"))

cleanProg(r2lBivTestKhi2)
cat(r2lBivTestKhi2(f1,f2,out="latex"))
cat(r2lBivTestKhi2(f1,f1b,out="latex"))

cleanProg(r2lBivTestFisherExact)
cat(r2lBivTestFisherExact(f1,f2,out="latex"))
cat(r2lBivTestFisherExact(f1,f1b,out="latex"))

cleanProg(r2lBivTestOddsRatio)
cat(r2lBivTestOddsRatio(f1,f1b,out="latex"))

cleanProg(r2lBivTestRelativeRisk)
cat(r2lBivTestRelativeRisk(f1,f1b,out="latex"))

cleanProg(r2lBivTest)
cat(r2lBivTest(n1,f1,test="Student",line=c(T,F),out="latex"))
cat(r2lBivTest(n1,f2,test="Anova",line=c(F,T),out="latex"))
cat(r2lBivTest(n1,o1,test=c("Anova","CorPearson"),line=c(T,F,F),out="latex"))
cat(r2lBivTest(n1,o1,test=c("Anova","CorPearson"),line=c(T,T,T),out="latex"))
cat(r2lBivTest(n1,n2,test="CorPearson",line=c(T,F),out="latex"))
cat(r2lBivTest(f1,f1b,test=c("Khi2","OddsRatio"),line=c(F,F,F),out="latex"))


cat("###################
#   Html output   #
###################\n")

cleanProg(r2lBivTestCorPearson)
cat(r2lBivTestCorPearson(n1,n2,out="html"))
cat(r2lBivTestCorPearson(n2,n3,out="html"))

cleanProg(r2lBivTestCorSpearman)
cat(r2lBivTestCorSpearman(n1,n2,out="html"))
cat(r2lBivTestCorSpearman(n2,n3,out="html"))

cleanProg(r2lBivTestStudent)
cat(r2lBivTestStudent(n1,f1,out="html"))

cleanProg(r2lBivTestWilcoxon)
cat(r2lBivTestWilcoxon(n1,f1,out="html"))

cleanProg(r2lBivTestAnova)
cat(r2lBivTestAnova(n1,f2,out="html"))

cleanProg(r2lBivTestKruskalWallis)
cat(r2lBivTestKruskalWallis(n1,f2,out="html"))

cleanProg(r2lBivTestKhi2)
cat(r2lBivTestKhi2(f1,f2,out="html"))
cat(r2lBivTestKhi2(f1,f1b,out="html"))

cleanProg(r2lBivTestOddsRatio)
cat(r2lBivTestOddsRatio(f1,f1b,out="html"))

cleanProg(r2lBivTest)
cat(r2lBivTest(n1,f1,test="Student",line=c(T,F),out="html"))
cat(r2lBivTest(n1,f2,test="Anova",line=c(F,T),out="html"))
cat(r2lBivTest(n1,o1,test=c("Anova","CorPearson"),line=c(T,F,F),out="html"))
cat(r2lBivTest(n1,o1,test=c("Anova","CorPearson"),line=c(T,T,T),out="html"))
cat(r2lBivTest(n1,n2,test="CorPearson",line=c(T,F),out="html"))
cat(r2lBivTest(f1,f1b,test=c("Khi2","OddsRatio"),line=c(F,F,F),out="html"))


cat("###############################################################
###                          graphs                         ###
###############################################################\n")

cleanProg(r2lDensities)
r2lDensities(y=n1,x=f1,main="", xlab="", ylab="")
r2lDensities(n1,f2)
r2lDensities(n2,f2)

cat("####################
#   Latex output   #
####################\n")

cleanProg(r2lMakePlot)
r2lMakePlot(kind="boxplot",arguments=list(x=n1),type="png",graphName="V1",graphDir="toto", out="latex")
r2lMakePlot(kind="boxplot",arguments=list(n1~f2),type="png",graphName="V2",graphDir="toto", out="latex")

cleanProg(r2lGraphQQPlot)
cat(r2lGraphQQPlot(x=n1,type="png",graphName="V3",graphDir="toto", out="latex"))
cat(r2lGraphQQPlot(x=n2,type="png",graphName="V4",graphDir="toto", out="latex"))

cleanProg(r2lGraphHist)
cat(r2lGraphHist(x=n1,type="png",graphName="V5",graphDir="toto", out="latex"))
cat(r2lGraphHist(x=n2,type="png",graphName="V6",graphDir="toto", out="latex"))

cleanProg(r2lGraphBoxplot)
cat(r2lGraphBoxplot(n1,type="png",graphName="V7",graphDir="toto", out="latex"))
cat(r2lGraphBoxplot(n1,f1,type="png",graphName="V8",graphDir="toto", out="latex"))
cat(r2lGraphBoxplot(n1,f2,type="png",graphName="V9",graphDir="toto", out="latex"))

cleanProg(r2lGraphMosaicPlot)
cat(r2lGraphMosaicPlot(f1,f2,type="png",graphName="V10",graphDir="toto", out="latex"))
cat(r2lGraphMosaicPlot(f2,f1b,type="png",graphName="V11",graphDir="toto", out="latex"))

cleanProg(r2lGraphScatterPlot)
cat(r2lGraphScatterPlot(n1,n2,type="png",graphName="V12",graphDir="toto", out="latex"))
cat(r2lGraphScatterPlot(n2,n3,type="png",graphName="V13",graphDir="toto", out="latex"))

cleanProg(r2lGraphBarplot)
cat(r2lGraphBarplot(f1,type="png",graphName="V14",graphDir="toto", out="latex"))
cat(r2lGraphBarplot(f2,type="png",graphName="V15",graphDir="toto", out="latex"))
cat(r2lGraphBarplot(f2,f1,type="png",graphName="V16",graphDir="toto", out="latex"))

cleanProg(r2lGraphDensity)
cat(r2lGraphDensity(n1,f1,type="png",graphName="V",graphDir="toto", out="latex"))
cat(r2lGraphDensity(n1,f2,type="png",graphName="V",graphDir="toto", out="latex"))
cat(r2lGraphDensity(n2,f2,type="png",graphName="V",graphDir="toto", out="latex"))


cat("###################
#   Html output   #
###################\n")

cleanProg(r2lMakePlot)
r2lMakePlot(kind="boxplot",arguments=list(x=n1),type="png",graphName="V1",graphDir="toto", out="html")
r2lMakePlot(kind="boxplot",arguments=list(n1~f2),type="png",graphName="V2",graphDir="toto", out="html")

cleanProg(r2lGraphQQPlot)
cat(r2lGraphQQPlot(x=n1,type="png",graphName="V3",graphDir="toto", out="html"))
cat(r2lGraphQQPlot(x=n2,type="png",graphName="V4",graphDir="toto", out="html"))

cleanProg(r2lGraphHist)
cat(r2lGraphHist(x=n1,type="png",graphName="V5",graphDir="toto", out="html"))
cat(r2lGraphHist(x=n2,type="png",graphName="V6",graphDir="toto", out="html"))

cleanProg(r2lGraphBoxplot)
cat(r2lGraphBoxplot(n1,type="png",graphName="V7",graphDir="toto", out="html"))
cat(r2lGraphBoxplot(n1,f1,type="png",graphName="V8",graphDir="toto", out="html"))
cat(r2lGraphBoxplot(n1,f2,type="png",graphName="V9",graphDir="toto", out="html"))

cleanProg(r2lGraphMosaicPlot)
cat(r2lGraphMosaicPlot(f1,f2,type="png",graphName="V10",graphDir="toto", out="html"))
cat(r2lGraphMosaicPlot(f2,f1b,type="png",graphName="V11",graphDir="toto", out="html"))

cleanProg(r2lGraphScatterPlot)
cat(r2lGraphScatterPlot(n1,n2,type="png",graphName="V12",graphDir="toto", out="html"))
cat(r2lGraphScatterPlot(n2,n3,type="png",graphName="V13",graphDir="toto", out="html"))

cleanProg(r2lGraphBarplot)
cat(r2lGraphBarplot(f1,type="png",graphName="V14",graphDir="toto", out="html"))
cat(r2lGraphBarplot(f2,type="png",graphName="V15",graphDir="toto", out="html"))
cat(r2lGraphBarplot(f2,f1,type="png",graphName="V16",graphDir="toto", out="html"))

cleanProg(r2lGraphDensity)
cat(r2lGraphDensity(n1,f1,type="png",graphName="V",graphDir="toto", out="html"))
cat(r2lGraphDensity(n1,f2,type="png",graphName="V",graphDir="toto", out="html"))
cat(r2lGraphDensity(n2,f2,type="png",graphName="V",graphDir="toto", out="html"))



cat("---------------------------------------------------------------
---                    Fin testFunctions                    ---
---------------------------------------------------------------\n")
