############ 4.1
### Functions for Ordered ~ Logical

r2lBivOrderedLogical <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    cat(r2lComment("r2lBivOrderedLogical",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out)),out=out,span=2))

    # Second line : Quartiles
    cat(r2lBuildRow(c(r2lBold("Quartiles",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivQuartilesTable(y,x,out=out)),out=out,span=2))

    # Third line : Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Fourth line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","Anova"),line=c(T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}
#  r2lBivOrderedLogical(o1,f1,tabTitle,graphDir="graphBiv",graphName="V1",out="latex")


############ 4.2
### Functions for Ordered ~ Factor

r2lBivOrderedFactorWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedFactorWide",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|cccc|",out=out))

    # First line : Table, Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Table",out),r2lBold("Quartiles",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
					  r2lBivQuartilesTable(y,x,out=out),
                      r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Second line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=4))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","Anova"),line=c(T,F,F),out=out), span=4, out=out))

    cat(r2lEndStruct(out))
}


r2lBivOrderedFactorMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedFactorMixed",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out)),
                    span=3, out=out))

    # Second line : Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Quartiles",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivQuartilesTable(y,x,out=out),
                      r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
	cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","Anova"),line=c(T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}


r2lBivOrderedFactorLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedFactorLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out)),
                    span=2, out=out))

   # Second line : Quartiles
    cat(r2lBuildRow(c(r2lBold("Quartiles",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivQuartilesTable(y,x,out=out)),
                    span=2, out=out))

    # Third line : Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Fourth line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","Anova"),line=c(T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}


############ 4.3
### Functions for Ordered ~ Ordered


r2lBivOrderedOrderedWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedOrderedWide",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|cccc|",out=out))

    # First line : Table, Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Table",out),r2lBold("Quartiles",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out),
					  r2lBivQuartilesTable(y,as.factor(x),out=out),
                      r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                      out=out))

    # Second line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,T,F,F),out=out), span=4, out=out))

    cat(r2lEndStruct(out))
}


r2lBivOrderedOrderedMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedOrderedMixed",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out)),
                    span=3, out=out))

    # Second line : Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Quartiles",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivQuartilesTable(y,as.factor(x),out=out),
                      r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}


r2lBivOrderedOrderedLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedOrderedLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out)),
                    span=2, out=out))

   # Second line : Quartiles
    cat(r2lBuildRow(c(r2lBold("Quartiles",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivQuartilesTable(y,as.factor(x),out=out)),
                    span=2, out=out))

    # Third line : Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                    out=out))

    # Fourth line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}




############ 4.4
### Functions for Ordered ~ Discrete

r2lBivOrderedDiscreteWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedDiscreteWide",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|cccc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out),r2lBold("Quartiles",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out),
					  r2lBivQuartilesTable(y,as.factor(x),out=out),
                      r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                      out=out))

    # Second line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,F,T,F),out=out), span=4, out=out))

    cat(r2lEndStruct(out))
}


r2lBivOrderedDiscreteMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedDiscreteMixed",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out)),
                    span=3, out=out))

    # Second line : Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Quartiles",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivQuartilesTable(y,as.factor(x),out=out),
                      r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,F,T,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}


r2lBivOrderedDiscreteLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedDiscreteLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out)),
                    span=2, out=out))

   # Second line : Quartiles
    cat(r2lBuildRow(c(r2lBold("Quartiles",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivQuartilesTable(y,as.factor(x),out=out)),
                    span=2, out=out))

    # Third line : Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                    out=out))

    # Fourth line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,F,T,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}


############ 4.5
### Functions for Ordered ~ Numeric

r2lBivOrderedContinuousWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedContinuousWide",out=out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out=out))

    # First line : Summary, Boxplot, Density
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Density",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y=x,x=y,out=out),
                         r2lGraphBoxplot(y=x,x=y,graphDir,graphName,type,out=out),
                         r2lGraphDensity(y=x,x=y,graphDir,graphName,type,out=out)),
                    out=out))

    # Second line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y=x,x=y,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out=out))
}


r2lBivOrderedContinuousMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedContinuousMixed",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=5,tabSpec="|ccccc|",out))

    # First line : Summay, Boxplot, ScatterPlot
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Scatter plot",out)),span=c(3,1,1),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y=x,x=y,out),
                      r2lGraphBoxplot(y=x,x=y,graphDir,graphName,type,out),
                      r2lGraphScatterPlot(x,jitter(as.numeric(y)),graphDir,graphName,type,out)),
                    span=c(3,1,1), out=out))

    # Second line : density, qqplot qqplot, test
    cat(r2lBuildRow(c(r2lBold("Density",out),r2lBold("QQplot (Ord.)",out),r2lBold("QQplot (Cont.)",out),r2lBold("Tests",out)),span=c(1,1,1,2),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphDensity(y=x,x=y,graphDir,graphName,type,out),
                      r2lGraphQQPlot(as.numeric(y),graphDir,paste(graphName,"y",sep="-"),type,out),
                      r2lGraphQQPlot(x,graphDir,paste(graphName,"x",sep="-"),type,out),
                      r2lBivTest(y=x,x=y,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out)),
                    span=c(1,1,1,2), out=out))

    cat(r2lEndStruct(out))
}


r2lBivOrderedContinuousLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivOrderedContinuousLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|cccc|",out))

    # First line : Summary
    cat(r2lBuildRow(r2lBold("Summary",out),span=4,hline=FALSE,out=out))
    cat(r2lBuildRow(r2lBivSummary(y=x,x=y,out), span=4, out=out))

    # Second line : Boxplot, ScatterPlot
    cat(r2lBuildRow(c(r2lBold("Boxplot",out),r2lBold("Scatter plot",out)),span=c(3,1),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBoxplot(y=x,x=y,graphDir,graphName,type,out),
                      r2lGraphScatterPlot(x,jitter(as.numeric(y)),graphDir,graphName,type,out)),
                    span=c(3,1), out=out))

    # Third line : density, qqplot qqplot, test
    cat(r2lBuildRow(c(r2lBold("Density",out),r2lBold("QQplot (Ord.)",out),r2lBold("QQplot (Cont.)",out),r2lBold("Tests",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphDensity(y=x,x=y,graphDir,graphName,type,out),
                      r2lGraphQQPlot(as.numeric(y),graphDir,paste(graphName,"y",sep="-"),type,out),
                      r2lGraphQQPlot(x,graphDir,paste(graphName,"x",sep="-"),type,out),
                      r2lBivTest(y=x,x=y,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out)), out=out))

    cat(r2lEndStruct(out))
}

