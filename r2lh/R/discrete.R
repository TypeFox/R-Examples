############ 5.1
### Functions for Discrete ~ Logical

r2lBivDiscreteLogical <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    cat(r2lComment("r2lBivDiscreteLogical",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    ## First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out)),
                    span=2, out=out))

    ## Second line : Quartiles
    cat(r2lBuildRow(c(r2lBold("Summary",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivSummary(y,x,out=out)),
                    span=2, out=out))

    ## Third line : Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    ## Fourth line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","Student","Wilcoxon"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}
#  r2lBivDiscreteLogical(o1,f1,tabTitle,graphDir="graphBiv",graphName="V1",out="latex")


############ 5.2
### Functions for Discrete ~ Factor

r2lBivDiscreteFactorWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteFactorWide",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|cc|cc|",out=out))

    # First line : Table, Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Table",out),r2lBold("Summary",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
					  r2lBivSummary(y,x,out=out),
                      r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Second line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=4))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","Anova","KruskalWallis"),line=c(T,F,T,F,F),out=out), span=4, out=out))

    cat(r2lEndStruct(out))
}


r2lBivDiscreteFactorMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteFactorMixed",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|c|cc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out)),
                    span=3, out=out))

    # Second line : Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y,x,out=out),
                      r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
	cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","Anova","KruskalWallis"),line=c(T,F,T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}


r2lBivDiscreteFactorLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteFactorLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out)),
                    span=2, out=out))

   # Second line : Quartiles
    cat(r2lBuildRow(c(r2lBold("Summary",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivSummary(y,x,out=out)),
                    span=2, out=out))

    # Third line : Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Fourth line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","Anova","KruskalWallis"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}


############ 5.3
### Functions for Discrete ~ Ordered


r2lBivDiscreteOrderedWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteOrderedWide",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|cc|cc|",out=out))

    # First line : Table, Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Table",out),r2lBold("Summary",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out),
					  r2lBivSummary(y,as.factor(x),out=out),
                      r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                      out=out))

    # Second line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","Anova","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,F,T,F),out=out), span=4, out=out))

    cat(r2lEndStruct(out))
}


r2lBivDiscreteOrderedMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteOrderedMixed",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|c|cc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out)),
                    span=3, out=out))

    # Second line : Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y,as.factor(x),out=out),
                      r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","Anova","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,F,T,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}


r2lBivDiscreteOrderedLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteOrderedLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|c|c|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out)),
                    span=2, out=out))

   # Second line : Quartiles
    cat(r2lBuildRow(c(r2lBold("Summary",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivSummary(y,as.factor(x),out=out)),
                    span=2, out=out))

    # Third line : Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                    out=out))

    # Fourth line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","Anova","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,F,T,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}




############ 5.4
### Functions for Discrete ~ Discrete

r2lBivDiscreteDiscreteWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteDiscreteWide",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|cc|cc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out),r2lBold("Summary",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out),
					  r2lBivSummary(y,as.factor(x),out=out),
                      r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                      out=out))

    # Second line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,F,T,F),out=out), span=4, out=out))

    cat(r2lEndStruct(out))
}


r2lBivDiscreteDiscreteMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteDiscreteMixed",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|c|cc|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out)),
                    span=3, out=out))

    # Second line : Quartiles, Barplot, Mosaic
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Barplot",out),r2lBold("Mosaic",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y,as.factor(x),out=out),
                      r2lGraphBarplot(y,as.factor(x),graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,as.factor(x),graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis","KruskalWallisInv","CorPearson","CorSpearman"),line=c(T,F,T,F,F,T,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}


r2lBivDiscreteDiscreteLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteDiscreteLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|c|c|",out=out))

    # First line : Table
    cat(r2lBuildRow(c(r2lBold("Table",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,as.factor(x),out=out)),
                    span=2, out=out))

   # Second line : Quartiles
    cat(r2lBuildRow(c(r2lBold("Summary",out)),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(c(r2lBivSummary(y,as.factor(x),out=out)),
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


############ 5.5
### Functions for Discrete ~ Numeric

r2lBivDiscreteContinuousWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteContinuousWide",out=out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|c|cc|",out=out))

    # First line : Summary, Boxplot, Density
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Density",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y=x,x=as.factor(y),out=out),
                         r2lGraphBoxplot(y=x,x=y,graphDir,graphName,type,out=out),
                         r2lGraphDensity(y=x,x=as.factor(y),graphDir,graphName,type,out=out)),
                    out=out))

    # Second line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y=x,x=y,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out=out))
}


r2lBivDiscreteContinuousMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteContinuousMixed",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=5,tabSpec="|ccccc|",out))

    # First line : Summay, Boxplot, ScatterPlot
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Scatter plot",out)),span=c(3,1,1),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y=x,x=as.factor(y),out),
                      r2lGraphBoxplot(y=x,x=y,graphDir,graphName,type,out),
                      r2lGraphScatterPlot(x,jitter(as.numeric(y)),graphDir,graphName,type,out)),
                    span=c(3,1,1), out=out))

    # Second line : density, qqplot qqplot, test
    cat(r2lBuildRow(c(r2lBold("Density",out),r2lBold("QQplot (Ord.)",out),r2lBold("QQplot (Cont.)",out),r2lBold("Tests",out)),span=c(1,1,1,2),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphDensity(y=x,x=as.factor(y),graphDir,graphName,type,out),
                      r2lGraphQQPlot(as.numeric(y),graphDir,paste(graphName,"y",sep="-"),type,out),
                      r2lGraphQQPlot(x,graphDir,paste(graphName,"x",sep="-"),type,out),
                      r2lBivTest(y=x,x=y,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out)),
                    span=c(1,1,1,2), out=out))

    cat(r2lEndStruct(out))
}


r2lBivDiscreteContinuousLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivDiscreteContinuousLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|ccc|c|",out))

    # First line : Summary
    cat(r2lBuildRow(r2lBold("Summary",out),span=4,hline=FALSE,out=out))
    cat(r2lBuildRow(r2lBivSummary(y=x,x=as.factor(y),out), span=4, out=out))

    # Second line : Boxplot, ScatterPlot
    cat(r2lBuildRow(c(r2lBold("Boxplot",out),r2lBold("Scatter plot",out)),span=c(3,1),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBoxplot(y=x,x=y,graphDir,graphName,type,out),
                      r2lGraphScatterPlot(x,jitter(as.numeric(y)),graphDir,graphName,type,out)),
                    span=c(3,1), out=out))

    # Third line : density, qqplot qqplot, test
    cat(r2lBuildRow(c(r2lBold("Density",out),r2lBold("QQplot (Ord.)",out),r2lBold("QQplot (Cont.)",out),r2lBold("Tests",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphDensity(y=x,x=as.factor(y),graphDir,graphName,type,out),
                      r2lGraphQQPlot(as.numeric(y),graphDir,paste(graphName,"y",sep="-"),type,out),
                      r2lGraphQQPlot(x,graphDir,paste(graphName,"x",sep="-"),type,out),
                      r2lBivTest(y=x,x=y,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out)), out=out))

    cat(r2lEndStruct(out))
}

