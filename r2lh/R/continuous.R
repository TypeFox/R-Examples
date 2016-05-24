############ 1.1
### Functions for Numeric ~ Logical

### 2 modalities
r2lBivContinuousLogical <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    cat(r2lComment("r2lBivContinuousLogical",out=out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|c|cc|",out=out))

    # First line : Summary, Boxplot, Density
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Density",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y,x,out=out),
                         r2lGraphBoxplot(y,x,graphDir,graphName,type,out=out),
                         r2lGraphDensity(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Second line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Student","Wilcoxon"),line=c(T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out=out))
}



############ 1.2
### Functions for Numeric ~ Factor

### From 3 to 8 modalities
r2lBivContinuousFactorWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivContinuousFactorWide",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|c|cc|",out))

    # First line : Summary, Boxplot, Density
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Density",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y,x,out),
                         r2lGraphBoxplot(y,x,graphDir,graphName,type,out),
                         r2lGraphDensity(y,x,graphDir,graphName,type,out)),
                    out=out))

    # Second line : test
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Anova","KruskalWallis"),line=c(T,F,F),out), span=3, out=out))

    cat(r2lEndStruct(out))
}


### 9 modalities or more
r2lBivContinuousFactorLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivContinuousFactorLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out))

    # First line : Summary
    cat(r2lBuildRow(r2lBold("Summary",out),span=2,hline=FALSE,out=out))
    cat(r2lBuildRow(r2lBivSummary(y,x,out), span=2, out=out))

    # Second line : Boxplot, Density
    cat(r2lBuildRow(c(r2lBold("Boxplot",out),r2lBold("Density",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBoxplot(y,x,graphDir,graphName,type,out),
                         r2lGraphDensity(y,x,graphDir,graphName,type,out)), out=out))

    # Third line : tests
    cat(r2lBuildRow(r2lBold("Tests",out),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Anova","KruskalWallis"),line=c(T,F,F),out), span=2, out=out))

    cat(r2lEndStruct(out))
}



############ 1.3
### Functions for Numeric ~ Ordered

r2lBivContinuousOrderedWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivContinuousOrderedWide",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=5,tabSpec="|ccccc|",out))

    # First line : Summay, Boxplot, ScatterPlot
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Scatter plot",out)),span=c(3,1,1),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y,x,out),
                      r2lGraphBoxplot(y,x,graphDir,graphName,type,out),
                      r2lGraphScatterPlot(y,jitter(as.numeric(x)),graphDir,graphName,type,out)),
                    span=c(3,1,1), out=out))

    # Second line : density, qqplot qqplot, test
    cat(r2lBuildRow(c(r2lBold("Density",out),r2lBold("QQplot (Cont.)",out),r2lBold("QQplot (Ord.)",out),r2lBold("Tests",out)),span=c(1,1,1,2),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphDensity(y,x,graphDir,graphName,type,out),
                      r2lGraphQQPlot(y,graphDir,paste(graphName,"y",sep="-"),type,out),
                      r2lGraphQQPlot(as.numeric(x),graphDir,paste(graphName,"x",sep="-"),type,out),
                      r2lBivTest(y,x,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out)),
                    span=c(1,1,1,2), out=out))

    cat(r2lEndStruct(out))
}
#r2lBivContinuousOrderedMixed(y1,o2,tabTitle,graphDir="graphBiv",graphName="V5",out="latex")

r2lBivContinuousOrderedLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivContinuousOrderedLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|cccc|",out))

    # First line : Summary
    cat(r2lBuildRow(r2lBold("Summary",out),span=4,hline=FALSE,out=out))
    cat(r2lBuildRow(r2lBivSummary(y,x,out), span=4, out=out))

    # Second line : Boxplot, ScatterPlot
    cat(r2lBuildRow(c(r2lBold("Boxplot",out),r2lBold("Scatter plot",out)),span=c(3,1),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBoxplot(y,x,graphDir,graphName,type,out),
                      r2lGraphScatterPlot(y,jitter(as.numeric(x)),graphDir,graphName,type,out)),
                    span=c(3,1), out=out))

    # Third line : density, qqplot qqplot, test
    cat(r2lBuildRow(c(r2lBold("Density",out),r2lBold("QQplot (Cont.)",out),r2lBold("QQplot (Ord.)",out),r2lBold("Tests",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphDensity(y,x,graphDir,graphName,type,out),
                      r2lGraphQQPlot(y,graphDir,paste(graphName,"y",sep="-"),type,out),
                      r2lGraphQQPlot(as.numeric(x),graphDir,paste(graphName,"x",sep="-"),type,out),
                      r2lBivTest(y,x,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out)), out=out))

    cat(r2lEndStruct(out))
}
#r2lBivContinuousOrderedLong(y1,o3,tabTitle,graphDir="graphBiv",graphName="V6",out="latex")




############ 1.4
### Functions for Numeric ~ Discrete

r2lBivContinuousDiscreteWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivContinuousDiscreteWide",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=5,tabSpec="|ccccc|",out))

    # First line : Summay, Boxplot, ScatterPlot
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Scatter plot",out)),span=c(3,1,1),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y,as.factor(x),out),
                      r2lGraphBoxplot(y,as.factor(x),graphDir,graphName,type,out),
                      r2lGraphScatterPlot(y,jitter(x),graphDir,graphName,type,out)),
                    span=c(3,1,1), out=out))

    # Second line : density, qqplot qqplot, test
    cat(r2lBuildRow(c(r2lBold("Density",out),r2lBold("QQplot (Cont.)",out),r2lBold("QQplot (Disc.)",out),r2lBold("Tests",out)),span=c(1,1,1,2),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphDensity(y,as.factor(x),graphDir,graphName,type,out),
                      r2lGraphQQPlot(y,graphDir,paste(graphName,"y",sep="-"),type,out),
                      r2lGraphQQPlot(as.numeric(x),graphDir,paste(graphName,"x",sep="-"),type,out),
                      r2lBivTest(y,x,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out)), span=c(1,1,1,2), out=out))

    cat(r2lEndStruct(out))
}


r2lBivContinuousDiscreteLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex") {
    cat(r2lComment("r2lBivContinuousDiscreteLong",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=4,tabSpec="|cccc|",out))

    # First line : Summary
    cat(r2lBuildRow(r2lBold("Summary",out),span=4,hline=FALSE,out=out))
    cat(r2lBuildRow(r2lBivSummary(y,as.factor(x),out), span=4, out=out))

    # Second line : Boxplot, ScatterPlot
    cat(r2lBuildRow(c(r2lBold("Boxplot",out),r2lBold("Scatter plot",out)),span=c(3,1),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBoxplot(y,x,graphDir,graphName,type,out),
                      r2lGraphScatterPlot(y,jitter(x),graphDir,graphName,type,out)),
                    span=c(3,1), out=out))

    # Third line : density, qqplot qqplot, test
    cat(r2lBuildRow(c(r2lBold("Density",out),r2lBold("QQplot (Cont.)",out),r2lBold("QQplot (Disc.)",out),r2lBold("Tests",out)),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphDensity(y,as.factor(x),graphDir,graphName,type,out),
                      r2lGraphQQPlot(y,graphDir,paste(graphName,"y",sep="-"),type,out),
                      r2lGraphQQPlot(as.numeric(x),graphDir,paste(graphName,"x",sep="-"),type,out),
                      r2lBivTest(y,x,test=c("Anova","KruskalWallis","CorPearson","CorSpearman"),line=c(T,F,T,F,F),out)), out=out))

    cat(r2lEndStruct(out))
}


############ 1.5
### Functions for Numeric ~ Numeric


r2lBivContinuousContinuous <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    cat(r2lComment("r2lBivContinuousContinuous",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=5,tabSpec="|ccccc|",out))

    # First line : Summary boxplot scatterplot
    cat(r2lBuildRow(c(r2lBold("Summary",out),r2lBold("Boxplot",out),r2lBold("Boxplot",out),r2lBold("Scatter plot",out)),span=c(2,1,1,1),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(y,x,out),
                         r2lGraphBoxplot(y=y,graphDir=graphDir,graphName=paste(graphName,"y",sep="-"),type=type,out=out),
                         r2lGraphBoxplot(y=x,graphDir=graphDir,graphName=paste(graphName,"x",sep="-"),type=type,out=out),
                         r2lGraphScatterPlot(y,x,graphDir,graphName,type,out)),
                       span=c(2,1,1,1), out=out))

    # Second line : density, qqplot qqplot, test
    cat(r2lBuildRow(c(r2lBold("Density",out),r2lBold("QQplot (y)",out),r2lBold("QQplot (x)",out),r2lBold("Tests",out)),span=c(1,1,1,2),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphDensity(y,x,graphDir,graphName,type,out),
                      r2lGraphQQPlot(y,graphDir,paste(graphName,"y",sep="-"),type,out),
                      r2lGraphQQPlot(x,graphDir,paste(graphName,"x",sep="-"),type,out),
                      r2lBivTest(y,x,test=c("CorPearson","CorSpearman"),line=c(T,F,F),out)), span=c(1,1,1,2), out=out))

    cat(r2lEndStruct(out))
}


