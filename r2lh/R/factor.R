############ 1.1
### Functions for Logical ~ Logical

r2lBivFactorLogical <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide"){
    cat(r2lComment("rtlb Factor~Logical",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out))


    ## First line : Table
    cat(r2lBuildColumnTitle(c("Table"),span=2,hline=FALSE,out=out))
    cat(r2lBuildRow(r2lBivContingencyTable(y,x,out=out),span=2,out=out))

    ## Second line : Barplot, Mosaic
    cat(r2lBuildColumnTitle(c("Barplot","Mosaic"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : Tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact"),line=c(T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}
#r2lBivFactorFactor22(f2,f3)


############ 1.2
### Functions for Logical ~ Factor

r2lBivFactorFactorWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Factor (Wide)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out))

    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Table","Barplot","Mosaic"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
                      r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Second line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact"),line=c(T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}


r2lBivFactorFactorMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Factor (Mixed)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    # For the second column of the first line
	attrs <- switch(out,
        "html" = " align='center' cellpadding=5",
        "latex" = "c"
	)
    twoGraph <- r2lStartTable(attrs=attrs,hline=FALSE,out=out)
    twoGraph <- paste(twoGraph,
                      r2lBuildColumnTitle("Barplot",hline=FALSE,out=out),
                      r2lBuildRow(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),out=out),
                      r2lBuildColumnTitle("Mosaic",hline=FALSE,out=out),
                      r2lBuildRow(r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out),hline=FALSE,out=out),
                      r2lEndTable(out))


    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Table",""),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out), twoGraph),out=out))

    # Second line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact"),line=c(T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}
# r2lBivFactorFactorMixed(f2,f3,graphDir="graphBiv",graphName="V4",out="latex")


r2lBivFactorFactorLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor ~ Factor (long)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out))

    # First line : Table
    cat(r2lBuildColumnTitle(c("Table"),span=2,hline=FALSE,out=out))
    cat(r2lBuildRow(r2lBivContingencyTable(y,x,out=out),span=2,out=out))

    # Second line : Barplot, Mosaic
    cat(r2lBuildColumnTitle(c("Barplot","Mosaic"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact"),line=c(T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}





############ 1.3
### Functions for Logical ~ Ordered

r2lBivFactorOrderedWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Ordered (Wide)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out=out))

    # For the second column of the first line
	attrs <- switch(out,
        "html" = " align='center' cellpadding=5",
        "latex" = "c"
	)
    twoGraph <- r2lStartTable(attrs=attrs,hline=FALSE,out=out)
    twoGraph <- paste(twoGraph,
                      r2lBuildRow(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),out=out),
                      r2lBuildColumnTitle("Mosaic",hline=FALSE,out=out),
                      r2lBuildRow(r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out),hline=FALSE,out=out),
                      r2lEndTable(out))


    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Table","Summary","Barplot"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
                      r2lBivQuartilesTable(x,y,out=out),
                      twoGraph
                      ),out=out))

    # Second line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact"),line=c(T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}

r2lBivFactorOrderedMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Ordered (Mixed)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    # For the second column of the first line
	attrs <- switch(out,
        "html" = " align='center' cellpadding=5",
        "latex" = "cc"
	)
    columnTwo <- r2lStartTable(attrs=attrs,hline=FALSE,out=out)
    columnTwo <- paste(columnTwo,
                       r2lBuildRow(c(r2lBivQuartilesTable(x,y,out=out),
                                     r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),hline=TRUE,out=out),
                       r2lBuildColumnTitle("Barplot",span=2,hline=FALSE,out=out,border=""),
                       r2lBuildRow(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),span=2,out=out,border="",hline=FALSE),
                       r2lEndTable(out))

    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Table","Summary - Mosaic"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
                      columnTwo),out=out))

    # Second line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}



r2lBivFactorOrderedLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Ordered (Wide)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out))

    # First line : Table Summary
    cat(r2lBuildColumnTitle(c("Table","Summary"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
                      r2lBivQuartilesTable(x,y,out=out)),out=out))

    # Second line : Barplot, Mosaic
    cat(r2lBuildColumnTitle(c("Barplot","Mosaic"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}



############ 1.4
### Functions for logical ~ discrete

r2lBivFactorDiscreteWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Discrete (Wide)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out=out))

    # For the second column of the first line
	attrs <- switch(out,
        "html" = " align='center' cellpadding=5",
        "latex" = "c"
	)
    twoGraph <- r2lStartTable(attrs=attrs,hline=FALSE,out=out)
    twoGraph <- paste(twoGraph,
                      r2lBuildRow(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),out=out),
                      r2lBuildColumnTitle("Mosaic",hline=FALSE,out=out),
                      r2lBuildRow(r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out),hline=FALSE,out=out),
                      r2lEndTable(out))


    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Table","Summary","Barplot"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
                      r2lBivSummary(x,as.factor(y),out=out),
                      twoGraph
                      ),out=out))

    # Second line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis"),line=c(T,F,T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}
#r2lBivFactorDiscreteWide(concours,nbRedoublement,graphDir="graphBiv",graphName="Vc1",out="latex")


r2lBivFactorDiscreteMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Discrete (Mixed)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    # For the second column of the first line
	attrs <- switch(out,
        "html" = " align='center' cellpadding=5",
        "latex" = "cc"
	)
    columnTwo <- r2lStartTable(attrs=attrs,hline=FALSE,out=out)
    columnTwo <- paste(columnTwo,
                       r2lBuildRow(c(r2lBivSummary(x,as.factor(y),out=out),
                                     r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),hline=TRUE,out=out),
                       r2lBuildColumnTitle("Barplot",span=2,hline=FALSE,out=out,border=""),
                       r2lBuildRow(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),span=2,out=out,border="",hline=FALSE),
                       r2lEndTable(out))

    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Table","Summary - Mosaic"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
                      columnTwo),out=out))

    # Second line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}


r2lBivFactorDiscreteLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Discrete (Wide)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out))

    # First line : Table Summary
    cat(r2lBuildColumnTitle(c("Table","Summary"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
                      r2lBivSummary(x,as.factor(y),out=out)),out=out))

    # Second line : Barplot, Mosaic
    cat(r2lBuildColumnTitle(c("Barplot","Mosaic"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Third line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Anova","KruskalWallis"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}



############ 1.5
### Functions for logical ~ continuous

r2lBivFactorContinuousWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Continuous (Wide)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out=out))

    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Summary","Boxplot","Density"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(x,as.factor(y),out=out),
                      r2lGraphBoxplot(x,y,graphDir,graphName,type,out=out),
                      r2lGraphDensity(x,as.factor(y),graphDir,graphName,type,out=out)
                      ),out=out))

    # Second line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Anova","KruskalWallis"),line=c(T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}


r2lBivFactorContinuousLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Factor~Continuous (Long)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))



    ## # For the second column of the first line
    ##     attrs <- switch(out,
    ##     "html" = " align='center' cellpadding=5",
    ##     "latex" = "c"
    ##     )
    ## twoGraph <- r2lStartTable(attrs=attrs,hline=FALSE,out=out)
    ## twoGraph <- paste(twoGraph,
    ##                   r2lBuildRow(r2lGraphBoxplot(x,y,graphDir,graphName,type,out=out),out=out),
    ##                   r2lBuildColumnTitle("Densities",hline=FALSE,out=out),
    ##                   r2lBuildRow(r2lGraphDensity(x,as.factor(y),graphDir,graphName,type,out=out),hline=FALSE,out=out),
    ##                   r2lEndTable(out))


    ## # First line : Table Barplot Mosaic
    ## cat(r2lBuildColumnTitle(c("Summary","Boxplot"),hline=FALSE,out=out))
    ## cat(r2lBuildRow(c(r2lBivSummary(x,as.factor(y),out=out),
    ##                   twoGraph
    ##                   ),out=out))

    # First line : Summary
    cat(r2lBuildColumnTitle(c("Summary"),hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivSummary(x,as.factor(y),out=out),out=out,span=2))

    # Second line : Boxplot, Densities
    cat(r2lBuildColumnTitle(c("Boxplot","Densities"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lGraphBoxplot(x,y,graphDir,graphName,type,out=out),
                      r2lGraphDensity(x,as.factor(y),graphDir,graphName,type,out=out)),
                    out=out))

    # Thrid line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Anova","KruskalWallis"),line=c(T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}
# sink(file="univLogicalContinuous.tex")
# r2lBivFactorContinuousWide(concours,taille,graphDir="graphBiv",graphName="Vc1",out="latex")
# r2lBivFactorContinuousWide(concours,poids,graphDir="graphBiv",graphName="Vc1",out="latex")
# r2lBivFactorContinuousLong(concours,taille,graphDir="graphBiv",graphName="Vc2",out="latex")
# r2lBivFactorContinuousLong(concours,poids,graphDir="graphBiv",graphName="Vc2",out="latex")
# sink()

#r2lBivFactorDiscreteWide(concours,nbRedoublement,graphDir="graphBiv",graphName="Vc1",out="latex")
