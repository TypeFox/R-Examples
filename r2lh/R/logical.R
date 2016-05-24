#######################################################
########### Functions for Logical ~ Logical ###########
#######################################################

r2lBivLogicalLogical <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide"){
    cat(r2lComment("rtlb Logical~Logical",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out))

    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Table","Barplot","Mosaic"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
                      r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Second line : Tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact","OddsRatio","RelativeRisk"),line=c(T,F,T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}




#######################################################
########### Functions for Logical ~ Factor ############
#######################################################

r2lBivLogicalFactorWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Logical~Factor (Wide)",out))
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

### Non utilisé
r2lBivLogicalFactorMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Logical~Factor (Mixed)",out))
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
# r2lBivLogicalFactorMixed(f2,f3,graphDir="graphBiv",graphName="V4",out="latex")


r2lBivLogicalFactorLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Logical ~ Factor (long)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out))

    # First line : Table and Mosaic
    cat(r2lBuildColumnTitle(c("Table","Mosaic"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
                      r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out)),
                    out=out))

    # Second line : Barplot, Mosaic
    cat(r2lBuildColumnTitle("Barplot",span=2,hline=FALSE,out=out))
    cat(r2lBuildRow(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),span=2,out=out))

    # Third line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(y,x,test=c("Khi2","FisherExact"),line=c(T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}


r2lBivLogicalFactor <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivLogicalFactorWide(y,x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivLogicalFactorLong(y,x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}



#######################################################
########### Functions for Logical ~ Ordered ###########
#######################################################

r2lBivLogicalOrderedWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Logical~Ordered (Wide)",out))
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

r2lBivLogicalOrderedMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Logical~Ordered (Mixed)",out))
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
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Student","Wilcoxon"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}



r2lBivLogicalOrderedLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Logical~Ordered (Wide)",out))
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
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Student","Wilcoxon"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}


r2lBivLogicalOrdered <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivLogicalOrderedWide(y,x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivLogicalOrderedLong(y,x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


#######################################################
########### Functions for logical ~ discrete ##########
#######################################################

r2lBivLogicalDiscreteWide <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Logical~Discrete (Wide)",out))
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
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Student","Wilcoxon"),line=c(T,F,T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}
#r2lBivLogicalDiscreteWide(concours,nbRedoublement,graphDir="graphBiv",graphName="Vc1",out="latex")


r2lBivLogicalDiscreteMixed <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Logical~Discrete (Mixed)",out))
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
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Student","Wilcoxon"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}


r2lBivLogicalDiscreteLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex"){
    cat(r2lComment("rtlb Logical~Discrete (Wide)",out))
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

   # For the second column of the first line
    ## attrs <- switch(out,
    ##     "html" = " align='center' cellpadding=5",
    ##     "latex" = "c"
    ## )
    ## columnTwo <- r2lStartTable(attrs=attrs,hline=FALSE,out=out)


    ## columnTwo <- paste(columnTwo,
    ##                    r2lBuildRow(r2lBivSummary(x,as.factor(y),out=out),hline=TRUE,out=out),
    ##                    r2lBuildColumnTitle("Mosaic",hline=FALSE,out=out,border=""),
    ##                    r2lBuildRow(r2lGraphMosaicPlot(y,x,graphDir,graphName,type,out=out),out=out,border="",hline=FALSE),
    ##                    r2lEndTable(out))

    ## # First line : Table Summary / Mosaic
    ## cat(r2lBuildColumnTitle(c("Table","Summary"),hline=FALSE,out=out))
    ## cat(r2lBuildRow(c(r2lBivContingencyTable(y,x,out=out),
    ##                   columnTwo),out=out))

    ## # Second line : Barplot, Mosaic
    ## cat(r2lBuildColumnTitle(c("Barplot"),span=2,hline=FALSE,out=out))
    ## cat(r2lBuildRow(r2lGraphBarplot(y,x,graphDir,graphName,type,out=out),span=2,out=out))

    # Third line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Khi2","FisherExact","Student","Wilcoxon"),line=c(T,F,T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}



r2lBivLogicalDiscrete <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
    if(displayStyle=="wide"){
        r2lBivLogicalDiscreteWide(y,x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }else{
        r2lBivLogicalDiscreteLong(y,x,tabTitle,graphDir=graphDir,graphName=graphName,type=type,out=out)
    }
}


#######################################################
########## Functions for logical ~ continuous #########
#######################################################

r2lBivLogicalContinuous <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide"){
    cat(r2lComment("rtlb Logical~Continuous (Wide)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=3,tabSpec="|ccc|",out=out))

    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Summary","Boxplot","Density"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(x,as.factor(y),out=out),
                      r2lGraphBoxplot(x,y,graphDir,graphName,type,out=out),
                      r2lGraphDensity(x,as.factor(y),graphDir,graphName,type,out=out)
                      ),out=out))

    # Second line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=3))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Student","Wilcoxon"),line=c(T,F,F),out=out), span=3, out=out))

    cat(r2lEndStruct(out))
}

### Non utilisé
r2lBivLogicalContinuousLong <- function(y,x,tabTitle,graphDir="graphBiv",graphName="V",type="png",out="latex",dsplayStyle="wide"){
    cat(r2lComment("rtlb Logical~Continuous (Long)",out))
    cat(r2lBivBeginStruct(y,x,tabTitle,nbColumn=2,tabSpec="|cc|",out=out))

    # For the second column of the first line
	attrs <- switch(out,
        "html" = " align='center' cellpadding=5",
        "latex" = "c"
	)
    twoGraph <- r2lStartTable(attrs=attrs,hline=FALSE,out=out)
    twoGraph <- paste(twoGraph,
                      r2lBuildRow(r2lGraphBoxplot(x,y,graphDir,graphName,type,out=out),out=out),
                      r2lBuildColumnTitle("Densities",hline=FALSE,out=out),
                      r2lBuildRow(r2lGraphDensity(x,as.factor(y),graphDir,graphName,type,out=out),hline=FALSE,out=out),
                      r2lEndTable(out))


    # First line : Table Barplot Mosaic
    cat(r2lBuildColumnTitle(c("Summary","Boxplot"),hline=FALSE,out=out))
    cat(r2lBuildRow(c(r2lBivSummary(x,as.factor(y),out=out),
                      twoGraph
                      ),out=out))

    # Second line : tests
    cat(r2lBuildColumnTitle("Tests",hline=FALSE,out=out,span=2))
    cat(r2lBuildRow(r2lBivTest(x,y,test=c("Student","Wilcoxon"),line=c(T,F,F),out=out), span=2, out=out))

    cat(r2lEndStruct(out))
}

#r2lBivLogicalContinuous <- function(y,x,graphDir="graphBiv",graphName="V",type="png",out="latex",displayStyle="wide") {
#    if(displayStyle=="wide"){
#        r2lBivLogicalContinuousWide(y=y,x=x,graphDir=graphDir,graphName=graphName,type=type,out=out)
#    }else{
#        r2lBivLogicalContinuousLong(y=y,x=x,graphDir=graphDir,graphName=graphName,type=type,out=out)
#    }
#}
