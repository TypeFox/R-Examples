# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.
randtest.hist <-
    function(results.file="raw_results_run1.csv",
             df=1,
             p.val=0.05,
             main="Default",
             xlim=NULL,
             PCTSlcol = "black",
             vlcol=c("red","orange"),
             ...)
{
    ##rm(list=ls())
    ##p.val <- 0.05
    ##df <- 2
    ##main="Default"
    ##results <- read.csv("randtest_dir1/raw_results_run1.csv")
    
    crit.val.nom <- qchisq(1-p.val, df=df)

    ##Read in all results
    results <- read.csv(results.file)

    ## check that classes are present
    createXposeClasses()
    
    ## Create the object
    xpobj       <- new("xpose.data",
                       Runno="PsN Randomization Test",
                       Data = NULL)
    
    ## read local options
    if (is.readable.file("xpose.ini")) {
        xpobj <- xpose.read(xpobj, file="xpose.ini")
    } else {
        ## read global options
        rhome   <- R.home()
        xdefini <- paste(rhome, "\\library\\xpose4\\xpose.ini", sep="")
        if (is.readable.file(xdefini)) {
            xpobj <- xpose.read(xpobj, file=xdefini)
        }else{
            xdefini2 <- paste(rhome, "\\library\\xpose4\\xpose.ini", sep="")
            if (is.readable.file(xdefini2)) {
                xpobj <- xpose.read(xpobj, file=xdefini2)
            } 
        }
    }

    results$ID <-1
    results$WRES <- 1
    Data(xpobj) <- results[-c(1,2),]
    
    crit.val.emp <- quantile(xpobj@Data$deltaofv,probs=p.val)    

    orig = results$deltaofv[2]       #dOFV for original dataset

    xpose.plot.histogram("deltaofv",
                         xpobj,
                         bins.per.panel.equal=FALSE,
                                        #layout=layout,
                                        #vdline=if(showOriginal){c(o1[all.index])} else {NULL},
                                        #showMean=showMean,
                                        #showMedian=showMedian,
                         xlim = if(is.null(xlim)){c(min(c(orig,crit.val.emp,xpobj@Data$deltaofv))-1,
                             max(c(orig,crit.val.emp,xpobj@Data$deltaofv,0))+0.2)},
                         showPCTS=TRUE,
                         PCTS=c(p.val),
                         PCTSlcol = PCTSlcol,
                         vline=c(orig,-crit.val.nom),
                         vlcol=vlcol,
                         main=if(main=="Default"){"Change in OFV for Randomization Test"}else{main},
                         key=list(#title = "Critical value lines",
                             columns = 1,
                             lines = list(type="l",col =c(vlcol,PCTSlcol)),#,lty=c(vlty,PCTSlty)),
                             #lines = list(type="l",col =c("red","orange","black")),
                             text = list(c("Original data", "Nominal", "Empirical (rand. test)")),
                             ##space="right",
                             corner=c(0.05,0.95),
                             border=T,
                             #transparent=FALSE,
                             alpha.background=1,
                             background = "white"
                             ),
                         ...
                         )

}
                     



