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

dOFV1.vs.dOFV2 <-
    function(xpdb1,
             xpdb2,
             xpdb3,
             ##covariates=xvardef("covariates",xpdb1),
             #sig.drop=-3.84,
             ##decrease.label.number=3,
             ##increase.label.number=3,
             ##id.lab.cex=0.6,
             ##id.lab.pos=2,
             #type="p",
             #xlb="Covariate",
             ylb=expression(paste(Delta, OFV1[i])),
             xlb=expression(paste(Delta, OFV2[i])),
             main="Default",
             #sig.line.col = "red",
             #sig.line.lty = "dotted",
             ##tot.line.col = "grey",
             ##tot.line.lty = "dashed",
             #key=list(#title = expression(paste("Individual influence on ",Delta,"OFV")),
             #       columns = 1,
             #       lines = list(pch = c(super.sym$pch[1:2],NA,NA),
             #       type=list("o","o","l","l"),
             #       col = c(super.sym$col[1:2],sig.line.col,tot.line.col),
             #       lty = c(super.sym$lty[1:2],sig.line.lty,tot.line.lty)
             #       ),
                    ## points = list(pch = c(super.sym$pch[1:2],NA),
                    ## col = c(super.sym$col[1:2],"red")),
             #       text = list(c(
             #       expression(paste(Delta, OFV[i] < 0)),
             #       expression(paste(Delta, OFV[i] > 0)),
             #       expression(paste("Significant  ",Delta, OFV)),
             #       expression(paste("Total  ",Delta, OFV))
                    ##"Individuals with OFV drop",
                    ##"Individuals with OFV increase",
                    ##"Significant drop in ")
             #       )),
                    ##space="right",
             #       corner=c(0.95,0.5),border=T
             #),
             smooth=NULL,
             abline=c(0,1),
             ablcol="grey",
             abllwd=2,
             abllty="dashed",
             lmline=TRUE,
             ##max.plots.per.page=1,
             ...)
{


  if(is.null(xpdb2)){
    cat("1st Comparison database needed for this plot!")
    return(NULL)
  }
  if(is.null(xpdb3)){
    cat("2nd Comparison database needed for this plot!")
    return(NULL)
  }

  

  iv1 <- xpdb1@Data.firstonly
  ##str(iv2)

  iv2 <- xpdb2@Data.firstonly
  ##str(iv2)

  iv3 <- xpdb3@Data.firstonly

  if(!all(iv1$ID == iv2$ID)){
    cat("All ID labels for 1st and 2nd databases must match\n")
    return(NULL)
  }
  if(!all(iv1$ID == iv3$ID)){
    cat("All ID labels for 1st and 3rd databases must match\n")
    return(NULL)
  }

  comp.frame <- data.frame(id=iv1$ID,obj1=iv1$OBJ,obj2=iv2$OBJ,obj3=iv3$OBJ)
  comp.frame$d.obj1 <- comp.frame$obj2 - comp.frame$obj1
  comp.frame$d.obj2 <- comp.frame$obj3 - comp.frame$obj1
  ##str(comp.frame)
  ##hist(comp.frame$d.obj)
  
  ## add to dataframe


  xpdb1@Data[!duplicated(xpdb1@Data[,xvardef("id",xpdb1)]),"d.obj1"] <- comp.frame$d.obj1
  xpdb1@Data[!duplicated(xpdb1@Data[,xvardef("id",xpdb1)]),"d.obj2"] <- comp.frame$d.obj2
  xpdb1@Prefs@Labels$d.obj1 <- "Change in individual OFV 1"
  xpdb1@Prefs@Labels$d.obj2 <- "Change in individual OFV 2"
    ##    str(xpdb1)
    ##xpdb1@Data[1:10,]

    ##covariate


  ## create list for plots
  
  ## cat.covs <- c()
  ## cont.covs <- c()
  ## for (i in covariates) {
  ##   if(is.factor(xpdb1@Data[[i]])){
  ##     cat.covs <- c(cat.covs,i)
  ##   }else{
  ##     cont.covs <- c(cont.covs,i)
  ##   }
  ## }
  
  ##Number.of.plots <- 0
  ##if(length(cat.covs)>0)   number.of.plots <- number.of.plots +1
  ##if(length(cont.covs)>0)   number.of.plots <- number.of.plots +1

  #cov.list <- list(cont.covs=cont.covs,cat.covs=cat.covs)
  
  #plotList <- vector("list",length(cov.list))
  #plot.num <- 0 # initialize plot number

  #for (j in 1:length(cov.list)) {

  default.plot.title <- paste("Individual change in OFV1 vs. OFV2\n",
                              "(OFV 1 = Run",xpdb2@Runno," - Run",xpdb1@Runno,")\n",
                              "(OFV 2 = Run",xpdb3@Runno," - Run",xpdb1@Runno,")",
                              sep="" )
  
  plotTitle <- xpose.multiple.plot.title(object=xpdb1,
                                         plot.text = default.plot.title,
                                         main=main,
                                         no.runno=T,
                                         ...)

  xplot <- xpose.plot.default("d.obj2",
                              "d.obj1",
                              xpdb1,
                              onlyfirst=T,
                              inclZeroWRES=T,
                              smooth=NULL,
                              lmline=lmline,
                              abline=abline,
                              ablcol=ablcol,
                              abllwd=abllwd,
                              abllty=abllty,
                              ylb=ylb,
                              xlb=xlb,
                              main=plotTitle,
                                        #pass.plot.list=TRUE,
                              ...)

  return(xplot)
}
