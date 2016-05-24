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

dOFV.vs.cov <-
    function(xpdb1,
             xpdb2,
             covariates=xvardef("covariates",xpdb1),
             #sig.drop=-3.84,
             ##decrease.label.number=3,
             ##increase.label.number=3,
             ##id.lab.cex=0.6,
             ##id.lab.pos=2,
             #type="p",
             #xlb="Covariate",
             ylb=expression(paste(Delta, OFV[i])),
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
             smooth=TRUE,
             abline=c(0,0),
             ablcol="grey",
             abllwd=2,
             abllty="dashed",
             max.plots.per.page=1,
             ...)
{


  if(is.null(xpdb2)){
    cat("Comparison database needed for this plot!")
    return(NULL)
  }

    #require(graphics)

    iv1 <- xpdb1@Data.firstonly
    ##str(iv2)

    iv2 <- xpdb2@Data.firstonly
    ##str(iv2)

    if(!all(iv1$ID == iv2$ID)){
        cat("All ID labels for both databases must match\n")
        return(NULL)
    }

    comp.frame <- data.frame(id=iv1$ID,obj1=iv1$OBJ,obj2=iv2$OBJ)
    comp.frame$d.obj <- comp.frame$obj2 - comp.frame$obj1
    ##str(comp.frame)
    ##hist(comp.frame$d.obj)

    ## add to dataframe


    xpdb1@Data[!duplicated(xpdb1@Data[,xvardef("id",xpdb1)]),"d.obj"] <- comp.frame$d.obj
    xpdb1@Prefs@Labels$d.obj <- "Change in individual OFV"
    ##    str(xpdb1)
    ##xpdb1@Data[1:10,]

    ##covariate


  ## create list for plots
  
  cat.covs <- c()
  cont.covs <- c()
  for (i in covariates) {
    if(is.factor(xpdb1@Data[[i]])){
      cat.covs <- c(cat.covs,i)
    }else{
      cont.covs <- c(cont.covs,i)
    }
  }

  #number.of.plots <- 0
  #if(length(cat.covs)>0)   number.of.plots <- number.of.plots +1
  #if(length(cont.covs)>0)   number.of.plots <- number.of.plots +1


  cov.list <- list(cont.covs=cont.covs,cat.covs=cat.covs)
  if(is.null(cat.covs))   cov.list$cat.covs <- NULL
  if(is.null(cont.covs))   cov.list$cont.covs <- NULL
  
  plotList <- vector("list",length(cov.list))
  plot.num <- 0 # initialize plot number

  for (j in 1:length(cov.list)) {

    xplot <- xpose.plot.default(cov.list[[j]],
                                "d.obj",
                                xpdb1,
                                onlyfirst=T,
                                inclZeroWRES=T,
                                smooth=smooth,
                                abline=c(0,0),
                                ablcol=ablcol,
                                abllwd=abllwd,
                                abllty=abllty,
                                ylb=ylb,
                                main=NULL,
                                pass.plot.list=TRUE,
                                ...)
    plot.num <- plot.num+1
    plotList[[plot.num]] <- xplot
  }
  
  default.plot.title <- paste("Individual change in OFV vs. Covariate(s)\n",
                              "(Run",xpdb2@Runno," - Run",xpdb1@Runno,")",sep="" )
  
  plotTitle <- xpose.multiple.plot.title(object=xpdb1,
                                         plot.text = default.plot.title,
                                         main=main,
                                         no.runno=T,
                                         ...)
  obj <- xpose.multiple.plot(plotList,plotTitle,max.plots.per.page=max.plots.per.page,
                             ...)
  return(obj)
}
