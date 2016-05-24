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

"xpose.VPC.both" <-
  function(vpc.info="vpc_results.csv",  #name of PSN file to use
           vpctab = dir(pattern="^vpctab")[1],
           object = NULL,
           #ids=NULL,
           #type="p",
           #by=NULL,
           #PI="area",
           subset=NULL,
           main="Default",
           main.sub=NULL,               #used for names above each plot
                                        #when using multiple plots
                                        #Should be a vector c("","")
           ##main.sub.cex=0.85, # size of main.sub 
           inclZeroWRES=FALSE,
           ##plot.cont.table=TRUE,
                                        #plot.cat.table=TRUE,
           #force.x.continuous=FALSE,
           #real.col=4,
           #median.line=FALSE,
           #median.col="darkgrey",
           #ci.lines=FALSE,
           #ci.col="blue",
           #ci.lines.col="darkblue",
           #xlb="Default",
           #ylb="Proportion of Total",
           #force.x.continuous=FALSE,
           #level.to.plot=NULL,
           #max.plots.per.page=1,
           #strip="Default",
           #rug=TRUE,
           #rug.col="orange",
           cont.logy=F,
           hline="default",
           add.args.cont=list(),
           add.args.cat=list(),
           ...) {


    ## Make sure we have the necessary variables defined
    if(is.null(object) & is.null(vpctab)){
      cat(paste("Both the arguments object and vpctab are NULL\n"))
      cat(paste("At least one of these must be defined\n"))
      return(NULL)
    }

    if(!is.null(vpctab)){
      tmp <- FALSE
      if(is.null(object)) tmp <- TRUE
      object.2 <- read.vpctab(vpctab=vpctab,
                              object=object,
                              inclZeroWRES=inclZeroWRES,
                              ...)
      if(tmp==TRUE) inclZeroWRES=TRUE
    }
    
    
    file.info <- read.npc.vpc.results(vpc.results=vpc.info,...)
    num.tables <- file.info$num.tables
    dv.var <- file.info$dv.var
    idv.var <- file.info$idv.var

    arg.list.1 <- list(vpc.info=vpc.info,
                       vpctab = vpctab,
                       object = object,
                       subset=subset,
                       main=NULL,
                       main.sub=main.sub,
                       #aspect="fill",
                       #xlb=NULL,
                       censored=T,
                       ...)
    arg.list <- c(arg.list.1,add.args.cat)
    cat.plots.tmp <- do.call(xpose.VPC.categorical,arg.list,quote=F)

    
    ## cat.plots.tmp <-
    ##   xpose.VPC.categorical(vpc.info=vpc.info,
    ##                         vpctab = vpctab,
    ##                         object = object,
    ##                         subset=subset,
    ##                         main=NULL,
    ##                         main.sub=NULL,
    ##                         #xlb=NULL,
    ##                         censored=T,
    ##                         ...)
    
    cat.plots <- cat.plots.tmp@plotList

    if(!is.na(match(hline,"default"))) {
      hline=c(file.info$lloq,file.info$uloq)
      if(cont.logy) hline=log10(hline)
    }
    
    if(num.tables==1) {
    cont.plots <- vector("list",1)
    arg.list.1 <- list(vpc.info=vpc.info,
                       vpctab = vpctab,
                       object = object,
                       subset=subset,
                       main=NULL,
                       aspect="fill",
                       xlb=NULL,
                       hline=hline,
                       logy=cont.logy,
                       ...)
    arg.list <- c(arg.list.1,add.args.cont)
    cont.plots[[1]] <- do.call(xpose.VPC,arg.list,quote=F)

    ## cont.plots[[1]] <-
    ##   xpose.VPC(vpc.info=vpc.info,
    ##             vpctab = vpctab,
    ##             object = object,
    ##             subset=subset,
    ##             main=NULL,
    ##             aspect="fill",
    ##             xlb=NULL,
    ##             hline=c(file.info$lloq,file.info$uloq),
    ##             add.args.cont,
    ##             ...)
    } else {
      arg.list.1 <- list(vpc.info=vpc.info,
                         vpctab = vpctab,
                         object = object,
                         subset=subset,
                         ##main=NULL,
                         aspect="fill",
                         xlb=NULL,
                         hline=hline,
                         logy=cont.logy,
                         ...)
      arg.list <- c(arg.list.1,add.args.cont)
      cont.plots.tmp <- do.call(xpose.VPC,arg.list,quote=F)
      
      ## cont.plots.tmp <-
      ##   xpose.VPC(vpc.info=vpc.info,
      ##             vpctab = vpctab,
      ##             object = object,
      ##             subset=subset,
      ##             aspect="fill",
      ##             xlb=NULL,
      ##             hline=c(file.info$lloq,file.info$uloq),
      ##             ...)
      
      cont.plots <- cont.plots.tmp@plotList
    }


    plotList <- vector("list",num.tables*2)
    j=0
    for(i in 1:num.tables){
      plotList[[i+j]] <- cont.plots[[i]]
      plotList[[i+j+1]] <- cat.plots[[i]]
      j=j+1
    }
    if(main=="Default"){
      no.runno <- FALSE
      text2 <- NULL
      if(object.2@Runno=="0"){
        no.runno <- TRUE
        text2 <- paste("\n(",file.info$model.file,")",sep="")
      }
      main <- xpose.create.title.text(NULL,dv.var,
                                      "VPC for",
                                        text2=text2,
                                      no.runno=no.runno,
                                      object.2,subset=subset,...)
    }
    
                                        #print(cont.plots[[i]],position=c(0,0.2,1,1),more=TRUE)
                                        #  print(cat.plots[[i]],position=c(0,0,1,0.33),more=TRUE)
    obj <- xpose.multiple.plot(plotList,plotTitle=main,bql.layout=T,...)
    return(obj)


  }
