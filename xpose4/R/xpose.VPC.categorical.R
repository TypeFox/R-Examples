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

"xpose.VPC.categorical" <-
  function(vpc.info="vpc_results.csv",  #name of PSN file to use
           vpctab = dir(pattern="^vpctab")[1],
           object = NULL,
           #ids=NULL,
           #type="p",
           #by=NULL,
           #PI="area",
           subset=NULL,
           main="Default",
           main.sub="Default",  # used for names above each plot when using multiple plots
                                        #Should be a vector c("","")
           main.sub.cex=0.85, # size of main.sub
           #inclZeroWRES=FALSE,
           #plot.cont.table=TRUE,
           #plot.cat.table=TRUE,
           #force.x.continuous=FALSE,
           real.col=4,
           real.lty="b",
           real.cex=1,
           real.lwd=1,
           median.line=FALSE,
           median.col="darkgrey",
           median.lty=1,
           ci.lines=FALSE,
           ci.col="blue",
           ci.lines.col="darkblue",
           ci.lines.lty=3,
           xlb="Default",
           ylb="Proportion of Total",
           force.x.continuous=FALSE,
           level.to.plot=NULL,
           max.plots.per.page=1,
           #strip="Default",
           rug=TRUE,
           rug.col="orange",
           censored=FALSE,
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
      object <- read.vpctab(vpctab=vpctab,
                            object=object,
                            inclZeroWRES=inclZeroWRES,
                            ...)
      if(tmp==TRUE) inclZeroWRES=TRUE
    }

    file.info <- read.npc.vpc.results(vpc.results=vpc.info,...)
    ##num.tables <- file.info$num.tables
    dv.var <- file.info$dv.var
    idv.var <- file.info$idv.var
    ##bin.table <- file.info$result.tables

    tmp <- c()
    if(is.null(object@Data[[dv.var]])) tmp <- c(tmp,dv.var)
    if(is.null(object@Data[[idv.var]])) tmp <- c(tmp,idv.var)
    if (!is.null(tmp)){
      cat("\n-----------Variable(s) not defined!-------------\n",
          tmp, "is/are not defined in the current database\n",
          "and must be defined for this command to work!\n",
          "------------------------------------------------\n")
      return(NULL)
    }

    ## check if the tables are present
    tables.exist <- FALSE
    if(censored){
        if(!is.null(file.info$num.tables.cen)) tables.exist <- TRUE
    } else {
        if(!is.null(file.info$num.tables.cat)) tables.exist <- TRUE
    }
    if(!tables.exist){
      cat("\n-----------Tables do not exist!-------------\n",
          "There are no tables that correspond to the categorical plots\n",
          "you would like to make.  Please rerun the vpc calculations.\n",
          "------------------------------------------------\n")
      return(NULL)
    }

    ## make the plots
    level.names <- c()
    if(censored){
      num.tables <- file.info$num.tables.cen
      plotList <- vector("list",num.tables)
      result.tables <- file.info$result.tables.cen
      n.levs <- 0
      if(!is.na(file.info$lloq)){
        level.names <- c(level.names,"LLOQ")
        n.levs <- n.levs+1
      }
      if(!is.na(file.info$uloq)) {
        level.names <- c(level.names,"ULOQ")
        n.levs <- n.levs+1
      }
    } else {
      num.tables <- file.info$num.tables.cat
      plotList <- vector("list",num.tables)
      result.tables <- file.info$result.tables.cat
      n.levs <- length(file.info$cat.boundaries)+1
      for(LEVS in 1:n.levs){
        if(LEVS==1){
          tmp.lev=paste(dv.var,"<=",file.info$cat.boundaries[[LEVS]],sep=" ")
        } else {
          if(LEVS==n.levs){
            tmp.lev=paste(file.info$cat.boundaries[[LEVS-1]],"<",dv.var,sep=" ")
          } else {
            tmp.lev=paste(file.info$cat.boundaries[[LEVS-1]],"<",
              dv.var,"<=",file.info$cat.boundaries[[LEVS]],sep=" ")
          }
        }
        level.names <- c(level.names,tmp.lev)
      }
    }



    plot.num <- 0 # initialize plot number
    for (i in 1:num.tables){
      if(num.tables==1){
        tmp.table <- result.tables
      } else {
        tmp.table <- result.tables[[i]]
      }



      ## Set up the data frame for a VPC
      tmp.table.2 <- rbind(tmp.table)
      tmp.table.2$idv <- rowMeans(tmp.table.2[c("upper","lower")], na.rm = TRUE, dims = 1)
      num.col.new <- 6
      n.idv.levs <- length(tmp.table.2[,"idv"])
      num.row.new <- n.levs*n.idv.levs
      ret.new   <- data.frame(matrix(nrow = num.row.new,
                                     ncol = num.col.new))
      names(ret.new) <- c("idv","real","lower","median","upper","by.var")

      tab.names <- names(tmp.table.2)
      real.index <- grep("Real.",tab.names)
      lower.index <- grep(".from",tab.names)
      upper.index <- grep(".to",tab.names)
      median.index <- grep("Sim.",tab.names)
      idv.index <- grep("idv",tab.names)

      ## Here is the problem when we have intervals of IDV
      ##         tab.names <- names(tmp.table)
      ##         real.index <- grep("Real.",tab.names)
      ##         lower.index <- grep(".from",tab.names)
      ##         upper.index <- grep(".to",tab.names)
      ##         median.index <- grep("Sim.",tab.names)

      ##         PPI <- tmp.table[c(lower.index[[1]],
      ##                            upper.index[[i]],
      ##                            median.index[[1]])]
      ##         names(PPI) <- c("upper","lower","median")
      ##         PPI$Xupper <- tmp.table$upper
      ##         PPI$Xlower <- tmp.table$lower

      ##         get.polygon.regions(PPI,NULL)

      for(LEVS in 1:n.levs){
        ret.new[(1+(LEVS-1)*n.idv.levs):(n.idv.levs*LEVS),1:5]<-
          tmp.table.2[,c(idv.index,
                         real.index[[LEVS]],
                         lower.index[[LEVS]],
                         median.index[[LEVS]],
                         upper.index[[LEVS]]
                         )]
        ret.new[(1+(LEVS-1)*n.idv.levs):(n.idv.levs*LEVS),"by.var"] <-
          level.names[LEVS]
      }

      ## check if x should be categorical
      if(length(unique(ret.new$idv))<= object@Prefs@Cat.levels) {
        if(!is.factor(ret.new$idv)){
          cat("\n  Inferring that ",idv.var," is categorical\n",sep="")
          ##cat("  Transforming ",idv.var," from continuous to categorical\n",sep="")
          tmp.levs <- unique(ret.new[,"idv"])
          tmp.levs <- tmp.levs[order(tmp.levs)]
          ret.new$idv <- factor(ret.new$idv,levels=tmp.levs,ordered=T)
        }
      }

      ## categorize the by.var
      ret.new$by.var <- factor(ret.new$by.var,levels=level.names,ordered=T)

      ## subset the by.var levels
      if(!is.null(level.to.plot)){
        ret.new <- ret.new[ret.new["by.var"]==level.names[level.to.plot],]
        panel.level <- level.names[level.to.plot]
      }

      ## set up formula
      if (n.levs==1){
        formel <- paste("real~idv|by.var",sep="")
      } else {
        formel <- paste("real~idv|by.var",sep="")
      }

      ## set up labels
      if(xlb[1]=="Default"){
        xlb <- idv.var
      }
      strata <- file.info$strata.names[i] # this can be fixed to aviod overwriting subsets
      if(is.null(main.sub)){
        sub.main=NULL
      } else {
        if(main.sub[1]=="Default"){
          sub.main=strata
        } else {
          sub.main=main.sub[i]
        }
      }
      ## if(!is.null(main.sub)){
      ##   sub.main=main.sub[i]
      ## } else {
      ##   sub.main=strata
      ## }


      ## make plot
      plot.num <- plot.num+1
      plotList[[plot.num]] <-
        xyplot(formula(formel),
               data=ret.new,
               type=real.lty,
               data2=ret.new,
               prepanel = function(x,y,subscripts,data2=data2,...) {
                 ## if(!is.null(cat.level.to.plot)){
                 ##                       tmp.data <- ret.new[ret.new$by.var==panel.level,]
                 ##                     } else {
                 ##                       tmp.data <- ret.new[ret.new$by.var==level.names[panel.number()],]
                 ##                     }
                 tmp.data <- data2[subscripts,]
                 if(is.factor(x)){#length(levs <- unique(x)) < object@Prefs@Cat.levels) {
                   xlim <- levels(x)
                 } else {
                   xlim <- range(c(x,tmp.table.2$lower,tmp.table.2$upper),na.rm=T)
                                        #xlim <- range(x)
                 }
                 ylim <- range(c(y,tmp.data$lower,tmp.data$upper))
                 list(xlim=xlim,ylim=ylim)
               },
               xlab=xlb,ylab=ylb,
               main=sub.main,
                                        #strip=strip,
               ...,
               panel=function(x,y,subscripts,data2=data2,...){
                 if(!is.null(level.to.plot)){
                   ##tmp.data <- ret.new[ret.new$by.var==panel.level,]
                   tmp.data <- data2[data2$by.var==panel.level,]
                 } else {
                   ##tmp.data <- ret.new[ret.new$by.var==level.names[panel.number()],]
                   tmp.data <- data2[data2$by.var==level.names[panel.number()],]
                 }
                 grid.polygon(c(tmp.data$idv,rev(tmp.data$idv)),
                              c(tmp.data$upper,rev(tmp.data$lower)),
                              default.units="native",
                              gp=gpar(fill=ci.col,alpha=0.3,col=NULL,lty=0)
                              )
                 panel.xyplot(x,y,col=real.col,cex=real.cex,lwd=real.lwd,...)
                 if(median.line){
                   panel.lines(tmp.data$idv,tmp.data$median,type="b",
                               col=median.col,
                               lty=median.lty)
                 }
                 if(ci.lines){
                   panel.lines(tmp.data$idv,tmp.data$lower,type="b",
                               col=ci.lines.col,
                               lty=ci.lines.lty)
                   panel.lines(tmp.data$idv,tmp.data$upper,type="b",
                               col=ci.lines.col,
                               lty=ci.lines.lty)
                 }
                 if(rug){
                   panel.rug(x=c(tmp.table.2$lower,tmp.table.2$upper),y=NULL,
                             col=rug.col, lwd=3)
                 }
               }
               )
    }

    if(!is.null(main)){
      if(main=="Default"){
        no.runno <- FALSE
        text2 <- NULL
        if(object@Runno=="0"){
          no.runno <- TRUE
          text2 <- paste("\n(",file.info$model.file,")",sep="")
        }
        if(censored) {
          main <- xpose.create.title.text(NULL,dv.var,
                                          "Categorical VPC for Censored",
                                          text2=text2,
                                          no.runno=no.runno,
                                          object,subset=subset,...)
        } else {
          main <- xpose.create.title.text(NULL,dv.var,
                                          "Categorical VPC for",
                                          text2=text2,
                                          no.runno=no.runno,
                                          object,subset=subset,...)
        }
      }
    }
    obj <- xpose.multiple.plot(plotList,plotTitle=main,
                               max.plots.per.page=max.plots.per.page,...)
    return(obj)
    ## if(!dont.plot){
    ##   xpose.multiple.plot.default(plotList,plotTitle=main,
    ##                               ...)
    ## }
    ## return(invisible(plotList))

  }
