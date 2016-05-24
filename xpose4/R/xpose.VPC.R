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

"xpose.VPC" <-
  function(vpc.info="vpc_results.csv",  #name of PSN file to use
           vpctab = dir(pattern="^vpctab")[1],
           object = NULL,
           ids=FALSE,
           type="p",
           by=NULL,
           PI=NULL,#"area",
           PI.ci="area",
           PI.real=T,
           PI.ci.med.arcol="red",
           subset=NULL,
           main="Default",
           main.sub=NULL,  # used for names above each plot when using multiple plots
                                        #Should be a vector c("","")
           main.sub.cex=0.85, # size of main.sub 
           inclZeroWRES=FALSE,
           force.x.continuous=FALSE,
           #strip="Default",
           #dont.plot=F,
           funy=NULL,
           logy=FALSE,
           ylb = "Default",
           verbose=FALSE,
           ...) {

    ## for testing
    ##vpctab="./vpc_strat_WT_4_mirror_5/vpctab5"
    ##vpctab="./vpc_strat_SEX_mirror_5/vpctab5"
    ##object <- xpdb
    ##inclZeroWRES <- FALSE
    
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
                            verbose=verbose,
                            ...)
      if(tmp==TRUE) inclZeroWRES=TRUE
    }

    file.info <- read.npc.vpc.results(vpc.results=vpc.info,verbose=verbose,...)
    num.tables <- file.info$num.tables
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
  
    if(is.factor(object@Data[[dv.var]])){
      change.cat.cont(object) <- c(dv.var)
    }
    
    if(force.x.continuous){
      if(is.factor(object@Data[[idv.var]])){
        change.cat.cont(object) <- c(idv.var)
      }
    }

    ## decide on the conditioning
    if (is.null(by) && num.tables!=1){
      ## get conditioning veriable name
      # for future use to automatically start conditioning
      #for (i in 1:num.tables){
      #  tmp.strata <- strata.names[i]
      #  strata.loc <- regexpr(strata.start.pat,strata.line)+7
      #  strata.names <- c(strata.names,substring(strata.line,strata.loc))
      #}

      ## use subsetting to get things working
      if(!is.null(subset)){ # this can be fixed below
        if(verbose) cat(paste("Overwriting the subset expression to handle multiple STRATA\n"))
      }

      plotList <- vector("list",num.tables)
      plot.num <- 0 # initialize plot number
      

      ## this can be updated as in npc.coverage.R
      for (i in 1:num.tables){ 
        ##subset <- file.info$result.tables[[num.tables+1]][i] # this can be fixed to aviod overwriting subsets
        subset <- file.info$strata.names[i] # this can be fixed to aviod overwriting subsets 
        final.bin.table <- file.info$result.tables[[i]]
        if(!is.null(main.sub)){
          sub.main=main.sub[i]
        } else {
          sub.main=subset
        }

        if(!is.character(ylb)){
        } else if(ylb != "Default"){
        } else {
          tmp.label <- xpose.create.label(dv.var,
                                          object,
                                          funy,
                                          logy,...)
        
          if(file.info$pred.corr && !file.info$var.corr){
            tmp.label <- paste(tmp.label,"\n(Pred Corr)")
          }
          if(file.info$pred.corr && file.info$var.corr){
            tmp.label <- paste(tmp.label,"\n(Pred and Var Corr)")
          }
          ylb=tmp.label
        }

        ## make the VPC
        xplot <- xpose.plot.default(idv.var,#xvardef("idv",object),
                                    dv.var,#xvardef("dv",object),
                                    object,
                                    ids=ids,
                                    type=type,
                                    subset=subset,
                                    PI=PI,
                                    PI.ci=PI.ci,
                                    PI.real=PI.real,
                                    PI.ci.med.arcol=PI.ci.med.arcol,
                                    PI.bin.table=final.bin.table,
                                    pass.plot.list=TRUE,
                                    main=sub.main,
                                    main.cex=main.sub.cex,
                                    inclZeroWRES=inclZeroWRES,
                                    ylb = ylb,
                                    funy=funy,
                                    logy=logy,
                                    ...)
        plot.num <- plot.num+1
        plotList[[plot.num]] <- xplot
      }
        
      default.plot.title <- "Visual Predictive Check\n"
      if(file.info$pred.corr && !file.info$var.corr){
        default.plot.title <- "Visual Predictive Check\n (Prediction Corrected)\n"
      }
      if(file.info$pred.corr && file.info$var.corr){
        default.plot.title <- "Visual Predictive Check\n (Prediction and Variance Corrected)\n"
      }

      default.plot.title <- paste(default.plot.title,
                                  xpose.create.title(idv.var,dv.var,object,
                                                     no.runno=T,...),sep="")
      plotTitle <- xpose.multiple.plot.title(object=object,
                                             plot.text = default.plot.title,
                                             main=main,
                                             #subset=subset,
                                             ...)

#      if(!dont.plot){
#        xpose.multiple.plot.default(plotList,plotTitle=plotTitle,...)
#      }
      obj <- xpose.multiple.plot(plotList,plotTitle,...)
#      return(invisible(plotList))
      return(obj)
      
    } else { ## either plot stratification with by or only one strata 
      ## check structure of stratification variable
      if(!is.null(by) && num.tables!=1){
        if(all(is.null(file.info$by.interval))){
          ## categorical variable
          if(!is.factor(object@Data[[by]])) change.cat.cont(object) <- by
        } else {
          ## continuous variable
          if(is.factor(object@Data[[by]])) change.cat.cont(object) <- by
        }
      }

      default.plot.title <- "Visual Predictive Check\n"
      if(file.info$pred.corr && !file.info$var.corr){
        default.plot.title <- "Visual Predictive Check\n (Prediction Corrected)\n"
      }
      if(file.info$pred.corr && file.info$var.corr){
        default.plot.title <- "Visual Predictive Check\n (Prediction and Variance Corrected)\n"
      }
      default.plot.title <- paste(default.plot.title,
                                  xpose.create.title(idv.var,dv.var,object,
                                                     no.runno=T,subset=subset,...),sep="")
      plotTitle <- xpose.multiple.plot.title(object=object,
                                             plot.text = default.plot.title,
                                             main=main,
                                             subset=subset,
                                             ...)

      if(!is.character(ylb)){
      } else if(ylb != "Default"){
      } else {
        tmp.label <- xpose.create.label(dv.var,
                                        object,
                                        funy,
                                        logy,...)
        
        if(file.info$pred.corr && !file.info$var.corr){
          tmp.label <- paste(tmp.label,"\n(Pred Corr)")
        }
        if(file.info$pred.corr && file.info$var.corr){
          tmp.label <- paste(tmp.label,"\n(Pred and Var Corr)")
        }
        ylb=tmp.label
      }

      ## make the VPC
      xplot <- xpose.plot.default(idv.var,#xvardef("idv",object),
                                  dv.var,#xvardef("dv",object),
                                  object,
                                  ids=ids,
                                  type=type,
                                  by=by,
                                  subset=subset,
                                  PI=PI,
                                  PI.ci=PI.ci,
                                  PI.real=PI.real,
                                  PI.ci.med.arcol=PI.ci.med.arcol,
                                  PI.bin.table=file.info$result.tables,
                                  #force.by.factor=TRUE,
                                  main=plotTitle,
                                  by.interval=file.info$by.interval,
                                  inclZeroWRES=inclZeroWRES,
                                  ylb = ylb,#tmp.label,
                                  funy=funy,
                                  logy=logy,
                                  ...)
      return(xplot)
    }
  }
