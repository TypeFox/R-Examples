### $Id: plotOptim.R 193 2012-06-24 21:13:42Z kristl $

plotOptim <- function(results){
  ## Plot optimisation through GUI
  
  if(requireNamespace("gWidgets", quietly = TRUE)){
    if(requireNamespace("lattice", quietly = TRUE)){
      
      # Set up main window with internal container
      win <- gWidgets::gwindow("Optimisation results", width=550)
      main <- gWidgets::ggroup(horizontal=FALSE)
      gWidgets::add(win,main)
      
      # Initialize parameters and paramter lists
      nAlgs <- 0
      progress <- toPlot <- groupPlots <- parameterList <- parameterPlots <- rGroups <- groups <- labels <- plotOne <- plotTwo <- plotFlip <- param <- list()
      oneDimNames <- character(length(results$results))
      if(exists("baselineAlgorithmsGUI",envir=.GlobalEnv)){
        bAGUI <- get("baselineAlgorithmsGUI",envir=.GlobalEnv)
      } else {
        bAGUI <- baselineAlgorithmsGUI
      }
      if(exists("baselineAlgorithms",envir=.GlobalEnv)){
        bA <- get("baselineAlgorithms",envir=.GlobalEnv)
      } else {
        bA <- baselineAlgorithms
      }
      
      # Functions for parsing beween method numbers and mehod names
      methodParse2 <- character(length(bA))
      for(i in 1:length(bA)){
        methodParse2[i] <- baselineAlgorithms[[i]]@funcName
      }
      
      # Function for adding a baseline correction method, including sub-functions and variables
      addAlg <- function(nm, result){
        gWidgets::delete(groups[[i]], progress[[i]])
        # Initialize parameters and groups
        nAlgs <<- nAlgs + 1
        name <- names(bA)[nm]
        nameLong <- bA[[nm]]@description
        method <- name
        
        # Add some space and the name of the baseline correction algorithm
        gWidgets::addSpace(groups[[nAlgs]],10, horizontal=FALSE)
        gWidgets::add(groups[[nAlgs]],gWidgets::glabel(nameLong))
        
        # GUI for plotting
        groupPlots[[nAlgs]] <<- gWidgets::ggroup(horizontal=FALSE)
        
        # Function for setting up parameter plotting array in GUI
        addparameterPlots <- function(nameStr, lineNo){
          if(lineNo>1){ # Baseline parameter
            if(is.null(result@param[[nameStr]])){ # Used default value during optimization
              parameterList[[nAlgs]][[lineNo]]  <<- c(gWidgets::gradio(c('Overall min.','Min.','Avg.','All','Chosen'),horizontal=TRUE), gWidgets::gcheckbox('default', checked=TRUE))
              gWidgets::enabled(parameterList[[nAlgs]][[lineNo]][[1]]) <- FALSE
              gWidgets::tag(parameterList[[nAlgs]][[lineNo]][[1]], "default") <- TRUE
            } else {
              if(max(nchar(as.character(result@param[[nameStr]],scientific=TRUE)))>8){ # More than 8 digits => use scientific format
                parameterList[[nAlgs]][[lineNo]]  <<- c(gWidgets::gradio(c('Overall min.','Min.','Avg.','All','Chosen'),horizontal=TRUE), gWidgets::gcheckboxgroup(format(result@param[[nameStr]], scientific=TRUE, digits=3), horizontal=TRUE))
              } else {
                parameterList[[nAlgs]][[lineNo]]  <<- c(gWidgets::gradio(c('Overall min.','Min.','Avg.','All','Chosen'),horizontal=TRUE), gWidgets::gcheckboxgroup(format(result@param[[nameStr]], scientific=FALSE), horizontal=TRUE))
              }
            }
          } else { # Regression parameter
            if(is.null(result@param[[1]])){
              parameterList[[nAlgs]][[lineNo]]  <<- c(gWidgets::gradio(c('Overall min.','Min.','Avg.','All','Chosen'),horizontal=TRUE), gWidgets::gcheckbox('default', checked=TRUE))
              gWidgets::enabled(parameterList[[nAlgs]][[lineNo]][[1]]) <- FALSE
              gWidgets::tag(parameterList[[nAlgs]][[lineNo]][[1]], "default") <- TRUE
            } else {
              if(max(nchar(as.character(result@param[[1]])))){
                parameterList[[nAlgs]][[lineNo]]  <<- c(gWidgets::gradio(c('Overall min.','Min.','Avg.','All','Chosen'),horizontal=TRUE, selected=4), gWidgets::gcheckboxgroup(format(result@param[[1]], scientific=TRUE, digits=3), horizontal=TRUE))
              } else {
                parameterList[[nAlgs]][[lineNo]]  <<- c(gWidgets::gradio(c('Overall min.','Min.','Avg.','All','Chosen'),horizontal=TRUE, selected=4), gWidgets::gcheckboxgroup(format(result@param[[1]], scientific=FALSE), horizontal=TRUE))
              }
            }
          }
          
          # Check which buttons should be available for current choices
          gWidgets::tag(parameterList[[nAlgs]][[lineNo]][[1]], "no") <- lineNo
          gWidgets::tag(parameterList[[nAlgs]][[lineNo]][[1]], "alg") <- nAlgs
          gWidgets::tag(parameterList[[nAlgs]][[lineNo]][[1]], "name") <- nameStr
          gWidgets::addhandlerchanged(parameterList[[nAlgs]][[lineNo]][[1]], handler = function(h,...){
            linNo <- gWidgets::tag(h$obj)$no
            nAlg <- gWidgets::tag(h$obj)$alg
            if(gWidgets::svalue(h$obj,index=TRUE) == 1){
              gWidgets::svalue(parameterList[[nAlg]][[linNo]][[2]],index=TRUE) <- NULL}
            else {
              gWidgets::svalue(parameterList[[nAlg]][[linNo]][[2]],index=TRUE) <- 1:length(result@param[[nameStr]])}
            if(gWidgets::svalue(h$obj,index=TRUE) == 1 || gWidgets::svalue(h$obj,index=TRUE) == 4){
              gWidgets::enabled(parameterList[[nAlg]][[linNo]][[2]]) <- FALSE}
            else{
              gWidgets::enabled(parameterList[[nAlg]][[linNo]][[2]]) <- TRUE}
            checkPlot(nAlg)
          })
          gWidgets::tag(parameterList[[nAlgs]][[lineNo]][[2]], "no") <- lineNo
          gWidgets::tag(parameterList[[nAlgs]][[lineNo]][[2]], "alg") <- nAlgs
          gWidgets::tag(parameterList[[nAlgs]][[lineNo]][[2]], "name") <- nameStr
          if(lineNo!=1)
            gWidgets::enabled(parameterList[[nAlgs]][[lineNo]][[2]]) <- FALSE
          if(lineNo==1)
            gWidgets::svalue(parameterList[[nAlgs]][[lineNo]][[2]],index=TRUE) <- 1:length(result@param[[nameStr]])
          gWidgets::addhandlerchanged(parameterList[[nAlgs]][[lineNo]][[2]], handler = function(h,...){
            linNo <- gWidgets::tag(h$obj)$no
            nAlg <- gWidgets::tag(h$obj)$alg
            checkPlot(nAlg)
          })
          
          # Set up visual optimization array for parameters
          parameterPlots[[nAlgs]][lineNo+1,1] <<- gWidgets::glabel(text=paste(nameStr,":",sep=""))
          parameterPlots[[nAlgs]][lineNo+1,2] <<- parameterList[[nAlgs]][[lineNo]][[1]]
          parameterPlots[[nAlgs]][lineNo+1,3] <<- parameterList[[nAlgs]][[lineNo]][[2]]
        }
        # Set up visual optimization array for parameters
        parameterList[[nAlgs]] <<- list()
        parameterPlots[[nAlgs]] <<- glayout(homogeneous = FALSE, spacing = 5, container=groupPlots[[nAlgs]])
        parameterPlots[[nAlgs]][1,1] <<- gWidgets::glabel("")
        parameterPlots[[nAlgs]][1,2] <<- gWidgets::glabel("Which:")
        parameterPlots[[nAlgs]][1,3] <<- gWidgets::glabel("")
        nameStrs <- c(names(result@param[1]), rownames(bA[[method]])) ##### ##### #### ####### ##### ####
        lns <- length(nameStrs)
        for(i in 1:lns)
          addparameterPlots(nameStrs[i],i)
        
        gWidgets::addSpace(groups[[nAlgs]], 10, horizontal=FALSE)
        gWidgets::add(groups[[nAlgs]], groupPlots[[nAlgs]], expand=FALSE)
        
        # Plot buttons
        plotOne[[nAlgs]] <<- gWidgets::gbutton("Curve plot")
        plotTwo[[nAlgs]] <<- gWidgets::gbutton("Level plot")
        gWidgets::tag(plotOne[[nAlgs]], "alg") <- nAlgs
        gWidgets::tag(plotTwo[[nAlgs]], "alg") <- nAlgs
        
        # Plot curves
        gWidgets::addhandlerchanged(plotOne[[nAlgs]], handler = function(h,...){
          nAlg <- gWidgets::tag(h$obj)$alg
          if(is.vector(toPlot[[nAlg]])){ # Single curve
            plot(names(toPlot[[nAlg]]),toPlot[[nAlg]], type='l', ylab=results$results[[nAlg]]@qualMeasName, xlab=oneDimNames[nAlg])
          } else { # Multiple curves
            if(gWidgets::svalue(plotFlip[[nAlg]])==FALSE){
              plot(rownames(toPlot[[nAlg]]),toPlot[[nAlg]][,1], type='l', xlab=names(dimnames(toPlot[[nAlg]]))[1], ylim=c(min(toPlot[[nAlg]]),max(toPlot[[nAlg]])), ylab=results$results[[nAlg]]@qualMeasName, axes=FALSE)
              axis(1, at=rownames(toPlot[[nAlg]]))
              axis(2)
              box()
              for(i in 2:dim(toPlot[[nAlg]])[2]){
                lines(rownames(toPlot[[nAlg]]),toPlot[[nAlg]][,i], col=i)
              }
              legend(x="topright", legend=colnames(toPlot[[nAlg]]), col=1:dim(toPlot[[nAlg]])[2], lty=1, title=names(dimnames(toPlot[[nAlg]]))[2])
            } else { # Transpose matrix before plotting
              plot(colnames(toPlot[[nAlg]]),toPlot[[nAlg]][1,], type='l', xlab=names(dimnames(toPlot[[nAlg]]))[2], ylim=c(min(toPlot[[nAlg]]),max(toPlot[[nAlg]])), ylab=results$results[[nAlg]]@qualMeasName, axes=FALSE)
              axis(1, at=colnames(toPlot[[nAlg]]))
              axis(2)
              box()
              for(i in 2:dim(toPlot[[nAlg]])[1]){
                lines(colnames(toPlot[[nAlg]]),toPlot[[nAlg]][i,], col=i)
              }
              legend(x="topright", legend=rownames(toPlot[[nAlg]]), col=1:dim(toPlot[[nAlg]])[1], lty=1, title=names(dimnames(toPlot[[nAlg]]))[1])
            }
          }
        })
        
        # Level plot
        gWidgets::addhandlerchanged(plotTwo[[nAlgs]], handler = function(h,...){
          nAlg <- gWidgets::tag(h$obj)$alg
          Gray <- function(n){gray(seq(0,1, length.out=n))}
          if(gWidgets::svalue(plotFlip[[nAlg]])==FALSE){
            plot(lattice::levelplot(toPlot[[nAlg]], xlab=names(dimnames(toPlot[[nAlg]]))[1], ylab=names(dimnames(toPlot[[nAlg]]))[2], col.regions = Gray))
          } else { # Transpose matrix before plotting
            plot(lattice::levelplot(t(toPlot[[nAlg]]), xlab=names(dimnames(toPlot[[nAlg]]))[2], ylab=names(dimnames(toPlot[[nAlg]]))[1], col.regions = Gray))
          }
        })
        plotFlip[[nAlgs]] <<- gWidgets::gcheckbox("flip")
        gWidgets::enabled(plotTwo[[nAlgs]]) <<- FALSE
        gWidgets::enabled(plotFlip[[nAlgs]]) <<- FALSE
        rGroups[[nAlgs]] <<- gWidgets::ggroup(horizontal=TRUE)
        gWidgets::add(rGroups[[nAlgs]], plotOne[[nAlgs]],expand=FALSE)
        gWidgets::add(rGroups[[nAlgs]], plotTwo[[nAlgs]],expand=FALSE)
        gWidgets::add(rGroups[[nAlgs]], plotFlip[[nAlgs]],expand=FALSE)
        gWidgets::addSpace(groups[[nAlgs]], 10, horizontal=FALSE)
        gWidgets::add(groups[[nAlgs]], rGroups[[nAlgs]],expand=FALSE)
        
        # gWidgets::add(nb,groups[[nAlgs]],label=name)
        visible(parameterPlots[[nAlgs]]) <<- TRUE
        checkPlot(nAlgs)
      }
      
      # Check what can be plotted
      checkPlot <- function(nAlg){
        m <- list(results$results[[nAlg]])
        nam <- ""
        mins <- avg <- character(0)
        j <- 2
        k <- 0
        l <- 0
        # Prepare call to function 'qualMeas'
        for(i in 1:length(parameterList[[nAlg]])){
          def <- gWidgets::tag(parameterList[[nAlg]][[i]][[1]])$default
          if(is.null(def)){
            if(gWidgets::svalue(parameterList[[nAlg]][[i]][[1]])=="Overall min."){
              m[[j]] <- "overall"
              nam[j] <- gWidgets::tag(parameterList[[nAlg]][[i]][[1]])$name
              j <- j+1
            } else
              if(gWidgets::svalue(parameterList[[nAlg]][[i]][[1]])=="All"){
                m[[j]] <- "all"
                nam[j] <- gWidgets::tag(parameterList[[nAlg]][[i]][[1]])$name
                j <- j+1
              } else
                if(gWidgets::svalue(parameterList[[nAlg]][[i]][[1]])=="Chosen"){
                  m[[j]] <- gWidgets::svalue(parameterList[[nAlg]][[i]][[2]],index=TRUE)
                  nam[j] <- gWidgets::tag(parameterList[[nAlg]][[i]][[1]])$name
                  j <- j+1
                } else
                  if(gWidgets::svalue(parameterList[[nAlg]][[i]][[1]])=="Min."){
                    k <- k+1
                    mins[k] <- gWidgets::tag(parameterList[[nAlg]][[i]][[1]])$name
                    m[[j]] <- gWidgets::svalue(parameterList[[nAlg]][[i]][[2]],index=TRUE)
                    nam[j] <- gWidgets::tag(parameterList[[nAlg]][[i]][[1]])$name
                    j <- j+1
                  }
            if(gWidgets::svalue(parameterList[[nAlg]][[i]][[1]])=="Avg."){
              l <- l+1
              avg[l] <- gWidgets::tag(parameterList[[nAlg]][[i]][[1]])$name
              m[[j]] <- gWidgets::svalue(parameterList[[nAlg]][[i]][[2]],index=TRUE)
              nam[j] <- gWidgets::tag(parameterList[[nAlg]][[i]][[1]])$name
              j <- j+1
            }
          }
        }            # End for loop
        if(k>0){
          m[[j]] <- mins
          nam[j] <- "MIN"
          j <- j+1
        }
        if(l>0){
          m[[j]] <- avg
          nam[j] <- "AVG"
        }
        
        names(m) <- nam
        toPlot[[nAlg]] <<- drop(do.call(qualMeas,m)) # What will be plotted?
        if(length(dim(toPlot[[nAlg]]))>2){ # Too many dimensions to plot
          gWidgets::enabled(plotOne[[nAlg]]) <- FALSE
          gWidgets::enabled(plotTwo[[nAlg]]) <- FALSE
          gWidgets::enabled(plotFlip[[nAlg]]) <- FALSE
        } else if((length(dim(toPlot[[nAlg]]))==2) && (prod(dim(toPlot[[nAlg]]))>0)){ # Two dimensions to plot
          gWidgets::enabled(plotOne[[nAlg]]) <- TRUE
          gWidgets::enabled(plotTwo[[nAlg]]) <- TRUE
          gWidgets::enabled(plotFlip[[nAlg]]) <- TRUE
        } else if(is.null(dim(toPlot[[nAlg]])) && (length(toPlot[[nAlg]])>1)){ # One dimension to plot
          gWidgets::enabled(plotOne[[nAlg]]) <- TRUE
          gWidgets::enabled(plotTwo[[nAlg]]) <- FALSE
          gWidgets::enabled(plotFlip[[nAlg]]) <- FALSE
          for(i in 1:length(parameterList[[nAlg]])){
            if((gWidgets::svalue(parameterList[[nAlg]][[i]][[1]])=="Chosen" || gWidgets::svalue(parameterList[[nAlg]][[i]][[1]])=="Avg." || gWidgets::svalue(parameterList[[nAlg]][[i]][[1]])=="All" || gWidgets::svalue(parameterList[[nAlg]][[i]][[1]])=="Avg.") && length(gWidgets::svalue(parameterList[[nAlg]][[i]][[2]]))>1){
              oneDimNames[nAlg] <<- gWidgets::tag(parameterList[[nAlg]][[i]][[1]])$name
            }
          }
        } else { # No dimensions to plot
          gWidgets::enabled(plotOne[[nAlg]]) <- FALSE
          gWidgets::enabled(plotTwo[[nAlg]]) <- FALSE
          gWidgets::enabled(plotFlip[[nAlg]]) <- FALSE
        }
      }
      
      # ############# #
      # Main notebook #
      # ############# #
      nb <- gWidgets::gnotebook()
      gWidgets::add(main,nb,expand=TRUE)
      
      # Prepare notebook for algorithms
      for(i in 1:length(results$results)){
        groups[[i]] <- gWidgets::ggroup(horizontal=FALSE)
        name <- names(bA)[which(methodParse2==results$baselineTests[[i]]@algorithm@funcName)]
        gWidgets::add(nb,groups[[i]],label=name)
        progress[[i]] <- gWidgets::glabel('Setting up GUI...')
        gWidgets::add(groups[[i]], progress[[i]], expand=FALSE)
      }
      # ######### #
      # Comparing #
      # ######### #
      algNames <- c("-> Choose first result set")
      for(i in 1:length(results$results)){ # Prepare choices for result set chooser
        name <- names(bA)[which(methodParse2==results$baselineTests[[i]]@algorithm@funcName)]
        algNames[i+1] <- name
      }
      
      # Label, button, radio buttons and groups
      compareLabel <- gWidgets::glabel('Compare algorithms')
      curveButton <- gWidgets::gbutton(paste('Plot',results$results[[1]]@qualMeasName,'against',names(results$results[[1]]@param)[1]), handler = function(h,...){
        ns <- list()
        if(gWidgets::svalue(minBest,index=TRUE)==3){
          nams <- character(length(results$results))
          for(i in 1:length(results$results)){
            nam <- c("", names(results$results[[1]]@param)[1], "AVG")
            m <- list(results$results[[i]])
            minOver <- gWidgets::svalue(minWhich,index=TRUE)
            m[[2]] <- minOver
            m[[3]] <- names(results$results[[1]]@param)[1]
            names(m) <- nam
            ns.tmp <- do.call(qualMeas,m)
            the.min <- which(ns.tmp==min(ns.tmp), arr.ind = TRUE)
            the.dim <- dimnames(ns.tmp)
            
            nam <- c("")
            m <- list(results$results[[i]])
            nam <- append(nam, names(results$results[[1]]@param)[1])
            m[[2]]   <- "all"
            for(j in 2:length(results$results[[i]]@param)){
              nam <- append(nam, names(results$results[[i]]@param)[j])
              m[[j+1]] <- the.min[j]
            }
            names(m) <- nam
            ns[[i]] <- drop(do.call(qualMeas,m))
            nams[i] <- names(bA)[which(methodParse2==results$baselineTests[[i]]@algorithm@funcName)]
          }
        } else {
          nams <- character(length(results$results))
          nam <- c("", names(results$results[[1]]@param)[1], "DEFAULT")
          for(i in 1:length(results$results)){
            m <- list(results$results[[i]])
            m[[2]] <- "all"
            if(gWidgets::svalue(minBest)=="Overall min.")
              m[[3]] <- "overall.min"
            else
              m[[3]] <- "cond.min"
            names(m) <- nam
            ns[[i]] <- drop(do.call(qualMeas,m))
            nams[i] <- names(bA)[which(methodParse2==results$baselineTests[[i]]@algorithm@funcName)]
          }
        }
        plot(names(ns[[1]]),ns[[1]], type='l', xlab=nam[2], ylab=results$results[[1]]@qualMeasName)
        legend(x="topright", legend=nams, col=1:length(ns), lty=1)
        if(length(ns)>1){
          for(i in 2:length(ns)){
            lines(names(ns[[i]]),ns[[i]], col=i)
          }
        }
      })
      #	minBest <- gWidgets::gradio(c("overall.min","cond.min"), horizontal=TRUE)
      result <- results$results[[1]]
      method <- names(bA)[which(methodParse2==results$baselineTests[[1]]@algorithm@funcName)]
      nameStrs <- c(names(result@param[1]), rownames(bA[[method]]))
      if(max(nchar(as.character(result@param[[nameStrs[1]]],scientific=TRUE)))>8){ # More than 8 digits => use scientific format
        minBest  <- gWidgets::gradio(c("Overall min.","Min.", "Min. avg."))
        minWhich <- gWidgets::gcheckboxgroup(format(result@param[[nameStrs[1]]], checked=!logical(length(result@param[[nameStrs[1]]])), scientific=TRUE, digits=3), horizontal=TRUE)
        gWidgets::enabled(minWhich) <- FALSE
      } else {
        minBest  <- gWidgets::gradio(c("Overall min.","Min.", "Min. avg."))
        minWhich <- gWidgets::gcheckboxgroup(result@param[[nameStrs[1]]], checked=logical(length(result@param[[nameStrs[1]]])), horizontal=TRUE)
        gWidgets::enabled(minWhich) <- FALSE}
      gWidgets::addhandlerchanged(minBest, handler = function(h,...){
        if(gWidgets::svalue(minBest,index=TRUE)==3){
          gWidgets::enabled(minWhich) <- TRUE
          gWidgets::svalue(minWhich, index=TRUE) <- 1:length(result@param[[nameStrs[1]]])
        } else {
          gWidgets::enabled(minWhich) <- FALSE
          gWidgets::svalue(minWhich, index=TRUE) <- integer()}
      })
      gWidgets::addhandlerchanged(minWhich, handler = function(h,...){
        if(length(gWidgets::svalue(minWhich, index=TRUE)) == 0 && gWidgets::svalue(minBest,index=TRUE)==3){
          gWidgets::enabled(curveButton) <- FALSE
        } else {
          gWidgets::enabled(curveButton) <- TRUE}
      })
      
      sGroup <- gWidgets::ggroup(horizontal=FALSE)
      mGroup <- gWidgets::ggroup(horizontal=TRUE)
      
      # Combine elements
      gWidgets::addSpace(sGroup,10,horizontal=FALSE)
      gWidgets::add(sGroup,compareLabel,expand=FALSE)
      gWidgets::addSpace(sGroup,10,horizontal=FALSE)
      gWidgets::add(mGroup,minBest,expand=FALSE)
      gWidgets::add(mGroup,minWhich,expand=FALSE)
      gWidgets::addSpace(mGroup,10,horizontal=TRUE)
      gWidgets::add(mGroup,curveButton,expand=FALSE)
      gWidgets::add(sGroup,mGroup,expand=FALSE)
      gWidgets::add(nb, sGroup, label="Compare")
      
      # Build GUI for algorithms
      for(i in 1:length(results$results)){
        nm <- which(methodParse2==results$baselineTests[[i]]@algorithm@funcName)
        addAlg(nm, results$results[[i]]);
      }
    } else {
      warning('Package lattice not installed')
      return(list())
    }
  } else {
    warning('Package gWidgets not installed')
    return(list())
  }
}