### $Id: optimWizard.R 193 2012-06-24 21:13:42Z kristl $

optimWizard <- function(X, y, postproc, predictionTest, cvsegments){
  ## Organize optimization through GUI
  
  if(missing(X))
    stop('No data specified')
  if(missing(y))
    stop('No response specified')
  if(missing(predictionTest)){
    rr <- FALSE
    predictionTest <- NULL}
  else
    rr <- TRUE
  if(missing(postproc)){
    pp <- FALSE
    postproc <- NULL}
  else
    pp <- TRUE
  if(missing(cvsegments)){
    cc <- FALSE
    cvsegments <- NULL}
  else
    cc <- TRUE
  bltest <- NULL
  
  if(requireNamespace("gWidgets", quietly = TRUE)){
    if(requireNamespace("pls", quietly = TRUE)){
      
      # Set up main window with internal container
      win <- gWidgets::gwindow("Optimisation wizard", width=450, height=300)
      main <- gWidgets::ggroup(horizontal=FALSE)
      gWidgets::add(win,main)
      
      # Initialize parameters and paramter lists
      nAlgs <- 0
      used <- numeric(0)
      if(exists("baselineAlgorithmsGUI",envir=.GlobalEnv)){
        bAGUI <- get("baselineAlgorithmsGUI",envir=.GlobalEnv)
      } else {
        bAGUI <- baselineAlgorithmsGUI
      }
      GUI.names <- sort(names(bAGUI))
      if(exists("baselineAlgorithms",envir=.GlobalEnv)){
        bA <- get("baselineAlgorithms",envir=.GlobalEnv)
      } else {
        bA <- baselineAlgorithms
      }
      genChoosers <- rGroups <- groups <- removes <- parameterGroup <- parameterList <- method <- list()
      
      
      # Function for adding a baseline correction method, including sub-functions and variables
      addAlg <- function(nm){
        # Initialize parameters and groups
        nAlgs <<- nAlgs + 1
        used[nAlgs] <<- 1
        groups[[nAlgs]] <<- gWidgets::ggroup(horizontal=FALSE)
        rGroups[[nAlgs]] <<- gWidgets::ggroup(horizontal=TRUE)
        #    name <- names(bA)[nm]
        name <- nm
        nameLong <- bA[[name]]@description
        method[[nAlgs]] <<- name
        
        # Add some space and the name of the baseline correction algorithm
        gWidgets::addSpace(groups[[nAlgs]],20)
        gWidgets::add(groups[[nAlgs]],gWidgets::glabel(nameLong))
        gWidgets::addSpace(groups[[nAlgs]],10)
        
        # Function for setting up parameter array in GUI together with buttons and functions
        addParameterGroup <- function(nameStr, lineNo){
          genChoosers[[nAlgs]][[lineNo]] <<- gdroplist(c("-> Generate", "Linear","Exponential"), selected=1)
          gWidgets::tag(genChoosers[[nAlgs]][[lineNo]],"no") <<- lineNo
          gWidgets::tag(genChoosers[[nAlgs]][[lineNo]],"name") <<- nameStr
          parameterList[[nAlgs]][[lineNo]]  <<- c(gWidgets::gedit(text = "", width=1,coerce.with=as.numeric),gWidgets::gedit(width=5,coerce.with=as.numeric),gWidgets::gedit(width=10,coerce.with=as.numeric),gWidgets::gedit(width=15),genChoosers[[nAlgs]][[lineNo]])
          # Generate sequence based on 'From', 'To', 'Steps' and choice from droplist
          gWidgets::addhandlerchanged(parameterList[[nAlgs]][[lineNo]][[5]], handler = function(h,...){
            if(gWidgets::svalue(genChoosers[[nAlgs]][[lineNo]],index=TRUE)>1){
              type <- gWidgets::svalue(h$obj,index=TRUE)-1
              gWidgets::svalue(genChoosers[[nAlgs]][[lineNo]],index=TRUE) <<- 1
              linNo <- gWidgets::tag(h$obj)$no
              if(is.finite(gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[1]])) && is.finite(gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[2]])) && is.finite(gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[3]]))){
                if(type==1){
                  linSeq <- seq(gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[1]]),gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[2]]),length.out=gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[3]]))
                  linOut <- linSeq[1]
                  if(length(linSeq)>1)
                    for(i in 2:length(linSeq))
                      linOut <- paste(linOut, linSeq[i], sep=", ")
                  gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[4]]) <- linOut
                }
                if(type==2){
                  linSeq <- exp(seq(log(gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[1]])),log(gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[2]])),length.out=gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[3]])))
                  linOut <- linSeq[1]
                  if(length(linSeq)>1)
                    for(i in 2:length(linSeq))
                      linOut <- paste(linOut, linSeq[i], sep=", ")
                  gWidgets::svalue(parameterList[[nAlgs]][[linNo]][[4]]) <- linOut
                }
              } else {
                gmessage("Missing value(s) in 'From', 'To' or 'Steps'", title="Sequence", icon = "warning")
              }
            }
          })
          # Set up visual optimization array for parameters
          parameterGroup[[nAlgs]][lineNo+1,1] <<- gWidgets::glabel(text=nameStr)
          size(parameterList[[nAlgs]][[lineNo]][[1]]) <<- c(60,28)
          size(parameterList[[nAlgs]][[lineNo]][[2]]) <<- c(60,28)
          size(parameterList[[nAlgs]][[lineNo]][[3]]) <<- c(20,28)
          size(parameterList[[nAlgs]][[lineNo]][[4]]) <<- c(220,28)
          size(parameterList[[nAlgs]][[lineNo]][[5]]) <<- c(120,28)
          parameterGroup[[nAlgs]][lineNo+1,2] <<- parameterList[[nAlgs]][[lineNo]][[1]]
          parameterGroup[[nAlgs]][lineNo+1,3] <<- parameterList[[nAlgs]][[lineNo]][[2]]
          parameterGroup[[nAlgs]][lineNo+1,4] <<- parameterList[[nAlgs]][[lineNo]][[3]]
          parameterGroup[[nAlgs]][lineNo+1,5] <<- parameterList[[nAlgs]][[lineNo]][[4]]
          parameterGroup[[nAlgs]][lineNo+1,6] <<- parameterList[[nAlgs]][[lineNo]][[5]]
        }
        
        # Set up visual optimization array for parameters
        parameterList[[nAlgs]] <<- genChoosers[[nAlgs]] <<- list()
        parameterGroup[[nAlgs]] <<- gWidgets::glayout(homogeneous = FALSE, spacing = 5, container=groups[[nAlgs]])
        parameterGroup[[nAlgs]][1,1] <<- gWidgets::glabel("")
        parameterGroup[[nAlgs]][1,2] <<- gWidgets::glabel("From:")
        parameterGroup[[nAlgs]][1,3] <<- gWidgets::glabel("To:")
        parameterGroup[[nAlgs]][1,4] <<- gWidgets::glabel("Steps:")
        parameterGroup[[nAlgs]][1,5] <<- gWidgets::glabel("Sequence:")
        nameStrs <- rownames(bAGUI[[method[[nAlgs]]]])
        lns <- length(nameStrs)
        for(i in 1:lns)
          addParameterGroup(nameStrs[i],i)
        parameterGroup[[nAlgs]][lns+2,2] <- gWidgets::gbutton("Collect", handler = function(h,...){
          #      if(exists("baseline.current")){
          if(getBaselineEnv("baseline.current")$method == name){
            for(i in 1:lns){
              gWidgets::svalue(parameterList[[nAlgs]][[i]][[1]]) <- getBaselineEnv("baseline.current")$parValues[i]
            }
          } else {
            gmessage(paste("'baseline.current$method' is not equal to '", name, "'", sep=""), title="Sequence", icon = "warning")
          }
          #      } else {
          #        gmessage("'baseline.current' not found", title="Sequence", icon = "warning")
          #      }
        })
        parameterGroup[[nAlgs]][lns+2,3] <- gWidgets::gbutton("Collect", handler = function(h,...){
          if(getBaselineEnv("baseline.current")$method == name){
            for(i in 1:lns){
              gWidgets::svalue(parameterList[[nAlgs]][[i]][[2]]) <- getBaselineEnv("baseline.current")$parValues[i]
            }
          } else {
            gmessage(paste("'baseline.current$method' is not equal to '", name, "'", sep=""), title="Sequence", icon = "warning")
          }
        })
        parameterGroup[[nAlgs]][lns+2,4] <- gWidgets::glabel("<- from current algorithm in baselineGUI")
        
        # Button for removal of method and adding the method to the main window's notebook
        removes[[nAlgs]] <<- gWidgets::gbutton("Remove baseline correction method")
        gWidgets::tag(removes[[nAlgs]],"nAlg") <<- nAlgs
        addhandlerclicked(removes[[nAlgs]], handler = function(h,...){
          nAlg <- gWidgets::tag(h$obj)$nAlg
          gWidgets::dispose(nb)
          used[nAlg] <<- 0
        })
        gWidgets::add(rGroups[[nAlgs]],removes[[nAlgs]],expand=FALSE)
        gWidgets::addSpace(groups[[nAlgs]],20)
        gWidgets::add(groups[[nAlgs]],rGroups[[nAlgs]],expand=FALSE)
        gWidgets::add(nb,groups[[nAlgs]],label=name)
        gWidgets::visible(parameterGroup[[nAlgs]]) <- TRUE
      }
      
      # Droplists for adding a baseline correction method and choosing post processing
      namesG <- character(length(GUI.names)+1)
      namesG[1] <- '-> Choose method for optimisation'
      for(i in 1:length(GUI.names)){ # Let bAGUI control, and bA have descriptions -------------
                                     namesG[i+1] <- paste("'", ifelse(is.null(bA[[GUI.names[i]]]@description),"",bA[[GUI.names[i]]]@description), " (", GUI.names[i], ")'", sep="")
      }
      methodChooser <- gdroplist(namesG,
                                 selected=1, handler = function(h,...){if(gWidgets::svalue(methodChooser,index=TRUE)>1) addAlg(GUI.names[gWidgets::svalue(methodChooser,index=TRUE)-1]); gWidgets::svalue(methodChooser,index=TRUE)<-1})
      postChooser <- gdroplist(c("None","Norm (L2)", "Mean", "Median",
                                 "Sum", "Sum of squares", "L1 postproc", "Maximum"),	selected=1)
      regChosen <- FALSE
      regOutGroup <- segChooser <- segNumber <- regParam <- numeric(0)
      regFrom <- regTo <- regSteps <- lambdaSequence <- numeric(0)
      # Analysis droplist
      regChooser <- gdroplist(c("-> Choose an analysis","PLSR / RMSEP","Ridge Regression / RMSEP"),	selected=1, handler = function(h,...){
        if(regChosen == TRUE){
          gWidgets::delete(regframe, regOutGroup) # Remove old analysis if chosen
        }
        if(gWidgets::svalue(regChooser,index=TRUE)>1){ # Analysis chosen
          regOutGroup <<- gWidgets::ggroup(horizontal=TRUE)
          regIntGroup1 <- gWidgets::ggroup(horizontal=FALSE)
          regIntGroup2 <- gWidgets::ggroup(horizontal=FALSE)
          regIntGroup3 <- gWidgets::ggroup(horizontal=FALSE)
          regIntGroup4 <- gWidgets::ggroup(horizontal=TRUE)
          regIntGroup5 <- gWidgets::ggroup(horizontal=FALSE)
          regIntGroupA <- gWidgets::ggroup(horizontal=TRUE)
          regIntGroupB <- gWidgets::ggroup(horizontal=FALSE)
          if(gWidgets::svalue(regChooser,index=TRUE)==2){ # Display extra parameters for PLSR
            if(cc)
              segChooser <<- gdroplist(c("random", "consecutive", "interleaved","custom"), selected=4, handler = function(h,...){
                if(gWidgets::svalue(segChooser, index=TRUE)==4){
                  gWidgets::svalue(segNumber) <- length(cvsegments)
                  gWidgets::enabled(segNumber) <- FALSE
                } else
                  gWidgets::enabled(segNumber) <- TRUE
              })
            else
              segChooser <<- gdroplist(c("random", "consecutive", "interleaved"), selected=1)
            segNumber <<- gWidgets::gedit("10")
            regParam <<- gWidgets::gedit()
            gWidgets::add(regIntGroup1, gWidgets::glabel("CV segment type"), expand=FALSE)
            gWidgets::add(regIntGroup1, segChooser, expand=FALSE)
            if(cc){
              # gWidgets::add(regIntGroup1, gWidgets::glabel(length(cvsegments)), expand=FALSE)
              gWidgets::svalue(segNumber) <- length(cvsegments)
              gWidgets::enabled(segNumber) <- FALSE
            }
            gWidgets::add(regIntGroup1, gWidgets::glabel("Number of segments"), expand=FALSE)
            gWidgets::add(regIntGroup1, segNumber, expand=FALSE)
            gWidgets::add(regIntGroup2, gWidgets::glabel("Number of components"), expand=FALSE)
            gWidgets::add(regIntGroup2, regParam, expand=FALSE)
            gWidgets::add(regOutGroup, regIntGroup1, expand=FALSE)
            gWidgets::add(regOutGroup, regIntGroup2, expand=FALSE)
          }
          if(gWidgets::svalue(regChooser,index=TRUE)==3){ # Display extra parameters for Ridge regression
            regFrom <<- gWidgets::gedit(coerce.with=as.numeric)
            size(regFrom) <<- c(60,24)
            regTo <<- gWidgets::gedit(coerce.with=as.numeric)
            size(regTo) <<- c(60,24)
            regSteps <<- gWidgets::gedit(coerce.with=as.numeric)
            size(regSteps) <<- c(60,24)
            lambdaSequence <<- gdroplist(c("-> Generate", "Linear","Exponential"))
            # Generate sequence based on 'From', 'To', 'Steps' and choice from droplist
            gWidgets::addhandlerchanged(lambdaSequence, handler = function(h,...){
              if(gWidgets::svalue(lambdaSequence,index=TRUE)>1){
                type <- gWidgets::svalue(lambdaSequence,index=TRUE)-1
                gWidgets::svalue(lambdaSequence,index=TRUE) <<- 1
                if(is.finite(gWidgets::svalue(regFrom)) && is.finite(gWidgets::svalue(regTo)) && is.finite(gWidgets::svalue(regSteps))){
                  if(type==1){
                    linSeq <- seq(gWidgets::svalue(regFrom),gWidgets::svalue(regTo),length.out=gWidgets::svalue(regSteps))
                    linOut <- linSeq[1]
                    if(length(linSeq)>1)
                      for(i in 2:length(linSeq))
                        linOut <- paste(linOut, linSeq[i], sep=", ")
                    gWidgets::svalue(regParam) <- linOut
                  }
                  if(type==2){
                    linSeq <- exp(seq(log(gWidgets::svalue(regFrom)),log(gWidgets::svalue(regTo)),length.out=gWidgets::svalue(regSteps)))
                    linOut <- linSeq[1]
                    if(length(linSeq)>1)
                      for(i in 2:length(linSeq))
                        linOut <- paste(linOut, linSeq[i], sep=", ")
                    gWidgets::svalue(regParam) <- linOut
                  }
                } else {
                  gmessage("Missing value(s) in 'From', 'To' or 'Steps'", title="Sequence", icon = "warning")
                }
              }
            })
            paramLabel <- gWidgets::glabel("Ridge parameter")
            regParam <<- gWidgets::gedit()
            gWidgets::add(regIntGroup1, gWidgets::glabel("From"), expand=FALSE)
            gWidgets::add(regIntGroup1, regFrom, expand=FALSE)
            gWidgets::add(regIntGroup2, gWidgets::glabel("To"), expand=FALSE)
            gWidgets::add(regIntGroup2, regTo, expand=FALSE)
            gWidgets::add(regIntGroup3, gWidgets::glabel("Steps"), expand=FALSE)
            gWidgets::add(regIntGroup3, regSteps, expand=FALSE)
            gWidgets::add(regIntGroup4, lambdaSequence, expand=FALSE)
            gWidgets::add(regIntGroup5, gWidgets::glabel("Lambda sequence"), expand=FALSE)
            gWidgets::add(regIntGroup5, regParam, expand=FALSE)
            gWidgets::add(regIntGroupA, regIntGroup1, expand=FALSE)
            gWidgets::add(regIntGroupA, regIntGroup2, expand=FALSE)
            gWidgets::add(regIntGroupA, regIntGroup3, expand=FALSE)
            gWidgets::add(regIntGroupA, regIntGroup4, expand=FALSE)
            gWidgets::add(regIntGroupB, regIntGroupA, expand=FALSE)
            gWidgets::add(regIntGroupB, regIntGroup5, expand=FALSE)
            gWidgets::add(regOutGroup, regIntGroupB, expand=FALSE)
          }
          gWidgets::add(regframe,regOutGroup,expand=FALSE)
          regChosen <<- TRUE
        }
      })
      
      # Notebook containing settings and baseline correction methods, group initialization
      nb <- gWidgets::gnotebook()
      sGroup <- gWidgets::ggroup(horizontal=TRUE)
      sGroup1 <- gWidgets::ggroup(horizontal=FALSE)
      # size(sGroup1) <- c(400,300)
      sGroup2 <- gWidgets::gframe("Optimisation", horizontal=FALSE)
      # size(sGroup2) <- c(170,300)
      methFrame <- gWidgets::gframe('Baseline correction', horizontal=FALSE)
      methGroup <- gWidgets::ggroup(horizontal=TRUE)
      normFrame <- gWidgets::gframe('Post processing', horizontal=FALSE)
      normGroup <- gWidgets::ggroup(horizontal=TRUE)
      regframe <- gWidgets::gframe('Analysis and quality measure', horizontal=FALSE)
      reggroup <- gWidgets::ggroup(horizontal=TRUE)
      verbCheck <- gWidgets::gcheckbox("Verbose", checked=TRUE)
      # Verification of optimization parameters
      verifyButton <- gWidgets::gbutton("Verify setup", handler = function(h,...){
        gWidgets::enabled(saveButton) <- FALSE
        gWidgets::enabled(startButton) <- FALSE
        # Check basic settings for optimisation
        options(warn=-1)
        if(gWidgets::svalue(regChooser, index=TRUE)>1)
          rp <- as.numeric(strsplit(gWidgets::svalue(regParam), c(","))[[1]])
        else
          rp <- NA
        options(warn=0)
        if(sum(used)==0){
          gmessage("No baseline correction algorithm chosen", title="Not ready", icon = "warning")
        } else if(gWidgets::svalue(regChooser, index=TRUE)==1){
          gmessage("No analysis chosen", title="Not ready", icon = "warning")
        } else if(!rr && nchar(gWidgets::svalue(regParam))==0){
          gmessage("Regression parameter not specified", title="Not ready", icon = "warning")
        } else if(!rr && (!is.finite(sum(rp)))){
          gmessage("Regression parameter incorrectly specified", title="Not ready", icon = "warning")
        } else {
          u <- 0
          faulty <- ""
          for(i in 1:nAlgs){
            if(used[i] == 1){
              for(j in 1:length(parameterList[[i]])){
                m <- as.numeric(strsplit(gWidgets::svalue(parameterList[[i]][[j]][[4]]), c(","))[[1]])
                if((length(m)==0 || sum(is.na(m))>0) && u==0){ # Where did error occur?
                  u <- u+1
                  faulty <- paste((rownames(bAGUI[[method[[i]]]])[j]), "of baseline correction algorithm", method[[i]])
                }
              }
            }
          }
          if(u==0){
            # Collect data for analysis
            if(!rr){
              if(gWidgets::svalue(regChooser, index=TRUE) == 2){
                if(cc)
                  predictionTest <<- new("PLSRTest", ncomp = as.numeric(gWidgets::svalue(regParam)), cvsegments = pls::cvsegments)
                else
                  predictionTest <<- new("PLSRTest", ncomp = as.numeric(gWidgets::svalue(regParam)), cvsegments = pls::cvsegments(dim(X)[1], as.numeric(gWidgets::svalue(segNumber)), type=gWidgets::svalue(segChooser)))
              } else if(gWidgets::svalue(regChooser, index=TRUE) == 3){
                predictionTest <<- new("ridgeRegressionTest", lambda = as.numeric(strsplit(gWidgets::svalue(regParam), c(","))[[1]]))
              }
            }
            bltest <<- list()
            q <- 0
            for(i in 1:nAlgs){
              if(used[i] == 1){
                q <- q+1
                params <- list()
                parNames <- character(length(parameterList[[i]]))
                for(j in 1:length(parameterList[[i]])){
                  params[[j]] <- as.numeric(strsplit(gWidgets::svalue(parameterList[[i]][[j]][[4]]), c(","))[[1]])
                  parNames[j] <- rownames(bAGUI[[method[[i]]]])[j]
                }
                names(params) <- parNames
                bltest[[q]] <<- new("baselineAlgTest", algorithm = bA[[method[[i]]]],
                                    param = params)
              }
            }
            # Choice of normalisation
            if(!pp){
              if(gWidgets::svalue(postChooser, index=TRUE)==1)
                postproc <<- NULL
              else {
                if(gWidgets::svalue(postChooser, index=TRUE)==2) # Norm (L2)
                  postproc <<- function(X){ for(i in 1:dim(X)[1]){ X[i,] <- X[i,]/sqrt(X[i,]%*%X[i,])};X}
                if(gWidgets::svalue(postChooser, index=TRUE)==3) # Mean
                  postproc <<- function(X){ for(i in 1:dim(X)[1]){ X[i,] <- X[i,]/mean(X[i,])};X}
                if(gWidgets::svalue(postChooser, index=TRUE)==4) # Median
                  postproc <<- function(X){ for(i in 1:dim(X)[1]){ X[i,] <- X[i,]/median(X[i,])};X}
                if(gWidgets::svalue(postChooser, index=TRUE)==5) # Sum
                  postproc <<- function(X){ for(i in 1:dim(X)[1]){ X[i,] <- X[i,]/sum(X[i,])};X}
                if(gWidgets::svalue(postChooser, index=TRUE)==6) # Sum of squares
                  postproc <<- function(X){ for(i in 1:dim(X)[1]){ X[i,] <- X[i,]/(X[i,]%*%X[i,])};X}
                if(gWidgets::svalue(postChooser, index=TRUE)==7) # L1 postproc
                  postproc <<- function(X){ for(i in 1:dim(X)[1]){ X[i,] <- X[i,]/sum(abs(X[i,]))};X}
                if(gWidgets::svalue(postChooser, index=TRUE)==8) # Maximum
                  postproc <<- function(X){ for(i in 1:dim(X)[1]){ X[i,] <- X[i,]/max(X[i,])};X}
              }
            }
            gWidgets::enabled(saveButton) <- TRUE
            gWidgets::enabled(startButton) <- TRUE
          } else
            gmessage(paste("Check parameter",faulty), title="Not ready", icon = "warning")
        }
      })
      saveButton <- gWidgets::gbutton("Save setup", handler = function(h,...){
        putBaselineEnv("bltests", bltest)
        putBaselineEnv("predictionTest", predictionTest)
        putBaselineEnv("postproc", postproc)
        cat(paste("\n# To run optimization later:\nopts <- getOptim()\noptimRes <- doOptim(opts$bltests, X, y, opts$predictionTest,\n        postproc = opts$postproc, verbose =", gWidgets::svalue(verbCheck), ", cleanTmp = TRUE)\n"))
      })
      startButton <- gWidgets::gbutton("START", handler = function(h,...){
        # Run optimisation
        putBaselineEnv("optimRes", doOptim(bltest, X, y, predictionTest,
                                           postproc = postproc, verbose = gWidgets::svalue(verbCheck),
                                           cleanTmp = TRUE))
        cat("# To retrieve optimisation results later:\nmyResults <- getOptimRes()")
      })
      gWidgets::enabled(saveButton) <- FALSE
      gWidgets::enabled(startButton) <- FALSE
      
      # Settings for correction, normalization and analysis
      gWidgets::addSpace(sGroup1, 10, horizontal=FALSE)
      gWidgets::add(methFrame, methodChooser, expand=FALSE)
      gWidgets::add(methGroup, methFrame, expand=FALSE)
      gWidgets::add(sGroup1, methGroup, expand=FALSE)
      gWidgets::addSpace(sGroup1, 15, horizontal=FALSE)
      if(!pp)
        gWidgets::add(normFrame, postChooser, expand=FALSE)
      else
        gWidgets::add(normFrame, gWidgets::glabel("User specified"), expand=FALSE)
      gWidgets::add(normGroup, normFrame, expand=FALSE)
      gWidgets::add(sGroup1, normGroup, expand=FALSE)
      gWidgets::addSpace(sGroup1, 15, horizontal=FALSE)
      if(!rr)
        gWidgets::add(regframe, regChooser, expand=FALSE)
      else
        gWidgets::add(normFrame, gWidgets::glabel("User specified"), expand=FALSE)
      gWidgets::add(reggroup, regframe, expand=FALSE)
      gWidgets::add(sGroup1, reggroup, expand=FALSE)
      
      # Settings for optimization
      gWidgets::addSpace(sGroup2, 10, horizontal=FALSE)
      gWidgets::add(sGroup2, verbCheck, expand=FALSE)
      gWidgets::addSpace(sGroup2, 30, horizontal=FALSE)
      gWidgets::add(sGroup2, verifyButton, expand=FALSE)
      gWidgets::addSpace(sGroup2, 30, horizontal=FALSE)
      gWidgets::add(sGroup2, saveButton, expand=FALSE)
      gWidgets::addSpace(sGroup2, 10, horizontal=FALSE)
      gWidgets::add(sGroup2, startButton, expand=FALSE)
      
      gWidgets::add(sGroup, sGroup1)
      gWidgets::addSpace(sGroup, 30, horizontal=TRUE)
      gWidgets::add(sGroup, sGroup2)
      gWidgets::add(nb, sGroup, label="Settings")
      gWidgets::add(main,nb,expand=TRUE)
      # gWidgets::add(main,gstatusbar("Tester statusbar"))
    } else {
      warning('Package pls not installed')
      return(list())
    }
  } else {
    warning('Package gWidgets not installed')
    return(list())
  }
}

getOptim <- function(){
  list(bltests        = getBaselineEnv("bltests"),
       predictionTest = getBaselineEnv("predictionTest"),
       postproc       = getBaselineEnv("postproc"))
}
getOptimRes <- function(){
  getBaselineEnv("optimRes")
}