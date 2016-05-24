###Plot

setMethod("plot",signature(x="Evaluation",y="missing"),
          function(x,
                   runs0=1:nrow(result(x)), dims0=1:ncol(result(x)),  ...
          ){
            dots <- list(...)
            ldim0 <- min(getdistrTEstOption("MaxNumberofPlottedEvaluationDims"),
                         length(dims0))
            if(ldim0<length(dims0))
                warning(paste("your evaluation is too big; only ", ldim0,
                              "evaluation dimensions are plotted"))
            boxplot(result(x)[runs0,dims0[1:ldim0]],...)
          }
          )

setMethod("plot",signature(x="EvaluationList",y="missing"),
           function(x,
                    runs0=1:nrow(result(Elist(x)[[1]])),
                    dims0= 1: ifelse(sum(grep("\\.id($|\\.)",
                                     colnames(result(Elist(x)[[1]]))))>0,
                                     ### contaminated data or not?
                                     ncol(result(Elist(x)[[1]]))/2,
                                     ncol(result(Elist(x)[[1]]))),
                    evals0=1:length(Elist(x)), ... )
{ dots <- list(...)

  ldim0 <- min(getdistrTEstOption("MaxNumberofPlottedEvaluationDims"),
               length(dims0))
  levals0 <- min(getdistrTEstOption("MaxNumberofPlottedEvaluations"),
                 length(evals0))

  if((ldim0<length(dims0))||(levals0<length(evals0)))
                warning(gettextf(
"your evaluation list is too big; only %i x%i evaluations x dimensions are plotted",
                        levals0, ldim0))

  evallist <- Elist(x)
  len <- length(evallist)
  resdim <- ncol(result(evallist[[1]]))
  nl <- nrow(result(evallist[[1]]))

  # the names for the different columns[1..resdim] of the different
  # Evaluations[1..len]
  resnames <- matrix(0,len,resdim)
  for(i in 1:len)
      resnames[i,] <- colnames(result(evallist[[i]]))

  ### is the data split into ideal/contaminated data?
  id0 <-grep("\\.id($|\\.)",resnames[1,])

  evallist0 <- Elist(x)[evals0[1:levals0]]
  len0 <- length(evallist0)

  if (sum(id0 > 0))

     {dims1 <- c(dims0[1:ldim0],dims0[1:ldim0]+resdim/2)
      resdim0 <- ncol(result(evallist0[[1]])[runs0,dims1])
      nl0 <- nrow(result(evallist0[[1]])[runs0,dims1])
     }

  else

     {dims1 <- dims0[1:ldim0]
      resdim0 <- ncol(result(evallist0[[1]])[runs0,dims1])
      nl0 <- nrow(result(evallist0[[1]])[runs0,dims1])
     }

  ma <- data.frame(matrix(0,nl0,len0*resdim0))

  ## reorganize the frames
  ### --- grouped by result-dimension (e.g. coordinates of an estimator)
  ###     and within this --- if possible --- grouped by ideal/real Data
  if( sum(id0 > 0))
     { for(i in 1:(resdim0/2))
           for(j in 1:len0)

              {ma[,(i-1)*2*len0+j] <- result(evallist0[[j]])[runs0,dims1[i]]
               ma[,(2*i-1)*len0+j] <-
                            result(evallist0[[j]])[runs0,dims1[i+resdim0/2]]
               colnames(ma)[(i-1)*2*len0+j] <-
                        colnames(result(evallist0[[j]]))[dims1[i]]
               colnames(ma)[(2*i-1)*len0+j] <-
                        colnames(result(evallist0[[j]]))[dims1[i+resdim0/2]]
               }
       resdim0 <- resdim0/2; len0 <- len0*2
     }
  else
     { for(i in 1:resdim0)
           for(j in 1:len0)
              {ma[,(i-1)*len0+j] <- result(evallist0[[j]])[runs0,i]
               colnames(ma)[(i-1)*len0+j] <- colnames(result(evallist0[[j]]))[i]
               }
     }
  
  o.warn <- getOption("warn")
  on.exit(options("warn"=o.warn))
        
  main0 <- character(resdim0)
  if("main" %in% names(dots))
      { options("warn" = -1)
        main0[1:resdim0] <- dots[["main"]]
        options("warn" = o.warn) }
  else
      main0 <- paste(gettextf("%d. coordinate",dims1[1:resdim0]))

  ylim0<-matrix(0,2,resdim0)
  if("ylim" %in% names(dots))
      { options("warn" = -1)
        ylim1 <- as.matrix(dots[["ylim"]])
        c1 <- ncol(ylim1); c2 <- resdim0%/%c1; c3 <- resdim0%%c1
        if(c2>0)
           ylim0[,1:(c2*c1)] <- ylim1
        if(c3>0)
           ylim0[,c2*c1+(1:c3)] <- ylim1[,1:c3]
        options("warn" = o.warn) }

  opar <- par(no.readonly = TRUE)
  # opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
  on.exit(par(opar))
  par(mfrow=c(resdim0,1))


  for(i in 1:resdim0)
      {if("main" %in% names(dots))
          dots[["main"]] <- main0[i]
       else {
          dots[[length(dots)+1]] <- main0[i]
          names(dots)[length(dots)] <- "main"
            }
       if("ylim" %in% names(dots))
          dots[["ylim"]] <- ylim0[,i]

       dots[["x"]]<- as.data.frame(ma[,(i-1)*len0+(1:len0)])
       do.call("boxplot", args = dots)

      }
   return(invisible())
})
