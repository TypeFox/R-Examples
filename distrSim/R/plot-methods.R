### Data-Class
## plot method

## changed w.r.t <1.8            
##setMethod("plot","Dataclass",
##          function(x,y=NULL,...){
##old:             y0<-1:runs(x)
##   matplot(y0,Data(x),xlab="Runindex",ylab="data",type="p",pch="*",col="blue")
##          })

setMethod("plot",signature(x = "Dataclass", y="missing"), 
           function(x, obs0=1:samplesize(x), dims0=1:obsDim(x), 
                    runs0=1:runs(x), ...){

            dots <- match.call(call = sys.call(sys.parent(1)), 
                               expand.dots = FALSE)$"..."
            doEnd <- FALSE
            if(!is.null(dots[["panel.first"]])) 
                {doEnd<- TRUE
                 dots[["panel.first"]] <- substitute(pf, 
                                         list(pf=dots[["panel.first"]]))
                }
            if(!is.null(dots[["panel.last"]])) 
                {doEnd<- TRUE
                 dots[["panel.last"]] <- substitute(pf, 
                                         list(pf=dots[["panel.last"]]))
                }
            lobs0 <- min(getdistrSimOption("MaxNumberofPlottedObs"), 
                         length(obs0))           
            lrun0 <- min(getdistrSimOption("MaxNumberofPlottedRuns"), 
                         length(runs0))           
            ldim0 <- min(getdistrSimOption("MaxNumberofPlottedObsDims"), 
                         length(dims0))           

            if( (lrun0 < length(runs0)) || (ldim0 < length(dims0)) || 
                (lobs0 < length(obs0)) )   
                warning(gettextf(
"Your data set is too big; only %i x %i x % i observations x dimensions x runs are plotted",
                        lobs0, ldim0, lrun0))
#            get(getOption("device"))()

            opar <- par(no.readonly = TRUE)
#            opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
            on.exit(par(opar))

            o.warn <- getOption("warn")
            on.exit(options("warn"=o.warn))
            
            options("warn" = -1)
            par(mfrow=c(1,lrun0))
            
            y0<-1:length(obs0)
            dots[["x"]] <- y0
#            if(lrun0==1)
#               matplot(y0,Data(x)[,dims0[1:ldim0]],
#      xlab=gettext("observation-index"),ylab=gettext("data"),type="p",pch="*")
#            else
#              {

            wylim <- FALSE ## with ylim-Argument
            if("ylim" %in% names(dots)) 
                { wylim <- TRUE
                  options("warn" = -1)
                  ylim1 <- as.matrix(dots[["ylim"]])
                  c1 <- ncol(ylim1); c2 <- ldim0%/%c1; c3 <- ldim0%%c1
                  if(c2>0)
                     ylim0[,1:(c2*c1)] <- ylim1
                  if(c3>0)
                     ylim0[,c2*c1+(1:c3)] <- ylim1[,1:c3]
                  options("warn" = o.warn) }  

            dots["xlab"] <- gettextf("observation-index")
            dots["ylab"] <- gettextf("data")
            dots["type"] <- "p"
            
            cex0 <- rep(1.3, ldim0, length = ldim0) 
            if("cex" %in% names(dots) ) 
                cex0 <- rep(unlist(dots["cex"]), ldim0, length = ldim0) 
            
            pch0 <- rep("*", ldim0, length = ldim0) 
            if("pch" %in% names(dots) ) 
                pch0 <- rep(unlist(dots["pch"]), ldim0, length = ldim0) 

            col0 <- rep((colors()[grep("blue",colors())])[65:1], 
                         ldim0, length = ldim0) 
            if("col" %in% names(dots) )
                col0 <- rep(unlist(dots["col"]), ldim0, 
                            length = ldim0) 

            dots[["cex"]] <- cex0
            dots[["pch"]] <- pch0
            dots[["col"]] <- col0
            
            for( i in 1: lrun0)
                   { if (wylim) dots[["ylim"]] <- ylim0[,i]
                     dots[["y"]] <- Data(x)[, dims0[1:ldim0], runs0[i]]                     
                     do.call("matplot", args = dots)
                    }                  
            #   }        
            if(doEnd)
               {dots[["add"]] <- TRUE;
                par(new=T)
                do.call("matplot", args = dots)}
            
            
          })


### Simulation-Class

##setMethod("plot","Simulation",
##          function(x,y=NULL,...){
##            if(is.null(Data(x)))
##             stop("No Data found -> simulate first")
##            
##            y0<-1:runs(x)
##  matplot(y0,Data(x),xlab="run-index",ylab="data",type="p",pch="*",col="blue")
##          })

## changed w.r.t <1.8            

setMethod("plot",signature(x="Simulation", y="missing"), 
           function(x, obs0=1:samplesize(x), dims0=1:obsDim(x), 
                    runs0 = 1:runs(x), ...){

            if(is.null(Data(x)))
               stop("No Data found -> simulate first")
  
           plot(as(x,"Dataclass"), y = NULL, obs0 = obs0, dims0 = dims0, 
                runs0 = runs0, ...)            
          })



### Contsimulation-Class

setMethod("plot",signature(x="Contsimulation", y="missing"), 
           function(x, obs0=1:samplesize(x), dims0=1:obsDim(x), 
                    runs0=1:runs(x), ...){

            dots <- match.call(call = sys.call(sys.parent(1)), 
                               expand.dots = FALSE)$"..."
            doEnd <- FALSE
            if(!is.null(dots[["panel.first"]])) 
                {doEnd<- TRUE
                 dots[["panel.first"]] <- substitute(pf, 
                                         list(pf=dots[["panel.first"]]))
                }
            if(!is.null(dots[["panel.last"]])) 
                {doEnd<- TRUE
                 dots[["panel.last"]] <- substitute(pf, 
                                         list(pf=dots[["panel.last"]]))
                }

            if(is.null(Data(x)))
               stop("No Data found -> simulate first")
            
            if(any(Data(x) == 0)) return("Warning: plot won't work properly")
            
            
            lobs0 <- min(getdistrSimOption("MaxNumberofPlottedObs"), 
                         length(obs0))           
            lrun0 <- min(getdistrSimOption("MaxNumberofPlottedRuns"), 
                         length(runs0))           
            ldim0 <- min(getdistrSimOption("MaxNumberofPlottedObsDims"), 
                         length(dims0))           
            if((lrun0 < length(runs0)) || (ldim0 < length(dims0))||
               (lobs0 < length(obs0)))   
                warning(gettextf(
"your data set is too big; only %i x %i x %i observations x dimensions x runs are plotted", 
                                 lobs0, ldim0,  lrun0)
                       )
            x.id <- array(aperm(aperm(Data(x), c(1,3,2)) * 
                    array(1-ind(x),
                          c(samplesize(x), runs(x), obsDim(x))), c(1,3,2)),
                          c(lobs0, ldim0, lrun0))
            x.id[x.id == 0] <- Inf
            
            x.c <-  array(aperm(aperm(Data(x), c(1,3,2)) * 
                    array(ind(x),
                          c(samplesize(x), runs(x), obsDim(x))), c(1,3,2)),
                          c(lobs0, ldim0, lrun0))
            x.c[x.c == 0] <- Inf
            
      #      get(getOption("device"))()
            o.warn <- getOption("warn")
            on.exit(options("warn"=o.warn))
            opar <- par(no.readonly = TRUE)
            opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
            on.exit(par(opar))
            
            options("warn" = -1)
            par(mfrow=c(1,lrun0))
            
            y0<-1:lobs0
            dots[["x"]] <- y0

            
            ## catch ylims given in ...
            ylim0 <- matrix(2*range(Data.id(x)), 2, lrun0)
            ##  wylim <- FALSE 
            ### is ylim specified? changed: ylim has to be set by default...
            if("ylim" %in% names(dots)) 
                { wylim <- TRUE
                  options("warn" = -1)
                  ylim1 <- as.matrix(dots[["ylim"]])
                  c1 <- ncol(ylim1); c2 <- ldim0%/%c1; c3 <- ldim0%%c1
                  if(c2>0)
                     ylim0[,1:(c2*c1)] <- ylim1
                  if(c3>0)
                     ylim0[,c2*c1+(1:c3)]<- ylim1[,1:c3]
                  options("warn" = o.warn) }  
            
            dots["xlab"] <- gettextf("observation-index")
            dots["ylab"] <- gettextf("data")
            dots["type"] <- "p"
            
            cex.id0 <- rep(1.3, ldim0, length = ldim0) 
            if("cex.id" %in% names(dots) )
                cex.id0 <- rep(unlist(dots["cex.id"]), ldim0, length = ldim0) 
            
            cex.c0 <- rep(0.8, ldim0, length = ldim0) 
            if("cex.c" %in% names(dots) )
                cex.c0 <- rep(unlist(dots["cex.c"]), ldim0, length = ldim0) 

            pch.id0 <- rep("*", ldim0, length = ldim0) 
            if("pch.id" %in% names(dots) )
                pch.id0 <- rep(unlist(dots["pch.id"]), ldim0, length = ldim0) 

            pch.c0 <- rep("x", ldim0, length = ldim0) 
            if("pch.c" %in% names(dots) )
                pch.c0 <- rep(unlist(dots["pch.c"]), ldim0, length = ldim0) 

            col.id0 <- rep((colors()[grep("blue",colors())])[65:1],
                           ldim0, length = ldim0) 
            if("col.id" %in% names(dots))
                col.id0 <- rep(unlist(dots["col.id"]), ldim0, length = ldim0) 
            
            col.c0 <- rep((colors()[grep("red",colors())]),
                          ldim0, length = ldim0) 
            if("col.c" %in% names(dots))
                col.c0 <- rep(unlist(dots["col.c"]), ldim0, length = ldim0) 

            if(!("add" %in% names(dots))) {
#                myadd <- dots["add"]; dots["add"] <- NULL
            } else dots[["add"]] <- TRUE
            
#            plot.new()
            for( i in 1: lrun0)
                   { ### if(wylim) 
                     dots[["ylim"]] <- ylim0[,i]
                     dots[["y"]] <- x.id[, dims0[1:ldim0], runs0[i]]
                     dots[["cex"]] <- cex.id0
                     dots[["pch"]] <- pch.id0
                     dots[["col"]] <- col.id0
                     do.call("matplot", args = dots)
                   
                    if(any(x.c[,dims0[1:ldim0],runs0[i]] != Inf)) 
                       { dots[["cex"]] <- cex.c0
                         dots[["pch"]] <- pch.c0
                         dots[["col"]] <- col.c0
                         dots[["y"]] <- x.c[, dims0[1:ldim0], runs0[i]]
                         do.call("matpoints", args = dots)                                              
                       }   
                   }                  
            #   }        
            if(doEnd)
               {dots[["add"]] <- TRUE;
                par(new=T)
                do.call("matplot", args = dots)}
            
          })

