"boa.plot" <-
function(type, dev = boa.par("dev"), mfdim = boa.par("plot.mfdim"),
                     newplot = boa.par("plot.new"),
                     onelink = boa.par("plot.onelink"),
                     title = boa.par("title"))
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   drawn <- FALSE
   switch(type,
      "acf"     = { text <- "Sampler Lag-Autocorrelations"
                    foo <- "boa.plot.acf"
                    foo.args <- expression(list(largs[[idx[i]]],
                                                pargs[[idx[i]]]))
                    oneparm <- TRUE
                    onelink <- TRUE },
      "bandg"   = { text <- "Brooks & Gelman Multivariate Shrink Factors"
                    foo <- "boa.plot.bandg"
                    foo.args <- list()
                    oneparm <- FALSE
                    onelink <- FALSE },
      "density" = { text <- "Estimated Posterior Density"
                    foo <- "boa.plot.density"
                    foo.args <- expression(list(largs[[idx[i]]],
                                                pargs[[idx[i]]]))
                    oneparm <- TRUE },
      "gandr"   = { text <- "Gelman & Rubin Shrink Factors"
                    foo <- "boa.plot.gandr"
                    foo.args <- expression(list(pargs[[idx[i]]]))
                    oneparm <- TRUE
                    onelink <- FALSE },
      "geweke"  = { text <- "Geweke Convergence Diagnostic"
                    foo <- "boa.plot.geweke"
                    foo.args <- expression(list(largs[[idx[i]]],
                                                pargs[[idx[i]]]))
                    oneparm <- TRUE
                    onelink <- TRUE },
      "history" = { text <- "Sampler Running Mean"
                    foo <- "boa.plot.history"
                    foo.args <- expression(list(largs[[idx[i]]],
                                                pargs[[idx[i]]]))
                    oneparm <- TRUE },
      "trace"   = { text <- "Sampler Trace"
                    foo <- "boa.plot.trace"
                    foo.args <- expression(list(largs[[idx[i]]],
                                                pargs[[idx[i]]]))
                    oneparm <- TRUE },
      { foo <- NULL
        cat("Warning: plot type does not exist\n") }
   )
   if(is.character(foo)) {
      work <- boa.chain("work")
      lnames <- names(work)
      largs <- pargs <- list(0)
      pidx <- NULL
      for(i in lnames) {
         pnames <- boa.pnames(work[[i]])
         if(oneparm) {
            for(j in pnames) {
                if(onelink) {
                   n <- length(largs)
                   largs[[n + 1]] <- i
                   pargs[[n + 1]] <- j
                   pidx <- c(pidx, paste(j, i))
                } else if(is.element(j, pidx)) {
                   largs[[j]] <- c(largs[[j]], i)
                } else {
                   largs[[j]] <- i
                   pargs[[j]] <- j
                   pidx <- c(pidx, j)
                }
            }
         } else if(length(pidx) > 0) {
            largs[[2]] <- c(largs[[2]], i)
            pargs[[2]] <- union(pargs[[2]], pnames)
         } else {
            largs[[2]] <- i
            pargs[[2]] <- pnames
            pidx <- 2
         }
      }
      largs[[1]] <- pargs[[1]] <- NULL

      if(!newplot)  boa.plot.close(boa.par("dev.list"))

      idx <- order(pidx)
      n <- length(idx)
      size <- prod(mfdim)
      newdim <- mfdim
      imin <- ifelse(mfdim[1] <= mfdim[2], 1, 2)
      imax <- imin %% 2 + 1
      ratio <- mfdim[imin] / mfdim[imax]
      for(i in 1:n) {
         if((size == 1) || ((i %% size) == 1)) {
            boa.plot.open()
            nleft <- n - i + 1
            while((nleft <= prod(newdim)) && (i == 1)) {
               if(newdim[1] > 1) {
                  mfdim <- c(ceiling(nleft / newdim[2]), newdim[2])
                  newdim[imax] <- newdim[imax] - 1
                  newdim[imin] <- round(ratio * newdim[imax])
               } else {
                  mfdim <- c(1, nleft)
                  newdim <- 0
               }
            }
            boa.plot.par(mfdim, title)
         }
         drawn <- do.call(foo, args = eval(foo.args)) || drawn
         if(title) boa.plot.title(text)
      }
   }

   return(drawn)
}
