"plot.fluct" <-
  function(x, col = NULL, signif = FALSE, p.cutoff = 0.005, 
           q.cutoff = 0.04, s.cutoff = 5, n.cutoff = 2, mean = FALSE, polygon = FALSE, 
           ncore = NULL, ...) {

    ## check input data
    if(is.vector(x))
      x = matrix(x, nrow=1)
    if(!is.matrix(x) || !is.numeric(x))
      stop("provide a numeric matrix or vector")
   
    ## check colors, which also define groups of input data if signif=TRUE
    if(is.null(col))
      col <- seq(1, nrow(x))
    if(length(col) != nrow(x)) 
       stop("length of col doesn't match dimension of x")
    if(any(is.na(col))) {
       x = x[!is.na(col), ]
       col = col[!is.na(col)]
    }

    ## check for significance calculation
    if(signif) {
       ns <- table(col)
       inds.signif <- which(ns >= s.cutoff)
       if(length(inds.signif) < 2) { 
         warning("Insufficient samples to calculate significance")
         signif = FALSE
      }
    }

    ## extract some values from '...' since we still do some plots here   
    ## These could be removed after merging this function with plot.bio3d()
    dots = list(...)
    if("rm.gaps" %in% names(dots)) rm.gaps = dots$rm.gaps
    else rm.gaps = formals(plotb3)$rm.gaps
    
    ## gaps positions
    gaps.pos <- gap.inspect(x)

    if(rm.gaps) yvals = x[, gaps.pos$f.inds, drop=FALSE]
    else yvals = x

    if("ylim2zero" %in% names(dots)) ylim2zero = dots$ylim2zero
    else ylim2zero = formals(plotb3)$ylim2zero

    if("xlim" %in% names(dots)) xlim = dots$xlim
    else xlim = c(1, ncol(yvals))

    if("ylim" %in% names(dots)) ylim = dots$ylim
    else ylim = range(yvals, na.rm = TRUE)

    if(ylim2zero) ylim[1] = 0

    if(! "ylab" %in% names(dots)) dots$ylab = "Fluctuation"
    dots$xlim = xlim
    dots$ylim = ylim
    #####################################################################

    if(signif) {

      ncore = setup.ncore(ncore)

#      op = par(no.readonly = TRUE)
      op = par()$new
      on.exit(par(new=op))
 
      pairs <- pairwise(length(inds.signif))

      ## get p-value and q-value for each non-gap position
      p.all <- mclapply(gaps.pos$f.inds, function(i) {
         p.i <- apply(pairs, 1, function(j) {
            inds1 <- which(col == names(ns)[inds.signif[j[1]]])
            inds2 <- which(col == names(ns)[inds.signif[j[2]]])
            p = t.test(x[inds1, i], x[inds2,i], alternative="two.sided")$p.value
            q <- abs(mean(x[inds1, i]) - mean(x[inds2, i]))
            c(p, q)
         }) 
         c(p=min(p.i[1, ]), q=p.i[2, which.min(p.i[1,])])
      })

      ## p-values with gaps inserted
      pvalue <- rep(NA, ncol(x))
      pvalue[gaps.pos$f.inds] <- sapply(p.all, "[", "p")

      ## q-values, i.e. difference of mean values, with gaps inserted
      qvalue <- rep(NA, ncol(x))
      qvalue[gaps.pos$f.inds] <- sapply(p.all, "[", "q")

      if(rm.gaps){
         pvalue = pvalue[gaps.pos$f.inds]
         qvalue = qvalue[gaps.pos$f.inds]
      }

      sig <- which(pvalue<=p.cutoff & qvalue >= q.cutoff)

      ## - start plotting
      if(length(sig) > 0) {
         ## Plot significance as shaded blocks
         bds <- bounds(sig)
         ii <- which(bds[, "length"] >= n.cutoff)

         if(length(ii) > 0) {
 
            plot.new()
            plot.window(xlim=xlim, ylim=ylim)
   
            ## to show bricks for single site significance
            adjust = 0.1
            rect(bds[ii,1]-adjust, rep(ylim[1], length(ii)), bds[ii,2]+adjust,
                 rep(ylim[2], length(ii)),
                 col=rep("lightblue", length(ii)), border=NA)
            
            ## add this for plot.bio3d on the same device
            par(new=TRUE)
         }
      }
    }

    if(mean) {
       # calculate mean values and replace
       yvals = apply(x, 2, tapply, col, mean, na.rm=TRUE)
       col = unique(col)

       if(!is.matrix(yvals))
          yvals = matrix(yvals, nrow=1)
       else 
          yvals = yvals[col, , drop=FALSE]  # correct order change due to tapply

       # still keep the same gaps in first row
       # this will help plot SSE in plot.bio3d()
       yvals[1, is.na(x[1, ])] = NA

       x = yvals
       if(rm.gaps) 
          yvals = x[, gaps.pos$f.inds, drop=FALSE]

       # trick to leave gap position unchanged.
       # Won't affect plot because plot.bio3d() only picks up the first row
       # All plots in this function should be done with yvals!!
       x = rbind(x, gap.mark=rep(0, ncol(x))) 
       x["gap.mark", gaps.pos$t.inds] = NA
    }

    ## Plot fluctuations
    if(polygon) {
       dots$type = "n" 
       do.call(plot.bio3d,  c(list(x=x), dots))
     
       xx = yvals[1, ] 

       ylim2 = range(xx, na.rm = TRUE)
       if(ylim2zero) ylim2[1] = 0

       n = bounds(which(is.na(xx)))
       if(length(n)>0) xx[n[, 1:2]] = ylim2[1]
       
       # color for polygon
       n.col = do.call(rgb, c(as.list(col2rgb(col[1])/255), list(alpha=0.4)))

       polygon(c(1, seq_along(xx), length(xx)), c(ylim2[1], xx, ylim2[1]),
          col = n.col, border=NA) 
    } else {

       do.call(plot.bio3d, c(list(x=x), dots))

    }
   
    ## Plot all lines 
    for(i in 1:nrow(yvals)) lines(yvals[i, ], col=col[i], lwd=2)

    if(signif)
       out <- list(signif=sig)
    else
       out <- NULL 
    invisible(out)
  }
