outlyingPlotIC <- function(data, 
                           IC.x, 
                           IC.y = IC.x,
                           dist.x = NormType(),
                           dist.y, 
                           cutoff.y = cutoff.chisq(), 
                           cutoff.x = cutoff.sememp(),
                           ...,
                           cutoff.quantile.x = 0.95,
                           cutoff.quantile.y = cutoff.quantile.x,
                           id.n,
                           lab.pts, 
                           adj, 
                           cex.idn,
                           col.idn, 
                           lty.cutoff, 
                           lwd.cutoff, 
                           col.cutoff,
                           robCov.x = TRUE,
                           robCov.y = TRUE,
                           tf.x = data,
                           tf.y = data,
                           jitt.fac=10,
                           main = gettext("Outlyingness \n by means of a distance-distance plot")
                           ){
     mc <- as.list(match.call(expand.dots = FALSE))[-1]
     dots <- mc$"..."
     if(is.null(dots$xlim)) dots$xlim <- TRUE
     if(is.null(dots$ylim)) dots$ylim <- TRUE
     if(is.null(mc$cutoff.quantile.x)) mc$cutoff.quantile.x <- 0.95
     if(is.null(mc$cutoff.quantile.y)) mc$cutoff.quantile.y <- cutoff.quantile.x
     if(is.null(mc$cutoff.x)) mc$cutoff.x <- cutoff.sememp()
     if(is.null(mc$cutoff.y)) mc$cutoff.y <- cutoff.chisq()
     if(missing(IC.x)) stop("Argument 'IC.x' must be given as argument to 'outlyingPlot'")
     if(missing(data)) stop("Argument 'data' must be given as argument to 'outlyingPlot'")
  
     if(missing(dist.y)){
       if(robCov.y){
          evIC = evalIC(IC.y,as.matrix(data))
          asVar = solve(getCov(CovMcd(data.frame(evIC[1,],evIC[2,]),alpha=0.5)))
          cat("\n", sep="", gettext("Robust asVar"), ":\n")
          print(asVar)
           }else{
          if("asCov" %in% names(Risks(IC.y)))
             if(is.matrix(Risks(IC.y)$asCov) || length(Risks(IC.y)$asCov) == 1)
                {asVar <- Risks(IC.y)$asCov
                  cat("\n", sep="", gettext("asVar"));# print("HHHH");
                  print(asVar)
                 }
              else{asVar <- Risks(IC.y)$asCov$value 
                   cat("\n", sep="", gettext("asVar"));#print("HHHH");
                   print(asVar)
                 }
              else{asVar <- getRiskIC(IC.y, risk = asCov())$asCov$value
                   cat("\n", sep="", gettext("Classic asVar"));
                 #  print("HHHH");
                 print(asVar)
                 }}
     
        asVar <- PosSemDefSymmMatrix(solve(asVar))
        mc$dist.y <- QFNorm(name = gettext("Mahalonobis-Norm"), QuadForm = asVar)
     }

  if(missing(dist.x)){
        #mc$dist.x <- NormType()
    if(robCov.x){
      evIC = evalIC(IC.x,as.matrix(data))
      asVar = getCov(CovMcd(data.frame(evIC[1,],evIC[2,]),alpha=0.5))
      #cat("\nRobust asVar:") ;print("KKKKK")
      #print(asVar)
      }
     else{
   if("asCov" %in% names(Risks(IC.y)))
      if(is.matrix(Risks(IC.x)$asCov) || length(Risks(IC.y)$asCov) == 1)
               {asVar <- Risks(IC.x)$asCov
                  cat("\n", sep="", gettext("asVar"));
                  #print("KKKKK2");
                  print(asVar)
                  }
         else
               {asVar <- Risks(IC.x)$asCov$value 
                  cat("\n", sep="", gettext("asVar"));
                  # print("KKKKK3");
                  print(asVar)
                  }
         else
            {asVar <- getRiskIC(IC.x, risk = asCov())$asCov$value
                   cat("\n", sep="", gettext("Classic asVar"));
                   #print("KKKKK4");
                   print(asVar)
                   }
       }
    
       asVar <- PosSemDefSymmMatrix(solve(asVar))
       mc$dist.x <- QFNorm(name = gettext("Mahalonobis-Norm"), QuadForm = asVar)
      }

    if(missing(tf.x)){
     tf.x <- function(x) apply(x,2,function(xx) evalIC(IC.x,xx))
     }else{tf.x <- mc$tf.x}
    if(missing(tf.y)){
     tf.y <- function(x) apply(x,2,function(xx) evalIC(IC.y,xx))
     }else{tf.y <- mc$tf.y}

     do.call(ddPlot,args=c(list(data=data),dots, 
       list(dist.x = mc$dist.x,
       dist.y = mc$dist.y, 
       cutoff.x = mc$cutoff.x, 
       cutoff.y = mc$cutoff.y,
       cutoff.quantile.x = mc$cutoff.quantile.x, 
       cutoff.quantile.y = mc$cutoff.quantile.y,
       transform.x = tf.x, 
       transform.y = tf.y,
       id.n = mc$id.n, 
       lab.pts = mc$lab.pts, 
       adj = mc$adj, 
       cex.idn = mc$cex.idn,
       col.idn = mc$col.idn, 
       lty.cutoff = mc$lty.cutoff,
       lwd.cutoff = mc$lwd.cutoff, 
       col.cutoff = mc$col.cutoff, 
       jitt.fac = mc$jitt.fac,
       main = main)))

     }

