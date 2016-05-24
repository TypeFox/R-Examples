#' Produce summary values for objects of class \code{clme}
#'
#' @description Summarizes the output of objects of class \code{clme}, such as those produced by \code{\link{clme}}.
#'
#' @param object an object of class \code{clme}.
#' @param nsim the number of bootstrap samples to use for inference.
#' @param seed the value for the seed of the random number generator.
#' @param verbose vector of logicals. First element will print progress for bootstrap test,
#'        second element is passed to the EM algorithm for every bootstrap sample.
#' @param ... additional arguments passed to other functions.
#'
#'
#' @return
#' The output of \code{summary.clme} is an object of the class \code{summary.clme}. This is a list
#' containing the input object (of class \code{clme}), along with elements:
#' \item{\code{p.value}}{ p-value for the global hypothesis}
#' \item{\code{p.value.ind}}{ p-values for each of the constraints}
#'
#' 
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' \dontrun{
#'   set.seed( 42 )
#'   data( rat.blood )
#'   cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#'   clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                    constraints = cons, seed = 42, nsim = 10)
#'   
#'   summary( clme.out )
#' }
#' 
#' 
#' @method summary clme
#' @export
#' 
summary.clme <- function( object, nsim=1000, seed=42, verbose=c(FALSE,FALSE), ... ){
  
  if( !is.clme(object) ){ stop("'object' is not of class clme")}
  
  ## Extract some values from the fitted object
  cust.const  <- object$cust.const
  search.grid <- object$search.grid 
  MNK         <- dim( search.grid )[1]  
  if( cust.const==TRUE ){
    loop.const <- object$constraints
    MNK        <- 1
  }
  
  mmat <- model_terms_clme( formula(object), data=object$dframe, ncon=object$ncon )
  P1 <- mmat$P1
  X1 <- mmat$X1
  X2 <- mmat$X2
  U  <- mmat$U
  
  Qs      <- object$gran
  tsf     <- object$tsf
  tsf.ind <- object$tsf.ind
  mq.phi  <- object$mq.phi
  
  est_const <- object$constraints
  est_order <- object$order$est_order
  
  if( length(verbose)<2 ){
    verbose <- c(verbose, rep(FALSE, 2-length(verbose) ) )
  }
  
  # Pick up some values from the input object if they are there
  if( !is.null(object$nsim) ){
    nsim2 <- eval( object$nsim )
  } else{
    nsim2       <- nsim
    object$nsim <- nsim
  }
  
  if( !is.null(object$seed) ){
    seed2 <- eval( object$seed )
  } else{
    seed2       <- seed
    object$seed <- seed
  }
  
  mr <- clme_resids( formula=formula(object), data=object$dframe, 
                     gfix=object$gfix_group, ncon=object$ncon )
    
  ## This is the loop for the bootstrap simulations
  # Eventually this should be (optionally) parallelized
  if( nsim2 > 0 ){
    ## Obtain bootstrap samples      
    Y_boot <- resid_boot( formula=formula(object), data=object$dframe, gfix=object$gfix_group, 
                          eps=mr$PA, xi=mr$xi, ssq=mr$ssq, tsq=mr$tsq, 
                          cov.theta=mr$cov.theta, nsim=nsim2, 
                          theta=object$theta.null, mySolver=object$mySolver,
                          seed=seed2, null.resids=FALSE, ncon=object$ncon  )
    
    ## EM for the bootstrap samples    
    p.value  <- rep( 0 , length(object$ts.glb) )
    pval.ind <- rep( 0 , nrow(est_const$A) )
    
    mprint <- round( seq( 1 , round(nsim2*0.9), length.out=10 ) )
    
    for( m in 1:nsim2 ){
      
      if( verbose[1]==TRUE & (m %in% mprint) ){
        print( paste( "Bootstrap Iteration " , m , " of " , nsim2 , sep=""))
      }
      
      ## Loop through the search grid
      ts.boot <- -Inf
      
      for( mnk in 1:MNK ){
        if( cust.const==FALSE ){
          grid.row <- list( order=search.grid[mnk,1], node=search.grid[mnk,3],
                            decreasing=search.grid[mnk,2] )
          loop.const <- create.constraints( P1=P1, constraints=grid.row )
        }
        
        clme.temp <- clme_em( Y=Y_boot[,m], X1=X1, X2=X2, U=U, Nks=object$gfix,
                              Qs=Qs, constraints=loop.const, mq.phi=mq.phi,
                              tsf=tsf, tsf.ind=tsf.ind, mySolver=object$mySolver,
                              verbose=verbose[2], ...)
        
        idx <- which(clme.temp$ts.glb > ts.boot)
        if( length(idx)>0 ){
          ts.boot[idx] <- clme.temp$ts.glb[idx]
        }
        
        update.ind <- (MNK==1) + (mnk == est_order)
        if( update.ind>0 ){
          ts.ind.boot <- clme.temp$ts.ind 
        }
      }
      p.value  <- p.value  + 1*( ts.boot    >= object$ts.glb )
      pval.ind <- pval.ind + 1*(ts.ind.boot >= object$ts.ind )
    }
    
    object$p.value     <- p.value/nsim2
    object$p.value.ind <- pval.ind/nsim2
    
    ## End of the SEQUENTIAL BOOTSTRAP LOOP
  } else{
    object$p.value     <- NA
    object$p.value.ind <- rep( NA, nrow(est_const$A) )
  }
  
  ## Collect the results and return the object
  
  class(object) <- "summary.clme"
  return(object)
}



#' S3 method to print a summary for objects of class \code{clme}
#'
#' @description Summarizes the output of objects of class \code{clme}, such as those produced by \code{\link{clme}}. Prints a tabulated display of global and individual tests, as well as parameter estimates.
#'
#' @param x an object of class \code{clme}.
#' @param alpha level of significance.
#' @param digits number of decimal digits to print.
#' @param ... additional arguments passed to other functions.
#'
#' @note 
#' The individual tests are performed on the specified order. If no specific order was specified, then the individual tests are performed on the estimated order.
#' 
#' @return
#' \code{NULL}, just prints results to the console.
#' 
#' @seealso
#' \code{\link{CLME-package}}
#' \code{\link{clme}}
#' 
#' @examples
#' \dontrun{
#'   set.seed( 42 )
#'   data( rat.blood )
#'   cons <- list(order = "simple", decreasing = FALSE, node = 1 )
#'   clme.out <- clme(mcv ~ time + temp + sex + (1|id), data = rat.blood , 
#'                    constraints = cons, seed = 42, nsim = 10)
#'   
#'   summary( clme.out )
#' }
#' 
#' @importFrom stringr str_pad
#' @importFrom stringr str_trim
#' @importFrom prettyR decimal.align
#' 
#' @method print summary.clme
#' @export
#' 
print.summary.clme <- function( x, alpha=0.05, digits=4, ...){
  
  object <- x
  
  if( class(object)=="summary.clme" ){
    class(object) <- "clme"  
  }
  
  ## Title and formula
  cat( "Linear mixed model subject to order restrictions\n" )
  cat( "Formula: ")
  print( object$formula )
  
  if( object$order$order=="unconstrained" ){
    cat( paste0("\nNo order restrictions (two-tailed alternatives)") )
  } else{
    ## Order statement
    if( object$order$order=="simple" ){
      order <- "simple order"
    }
    if( object$order$order=="umbrella" ){
      order <- paste0("umbrella order with node at ", object$order$node )
    }
    if( object$order$order=="simple.tree" ){
      order <- paste0("tree order with node at ", object$order$node )
    }
    
    if( object$order$order == "custom" ){
      cat( "\nCustom order constraints were provided" )
    } else{
      if( object$order$estimated ){
        ## Estimated the order
        cat( paste0("\nOrder estimated: " , object$order$inc.dec , " ", order ) )
      } else{
        cat( paste0("\nOrder specified: " , object$order$inc.dec , " ", order ) )
      }
    }
  }
  
    
  ## Diagnostic criterion
  crit <- c(logLik(object),
            AIC(object),
            BIC(object) )
  critc <- format( crit , digits=4)
  cat( "\n\nlog-likelihood:", critc[1] )
  cat( "\nAIC:           "  , critc[2] )
  cat( "\nBIC:           "  , critc[3] )
  cat( "\n(log-likelihood, AIC, BIC computed under normality)")  
  
  ## Tests
  est    <- fixef(object)
  tnames <- names(est)
  Amat   <- object$constraints$A
  Bmat   <- object$constraints$B

  ## Global tests
  if( object$order$order != "unconstrained" ){
  if( length(object$ts.glb)>1 ){
    glbs <- object$ts.glb
    grow <- matrix( "NA" , nrow=length(glbs), ncol=3 )
    
    for( ii in 1:length(glbs) ){
      #if( is.null(Bmat) ){
      #  glbn <- "Unknown"
      #  glbe <- "Unknown"
      #} else{
      #  glbn <- paste( tnames[Bmat[ii,2]] , "-", tnames[Bmat[ii,1]] , sep=" " )
      #  glbe <- round( est[Bmat[ii,2]] - est[Bmat[ii,1]], digits=3 )
      #}
      ## Estimate will generally be NA anyway (for LRT test), just drop it.
      # grow[ii,] <- c( glbn, glbe, round(object$ts.glb[ii],3) , sprintf("%.4f", object$p.value[ii]) )
      if( is.null(names(glbs)) ){
        glbn <- rep( "Unknown" , length(glbs) )
      } else{
        glbn <- names( object$ts.glb )
      }
      grow[ii,] <- c( glbn[ii], round(object$ts.glb[ii],3) , sprintf("%.4f", object$p.value[ii]) )
      
    }
    
    #colnames( grow ) <- c("Contrast", "Estimate", "Stat", "p-value")
    #grow <- .align_table.clme( grow )
    for( ii in 2:3){
      val1 <- str_trim( decimal.align( grow[,ii]), side="right"  )
      grow[,ii] <- str_pad(val1, width=max(nchar(val1)), side = "right", pad = "0")
    }  
    
    #colnames( grow ) <- c("Contrast", "Estimate", "Statistic", "p-value")
    colnames( grow ) <- c("Contrast", "Statistic", "p-value")
    grow1 <- c(colnames(grow)[1], grow[,1])
    grow1 <- str_pad( grow1, width=max(nchar(grow1)), side = "right", pad = " ")    
    grow2 <- .align_table.clme( grow[,2:3,drop=FALSE] )
    grow <- cbind( grow1[2:length(grow1)] , grow2)
    colnames(grow)[1] <- grow1[1]
    
    cat( "\n\nGlobal tests: ")
    cat( "\n", paste(colnames(grow) , collapse="  ") )
    for( ii in 1:length(glbs) ){
      cat( "\n", paste(grow[ii,] , collapse="  ")     )
    }
    
  } else{
    ## Single global tests
    #if( is.null(Bmat) ){
    #  glbn <- "Unknown"
    #  glbe <- "Unknown"
    #} else{
    #  glbn <- paste( tnames[Bmat[1,2]] , "-", tnames[Bmat[1,1]] )
    #  glbe <- round( est[Bmat[1,2]] - est[Bmat[1,1]], digits=3 )
    #}
    
    #grow <- cbind( glbn, glbe, round(object$ts.glb,3) , sprintf("%.4f", object$p.value) )
    #colnames( grow ) <- c("Contrast", "Estimate", "Statistic", "p-value")
    if( is.null(names(object$ts.glb)) ){
      glbn <- "Unknown"
    } else{
      glbn <- names( object$ts.glb )
    }
    grow <- cbind( glbn, round(object$ts.glb,3) , sprintf("%.4f", object$p.value) )
    
    colnames( grow ) <- c("Contrast", "Statistic", "p-value")
    grow1 <- c(colnames(grow)[1], grow[,1])
    grow1 <- str_pad( grow1, width=max(nchar(grow1)), side = "right", pad = " ")    
    grow2 <- .align_table.clme( grow[,2:3,drop=FALSE] )
    grow <- cbind( grow1[2:length(grow1)] , grow2)
    colnames(grow)[1] <- grow1[1]
    
    cat( "\n\nGlobal test: ")
    cat( "\n", paste0(colnames(grow) , collapse="  ") )   
    cat( "\n", paste0(grow , collapse="  ")     )
  }
  }
  
  ## Individual tests
  glbs <- object$ts.ind
  grow <- matrix( "NA" , nrow=length(glbs), ncol=4 )
  
  for( ii in 1:length(glbs) ){
    glbn <- paste( tnames[Amat[ii,2]] , "-", tnames[Amat[ii,1]] )
    glbe <- round( est[Amat[ii,2]] - est[Amat[ii,1]], digits=3 )
    grow[ii,] <- c( glbn, glbe, round(object$ts.ind[ii],3) , sprintf("%.4f", object$p.value.ind[ii]) )
  }
  
  #colnames( grow ) <- c("Contrast", "Estimate", "Statistic", "p-value")
  #grow <- .align_table.clme( grow[,2:4] )

  for( ii in 2:4){
    
    if( !any(grow[,ii]==rep("NA",nrow(grow))) ){
      val1 <- str_trim( decimal.align( grow[,ii]), side="right"  )
      grow[,ii] <- str_pad(val1, width=max(nchar(val1)), side = "right", pad = "0")   
    }
  }  
  
  colnames( grow ) <- c("Contrast", "Estimate", "Statistic", "p-value")
  grow1 <- c(colnames(grow)[1], grow[,1])
  grow1 <- str_pad( grow1, width=max(nchar(grow1)), side = "right", pad = " ")    
  grow2 <- .align_table.clme( grow[,2:4,drop=FALSE] )
  grow <- cbind( grow1[2:length(grow1)] , grow2)
  colnames(grow)[1] <- grow1[1]
    
  cat( "\n\nIndividual Tests (Williams' type tests): ")
  cat( "\n", paste(colnames(grow) , collapse="  ") )
  for( ii in 1:length(glbs) ){
    cat( "\n", paste(grow[ii,] , collapse="  ")     )
  }

  
  ## Random effects
  cat( "\n\nVariance components: \n")
  print( VarCorr.clme(object) )
  
  ## Fixed effects
  vars  <- diag( vcov.clme(object) )
  CIs   <- confint.clme( object , level=(1-alpha), ...)
  tvals <- cbind(tnames = str_pad(tnames, width=max(nchar(tnames)), side = "right", pad = " "), 
                 cest   = format(est       , digits=4),
                 cvars  = format(sqrt(vars), digits=4),
                 clcl   = format(CIs[,1]   , digits=4),
                 cucl   = format(CIs[,2]   , digits=4))
    
  
  cipct <- round(100*(1-alpha),2)
  colnames(tvals) <- c(" ", "Estimate", "Std. Err",
                       paste0(cipct, "% lower"), paste0(cipct, "% upper"))
  tvals <- .align_table.clme( tvals )
  
  cat( "\n\nFixed effect coefficients (theta): \n")
  cat( paste0(colnames(tvals), collapse="  ") )
  for( ii in 1:length(est) ){
    cat( "\n", paste0( c(tvals[ii,]),  collapse="  ") )
  }
  cat( "\nStd. Errors and confidence limits based on unconstrained covariance matrix")
  
  cat( "\n\nParameters are ordered according to the following factor levels:\n" )
  cat( paste(  names(fixef(object))[1:object$P1], collapse=", ") )
  cat( "\n\nModel based on", paste0(object$nsim), "bootstrap samples" )
  
}


