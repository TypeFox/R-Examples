#######################################################################
#######################################################################
## Method for combining the results of repeated-imputation analyses
#######################################################################
#######################################################################
miInference <- function( est.list, std.err.list, 
   method = "scalar", df.complete = NULL){
   #########################################
   # function for combining the results 
   # from a multiply-imputed analysis
   #########################################
   # check args
   if( !is.list( est.list ) )
      stop("Argument est.list must be a list.")
   if( !is.list( std.err.list ) )
      stop("Argument std.err.list must be a list.")
   if( length( est.list ) != length( std.err.list ) )
      stop("Arguments est.list and std.err.list have unequal length.")
   M <- length( est.list )
   if( M < 2 )
      stop("Length of est.list and std.err.list must be >= 2.")
   d <- length( est.list[[1]] )
   for( i in 1:M ){
      if( ( length( est.list[[i]] ) != d ) |
          ( length(  std.err.list[[i]] ) != d ) )
         stop("Elements of est.list or std.err.list not of equal length.")
      if( ( mode( est.list[[i]] ) != "numeric" ) |
          ( mode(  std.err.list[[i]] ) != "numeric" ) )
         stop("Non-numeric data in est.list or std.err.list.")
   }
   #
   if( method != "scalar" ) 
      stop( "Value for argument \"method\" not recognized." )
   if( method == "scalar" ){
      # reshape estimates into a matrix
      Q <- est.list[[1]]
      for( i in 2:M ){
         Q <- rbind( Q, est.list[[i]] )
      }
      colnames( Q ) <- names( est.list[[1]] )
      Qbar <- apply( Q, 2, mean)
      #
      # reshape standard errors^2 into a matrix
      U <- std.err.list[[1]]
      for( i in 2:M ){
         U <- rbind( U, std.err.list[[i]] )
      }
      U <- U^2
      colnames( U ) <- colnames( Q )
      # 
      # apply Rubin's rules
      Ubar <- apply( U, 2, mean )
      Bm <- apply( Q, 2, var)
      Tm <- Ubar + ( 1 + 1/M ) * Bm
      r <- ( 1 + 1/M )* Bm / Ubar
      mis.inf <- ( 1 + 1/M )* Bm / Tm
      df.rubin <- ( M - 1 ) / ( mis.inf^2 )
      if( is.null( df.complete )){
         df <- df.rubin
      }
      else{
         lambda <- ( df.complete + 1 ) / ( df.complete + 3 )
         dfhat <- lambda * df.complete * ( 1 - mis.inf )
         df <- 1 / ( 1/df.rubin + 1/dfhat )
         df[ df.complete == Inf ] <- df.rubin[ df.complete == Inf ]
      }
      p <- 2 * ( 1 - pt( abs( Qbar / sqrt(Tm) ), df ) )
      # 
      # return result
      result <- list(
         names = colnames(Q),
         est = Qbar,
         std.err = sqrt(Tm),
         df = df,
         p = p,
         rel.incr = r,
         mis.inf = mis.inf,
         method = method )
   }
   class( result ) <- "miInference"
   return( result )}
#######################################################################
print.miInference <- function( x, ... ){
   ############################################
   # check argument
   if( class(x) != "miInference" )
      stop("Argument should be of class \"miInference\".")
   ###############################################
   if( x$method == "scalar" ){
      coef.table <- cbind(
         Est = signif( x$est, 5 ),
         SE = signif( x$std.err, 5 ),
         "Est/SE" = round( x$est / x$std.err, 3 ),
         df = round( x$df, 1 ),
         p = round( x$p, 3 ),
         "Pct.mis" = round( 100*x$mis.inf, 1) )
      rownames( coef.table ) <- x$names
      print( coef.table )
   }
   else{
      stop( "x$method not recognized.")
   }
   return( invisible(x) )}
#######################################################################
