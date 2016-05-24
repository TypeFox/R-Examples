`predict.mex` <-
function(object, which, pqu = .99, nsim = 1000, trace=10, ...){
	theCall <- match.call()

  # Class can be either mex or bootmex
  theClass <- class(object)[1]
  if (! theClass %in% c("mex", "bootmex")){
      stop("object must have class 'mex' or 'bootmex'")
  }

	if (theClass == "bootmex" ){
      which <- object$which
      migpd <- object$simpleMar
      margins <- object$margins
      constrain <- object$constrain
      dall <- mexDependence( migpd , which=which , dqu=object$dqu, margins = margins, constrain=constrain )
  } else {
      which <- object$dependence$which
      migpd <- object$margins
      margins <- object$dependence$margins
      constrain <- object$dependence$constrain
      dall <- object
  }

	################################################################
  MakeThrowData <- function(dco,z,coxi,coxmi,data){
    ui <- runif( nsim , min=pqu )
    if( margins == "gumbel"){
      y <- -log( -log( ui ) )
      distFun <- function(x) exp(-exp(-x))
    } else if (margins == "laplace"){
      y <- ifelse(ui < .5,  log(2 * ui), -log(2 * (1 - ui) ))
      distFun <- function(x) ifelse(x<0, exp(x)/2, 1-exp(-x)/2)
    }

    z <- as.matrix(z[ sample( 1:( dim( z )[ 1 ] ), size=nsim, replace=TRUE ) ,])
    ymi <- sapply( 1:( dim( z )[[ 2 ]] ) , makeYsubMinusI, z=z, v=dco , y=y )

    xmi <- apply( ymi, 2, distFun )

    xi <- u2gpd( ui, p = 1 - migpd$mqu[ which ], th=migpd$mth[ which ], sigma=coxi[ 1 ], xi = coxi[ 2 ] )

  	for( i in 1:( dim( xmi )[[ 2 ]] ) ){
		  xmi[, i ] <- revTransform( xmi[ ,i ], as.matrix(data[,-which])[, i ],
								             th = migpd$mth[ -which ][ i ],
								             qu = migpd$mqu[ -which ][ i ],
								             sigma=coxmi[ 1,i ], xi=coxmi[ 2,i ] )
	  }
    sim <- data.frame( xi , xmi )
    names( sim ) <- c( colnames( migpd$data )[ which ], colnames( migpd$data )[ -which ])
    sim
  }

	################################################################
  makeYsubMinusI <- function( i, z, v , y ){
			v <- v[ , i ]
			z <- z[ , i ]
			if ( !is.na( v[ 1 ] ) ){
				if( v[ 1 ] < 10^(-5) & v[ 2 ] < 0 ){
					if( v[ 4 ] < 10^(-5 ) ) d <- 0
					else d <- v[ 4 ]
					a <- v[ 3 ] - d * log( y )
				}
				else a <- v[ 1 ] * y
			} # close if( !is.na...
			else a <- NA
			a + ( y^v[ 2 ] ) * z
		}

  ###############################################################
  if (theClass == "bootmex"){
	# The function lfun does most of the work
  lfun <- function( i , bo, pqu, nsim , migpd, which ){
	   if ( i %% trace == 0 ) cat( i, "sets done\n" )

     res <- MakeThrowData(dco=bo[[ i ]]$dependence,z=bo[[ i ]]$Z, coxi = bo[[i]]$GPD[,which],
                          coxmi = as.matrix(bo[[ i ]]$GPD[,-which]),
                          data = bo[[i]]$Y)
	   res
  }

	bootRes <- lapply( 1:length( object$boot ) , lfun ,
	        			    migpd=migpd, pqu=pqu, bo = object$boot, nsim=nsim,
	        			    which = which )
	    # bootRes contains the bootstrap simulated complete vectors X on the original
      # scale of the data, conditional on having the _which_ component above the pqu quantile.
	} else {
    bootRes <- NULL
  }

	##########################################################################
	# Get a sample using the point estimates of the parameters
	# that are suggested by the data

  cox <- coef(migpd)[3:4, which]
  coxmi <- as.matrix(coef(migpd)[3:4, -which])

  sim <- MakeThrowData(dco=dall$dependence$coefficients,z=dall$dependence$Z,coxi=cox,coxmi=coxmi,data=migpd$data)

  m <- 1 / ( 1 - pqu ) # Need to estimate pqu quantile
 	zeta <- 1 - migpd$mqu[ which ] # Coles, page 81
	pth <- migpd$mth[ which ] + cox[ 1 ] / cox[ 2 ] * ( ( m*zeta )^cox[ 2 ] - 1 )

	data <- list( real = data.frame( migpd$data[, which ], migpd$data[, -which] ), simulated = sim, pth=pth)
  names(data$real)[1] <- colnames(migpd$data)[which]

  res <- list( call = theCall , replicates = bootRes, data = data,
				       which = which, pqu = pqu,
				       mth=c( migpd$mth[ which ], migpd$mth[ -which ] ),
               gpd.coef = coef(migpd)[,c(which,(1:dim(data$real)[2])[-which])])

  oldClass( res ) <- "predict.mex"

  res
}

