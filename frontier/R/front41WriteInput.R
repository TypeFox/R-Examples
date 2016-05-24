front41WriteInput <- function( data, crossSectionName, timePeriodName = NULL,
   yName, xNames = NULL, qxNames = NULL, zNames = NULL, quadHalf = TRUE,
   modelType = ifelse( is.null( zNames ), 1, 2 ), functionType = 1,
   logDepVar = TRUE, mu = FALSE, eta = FALSE,
   insFile = "front41.ins", dtaFile = sub( "\\.ins$", ".dta", insFile ),
   outFile = sub( "\\.ins$", ".out", insFile ), startUpFile = "front41.000",
   iprint = 5, indic = 1, tol = 0.00001, tol2 = 0.001, bignum = 1.0E+16,
   step1 = 0.00001, igrid2 = 1, gridno = 0.1, maxit = 100, ite = 1 ) {

   if( qxNames == "all" && !is.null( qxNames ) ) {
      qxNames <- xNames
   }
   if( !is.character( qxNames ) && !is.null( qxNames ) ) {
      stop( "argument 'qxNames' must be either logical or a vector of strings" )
   }

   checkNames( c( crossSectionName, timePeriodName, yName, xNames, zNames,
      qxNames ), names( data ) )

   if( !modelType %in% c( 1, 2 ) ) {
      stop( "argument 'modelType' must be either 1 or 2" )
   }
   if( !functionType %in% c( 1, 2 ) ) {
      stop( "argument 'functionType' must be either 1 or 2" )
   }
   if( !is.logical( logDepVar ) ) {
      stop( "argument 'logDepVar' must be logical" )
   }
   if( !is.logical( mu ) ) {
      stop( "argument 'mu' must be logical" )
   }
   if( modelType == 1 ) {
      if( !is.logical( eta ) ) {
         stop( "argument 'eta' must be logical" )
      }
   }
   # iprint
   if( !is.numeric( iprint ) ) {
      stop( "argument 'iprint' must be numeric" )
   } else if( iprint != round( iprint ) ) {
      stop( "argument 'iprint' must be an iteger" )
   } else if( iprint < 0 ) {
      stop( "argument 'iprint' must be non-negative" )
   }
   iprint <- as.integer( iprint )
   # indic
   if( !is.numeric( indic ) ) {
      stop( "argument 'indic' must be numeric" )
   } else if( indic != round( indic ) ) {
      stop( "argument 'indic' must be an integer" )
   }
   indic <- as.integer( indic )
   # tol
   if( !is.numeric( tol ) ) {
      stop( "argument 'tol' must be numeric" )
   } else if( tol < 0 ) {
      stop( "argument 'tol' must be non-negative" )
   }
   # tol2
   if( !is.numeric( tol2 ) ) {
      stop( "argument 'tol2' must be numeric" )
   } else if( tol2 < 0 ) {
      stop( "argument 'tol2' must be non-negative" )
   }
   # bignum
   if( !is.numeric( bignum ) ) {
      stop( "argument 'bignum' must be numeric" )
   } else if( bignum <= 0 ) {
      stop( "argument 'bignum' must be positive" )
   }
   # step1
   if( !is.numeric( step1 ) ) {
      stop( "argument 'step1' must be numeric" )
   } else if( step1 <= 0 ) {
      stop( "argument 'step1' must be positive" )
   }
   # igrid2
   if( ! igrid2 %in% c( 1, 2 ) ) {
      stop( "argument 'igrid2' must be either '1' or '2'" )
   }
   # gridno
   if( !is.numeric( gridno ) ) {
      stop( "argument 'gridno' must be numeric" )
   } else if( gridno <= 0 ) {
      stop( "argument 'gridno' must be positive" )
   }
   # maxit
   if( !is.numeric( maxit ) ) {
      stop( "argument 'maxit' must be numeric" )
   } else if( maxit != round( maxit ) ) {
      stop( "argument 'maxit' must be an integer" )
   } else if( maxit <= 0 ) {
      stop( "argument 'maxit' must be positive" )
   }
   maxit <- as.integer( maxit )
   # ite
   if( ! ite %in% c( 0, 1 ) ) {
      stop( "argument 'ite' must be either '0' or '1'" )
   }


   nCrossSection <- length( unique( data[[ crossSectionName ]] ) )
   nTimePeriods  <- ifelse( is.null( timePeriodName ), 1,
      length( unique( data[[ timePeriodName ]] ) ) )
   nTotalObs     <- nrow( data )
   nXvars        <- length( xNames )
   nTLvars       <- length( qxNames )
   nXtotal       <- nXvars + nTLvars * ( nTLvars + 1 ) / 2
   nZvars        <- length( zNames )

   if( modelType == 2 ) {
      eta <- nZvars
   } else {
      eta <- ifelse( eta, "y", "n" )
   }

   commentRow <- max( 16, nchar( dtaFile ) + 1 )

   cat( modelType, rep( " ", commentRow - 1 ),
      "1=ERROR COMPONENTS MODEL, 2=TE EFFECTS MODEL\n",
      file = insFile, sep = "" )
   cat( dtaFile, rep( " ", commentRow - nchar( dtaFile ) ),
      "DATA FILE NAME\n", file = insFile, append = TRUE, sep = "" )
   cat( outFile, rep( " ", commentRow - nchar( outFile ) ),
      "OUTPUT FILE NAME\n", file = insFile, append = TRUE, sep = "" )
   cat( functionType, rep( " ", commentRow - 1 ),
      "1=PRODUCTION FUNCTION, 2=COST FUNCTION\n",
      file = insFile, append = TRUE, sep = "" )
   cat( ifelse( logDepVar, "y", "n" ), rep( " ", commentRow - 1 ),
      "LOGGED DEPENDENT VARIABLE (Y/N)\n",
      file = insFile, append = TRUE, sep = "" )
   cat( nCrossSection,
      rep( " ", commentRow - nchar( as.character( nCrossSection ) ) ),
      "NUMBER OF CROSS-SECTIONS\n",
      file = insFile, append = TRUE, sep = "" )
   cat( nTimePeriods,
      rep( " ", commentRow - nchar( as.character( nTimePeriods ) ) ),
      "NUMBER OF TIME PERIODS\n",
      file = insFile, append = TRUE, sep = "" )
   cat( nTotalObs,
      rep( " ", commentRow - nchar( as.character( nTotalObs ) ) ),
      "NUMBER OF OBSERVATIONS IN TOTAL\n",
      file = insFile, append = TRUE, sep = "" )
   cat( nXtotal,
      rep( " ", commentRow - nchar( as.character( nXtotal ) ) ),
      "NUMBER OF REGRESSOR VARIABLES (Xs)\n",
      file = insFile, append = TRUE, sep = "" )
   cat( ifelse( mu, "y", "n" ), rep( " ", commentRow - 1 ),
      "MU (Y/N) [OR DELTA0 (Y/N) IF USING TE EFFECTS MODEL]\n",
      file = insFile, append = TRUE, sep = "" )
   cat( eta, rep( " ", commentRow - nchar( as.character( eta ) ) ),
      "ETA (Y/N) [OR NUMBER OF TE EFFECTS REGRESSORS (Zs)]\n",
      file = insFile, append = TRUE, sep = "" )
   cat( "n", rep( " ", commentRow - 1 ),
      "STARTING VALUES (Y/N)\n",
      file = insFile, append = TRUE, sep = "" )

   ## create table for data
   # cross section identifier
   dataTable <- matrix( data[[ crossSectionName ]], ncol = 1 )

   # time period identifier
   if( is.null( timePeriodName ) ) {
      dataTable <- cbind( dataTable, rep( 1, nrow( dataTable ) ) )
   } else {
      dataTable <- cbind( dataTable, data[[ timePeriodName ]] )
   }

   # endogenous variable
   dataTable <- cbind( dataTable, data[[ yName ]] )

   # exogenous variables
   if( nXvars > 0 ) {
      for( i in 1:nXvars ) {
         dataTable <- cbind( dataTable, data[[ xNames[ i ] ]] )
      }
   }

   # exogenous variables: quadratic and interaction terms
   if( nTLvars > 0 ) {
      for( i in 1:nTLvars ) {
         for( j in i:nTLvars ) {
            dataTable <- cbind( dataTable,
               ifelse( i == j, 1 , 2 ) * ifelse( quadHalf, 0.5, 1 ) *
               data[[ qxNames[ i ] ]] * data[[ qxNames[ j ] ]] )
         }
      }
   }

   # variables explaining the efficiency level
   if( nZvars > 0 ) {
      for( i in 1:nZvars ) {
         dataTable <- cbind( dataTable, data[[ zNames[ i ] ]] )
      }
   }

   # write data file to disk
   write.table( dataTable, file = dtaFile, row.names = FALSE,
      col.names = FALSE, sep = "\t" )

   ## create start-up file
   if( !is.null( startUpFile ) ) {
      cat( "KEY VALUES USED IN FRONTIER PROGRAM (VERSION 4.1)\n",
         file = startUpFile )
      cat( "NUMBER:         DESCRIPTION:\n",
         file = startUpFile, append = TRUE )
      cat( iprint,
         rep( " ", 16 - nchar( as.character( iprint ) ) ),
         "IPRINT - PRINT INFO EVERY \"N\" ITERATIONS, 0=DO NOT PRINT\n",
         file = startUpFile, append = TRUE, sep = "" )
      cat( indic,
         rep( " ", 16 - nchar( as.character( indic ) ) ),
         "INDIC - USED IN UNIDIMENSIONAL SEARCH PROCEDURE - SEE BELOW\n",
         file = startUpFile, append = TRUE, sep = "" )
      tolString <- sub( "e", "D", format( tol, scientific = 2 ) )
      cat( tolString,
         rep( " ", 16 - nchar( tolString ) ),
         "TOL - CONVERGENCE TOLERANCE (PROPORTIONAL)\n",
         file = startUpFile, append = TRUE, sep = "" )
      tol2String <- sub( "e", "D", format( tol2, scientific = 2 ) )
      cat( tol2String,
         rep( " ", 16 - nchar( tol2String ) ),
         "TOL2 - TOLERANCE USED IN UNI-DIMENSIONAL SEARCH PROCEDURE\n",
         file = startUpFile, append = TRUE, sep = "" )
      bignumString <- sub( "e", "D", format( bignum, scientific = 2 ) )
      cat( bignumString,
         rep( " ", 16 - nchar( bignumString ) ),
         "BIGNUM - USED TO SET BOUNDS ON DEN & DIST\n",
         file = startUpFile, append = TRUE, sep = "" )
      step1String <- sub( "e", "D", format( step1, scientific = 2 ) )
      cat( step1String,
         rep( " ", 16 - nchar( step1String ) ),
         "STEP1 - SIZE OF 1ST STEP IN SEARCH PROCEDURE\n",
         file = startUpFile, append = TRUE, sep = "" )
      cat( igrid2,
         rep( " ", 16 - nchar( as.character( igrid2 ) ) ),
         "IGRID2 - 1=DOUBLE ACCURACY GRID SEARCH, 0=SINGLE\n",
         file = startUpFile, append = TRUE, sep = "" )
      cat( gridno,
         rep( " ", 16 - nchar( as.character( gridno ) ) ),
         "GRIDNO - STEPS TAKEN IN SINGLE ACCURACY GRID SEARCH ON GAMMA\n",
         file = startUpFile, append = TRUE, sep = "" )
      cat( maxit,
         rep( " ", 16 - nchar( as.character( maxit ) ) ),
         "MAXIT - MAXIMUM NUMBER OF ITERATIONS PERMITTED\n",
         file = startUpFile, append = TRUE, sep = "" )
      cat( ite,
         rep( " ", 16 - nchar( as.character( ite ) ) ),
         "ITE - 1=PRINT ALL TE ESTIMATES, 0=PRINT ONLY MEAN TE\n",
         file = startUpFile, append = TRUE, sep = "" )
      cat( "\n",
         file = startUpFile, append = TRUE )
      cat( "THE NUMBERS IN THIS FILE ARE READ BY THE FRONTIER PROGRAM WHEN IT BEGINS\n",
         file = startUpFile, append = TRUE )
      cat( "EXECUTION. YOU MAY CHANGE THE NUMBERS IN THIS FILE IF YOU WISH. IT IS\n",
         file = startUpFile, append = TRUE )
      cat( "ADVISED THAT A BACKUP OF THIS FILE IS MADE PRIOR TO ALTERATION.\n",
         file = startUpFile, append = TRUE )
      cat( "\n",
         file = startUpFile, append = TRUE )
      cat( "FOR MORE INFORMATION ON THESE VARIABLES SEE: COELLI (1996), CEPA WORKING\n",
         file = startUpFile, append = TRUE )
      cat( "PAPER 96/07, UNIVERSITY OF NEW ENGLAND, ARMIDALE, NSW, 2351, AUSTRALIA.\n",
         file = startUpFile, append = TRUE )
      cat( "\n",
         file = startUpFile, append = TRUE )
      cat( "INDIC VALUES:\n",
         file = startUpFile, append = TRUE )
      cat( "indic=2 says do not scale step length in unidimensional search\n",
         file = startUpFile, append = TRUE )
      cat( "indic=1 says scale (to length of last step) only if last step was smaller\n",
         file = startUpFile, append = TRUE )
      cat( "indic= any other number says scale (to length of last step) \n",
         file = startUpFile, append = TRUE )
   }

   returnList <- list( data = dataTable,
      crossSectionName = crossSectionName,
      timePeriodName = timePeriodName,
      yName = yName,
      xNames = xNames,
      qxNames = qxNames,
      zNames = zNames,
      quadHalf = quadHalf,
      functionType = functionType,
      logDepVar = logDepVar,
      mu = mu,
      eta = eta,
      insFile = insFile,
      dtaFile = dtaFile,
      outFile = outFile,
      startUpFile = startUpFile,
      iprint = iprint,
      indic = indic,
      tol = tol,
      tol2 = tol2,
      bignum = bignum,
      step1 = step1,
      igrid2 = igrid2,
      gridno = gridno,
      maxit = maxit,
      ite = ite,
      modelType = modelType,
      nCrossSection = nCrossSection,
      nTimePeriods = nTimePeriods,
      nTotalObs = nTotalObs,
      nXtotal = nXtotal,
      nZvars = nZvars )
   class( returnList ) <- "front41WriteInput"
   invisible( returnList )
}
