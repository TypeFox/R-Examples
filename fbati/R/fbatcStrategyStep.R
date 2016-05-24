rrange <- function( from=1, to=from, by=1 ) {
  if( from > to ) return( c() )
  return( seq( from=from, to=to, by=by ) )
}

## For strategies 3 and 4, returns marker names
univariateFilter <- function( ped, phe, markers, trait="trait", alpha=0.05, tempPrefix="temp_", sim=FALSE, FBATEXE="~/bin/fbat" ) {
  ## Either use the fbat package, or use FBAT... where is that code?
  res <- fbatShell( ped=ped, phe=phe, trait=trait, model="a", markers=markers, tempPrefix=tempPrefix, FBATEXE=FBATEXE )
  res$P <- as.numeric(res$P)
  markers <- intersect(  markers,  unique( res$Marker[res$P<alpha] )  )

  #print( res )
  #stop()

  ## If simulating, just return the markers
  if( sim )
    return( markers )

  ## Otherwise return the univariate results as well
  wh <- seq( from=1, to=nrow(res), by=2 )
  afreq <- as.numeric(res$afreq[wh])
  for( a in 1:length(afreq) ) afreq[a] <- min( afreq[a], 1-afreq[a] )
  univariateResults <- data.frame( marker=res$Marker[wh], afreq=afreq, numInf=as.numeric(res$"fam#")[wh], pvalue=as.numeric(res$P[wh]), stringsAsFactors=FALSE )
  return( list( filteredMarkers=markers, univariateResults=univariateResults ) )
}


correlationOrder <- function( m ) {
  #print( "nrow(m)")
  #print( nrow(m) )
  n <- nrow(m)
  if( n==1 ) return( 1 )  ## special case

  ## Zap lower triangle
  for( r in 1:n )
    m[r,1:r] <- NA

  ## Find the maximum
  w <- which( m == max(m,na.rm=TRUE) )[1]
  wc <- floor( (w-1)/n ) + 1
  wr <- w-(wc-1)*n
  #print( w )
  #print( wr )
  #print( wc )
  #print( m[wr,wc] )
  #stop()

  #print( "m" )
  #print( m )

  #print( "notchosen" )
  neworder <- c( wr, wc )
  #print( "neworder" )
  #print( neworder )
  #print( "n" )
  #print( n )
  while( length(neworder) != n ) {
    notchosen <- setdiff(1:n,neworder)
    #print( notchosen )

    ## Find maximum correlated with front
    maxFront <- NA
    maxFrontValue <- -999
    for( i in notchosen ) {
      newcorr <-  m[ min(neworder[1],i), max(neworder[1],i) ]
      if( newcorr > maxFrontValue ) {
        maxFrontValue <- newcorr
        maxFront <- i
      }
    }
    if( is.na(maxFront) ) stop("correlationOrder infinite loop!")
    neworder <- c( maxFront, neworder )
  }

  #print( "neworder final" )
  #print( neworder )

  return( neworder )
}

correlationOrder.debug <- function() {
  mat <- rbind( c(NA,0,0.1,0.5), c(0,NA,0.3,0.2), c(0.1,0.3,NA,0.15), c(0.5,0.2,0.15,NA) )
  mat-t(mat) ## should be zeros
  correlationOrder( mat )
}

fbatcStrategyStepUp <- function( ped, phe, markers=ped.markerNames(ped), trait="trait", traitType="auto", alphaMMarker=0.05, alphaStep=alphaMMarker, sortByCorrelation=TRUE, tempPrefix="temp_", sim=FALSE, debug=FALSE ) {
  fbat.install()
  sh.install()
  #FBAT <- fbat.exename()  ## 02/20/2014 codetools
  FBATEXE <- fbat.exename()

  ## First compute the multi-marker pvalue
  mmarkerPvalue <- fbatShellMM( ped=ped, phe=phe, markers=markers, trait=trait, FBATEXE=FBATEXE, tempPrefix=tempPrefix ) ########################

  if( mmarkerPvalue > alphaMMarker ) {
    ## No reason to go further
    if( sim )
      return( c( marker="", markerR="" ) )

    return( list( mmarkerPvalue=mmarkerPvalue ) )
  }

  ## Sort by the correlation?
  if( !sim & sortByCorrelation ) {
    ## Resort markers by correlation
    corr <- fbatShellCorrelation( ped=ped, markers=markers, FBATEXE=FBATEXE, tempPrefix=tempPrefix ) #####################
    neworder <- correlationOrder( corr )
    markers <- markers[ neworder ]
  }


  ### Compute the correlation of the markers
  #correlation <- fbatShellCorrelation( ped=ped, markers=markers, FBATEXE=FBATEXE, tempPrefix=tempPrefix )
  ## Compute the correlation of the markers (05.04.2009)
  correlation <- NULL
  if( !sim ) {
    try( {
      correlation <- fbatShellCorrelation( ped=ped, markers=markers, FBATEXE=FBATEXE, tempPrefix=tempPrefix ) #####################
    }, silent=TRUE )
  }

  ## Then compute the univariate results
  univariate <- univariateFilter( ped=ped, phe=phe, markers=markers, trait=trait, alpha=alphaStep, tempPrefix=tempPrefix, sim=FALSE, FBATEXE=FBATEXE )
  ## univariate$univariateResults is what we want
  univariate <- univariate$univariateResults

  if( debug )
    print( univariate ) ## DEBUG ONLY

  ## find the most significant marker (May also want to make sure enough informative families...)
  wh <- which( univariate$pvalue == min(univariate$pvalue,na.rm=TRUE) )
  if( length(wh) > 1 ) wh <- wh[1]

  markersChosen <- markersChosenR <- NULL
  step <- stepR <- list()
  if( length(wh)==1 && univariate$pvalue[wh] < alphaStep ) {
    markersChosen <- markersChosenR <- univariate$marker[wh]

    done <- doneR <- FALSE
    while( !done || !doneR ) {
      ## Are the robust and nonrobust the same, or do we need to do them differently?
      same <- length(markersChosen)==length(markersChosenR) && all( markersChosen==markersChosenR )

      ## Run the analysis on each of the markers
      markersA <- as.list(  setdiff( markers, markersChosen )  )
      markersC <- rep( list(markersChosen), length(markersA) )

      markersAR <- as.list(  setdiff( markers, markersChosenR )  )
      markersCR <- rep( list(markersChosenR), length(markersAR) )

      pvalue <- rep(1,length(markersA)); names(pvalue) <- markersA;
      pvalueR <- rep(1,length(markersAR)); names(pvalueR) <- markersAR;
      #numInf <- rep(0,length(markersA)); names(numInf) <- markersA;
      #numInfR <- rep(0,length(markersAR)); names(numInfR) <- markersAR;
      numInf <- rep("",length(markersA)); names(numInf) <- markersA;
      numInfR <- rep("",length(markersAR)); names(numInfR) <- markersAR;
      varExpl <- rep(NA,length(markersA)); names(varExpl) <- markersA; ## 01/03/09

      if( debug ) {
        print( "markersA" )
        print( markersA )
        print( "markersC" )
        print( markersC )
        print( "markersAR" )
        print( markersAR )
        print( "markersCR" )
        print( markersCR )
      }

      if( !done ) {
        for( m in 1:length(markersA) ) {
          try( {
                temp <- NULL
                temp <- fbatc( ped=ped, phe=phe, markerAnalyze=markersA[[m]], markerCondition=markersC[[m]], tempPrefix=tempPrefix, trait=trait )
                if( debug ) {
                  cat( "############ MODEL ############\n" )
                  cat( "markerAnalyze" ); print( markersA[[m]] );
                  cat( "markerCondition" ); print( markersC[[m]] );
                  print( temp )
                }
                pvalue[m] <- temp$pvalue
                ##numInf[m] <- as.numeric(as.character(temp$numInf))
                numInf[m] <- temp$numInf
                try( { varExpl[m] <- as.numeric(as.character(temp$varExpl)) }, silent=TRUE )  ## 01/03/09
                if( same && !doneR ) {
                  pvalueR[m] <- temp$pvalueR
                  ##numInfR[m] <- as.numeric(as.character(temp$numInfR))
                  numInfR[m] <- temp$numInfR
                }
              } )

        }
      }
      if( !doneR && !same ) {
        ## The robust needs to be done on different markers
        for( m in 1:length(markersAR) ) {
          try( {
                temp <- NULL
                temp <- fbatc( ped=ped, phe=phe, markerAnalyze=markersAR[[m]], markerCondition=markersCR[[m]], tempPrefix=tempPrefix, trait=trait )
                if( debug ) {
                  cat( "############ ROBUST ############\n" )
                  cat( "markerAnalyze" ); print( markersAR[[m]] );
                  cat( "markerCondition" ); print( markersCR[[m]] );
                  print( temp )
                }
                pvalueR[m] <- temp$pvalueR
                numInfR[m] <- temp$numInfR ##as.numeric(as.character(temp$numInfR))
              })
        }
      }

      ## Now see if we should add any other markers, or if we are done!
      if( !done ) {
        step[[length(step)+1]] <- list( pvalue=pvalue, numInf=numInf, markersAnalyze=markersA, markersCondition=markersC, varExpl=varExpl )
        if( any( pvalue < alphaStep, na.rm=TRUE ) ) {
          wh <- which( pvalue==min(pvalue) )[1]
          markersChosen <- c( markersChosen, markersA[[wh]] )
          if( length(pvalue) == 1 ) done <- TRUE  ## In case they all matter
        }else{
          done <- TRUE
        }
      }
      if( !doneR ){
        stepR[[length(stepR)+1]] <- list( pvalue=pvalueR, numInf=numInfR, markersAnalyze=markersAR, markersCondition=markersCR )
        if( any( pvalueR < alphaStep, na.rm=TRUE ) ) {
          wh <- which( pvalueR==min(pvalueR) )[1]
          markersChosenR <- c( markersChosenR, markersAR[[wh]] )
          if( length(pvalue) == 1 ) doneR <- TRUE  ## In case they all matter
        }else{
          doneR <- TRUE
        }
      }
    }
  }

  res <- list( mmarkerPvalue=mmarkerPvalue, correlation=correlation, univariate=univariate, step=step, markersChosen=markersChosen, stepR=stepR, markersChosenR=markersChosenR )
  class(res) <- c("fbatcSStep","list")
  return( res )
}

## file="" dumps it to screen (the default)
## preamble=TRUE means you can compile all the output
## build=TRUE means it will run pdflatex on the filename
fbatcStrategyStepLatex <- function( res, digits=4, ffile="", preamble=FALSE, build=preamble, pdf="" ) {
  fn <- function( num ) return( signif(num,digits) )

  ## Rogue input
  if( ( is.character(ffile) && ffile=="" ) || !is.character(ffile) ) {
    if( build ) build <- FALSE  ## Maybe print a message?
  }

  ## Allow us to dump to disk as well, stick this in every cat piece
  filename <- ffile  ## In case we wan't to compile the file
  print( filename )
  if( is.character(ffile) && nchar(ffile)>0 ) {
    ffile = file(ffile, 'w') ## Then open up a new file
    on.exit({ if(!is.null(ffile)) close(ffile) })
  }
  ccat <- function( ... )
    cat( ..., file=ffile )

  ## Output a header if necessary
  if( preamble ) {
    ccat( "\\documentclass[10pt,letterpaper]{article}\n")
    ccat( "\\usepackage{geometry}\n" )
    ccat( "\\geometry{verbose,letterpaper,tmargin=0.5in,bmargin=0.5in,lmargin=0.5in,rmargin=0.5in}\n" )
    ccat( "\\begin{document}\n\n" )
    if( is.character(filename) && filename!="" )
      ccat( "\\section*{", filename, "}\n\n", sep="" )
  }

  ## Print the multimarker p-value
  ccat( "Multimarker P-value = ", fn(res$mmarkerPvalue), "\n\n", sep="" )

  print( res )
  print( str(res) )

  #if( length(res) > 1 ) {
  if( !is.null(res$correlation) && !is.null(res$univariate) ) {
    ## format the numbers properly
    res$correlation <- fn( res$correlation )
    ##res$univariate[,2:3] <- fn( res$univariate[,2:3] )
    res$univariate[[2]] <- fn( res$univariate[[2]] )
    res$univariate[[4]] <- fn( res$univariate[[4]] )  ## 12/10/2008

    ## output the univariate results
    ccat( "\\begin{table}\n" )
    ccat( "  \\begin{tabular}{cccc}\n" )
    ccat( "    Marker & Allele Freq. & \\# Inf & pvalue \\\\\n" )
    ccat( "    \\hline\n" )
    for( r in 1:nrow(res$univariate) )
      ccat( "    ", res$univariate[r,1], " & ", res$univariate[r,2], " & ", res$univariate[r,3], " & ", res$univariate[r,4], "\\\\\n", sep="" )
    ccat( "  \\end{tabular}\n" )
    ccat( "  \\caption{Univariate results.}\n" )
    ccat( "\\end{table}\n\n\n" )

    ## output the correlation matrix
    n <- nrow(res$correlation)
    ccat( "\\begin{table}\n" )
    ccat( "  \\begin{tabular}{c|", paste(rep("c",n),sep=""), "}\n", sep="" )
    ccat( "    $m_{row}|m_{col}$ & ", paste(rownames(res$correlation),collapse=" & "), "\\\\\n" )
    ccat( "    \\hline\n" )
    for( r in 1:n ) {
      ccat( "    ", rownames(res$correlation)[r], " & ", sep="" )
      for( col in 1:n ) {
        if( col!=r ) ccat( res$correlation[r,col] ) else ccat(1) ##cat( "{\\it $\\rho^2$}" )
        if( col!=n ) ccat( " & ") else ccat( "\\\\\n")
      }
    }
    ccat( "  \\end{tabular}\n" )
    ccat( "  \\caption{Correlation matrix.}\n")
    ccat( "\\end{table}\n\n\n" )

    ## Helper function to output the results
    outstep <- function( step, pvalueHeader="pvalue", caption="" ) {
      if( length(step) == 0 ) return()

      useVarExpl <- FALSE
      ##if( !is.na(step[[1]]$varExpl[1]) )
      if( !is.null(step[[1]]$varExpl) && !is.na(step[[1]]$varExpl[1]) )
        useVarExpl <- TRUE

      ccat( "\\begin{table}\n" )
      if( is.null(step[[1]]$stepType) ) {
        if( !useVarExpl ) {
          ccat( "  \\begin{tabular}{cccc}\n" )
          ccat( "    Analyze & Condition & \\# Inf & ", pvalueHeader, "\\\\\n", sep="" )
        }else{
          ccat( "  \\begin{tabular}{ccccc}\n" )
          ccat( "    Analyze & Condition & \\# Inf & ", pvalueHeader, " & Var Expl\\\\\n", sep="" )
        }
      }else{
        if( !useVarExpl ) {
          ccat( "  \\begin{tabular}{ccccc}\n" )
          ccat( "    Step & Analyze & Condition & \\# Inf & ", pvalueHeader, "\\\\\n", sep="" )
        }else{
          ccat( "  \\begin{tabular}{cccccc}\n" )
          ccat( "    Step & Analyze & Condition & \\# Inf & ", pvalueHeader, " & Var Expl\\\\\n", sep="" )
        }
      }
      ccat( "    \\hline\n" )
      ccat( "    \\hline\n" )

      #print( str(step) )

      markerCondStrPrev <- ""
      markerAnalStrPrev <- ""
      for( s in 1:length(step) ) {
        for( i in 1:length(step[[s]]$pvalue) ) {
          if( !is.null(step[[s]]$stepType) ) {
            if( i==1 && ( s==1 || step[[s]]$stepType!=step[[s-1]]$stepType ) ) {
              if( s!= 1 )
                ccat( "    \\hline\n" )
              ccat( "    ", step[[s]]$stepType, " & ", sep="" )
            }else{
              ccat( "     & " )
            }
          }

          markerCondStr <- paste( step[[s]]$markersCondition[[i]], collapse="," )
          if( markerCondStrPrev=="" || markerCondStr!=markerCondStrPrev ) {
            markerCondStrPrev <- markerCondStr
          }else{
            markerCondStr <- ""
          }

          markerAnalStr <- paste( step[[s]]$markersAnalyze[[i]], collapse="," )
          if( markerAnalStrPrev=="" || markerAnalStr!=markerAnalStrPrev ) {
            markerAnalStrPrev <- markerAnalStr
          }else{
            markerAnalStr <- ""
          }

          if( !useVarExpl ) {
            ccat( "    ", markerAnalStr, " & ", markerCondStr, " & ", step[[s]]$numInf[i], " & ", step[[s]]$pvalue[i], "\\\\\n", sep="" )
          }else{
            ccat( "    ", markerAnalStr, " & ", markerCondStr, " & ", step[[s]]$numInf[i], " & ", step[[s]]$pvalue[i], " & ", step[[s]]$varExpl[i], "\\\\\n", sep="" )
          }
        }
        if( s!=length(step) ) ccat( "    \\hline\n")
      }

      ccat( "  \\end{tabular}\n" )
      ccat( "  \\caption{",caption,"}\n" )
      ccat( "\\end{table}\n\n" )
    }

    ## And output the results
    outstep( res$step, pvalueHeader="FBAT-C Model", caption="FBAT-C model-based stepwise result." );
    ccat( "Model-based markers chosen:" ); ccat( paste(res$markersChosen,collapse=", ") ); ccat( "\n\n" );

    outstep( res$stepR, pvalueHeader="FBAT-C Robust", caption="FBAT-C robust stepwise result." );
    ccat( "Model-free markers chosen:" ); ccat( paste(res$markersChosenR,collapse=", ") ); ccat( "\n\n" );
  }

  ## Close the file if necessary
  if( preamble ) {
    ccat( "\\end{document}\n")
    close(ffile)
    ffile <- NULL
  }

  ## And maybe build (only linux)
  if( build ) {
    system( paste( "pdflatex", filename ) )
    print( "filename" )
    print( filename )
    if( pdf != "" && is.character(filename) ) {
      pdfname <- paste( substring( filename, 1, nchar(filename)-4 ), ".pdf", sep="" )
      print( "pdfname" )
      print( pdfname )
      if( file.exists(pdfname) ) {
        system( paste( pdf, pdfname ) )
      }else{
        stop("Latex build failed.")
      }
    }
  }
}

## Written from the fbatcStrategyStepLatex
print.fbatcSStep <- function( x, ... ) {
  fn <- function( num ) return( signif(num,digits=4) )

  ## Output the correlation and the univariate results
  if( !is.null(x$correlation) ) {
    cat( "Correlation:\n" )
    print( fn(x$correlation) )
  }
  if( !is.null(x$univariate) ) {
    cat( "Univariate:\n" )
    x$univariate[[2]] <- fn( x$univariate[[2]] )
    x$univariate[[4]] <- fn( x$univariate[[4]] )
    print( x$univariate )
  }

  ## Helper function for each of the step-up/down
  outstep <- function( step, pvalueHeader="Pvalue" ) {
    Step <- c()
    Analyze <- c()
    Condition <- c()
    NumInf <- c()
    Pvalue <- c()

    if( !is.null(step) && length(step)>0 ) {
      for( s in 1:length(step) ) {
        for( i in 1:length(step[[s]]$pvalue) ) {
          if( !is.null(step[[s]]$stepType) ) {
            Step <- c( Step, step[[s]]$stepType )
          }else{
            Step <- c( Step, "" )
          }

          #Analyze <- c( Analyze, paste( step[[s]]$markersCondition[[i]], collapse="," ) )
          #Condition <- c( Condition, paste( step[[s]]$markersAnalyze[[i]], collapse="," ) )
          Analyze <- c( Analyze, paste( step[[s]]$markersAnalyze[[i]], collapse="," ) )
          Condition <- c( Condition, paste( step[[s]]$markersCondition[[i]], collapse="," ) )

          NumInf <- c( NumInf, step[[s]]$numInf[i] )
          Pvalue <- c( Pvalue, step[[s]]$pvalue[i] )
        }
      }
      df <- data.frame( Step=Step, Analyze=Analyze, Condition=Condition, NumInf=NumInf, Pvalue=Pvalue )
      names(df)[5] <- pvalueHeader

      print(df)
    }

    return(invisible())
  }

  ## And output the results
  cat( "FBAT-C model-based stepwise result:\n" )
  outstep( x$step, pvalueHeader="FBAT-C Model" )
  cat( "Model-based markers chosen:\n" )
  cat( paste(x$markersChosen,collapse=","), "\n\n" )

  cat( "FBAT-C robust stepwise result:\n" )
  outstep( x$stepR, pvalueHeader="FBAT-C Robust" )
  cat( "Model-free markers chosen:\n" )
  cat( paste(x$markersChosenR,collapse=","), "\n\n" )

  return(invisible())
}

fbatcStrategyStepDown <- function( ped, phe, markers=ped.markerNames(ped), markersChosen=ped.markerNames(ped), markersChosenR=markersChosen, trait="trait", traitType="auto", alphaMMarker=0.05, alphaStep=alphaMMarker, sortByCorrelation=TRUE, tempPrefix="temp_", sim=FALSE, debug=FALSE ) {
  fbat.install()
  sh.install()
  FBAT <- fbat.exename()
  FBATEXE <- FBAT

  ## First compute the multi-marker pvalue
  mmarkerPvalue <- fbatShellMM( ped=ped, phe=phe, markers=markers, trait=trait, FBATEXE=FBATEXE, tempPrefix=tempPrefix )

  if( mmarkerPvalue > alphaMMarker ) {
    ## No reason to go further
    if( sim )
      return( list( marker="", markerR="" ) )

    return( list( mmarkerPvalue=mmarkerPvalue ) )
  }

  ## Sort by the correlation?
  if( !sim & sortByCorrelation ) {
    ## Resort markers by correlation
    corr <- fbatShellCorrelation( ped=ped, markers=markers, FBATEXE=FBATEXE, tempPrefix=tempPrefix )
    neworder <- correlationOrder( corr )
    markers <- markers[ neworder ]
  }


  ### Compute the correlation of the markers
  #correlation <- fbatShellCorrelation( ped=ped, markers=markers, FBATEXE=FBATEXE, tempPrefix=tempPrefix )
  ## Compute the correlation of the markers (05.04.2009)
  correlation <- NULL
  if( !sim ) {
    try( {
      correlation <- fbatShellCorrelation( ped=ped, markers=markers, FBATEXE=FBATEXE, tempPrefix=tempPrefix ) #####################
    }, silent=TRUE )
  }

  ## Then compute the univariate results
  univariate <- univariateFilter( ped=ped, phe=phe, markers=markers, trait=trait, alpha=alphaStep, tempPrefix=tempPrefix, sim=FALSE, FBATEXE=FBATEXE )
  ## univariate$univariateResults is what we want
  univariate <- univariate$univariateResults

  if( debug )
    print( univariate ) ## DEBUG ONLY

  ## If there aren't any univariate significant results, then do we bother?
  #if( sum(univariate$pvalue<alphaStep) == 0 )
  #  return( list( mmarkerPvalue=mmarkerPvalue, univariate=univariate, correlation=correlation ) )

  step <- stepR <- list()

  done <- length(markersChosen)==0
  doneR <- length(markersChosenR)==0
  while( !done || !doneR ) {
   #print( "DOWN NOT DONE!" )
    ## Special cases if only one marker -- then just do the univariate test
    if( length(markersChosen) == 1 ) {
      done <- TRUE ## one way or the other, will be done
      unival <- univariate$pvalue[univariate$marker==markersChosen]
      ## Then we should only do the univariate test
      if( unival > alphaStep ) {
        markersChosen <- c()
        next
      }
    }
    if( length(markersChosenR) == 1 ) {
      doneR <- TRUE
      unival <- univariate$pvalue[univariate$marker==markersChosenR]
      if( unival > alphaStep ) {
        markersChosenR <- c()
        next
      }
    }

    ## Are the robust and nonrobust the same, or do we need to do them differently?
    same <- length(markersChosen)==length(markersChosenR) && all( markersChosen==markersChosenR )

    ## Run the analysis on each of the markers
    markersA <- as.list(markersChosen)
    markersC <- c()
    for( m in rrange(1,length(markersA)) )
      markersC[[m]] <- setdiff( markersChosen, markersA[[m]] )

    markersAR <- markersChosenR
    markersCR <- c()
    for( m in rrange(1,length(markersAR)) )
      markersCR[[m]] <- setdiff( markersChosenR, markersAR[[m]] )

    pvalue <- rep(1,length(markersA)); names(pvalue) <- markersA;
    pvalueR <- rep(1,length(markersAR)); names(pvalueR) <- markersAR;
    numInf <- rep("",length(markersA)); names(numInf) <- markersA;
    numInfR <- rep("",length(markersAR)); names(numInfR) <- markersAR;
    varExpl <- rep(NA,length(markersA)); names(varExpl) <- markersA;

    if( debug ) {
      print( "markersA" )
      print( markersA )
      print( "markersC" )
      print( markersC )
      print( "markersAR" )
      print( markersAR )
      print( "markersCR" )
      print( markersCR )
    }

    if( !done ) {
      for( m in 1:length(markersA) ) {
        try( {
              temp <- NULL
              temp <- fbatc( ped=ped, phe=phe, markerAnalyze=markersA[[m]], markerCondition=markersC[[m]], tempPrefix=tempPrefix, trait=trait )
              if( debug ) {
                cat( "############ MODEL ############\n" )
                cat( "markerAnalyze" ); print( markersA[[m]] );
                cat( "markerCondition" ); print( markersC[[m]] );
                print( temp )
              }
              pvalue[m] <- temp$pvalue
              numInf[m] <- temp$numInf #as.numeric(as.character(temp$numInf))
              try( { varExpl[m] <- as.numeric(as.character(temp$varExpl)) }, silent=TRUE )
              if( same && !doneR ) {
                pvalueR[m] <- temp$pvalueR
                numInfR[m] <- temp$numInfR #as.numeric(as.character(temp$numInfR))
              }
            } )

      }
    }
    if( !doneR && !same ) {
      ## The robust needs to be done on different markers
      for( m in 1:length(markersAR) ) {
        try( {
              temp <- NULL
              temp <- fbatc( ped=ped, phe=phe, markerAnalyze=markersAR[[m]], markerCondition=markersCR[[m]], tempPrefix=tempPrefix, trait=trait )
              if( debug ) {
                cat( "############ ROBUST ############\n" )
                cat( "markerAnalyze" ); print( markersAR[[m]] );
                cat( "markerCondition" ); print( markersCR[[m]] );
                print( temp )
              }
              pvalueR[m] <- temp$pvalueR
              numInfR[m] <- temp$numInfR #as.numeric(as.character(temp$numInfR))
            })
      }
    }

    ## Now see if we should add any other markers, or if we are done!
    if( !done ) {
      step[[length(step)+1]] <- list( pvalue=pvalue, numInf=numInf, markersAnalyze=markersA, markersCondition=markersC, varExpl=varExpl )
      if( !all( pvalue > alphaStep, na.rm=TRUE ) && any( pvalue > alphaStep, na.rm=TRUE ) ) { ## 05.01.2009
        wh <- which( pvalue==max(pvalue) )[1]
        markersChosen <- setdiff( markersChosen, markersA[[wh]] )
        ## No update to done given special case at the top of the loop for univariate
      }else{
        done <- TRUE
      }
    }
    if( !doneR ){
      stepR[[length(stepR)+1]] <- list( pvalue=pvalueR, numInf=numInfR, markersAnalyze=markersAR, markersCondition=markersCR )
      if( any( pvalueR > alphaStep, na.rm=TRUE ) ) {
        wh <- which( pvalueR==max(pvalueR) )[1]
        markersChosenR <- setdiff( markersChosenR, markersAR[[wh]] )
        ## No update to doneR given special case at the top of the loop for univariate
      }else{
        doneR <- TRUE
      }
    }
  }

  res <- list( mmarkerPvalue=mmarkerPvalue, correlation=correlation, univariate=univariate, step=step, markersChosen=markersChosen, stepR=stepR, markersChosenR=markersChosenR )
  class(res) <- c("fbatcSStep","list")
  return( res )
}

fbatcStrategyStep <- function( ped=NULL, phe=NULL, markers=ped.markerNames(ped), trait="trait", traitType="auto", alphaMMarker=0.05, alphaStep=alphaMMarker, sortByCorrelation=TRUE, tempPrefix="temp_", sim=FALSE, debug=FALSE ) {
  if( is.null(ped) )
    return( fbatcSSFuncGUI() )

  ## First do step-up routine
  if( debug ) cat( "*** fbatcStrategyStepUp ***\n")
  up <- fbatcStrategyStepUp( ped=ped, phe=phe, markers=markers, trait=trait, traitType=traitType, alphaMMarker=alphaMMarker, alphaStep=alphaStep, sortByCorrelation=sortByCorrelation, tempPrefix=tempPrefix, sim=sim, debug=debug )
  ## Then step-down, based on the step-up results
  if( debug ) cat( "*** fbatcStrategyStepDown ***\n")
  #print( up$markersChosen )
  #print( up$markersChosenR )
  #stop()
  #up$markersChosen <- markers[1:2]
  #up$markersChosenR <- markers[1:2]
  down <- fbatcStrategyStepDown( ped=ped, phe=phe, markersChosen=up$markersChosen, markersChosenR=up$markersChosenR, trait=trait, traitType=traitType, alphaMMarker=alphaMMarker, alphaStep=alphaStep, sortByCorrelation=sortByCorrelation, tempPrefix=tempPrefix, sim=sim, debug=debug )
  if( debug ) cat( "fbatcStrategy -- putting results together...\n")

  ## Now put those results together...
  step <- list()
  ########print( "up" )
  ########print( up )
  ########print( "step" ); print( step );
  for( i in rrange(1,length(up$step)) ) {
    step[[length(step)+1]] <- up$step[[i]]
    step[[length(step)]]$stepType <- "up"
  }
  for( i in rrange(1,length(down$step)) ) {
    step[[length(step)+1]] <- down$step[[i]]
    step[[length(step)]]$stepType <- "down"
  }

  stepR <- list()
  for( i in rrange(1,length(up$stepR)) ) {
    stepR[[length(stepR)+1]] <- up$stepR[[i]]
    stepR[[length(stepR)]]$stepType <- "up"
  }
  for( i in rrange(1,length(down$stepR)) ) {
    stepR[[length(stepR)+1]] <- down$stepR[[i]]
    stepR[[length(stepR)]]$stepType <- "down"
  }

  res <- list( mmarkerPvalue=up$mmarkerPvalue, correlation=up$correlation, univariate=up$univariate, step=step, markersChosen=down$markersChosen, stepR=stepR, markersChosenR=down$markersChosenR )
  class(res) <- c("fbatcSStep","list")
  return( res )
}

#######################
## GUI Functionality ##
fbatcSSFunc <- function( ped, phe, trait, traitType="auto", alphaMMarker=0.05, alphaStep=0.05, sortByCorrelation=TRUE, tempPrefix="temp_", print_results, latex_results ) {
  if( is.null(ped) || is.na(ped) || (is.character(ped) & nchar(ped)==0) )
    return( "A 'pedigree' file must be specified." )

  if( is.na(trait) || is.null(trait) )
    return( "A 'trait' must be specified." )

  fped <- fread.ped( ped, lowercase=FALSE )
  fphe <- NULL
  if( is.character(phe) & nchar(phe)!=0 )
    try( fphe <- fread.phe( phe, lowercase=FALSE ), silent=TRUE )

  res <- fbatcStrategyStep( ped=fped, phe=fphe, trait=trait, traitType=traitType, alphaMMarker=alphaMMarker, alphaStep=alphaStep, sortByCorrelation=sortByCorrelation, tempPrefix=tempPrefix )

  guiSet( "fbatcSStepRes", res )

  return( "Processed." )
}

updateFbatcSSFunc <- function( arg ) {
  ## Don't care about the pedigree -- does all markers in the pedigree, period.
  if( arg=="ped" ) {
    ## No point in setting the markers, like we would usually do...
    pedfile <- guiGetValue(1)
    file.strip.extension <- getFromNamespace( "file.strip.extension", "pbatR" )
    phename <- paste( file.strip.extension(pedfile), ".phe", sep="" )
    guiSetValue( 2, phename )
    arg <- "phe"
  }

  if( arg=="phe" ) {
    ##phefile <- guiGetValue(2)  ## 02/20/2014 codetools...
    phe <- NULL
    try( phe <- read.phe(guiGetValue(2)), silent=TRUE )
    if( !is.null(phe) ) {
      possibleEnv <- names(phe)[-c(1,2)]
      setListElements( "trait", c("AffectionStatus",possibleEnv) )
    }else{
      setListElements( "trait", c("AffectionStatus") )
    }
  }
}

printFbatcSSFunc <- function() {
  res <- guiGetSafe( "fbatcSStepRes" )

  if( !is.null(res) && !is.na(res) ) {
    print( res )
  }else{
    tkmessageBox( message="There are no results to print.", title="No Results" )
  }
}

latexFbatcSSFunc <- function() {
  res <- guiGetSafe( "fbatcSStepRes" )

  if( !is.null(res) && !is.na(res) ) {
    ## Get a filename to write the results to...
    defaultFile <- "results.tex"
    outStr <- tclvalue(tkgetSaveFile(title="Latex FBAT-C Results", filetypes="{{TEX (latex)} {.tex}}", initialfile=defaultFile))
    if( nchar(outStr) > 0 ) {
      fbatcStrategyStepLatex( res, ffile=outStr, preamble=TRUE, build=!isWindows() )
      cat( paste( "Check '", outStr, "' for the resulting LaTeX file that can be compiled (or may be done so automatically for you).\n", sep="" ) )
    }else{
      tkmessageBox( "Could not write file out to disk.", title="Write Failure" )
    }
  }else{
    tkmessageBox( message="There are no results to write.", title="No Results" )
  }
}

fbatcSSFuncGUI <- function() {
  ##require( fgui )

  gui( fbatcSSFunc,
      argFilename=list(ped=NULL,phe=NULL),
      argFilter=list(ped="{{Pedigree file} {.ped}}", phe="{{Phenotype file} {.phe}}"),
      argOption=list(traitType=c("auto","binary","continuous"),
        sortByCorrelation=c("TRUE","FALSE")),
      argList=list(trait=c("AffectionStatus")),
      argCommand=list(print_results=printFbatcSSFunc, latex_results=latexFbatcSSFunc),
      callback=updateFbatcSSFunc,
      helpsFunc=fbatcStrategyStep,
      title="FBAT-C Stepwise GUI",
      argText=list( ped="Pedigree File...", phe="Phenotype file...", tempPrefix="Temp Prefix", alphaMMarker="Multi-marker Alpha", alphaStep="Stepwise Alpha", trait="Trait", traitType="Trait Type", print_results="Print Results to R Console", latex_results="Print Results to LaTeX File")
      )
  return( guiGetSafe("fbatcSStepRes") )
}
