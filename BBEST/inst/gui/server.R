##########################################################################################
#                                                                                        #
#                                SHINY SERVER                                            #
#                                                                                        #
##########################################################################################
options(shiny.maxRequestSize=30*1024^2)
library(shiny)
##  MODULE STRUCTURE: (OLD...)
##
##  LOAD DATA
##    READ DATA              [[depends on: vals$nB; changes: vals$dat]]
##  SET ADDITIONAL PARAMETERS 
##    TRUNCATE DATA
##    SET LAMBDA
##    SET BASELINE
##    RETURN BASELINE
##    SET SIGMA
##    RENDERING:
##      TRUNCATE
##      BKG BOUNDS
##  SET R-SPACE LIKELIHOOD
##    SET GR
##    RENDERING:
##      SORRY GR (FOR BANKS)  [[depends on: vals$nB]]
##  SET PARAMETERS FOR DIFEV AND DO FIT
##    READ INPUTS AND DO FIT  [[depends on: vals$dat, vals$nB, vals$datGr; changes: vals$fitRes]]
##  FIT RESULTS
##    DOWNLOAD
##      DOWNLOAD RDATA BUTTON 
##      DOWNLOAD RDATA Handler
##      DOWNLOAD TEXT BUTTON 
##      DOWNLOAD TEXT Handler
##      DOWNLOAD FIX BUTTON 
##      DOWNLOAD FIX Handler
##      DOWNLOAD GR BUTTON 
##      DOWNLOAD GR Handler
##    CALCULATE
##      INPUTS FOR GR
##        HEADER
##        MIN(R)
##        MAX(R)
##        DR
##        CALC GR BUTTON
##      CALCULATE GR HANDLER
##    DO ITERATION
##      HEADER
##      EPS
##      N.ITER
##      DO BUTTON

##  RENDER OUTPUT
##    OUTPUT TABLE
##    OUTPUT DATA PLOT
##    DOWNLOAD DATA HANDLER
##    SHOWS PROGRESS
##    RENDER FIT RESULTS
##      FIT RESULTS PLOT SQ 
##      FIT RESULTS PLOT GR 

  
shinyServer(function(input, output, session) {
 
## initialization:
##   dat = main data variable
##   nB = number of banks 
##   Gr = calculated PDF 
##   datGr = data plus Gr likelihood info
##   fitRes = results of the fit
##   fitIter = results of the fit after iteration 
##   fitResFinal = fitIter, if exists; fitRes, if not.
##                 helps to leave fitRes untouched               

  vals <- reactiveValues(dat=list(list()), XInit=list(), nB=1, 
                         Gr=list(), estGr=list(), datGr=list(list()), 
                         fitRes=list(list()), fitResIter=list(list()),
                         fitResFinal=list(list()),
                         xlim=NA, ylim=NA, yRescale=c(0,1))
                         
                         
##########################################################################################
#                                                                                        #
#                                LOAD DATA                                               #
#                                                                                        #
##########################################################################################   
                      
########################   
##  ==  READ DATA  ==
  observe({
    inFile <- input$datafile
    isolate({
      if (is.null(inFile))
        return(NULL)
##    don't ask me why...      
      write(inFile$name, file="01x000.tmp")
      ext <- scan(file="01x000.tmp", what="list", sep='\n')  # get extension
      ext <- tail(strsplit(inFile$name, '[.]')[[1]], 1)
      file.remove("01x000.tmp")
      vals$nB <- 1
      if(ext=="sqa"){
        vals$dat  <- read.sqa(file=inFile$datapath)
        vals$nB <- length(vals$dat)
      }
      else if(ext=="sqb" ||  ext=="sq"){
        vals$dat[[1]] <- read.sqb(file=inFile$datapath)
      }
      else if(ext=="csv" || ext=="txt"){  # another don't ask me why...
        dat.tmp <- read.csv(inFile$datapath, header=input$headerCB, sep=input$separatorRB)  
        vals$dat[[1]] <- dat.tmp
      }
      else if(ext=="RData"){
        L <- sapply(inFile$datapath, function(x) mget(load(x)), simplify = FALSE)
        L <- L[[1]]$fit.results
        if(!is.null(L)){
          N <- length(L)
          vals$nB <- N
          for(i in 1:N){
            vals$fitRes[[i]] <- L[[i]]
            dat <- list(x=L[[i]]$x, y=L[[i]]$curves$y, SB=L[[i]]$curves$SB, 
                        sigma=L[[i]]$fit.details$sigma, lambda=L[[i]]$fit.details$lambda)
            vals$dat[[i]] <- dat      
          }
          if(N==1) vals$datGr[[1]] <- L[[1]]$fit.details$Gr
        }
      }
      vals$xlim <- vals$ylim <- matrix(NA, nrow=vals$nB, ncol=2)
      for(i in 1:vals$nB) vals$XInit[[i]] <- vals$dat[[i]]$x
      vals$yRescale<- c(0,1)
    })
  })
  
##########################################################################################
#                                                                                        #
#                         SET ADDITIONAL PARAMETERS                                      #
#                                                                                        #
##########################################################################################      

###########################
##  ==  TRUNCATE DATA  ==
  observe({
    input$truncLimits
    Sys.sleep(1)
    isolate({
      trunc <- input$truncLimits
      if(is.null(trunc))
        return(NULL)
      if (trunc != ""){
        tr <- as.numeric(unlist(strsplit(trunc, ","))) 
        if( (length(tr)==2) && !any(is.na(tr))){
          if((tr[1] == min(vals$dat[[1]]$x)) && (tr[2] == max(vals$dat[[1]]$x)))
            return(NULL)
          inFile <- input$datafile
          if (is.null(inFile))
            return(NULL)
          write(inFile$name, file="01x001.tmp")
          ext <- scan(file="01x001.tmp", what="list", sep='\n')  # get extension
          ext <- tail(strsplit(inFile$name, '[.]')[[1]], 1)
          file.remove("01x001.tmp")       
          vals$dat <- list(list())
          vals$nB <- 1
          if(ext=="sqa"){
            vals$dat <- read.sqa(file=inFile$datapath) 
          }
          else if(ext=="sqb" || ext=="sq"){
            vals$dat[[1]] <- read.sqb(file=inFile$datapath)
          }
          else if(ext=="csv" || ext=="txt"){  # another don't ask me why...
            dat.tmp <- read.csv(inFile$datapath, header=input$headerCB, sep=input$separatorRB)  
            vals$dat[[1]] <- dat.tmp
          }
          else if(ext=="RData"){
            L <- sapply(inFile$datapath, function(x) mget(load(x)), simplify = FALSE)
            L <- L[[1]]
            L <- L$fit.results
            if(!is.null(L)){
              if(is.null(L$x) && !is.null(L[[1]]$x)){ #number of banks
                N <- length(L)
                vals$nB <- N
                for(i in 1:N){
                  vals$fitRes[[i]] <- L[[i]]
                  dat <- list(x=L[[i]]$x, y=L[[i]]$curves$y, SB=L[[i]]$curves$SB, 
                              sigma=L[[i]]$fit.details$sigma, lambda=L[[i]]$fit.details$lambda)
                  vals$dat[[i]] <- dat      
                }
              }
              else{  #single function
                vals$fitRes[[1]] <- L
                dat <- list(x=L$x, y=L$curves$y, SB=L$curves$SB, sigma=L$fit.details$sigma, lambda=L$fit.details$lambda)
                datGr <- L$fit.details$Gr
                vals$dat[[1]] <- dat
                vals$datGr[[1]] <- datGr
                vals$nB <- 1         
              }
            }
          }
          vals$dat[[1]] <- trim.data(vals$dat[[1]], tr[1], tr[2])    
          lambda <- input$lambda
          if (lambda != ""){
            lam <- as.numeric(unlist(strsplit(lambda, ",")))
          if(length(lam)==5)
            vals$dat[[1]] <-  set.lambda(vals$dat[[1]], lambda=NA, lambda_1=lam[2], lambda_2=lam[4], 
                            lambda_0=lam[5], x_1=lam[1], x_2=lam[3]) 
          }
          if(input$setSB){
            n.atoms <- as.numeric(unlist(strsplit(input$SBNAtoms, ",")))
            f <- as.numeric(unlist(strsplit(input$SBScLen, ",")))
            oneADP <- input$oneADP
            if(!input$fitADP)
              ADP <- as.numeric(unlist(strsplit(input$ADP, ",")))
            else
              ADP <- NA
            if( (length(n.atoms)==length(f)) && ( (length(f)==length(ADP)) || ( (oneADP==TRUE) && (length(ADP)==1) ) || (input$fitADP==TRUE) )  && (length(f)>0) )
              vals$dat[[1]] <- set.SB(vals$dat[[1]], SB=NA, n.atoms=n.atoms, scatter.length=f, ADP=ADP, fit=input$fitADP, oneADP=oneADP)    
          }
          vals$dat[[1]]$sigma <- NULL  
          vals$xlim <- vals$ylim <- matrix(NA, nrow=vals$nB, ncol=2)
        }
      }
    })
  })

##########################  
##  ==  SET LAMBDA  ==
  observe({
    input$lambda
    isolate({    ## react on change
      if(is.null(input$lambda))                   
        return(NULL)
      lambda <- input$lambda
      if (lambda != ""){
        lam <- as.numeric(unlist(strsplit(lambda, ",")))
        if((length(lam)==5) && !any(is.na(lam))){
          for(i in 1:vals$nB){ 
            vals$dat[[i]] <- set.lambda(vals$dat[[i]], lambda=NA, lambda_1=lam[2], lambda_2=lam[4], 
                            lambda_0=lam[5], x_1=lam[1], x_2=lam[3])
          }
        }
      }
    }) 
  })

##########################  
##  ==  SET BASELINE  ==
  observe({
    input$SBNAtoms
    input$SBScLen
    input$ADP
    input$oneADP
    input$fitADP
    if(input$setSB){
      isolate({
        n.atoms <- as.numeric(unlist(strsplit(input$SBNAtoms, ",")))
        f <- as.numeric(unlist(strsplit(input$SBScLen, ",")))
        oneADP <- input$oneADP
        if(!input$fitADP)
          ADP <- as.numeric(unlist(strsplit(input$ADP, ",")))
        else
          ADP <- NA
        if( (length(n.atoms)==length(f)) &&                                       # numbers of atoms and sc lengths are ready
            ( (length(f)==length(ADP)) || ((oneADP==TRUE) && (length(ADP)==1)) || # ADP factor(s) is(are) ready
                (input$fitADP==TRUE) )  &&                                        # smth was indicated
            (length(f)>0) 
          ){                      
          for(i in 1:vals$nB)  
            vals$dat[[i]] <- set.SB(vals$dat[[i]], SB=NA, n.atoms=n.atoms, scatter.length=f, ADP=ADP, fit=input$fitADP, oneADP=oneADP)    
        }
        else{                                                                     # smth was indicated
          for(i in 1:vals$nB)  
            vals$dat[[i]]$SB <- rep(0, length(vals$dat[[i]]$x))           
        }
      })
    }
  })


############################  
##  ==  RETURN BASELINE  ==
## restores baseline to the value specified in datafile if 'set/recalculate baseline' was cancelled
  observe({
    input$setSB
    isolate({
      if(!input$setSB){
        inFile <- input$datafile
        if (is.null(inFile))
          return(NULL)
  ##    don't ask me why...      
        write(inFile$name, file="01x002.tmp")
        ext <- scan(file="01x002.tmp", what="list", sep='\n')  # get extension
        ext <- tail(strsplit(inFile$name, '[.]')[[1]], 1)
        file.remove("01x002.tmp")
        dat <- list(list())
        if(ext=="sqa")
          dat  <- read.sqa(file=inFile$datapath)
        else if(ext=="sqb" || ext=="sq" )
          dat[[1]] <- read.sqb(file=inFile$datapath)
        else{  # another don't ask me why...
          dat.tmp <- read.csv(inFile$datapath, header=input$headerCB, sep=input$separatorRB)  
          dat[[1]] <- dat.tmp
        }
        wis <-  whatIsSpecified(dat)
        if(wis[[1]]$SB==TRUE){
          tr <- as.numeric(unlist(strsplit(input$truncLimits, ","))) 
          if( !(length(tr)==2) || any(is.na(tr))){
            tr <- 0
            tr[1] = min(vals$dat[[1]]$x) 
            tr[2] = max(vals$dat[[1]]$x)
          }    
          for(i in 1:vals$nB){
            if(vals$nB==1) dat[[1]] <- trim.data(dat[[1]], tr[1], tr[2])       
            vals$dat[[i]]$SB <- dat[[i]]$SB         
          }
        }
        vals$xlim <- vals$ylim <- matrix(NA, nrow=vals$nB, ncol=2)
      }
    })   
  })

  
############################
##  ==  SET SIGMA  ==
  observe({
    input$calcSigmaButton
    isolate({
      sigPar <- as.numeric(unlist(strsplit(input$sigma, ",")))
      k <- as.numeric(unlist(strsplit(input$sigmaTS, ",")))
      progress <- Progress$new(session)
      mess <- "Calculating, please wait..."
      progress$set(message = mess, value = 0.1)
      if( length(sigPar)==1 && !is.na(sigPar) && !any(is.na(k)) ){
        for(i in 1:vals$nB){  
          vals$dat[[i]] <- set.sigma(vals$dat[[i]], n.regions=sigPar, thresh.scale=k)
          progress$set(message = mess, value = (i/vals$nB-0.01))
        }
      }
      if( length(sigPar)==2 && !any(is.na(sigPar))  && !is.na(k) ){
        for(i in 1:vals$nB){
          vals$dat[[i]] <- set.sigma(vals$dat[[i]], x.bkg.only=sigPar, thresh.scale=k)
          progress$set(message = mess, value = (i/vals$nB-0.01))
        }
      }
      progress$set(message = 'Calculating, please wait...', value = 0.999)
      progress$close()
    })
  })

    
    
############################
##  ==  SET R-SIGMA AND PLOT G(R) ==
  observe({
    input$plotPrelimGr
    isolate({
      gridparam <- as.numeric(unlist(strsplit(input$rGrid, ",")))

      progress <- Progress$new(session)
      mess <- "Calculating, please wait... \n\n"
      wis <- whatIsSpecified(vals$dat)
      if(!is.null(input$bankNo))
        bankNo <- as.numeric(input$bankNo)
      else
        bankNo <- 1
      if(wis[[bankNo]]$sigma && (length(gridparam)==3)){
        progress$set(message = mess, value = 0.1)
        minR =  gridparam[1]  
        maxR =  gridparam[2]  
        dr =  gridparam[3]  
        r <- seq(minR, maxR, dr)
        sigma.r <- 0
        delta <- c(diff(vals$dat[[bankNo]]$x)[1], diff(vals$dat[[bankNo]]$x))
        cat("Calculating r-space noise... \n\n")
        progress$set(message = mess, value = 0.25)
        for(j in 1:length(r)){
          sigma.r[j] <- sum((2/pi*delta*vals$dat[[bankNo]]$x*sin(vals$dat[[bankNo]]$x*r[j])*vals$dat[[bankNo]]$sigma)^2)
          sigma.r[j] <- sqrt(sigma.r[j])
        }
        # avoid dividing by zero  
        if(sigma.r[1]==0)
          sigma.r[1] <- sigma.r[2]
      
        progress$set(message = mess, value = 0.75)
        cat("Calculating FT of the experimental data... \n\n")
        gr <- sineFT(f.Q=vals$dat[[bankNo]]$y-1, Q=vals$dat[[bankNo]]$x, r=r)
        vals$estGr <- list(r=r, gr=gr, stdev=sigma.r)
        progress$set(message = mess, value = 0.999)
      }
      else if(!wis[[bankNo]]$sigma && wis[[bankNo]]$x){
        progress$set(message = "Estimate Q-space noise first!", value = 0.0)
        Sys.sleep(2)
      }      
      else if(length(gridparam)!=3 && wis[[bankNo]]$x){
        progress$set(message = "Set r-space grid!", value = 0.0)
        Sys.sleep(2)
      }         
      
      progress$close()
    })
  })

###################################
##                               ##
##        RENDERING OUTPUT       ##
##                               ##
###################################

####################################
## OUTPUT TRUNCATE DATA 
  output$truncLimitsR <- renderUI({
    if (vals$nB!=1)
      return(helpText("not available for banks..."))
      
    if (is.null(vals$dat[[1]]) || (length(vals$dat[[1]]$x)==0))
      return(textInput("truncLimits", label = c("Type minimum x, maximum x"), value =""))

    truncLim <- toString(c(min(vals$dat[[1]]$x), max(vals$dat[[1]]$x)))
    textInput("truncLimits", label = c("Type minimum x, maximum x"), value = truncLim)
  })
  
####################################
## OUTPUT SQA SPLIT DATA   
   output$sqaSplit <- renderUI({
    if (vals$nB==1)
      return(NULL)
    downloadButton('downloadSqaSplit', 'Split by banks and download')  
  }) 
  
  output$downloadSqaSplit <- downloadHandler(
    filename = function() { paste('banks', '.zip', sep='') }, 
      content = function(file) {
        inFile <- input$datafile
        sqa <- scan(file=inFile$datapath, what="list", sep="\n")
        N <- length(sqa)
        i.start <- 0
        nBanks <- 0
        for(i in 1:N){
          if(strsplit(sqa[i], split=" ")[[1]][1]=="#L"){
            i.start[nBanks+1] <- i+1 
            nBanks <- nBanks + 1 
          }
        }
        i.start[nBanks+1] <- length(sqa)+5
        
        name <- 0
        for(i in 1:nBanks){
          name[i] <- strsplit(inFile$name, '[.]')[[1]]
          name[i] <- paste(name[i], "_b", i, ".sqa", sep="") 
          writeLines(sqa[ (i.start[i]-4):(i.start[i+1]-5)], con = name[i], sep = "\n", useBytes = FALSE)  
        }
        zip(zipfile=file, files=name)
        if(file.exists(paste0(file, ".zip"))) {file.rename(paste0(file, ".zip"), file)}
  
      }
  )
  

       
       
####################################
## OUTPUT BKG BOUNDS 
  output$bkgBoundsR <- renderUI({   
    if (is.null(vals$dat[[1]]) || (length(vals$dat[[1]]$y)==0))
      return(textInput("bkgBounds", label = strong("Lower and upper bounds for background")))
    
    bkgBndsArray <- matrix(0, nrow=vals$nB, ncol=2)
    sbBndsArray <- matrix(0, nrow=vals$nB, ncol=2)
    for(i in 1:vals$nB)
      bkgBndsArray[i,] <- c(min(vals$dat[[i]]$y), max(vals$dat[[i]]$y))

    bkgBnds <- c(min(bkgBndsArray[,1]), max(bkgBndsArray[,2]))
    
    isSBAvail <- whatIsSpecified(vals$dat)[[1]]$SB
    if(isSBAvail==TRUE){
      for(i in 1:vals$nB)
        sbBndsArray[i,] <- c(min(vals$dat[[i]]$SB), max(vals$dat[[i]]$SB))
    }  
    sbBnds <- c(min(sbBndsArray[,1]), max(sbBndsArray[,2]))
    
    bkgBnds[1] <- signif(bkgBnds[1] - sbBnds[2] - 0.2*abs(bkgBnds[1]) - 0.2*abs(sbBnds[2]), 3)
    bkgBnds[2] <- signif(bkgBnds[2] - sbBnds[1] + 0.2*abs(bkgBnds[2]) + 0.2*abs(sbBnds[1]), 3)
    
    return(textInput("bkgBounds", label =  strong("Lower and upper bounds for background"), value = toString(bkgBnds)))
  })    

##########################################################################################
#                                                                                        #
#                            SET R-SPACE LIKELIHOOD                                      #
#                                                                                        #
##########################################################################################  

######################## 
##  ==  SET Gr  ==
  observe({
    input$setGrButton
    isolate({
      rmin <- input$rminInclGr
      rmax <- input$rmaxInclGr
      dr <- input$drInclGr
      rho <- input$rhoInclGr
      if(input$GrNoiseType=="gauss")
        type="gaussianNoise"
      if(input$GrNoiseType=="correlated")        
        type="correlatedNoise"
    
      sigmaIsAvail <- whatIsSpecified(vals$dat)[[1]]$sigma
      
      if(!is.na(rho) && !is.na(rmin) && 
          !is.na(rmax) && !is.na(dr) && sigmaIsAvail){
        progress <- Progress$new(session)
        mess <- "Calculating, please wait..."
        progress$set(message = mess, value = 0.1)
        r1 <- seq(rmin, rmax, dr)
        for(i in 1:vals$nB){      
          dat <- list(x=vals$dat[[i]]$x, y=vals$dat[[i]]$y, sigma=vals$dat[[i]]$sigma)
          dat <- set.Gr(dat, r1=r1, rho.0=rho, type1=type)               
          progress$set(message = mess, value = i/vals$nB-0.01)
          vals$datGr[[i]] <- dat$Gr
        }
        progress$set(message = mess, value = 0.999)
        progress$close()         
      }
    })
  })
  
###################################
##                               ##
##        RENDERING OUTPUT       ##
##                               ##
###################################
####################################
## SORRY GR
  output$GrNoteForBanks <- renderUI({
    if(vals$nB > 1)
      return(span(h4("We recommend not to use this option for individual data banks!"), style = "color:red"))
    else
      return(NULL)
  })   
  

##########################################################################################
#                                                                                        #
#                   SET PARAMETERS FOR DIFEV AND DO FIT                                  #
#                                                                                        #
##########################################################################################  

  output$fitWithR <- renderUI({
    if( (vals$nB>1))
      return(
        radioButtons('fitWith', strong('Fit background with'),
          choices=c("spline functions"='fitWith.splines',
                    "analytical function"='fitWith.analyt'),
          selected='fitWith.analyt')
      )
    else
      return(
        radioButtons('fitWith', strong('Fit background with'),
          choices=c("spline functions"='fitWith.splines',
                    "analytical function"='fitWith.analyt'),
          selected='fitWith.splines')
      )
  })
  
  
######################## 
##  ==  Do Fit  ==  
  observe({
    if(input$doFit==0)
      return(NULL)
    isolate({    ## react on change
      is.x <- is.y <- is.sigma <- is.lambda <-TRUE
      is.NP <- is.F <- is.CR <- is.itermax <- TRUE
      is.bounds <- is.knots <- is.scale  <- TRUE
    
      wis <- whatIsSpecified(vals$dat)
      for(i in 1:vals$nB){
        if(!wis[[i]]$x)is.x <- FALSE
        if(!wis[[i]]$y)is.y <- FALSE
        if(!wis[[i]]$SB)is.SB <- FALSE
        if(!wis[[i]]$sigma)is.sigma <- FALSE
        if(!wis[[i]]$lambda)is.lambda <- FALSE
      }      
      if( !( is.numeric(input$fitNP) && (input$fitNP>2) ) )
       is.NP <- FALSE      
      if( !( is.numeric(input$fitItermax) && (input$fitItermax>2) ) )
       is.itermax <- FALSE  
      if( !( is.numeric(input$fitCR) && (input$fitCR>0) && (input$fitCR<1) ) )
       is.CR <- FALSE  
      if( !( is.numeric(input$fitF) && (input$fitF>0) && (input$fitF<2) ) )
       is.F <- FALSE  
      if( !( !is.na(input$bkgBounds) && (length(as.numeric(unlist(strsplit(input$bkgBounds, ","))))==2) && 
             !any(is.na(as.numeric(unlist(strsplit(input$bkgBounds, ","))))) ) )
       is.bounds <- FALSE          
      if( !( !is.na(input$fitKnots) && !any(is.na(as.numeric(unlist(strsplit(input$fitKnots, ","))))) ) )
       is.knots <- FALSE        
      if( !( !is.na(input$fitScale) && (length(as.numeric(unlist(strsplit(input$fitScale, ","))))==2) && 
             !any(is.na(as.numeric(unlist(strsplit(input$fitScale, ","))))) ) )
       is.scale <- FALSE          
          
      if(!is.x || !is.y || !is.sigma || !is.lambda || 
         !is.NP || !is.itermax || !is.CR || !is.F || 
         !is.bounds || !is.knots || !is.scale)
        return(NULL)
      if(!is.null(input$fitADP) && input$fitADP==TRUE){  #if we want to fit ADP
        if(!((!is.null(vals$datGr[[1]])) && (length(vals$datGr[[1]])>1)))
          return(NULL)        
      }
      CR <- input$fitCR
      F <- input$fitF
      NP <- input$fitNP
      itermax <- input$fitItermax
      p.bkg <- input$pbkg
      ctrl <- set.control(CR=CR, F=F, NP=NP, itermax=itermax, parallelType=1)

      bounds <- as.numeric(unlist(strsplit(input$bkgBounds, ",")))
      scale <- as.numeric(unlist(strsplit(input$fitScale, ",")))
      knots <- as.numeric(unlist(strsplit(input$fitKnots, ",")))
      if(is.null(knots) || any(is.na(knots)))
        knots <- 20
      knots.n <- knots.x <- NA
      if(length(knots)==1)
        knots.n <- knots
      else
        knots.x <- knots      
             
      progress <- Progress$new(session)
      mess <- "Calculating, please wait. This may take a while..."
      progress$set(message = mess, value = 0.1)     

      if(length(vals$datGr[[1]])>1)
        Gr <- vals$datGr[[i]]
      else 
        Gr <- NULL
      for(i in 1:vals$nB){
        dat <- list(x=vals$dat[[i]]$x, y=vals$dat[[i]]$y, SB=vals$dat[[i]]$SB, 
                    sigma=vals$dat[[i]]$sigma, lambda=vals$dat[[i]]$lambda, 
                    Gr=Gr, fitADP=vals$dat[[i]]$fitADP, id=vals$dat[[i]]$id) 
        if(vals$nB>1) 
          progress$set(message = mess, value = (i/vals$nB-0.01))
        else
          progress$set(message = mess, value = 0.5)
        analyt <- {input$fitWith=='fitWith.analyt'}
        vals$fitRes[[i]] <- do.fit(dat, bounds.lower=bounds[1], bounds.upper=bounds[2], 
                        scale=scale, knots.x=knots.x, knots.n=knots.n, analytical=analyt, 
                        stdev=TRUE, control=ctrl, p.bkg=p.bkg, save.to="")
      }
      progress$set(message = 'Calculating, please wait...', value = 0.999)  
      vals$fitResFinal <- list(list())
      cat("\n Done! \n")
      progress$close()
    }) 
  })
  
##########################################################################################
#                                                                                        #
#                                  FIT RESULTS                                           #
#                                                                                        #
##########################################################################################  

#################################   
##                             ##
##         DOWNLOAD            ##
##                             ##
#################################

###############################
## DOWNLOAD RDATA BUTTON
  output$downloadRDataR <- renderUI({
    if( (length(vals$fitRes[[1]]) > 1))
      return(downloadButton('downloadRData', 'Download fit results as .RData file'))
    else
      return(NULL)
  }) 
####################################
## DOWNLOAD TO RDATA FILE!
  output$downloadRData <- downloadHandler(filename = function() { paste('fit.results', '.RData', sep='') }, content = function(file) {
    fit.results <- vals$fitResFinal
    save(fit.results, file=file)
  })
  

###############################
## DOWNLOAD TEXT BUTTON
  output$downloadFitResAsTxtR <- renderUI({
    if( (length(vals$fitRes[[1]]) > 1) && (vals$nB == 1) )
      return(downloadButton('downloadFitResAsTxt', HTML('Download fit results as text file&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;')))
    else
      return(NULL)
  }) 
####################################
## DOWNLOAD TO TEXT FILE!
  output$downloadFitResAsTxt <- downloadHandler(filename = function() { paste('fit.results', '.txt', sep='') }, content = function(file) {
    fit.res <- vals$fitResFinal[[1]]
    write.fit.results(fit.res, file = "fit.tmp")
    writeLines(readLines("fit.tmp"), file)
    file.remove("fit.tmp")
  })
  
###############################
## DOWNLOAD FIX BUTTON
  output$downloadFixR <- renderUI({
    if( (length(vals$fitRes[[1]]) > 1) && (vals$nB >= 1) )
      return(downloadButton('downloadFix', HTML(paste("Download .fix file for", em("PDFgetN")))))
    else
      return(NULL)
  }) 
    
####################################
## DOWNLOAD TO FIX FILE!
  output$downloadFix <- downloadHandler(filename = function() { paste('corrections', '.fix', sep='') }, content = function(file) {
      fit.res <- vals$fitRes
      for(i in 1:vals$nB){
        N <- length(fit.res[[i]]$x)
        NInit <- length(vals$XInit[[i]])
        if(N < NInit){
           fit.res[[i]]$x <- vals$XInit[[i]]
           fit.res[[i]]$curves$bkg <- c(fit.res[[i]]$curves$bkg, rep(0, NInit-N))
        }      
      }
      
      write.fix(fit.res, file = "fix.tmp")
      writeLines(readLines("fix.tmp"), file)
      file.remove("fix.tmp")
  })

####################################
## APPEND FIX BUTTON  
  output$messageFixR <- renderUI({
    if( (length(vals$fitRes[[1]]) > 1) && (vals$nB >= 1) )
      return((h4("Append to existing .fix file")))
    else
      return(NULL)
  }) 
  output$selectFixR <- renderUI({
    if( (length(vals$fitRes[[1]]) > 1) && (vals$nB >= 1) )
      return( fileInput('fixfile', strong('Append to existing .fix file'), accept=c('.fix')) )
    else
      return(NULL)
  })   
    
   output$appendFixR <- renderUI({
    if( (length(vals$fitRes[[1]]) > 1) && (vals$nB >= 1) && !is.null(input$fixfile))
      return(downloadButton('appendFix', HTML(paste("Download it here"))))
    else
      return(NULL)
  })               
                
####################################
## APPEND TO FIX FILE!
  output$appendFix <- downloadHandler(filename = function() { paste('corrections', '.fix', sep='') }, content = function(file) {
      inFile <- input$fixfile   
      fit.res <- vals$fitRes
      for(i in 1:vals$nB){
        N <- length(fit.res[[i]]$x)
        NInit <- length(vals$XInit[[i]])
        if(N < NInit){
           fit.res[[i]]$x <- vals$XInit[[i]]
           fit.res[[i]]$curves$bkg <- c(fit.res[[i]]$curves$bkg, rep(0, NInit-N))
        }      
      }
  
      writeLines(readLines(inFile$datapath), "01x001.tmp")
      N <- length(fit.res)
      if(!is.null(fit.res$fit.details)){
        fit.res <- list(fit.res)
        N <- 1
      }
      options(warn=-1)
      for(i in 1:N){
        write(c(paste("#S ",i," Correction File for Bank ",fit.res[[i]]$fit.details$id,sep=""), "#L Q MULT ADD"), file="01x001.tmp", append=TRUE)
        res <- cbind(fit.res[[i]]$x, rep(1,length(fit.res[[i]]$x)), -fit.res[[i]]$curves$bkg)
        write.table(res, file="01x001.tmp", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")  
      }
      options(warn=0)

      writeLines(readLines("01x001.tmp"), file)
      file.remove("01x001.tmp")
     
  })
####################################
## DOWNLOAD GR AS TEXT!       
  output$downloadGrR <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0))
      return(NULL)
    else
      return(downloadButton('downloadGr', 'Download G(r) as text file'))
  }) 
  output$downloadGr <- downloadHandler(filename = function() { paste('Gr', '.txt', sep='') }, content = function(file) {
    gr <- vals$Gr
    write.table(data.frame(gr), file, row.names=FALSE, quote=FALSE, sep="\t") 
  }) 
#################################   
##                             ##
##         CALCULATE           ##
##                             ##
#################################

####################################
##    ==INPUTS FOR GR==
  output$outHeaderGr <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0))
      return(NULL)
    else
      return(h4("Calculate and plot G(r)"))
  }) 
###  
  output$rminCalcGrR <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0))
      return(NULL)
    else
      return(numericInput("rminCalcGr", min=0, max=100, step=0.1,
                label = strong("min(r)"), 
                value = 0))
  })   
  output$rmaxCalcGrR <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0))
      return(NULL)
    else
      return(numericInput("rmaxCalcGr", min=2, max=100, step=0.1,
                label = strong("max(r)"),
                value = 10))
  })    
   output$drCalcGrR <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0))
      return(NULL)
    else
      return(numericInput("drCalcGr", min=0.001, max=0.5, step=0.001,
                label = div( span(strong("grid spacing")), span(strong(em("dr"))) ), 
                value = 0.01))
  })   
   output$calcGrButtonR <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0))
      return(NULL)
    else
      return(actionButton("calcGrButton", label = strong("Calculate")))
  })  
  
####################################
##    ==CALCULATE GR==
  observe({
    if(is.null(input$calcGrButton))
      return(NULL)
    if(input$calcGrButton==0)
      return(NULL)
    isolate({
      if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0))
        return(NULL)
      fit.res <- vals$fitResFinal[[1]]
      progress <- Progress$new(session)
      mess <- "Calculating G(r), please wait..."
      progress$set(message = mess, value = 0.5)
      if(is.numeric(input$rhoInclGr)) 
        rho.0 <- input$rhoInclGr
      else
        rho.0 <- 0
      if(is.numeric(input$rminCalcGr)) 
        minR <- input$rminCalcGr
      else
        minR <- 0
      if(is.numeric(input$rmaxCalcGr)) 
        maxR <- input$rmaxCalcGr
      else
        maxR <- 10        
      if(is.numeric(input$drCalcGr)) 
        dr <- input$drCalcGr
      else
        dr <- 0.01  
      vals$Gr <- calc.Gr(fit.results=fit.res, rho.0=rho.0, r.min=minR, r.max=maxR, dr=dr, plot=FALSE)
      progress$set(message = "Done!", value = 0.999)
      progress$close()      
    })  
  })

  
###############################
##                           ##
##      = ITERATIONS =       ##
##                           ## 
###############################            

####################################
##  ==  INPUTS ==      
  output$iterHeader <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0) || 
        length(vals$datGr[[1]])==0 || !is.na(vals$fitRes[[1]]$pars))
      return(NULL)
    else
      return(h4("Perform iterative Bayesian background estimation"))
  }) 
  output$iterTechniqueR <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0) || 
        length(vals$datGr[[1]])==0 || !is.na(vals$fitRes[[1]]$pars))
      return(NULL)
    else
      return(radioButtons('iterTechnique', '',
              choices=c("Local gradient descent algorithm"='local', "Global DifEv algorithm"='global'),
              selected='global'
            ))
  })
  output$iterEpsR <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0) || 
        length(vals$datGr[[1]])==0 || !is.na(vals$fitRes[[1]]$pars) )
      return(NULL)
    else
      return(numericInput("iterEps", label = strong("Convergence tolerance"),
                min=0, max=0.1, step=1e-4, value = 1e-3))
  })  
  output$iterNIterR <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0) || 
        length(vals$datGr[[1]])==0 || !is.na(vals$fitRes[[1]]$pars))
      return(NULL)
    else
      return(numericInput("iterNIter", label = strong("The maximum iteration for a gradient descent method"),
                min=0, max=1e6, step=1e5, value = 1e5))
  })
  output$doIterationR <- renderUI({
    if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0) || 
        length(vals$datGr[[1]])==0 || !is.na(vals$fitRes[[1]]$pars) )
      return(NULL)
    else
      return(actionButton("doIteration", label = strong("Try iteration")))
  })

####################################
##  ==  PERFOMING ITERATION == 
  observe({
    if(is.null(input$doIteration))
      return(NULL)
    if(input$doIteration==0)
      return(NULL)
    isolate({
      if( (vals$nB>1) || (length(vals$fitRes[[1]]) == 0) || (length(vals$datGr[[1]])==0)  )
        return(NULL)
      fit.res <- vals$fitRes[[1]]
      progress <- Progress$new(session)
      mess <- "Calculating, please wait..."
      progress$set(message = mess, value = 0.3)
      rho.0 <- fit.res$fit.details$Gr$rho.0     
      if(is.numeric(input$iterNIter)) 
        n.iter <- input$iterNIter
      else
        n.iter <- 100000
      if(is.numeric(input$iterEps)) 
        eps <- input$iterEps
      else
        eps <- 1e-3
      if(input$iterTechnique=="local")
        local=TRUE
      else   
        local=FALSE     
      fit.res <- do.iter(fit.results=fit.res, local=local, eps=eps, n.iter=n.iter, save.to="")  
      vals$fitResIter[[1]] <- fit.res
      mess <- "Recalculating G(r)..."
      progress$set(message = mess, value = 0.75)
      # Recalculating G(r)...
      if(is.numeric(input$rminCalcGr)) 
        minR <- input$rminCalcGr
      else
        minR <- 0
      if(is.numeric(input$rmaxCalcGr)) 
        maxR <- input$rmaxCalcGr
      else
        maxR <- 10        
      if(is.numeric(input$drCalcGr)) 
        dr <- input$drCalcGr
      else
        dr <- 0.01       
      vals$Gr <- calc.Gr(fit.results=fit.res, rho.0=rho.0, r.min=minR, r.max=maxR, dr=dr, plot=FALSE)
      progress$set(message = mess, value = 0.99)
      progress$close()      
    })  
  })  

#### 
observe({
  if(length(vals$fitResIter[[1]])==0)
    vals$fitResFinal <- vals$fitRes
  else
    vals$fitResFinal <- vals$fitResIter
})

##########################################################################################
#                                                                                        #
#                                RENDERING OUTPUT                                        #
#                                                                                        #
##########################################################################################  

####################################
##  ==  OUTPUT TABLE ==
  output$datatable <- renderTable({
    if (length(vals$dat[[1]])==0)
      return(data.frame()) 
      
    dat.table <- list()
    for(i in 1:vals$nB){
      dat.table[[i]] <- unclass(vals$dat[[i]])
      dat.table[[i]]$fitADP <- dat.table[[i]]$Gr <- NULL
      dat.table[[i]] <- data.frame(dat.table[[i]])
      for(j in 1:length(colnames(dat.table[[i]]))){
        colnames(dat.table[[i]])[j] <- if(vals$nB==1) paste(colnames(dat.table[[i]])[j],sep="") else paste(colnames(dat.table[[i]])[j],toString(i), sep="") 
      }
    }
    k <- 1
    while(k < vals$nB){
      k <- k + 1
      dat.table[[1]] <- cbind(dat.table[[1]], dat.table[[k]])
    }
    return(dat.table[[1]])
  })

####################################
##  DOWNLOAD DATA
  output$downloadData <- downloadHandler(filename = function() { paste('data', '.txt', sep='') }, content = function(file) {
 #   if(length(vals$dat[[1]])==0)
 #     return(NULL)
    dat.table <- list()
    for(i in 1:vals$nB){
      dat.table[[i]] <- unclass(vals$dat[[i]])
      dat.table[[i]]$fitADP <- dat.table[[i]]$Gr <- NULL
      dat.table[[i]] <- data.frame(dat.table[[i]])
      for(j in 1:length(colnames(dat.table[[i]]))){
        colnames(dat.table[[i]])[j] <- if(vals$nB==1) paste(colnames(dat.table[[i]])[j],sep="") else paste(colnames(dat.table[[i]])[j],toString(i), sep="") 
      }
    }
    k <- 1
    while(k < vals$nB){
      k <- k + 1
      dat.table[[1]] <- cbind(dat.table[[1]], dat.table[[k]])  
    }
    
    write.table(dat.table[[1]], file, row.names=FALSE, quote=FALSE, sep="\t") 
  })  

  

####################################
##  ==  OUTPUT DATA PLOT ==

###############
# SELECT BANK 
   output$selectBank <- renderUI({
    if (vals$nB==1)
       return(NULL)   
    choices <- list()
    for(i in 1:vals$nB){
      name <- paste("Showing: Bank #", i)
      id <- paste(i)
      choices[[name]] <- id      
    }
    return(
      selectInput("bankNo", label = "",
                  choices = choices,
                  selected = "1",
                  width='160px')
    )     
  })   
###############
# PLOT FUNCTION
  dataPlotFunc <- function(onHover=TRUE){
    dat <- vals$dat
    toPlot <- whatIsSpecified(dat)
    N <- vals$nB
    n.x <- n.y <- 1
    if(N>=2) n.y <- 2
    if(N>=3) n.x <- 2
    par(mfrow=c(1, 1), mar=c(5,4,1,1))
#    par(oma = c(2, 1, 1, 1))
    if(!is.null(input$bankNo))
      bankNo <- as.numeric(input$bankNo)
    else
      bankNo <- 1      
    if(N==1){
      xlab=paste("x")    
      ylab=paste("y")    
    }
    else{
      xlab=paste("x", bankNo, sep="")    
      ylab=paste("y", bankNo, sep="")  
    }
    xlim <- c(min(vals$dat[[bankNo]]$x), max(vals$dat[[bankNo]]$x))
    ylim <- c(min(vals$dat[[bankNo]]$y), max(vals$dat[[bankNo]]$y))
    if(!is.null(input$selectPlot) && input$selectPlot==paste("bank", bankNo, sep="")){
      xlim <- input$plotLimX
      ylim <- input$plotLimY
    }
    plot(x=dat[[bankNo]]$x, y=dat[[bankNo]]$y, t="l", xlab=xlab, ylab=ylab, 
         xlim=xlim, ylim=ylim, lwd=2)
    par(xpd=TRUE)
    if(onHover){      
      hover <- input$mainHover
      if(!is.null(hover)){
        abline(v=hover$x, lty=2)
        abline(h=hover$y, lty=2)
        legend(hover$x, hover$y, sprintf("x=%.4g   y=%.4g", hover$x, hover$y), bty="n", pt.lwd=0, text.col=2, cex=0.7)
      }
      click <- input$mainClick
      if(!is.null(click)){
        input$mainClick
        isolate({
          abline(v=click$x, lty=2)
          abline(h=click$y, lty=2)
          legend(click$x, click$y, sprintf("x=%.4g   y=%.4g", click$x, click$y), bty="n", pt.lwd=0, text.col=2, cex=0.7)
        })
      }
    }
    par(xpd=FALSE)
    if(toPlot[[bankNo]]$SB)  lines(dat[[bankNo]]$x, dat[[bankNo]]$SB, col=3, lwd=2)
    if(toPlot[[bankNo]]$sigma){
      if(toPlot[[bankNo]]$smoothed){
        lines(dat[[bankNo]]$x, dat[[bankNo]]$smoothed, col="cyan", lwd=2)
        lines(dat[[bankNo]]$x, dat[[bankNo]]$smoothed+2*dat[[bankNo]]$sigma, col=2)
        lines(dat[[bankNo]]$x, dat[[bankNo]]$smoothed-2*dat[[bankNo]]$sigma, col=2)       
      }
      else{
        lines(dat[[bankNo]]$x, dat[[bankNo]]$y+2*dat[[bankNo]]$sigma, col=2)
        lines(dat[[bankNo]]$x, dat[[bankNo]]$y-2*dat[[bankNo]]$sigma, col=2)       
      }     
    }
    if(toPlot[[bankNo]]$lambda)  lines(dat[[bankNo]]$x, dat[[bankNo]]$lambda, col=6, lwd=2)            
  }     
###############
# PLOT RENDER
  output$dataPlot <- renderPlot({
    dat <- vals$dat
    if (length(dat[[1]])==0)
       return(NA)   
    toPlot <- whatIsSpecified(dat)
    if (!toPlot[[1]]$x || !toPlot[[1]]$y)
      return(NA)
      
    dataPlotFunc()
  })
  
  legendPlotFunc <- function(){
    par(mfrow=c(1,1), mar=c(1, 2, 2, 2) + 0.1) 
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#    par(fig = c(0, 1, 0, 1), oma = c(3, 3, 3, 3), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")    
    legend("bottom", c("data", "baseline", "lambda", "smoothed", "+/-2*stdev"), xpd = TRUE, horiz = TRUE, 
           inset = c(0,0), bty = "n", lty=1, col = c(1,3,6,"cyan",2), lwd=2, cex = 1)#      par(xpd=FALSE)   
  }  
  output$legendPlot <- renderPlot({
    if (length(vals$dat[[1]])==0)
       return(NA)   
    legendPlotFunc()
  })    
  
###############
# DOWNLOAD BUTTON
  output$downloadMainPlotR <- renderUI({
    dat <- vals$dat
    if (length(dat[[1]])==0)
       return(NULL)   
    toPlot <- whatIsSpecified(dat)
    if (!toPlot[[1]]$x || !toPlot[[1]]$y)
      return(NULL)
    
    return(downloadButton('downloadMainPlot', 'Download plot'))     
   })  

####################################
## DOWNLOAD HADLER
  output$downloadMainPlot <- downloadHandler(
      filename = function() { 'data.png' }, 
      content = function(file) {
        plotToPng <- function(){
          dataPlotFunc(onHover=FALSE)
          legendPlotFunc()
        }
        png(file, width=12, height=8, units="in", res=600,  pointsize=12)
        print(plotToPng())
        dev.off()
     }
  )

###############
# DOWNLOAD BUTTON
  output$downloadestGrPlotR <- renderUI({
    dat <- vals$dat
    if (length(dat[[1]])==0)
       return(NULL)   
    if (length(vals$estGr)==0)
       return(NULL)  
    toPlot <- whatIsSpecified(dat)
    if (!toPlot[[1]]$x || !toPlot[[1]]$y)
      return(NULL)
    
    return(downloadButton('downloadestGrPlot', 'Download plot'))     
   })      
####################################
## DOWNLOAD HADLER
  output$downloadestGrPlot <- downloadHandler(
      filename = function() { 'estGr.png' }, 
      content = function(file) {
        PDF <- vals$estGr
        stdev <- PDF$stdev*2
        gr <- PDF$gr
        r <- PDF$r
        rho.0 <- 0
        xlim <- ylim <- NA
        if(!is.null(input$selectPlot) && input$selectPlot==paste("estgr")){
          xlim <- input$plotLimX
          ylim <- input$plotLimY
        }   
        png(file, width=12, height=8, units="in", res=600,  pointsize=12)
        print(fplot.Gr(r=r, gr=gr, stdev=stdev, rho.0=rho.0, xlim=xlim, ylim=ylim, title="Estimated G(r)"))                       
        dev.off()
     }
  )   
  
  
###############
# DOWNLOAD BUTTON
  output$downloadestGrDataR <- renderUI({
    dat <- vals$dat
    if (length(dat[[1]])==0)
       return(NULL)   
    if (length(vals$estGr)==0)
       return(NULL)  
    toPlot <- whatIsSpecified(dat)
    if (!toPlot[[1]]$x || !toPlot[[1]]$y)
      return(NULL)
    
    return(downloadButton('downloadestGrData', 'Download G(r) as text file'))     
   })      
####################################
## DOWNLOAD HADLER
  output$downloadestGrData <- downloadHandler(
      filename = function() { 'estGr.txt' }, 
      content = function(file) {
        PDF <- vals$estGr
        stdev <- PDF$stdev
        gr <- PDF$gr
        r <- PDF$r
        write.table(data.frame(r, gr, stdev, gr-2*stdev, gr+2*stdev), file, row.names=FALSE, quote=FALSE, sep="\t") 
     }
  )    
  
 #############################################
## SHOWS PROGRESS IN PARAMETER ESTIMATIONS
  output$progress <- renderUI({
#    if(length(vals$dat[[1]])==0)
#      return(h3(" "))
    
    turnGreen <- whatIsSpecified(vals$dat)
    
    x.pr <- span(" x ", style = "color:#33CC00") 
    y.pr <- span("y ", style = "color:#33CC00")
    SB.pr <- span("SB ", style = "color:#33CC00")
    sigma.pr <- span(HTML("&epsilon; "), style = "color:#33CC00")
    lambda.pr <- span(HTML("&lambda; "), style = "color:#33CC00")
    Gr.pr <- span("G(r) ", style = "color:#33CC00")
    DifEv.pr <- span("DifEv ", style = "color:#33CC00")
    
    
   # write(vals$ind, file="aaa.txt")      
    for(i in 1:vals$nB){
      if(!turnGreen[[i]]$x) x.pr <- span(" x ", style = "color:red")
      if(!turnGreen[[i]]$y) y.pr <- span("y ", style = "color:red")
      if(!turnGreen[[i]]$SB) SB.pr <- span("SB ", style = "color:#B8B8B8")
      if(!turnGreen[[i]]$sigma) sigma.pr <- span(HTML("&epsilon; "), style = "color:red")
      if(!turnGreen[[i]]$lambda) lambda.pr <- span(HTML("&lambda; "), style = "color:red")   
 #     write(turnGreen[[i]], file="aa.txt")      
    }
    
    if(!((!is.null(vals$datGr[[1]])) && (length(vals$datGr[[1]])>1)))
      if(is.null(input$fitADP) || input$fitADP==FALSE)
        Gr.pr <- span("G(r) ", style = "color:#B8B8B8 ")
      else
        Gr.pr <- span("G(r) ", style = "color:red")
    
    DifEv <- TRUE
    if( !( is.numeric(input$fitNP) && (input$fitNP>2) ) )
      DifEv <- FALSE      
    if( !( is.numeric(input$fitItermax) && (input$fitItermax>2) ) )
      DifEv <- FALSE  
    if( !( is.numeric(input$fitCR) && (input$fitCR>0) && (input$fitCR<1) ) )
      DifEv <- FALSE  
    if( !( is.numeric(input$fitF) && (input$fitF>0) && (input$fitF<2) ) )
      DifEv <- FALSE  
    if( is.null(input$bkgBounds) || is.na(input$bkgBounds) || (length(as.numeric(unlist(strsplit(input$bkgBounds, ","))))!=2) || 
           any(is.na(as.numeric(unlist(strsplit(input$bkgBounds, ","))))) )
      DifEv <- FALSE          
    if( is.na(input$fitKnots) || any(is.na(as.numeric(unlist(strsplit(input$fitKnots, ","))))) ) 
      DifEv <- FALSE        
    if( !( !is.na(input$fitScale) && (length(as.numeric(unlist(strsplit(input$fitScale, ","))))==2) && 
           !any(is.na(as.numeric(unlist(strsplit(input$fitScale, ","))))) ) )
      DifEv <- FALSE          

        
    if(DifEv==FALSE) DifEv.pr <- span("DifEv ", style = "color:red")     
                
##  returns
    h3(x.pr, y.pr, lambda.pr, SB.pr, sigma.pr,  Gr.pr, DifEv.pr, align="left")   
  })  

##################################################   
##                                              ##
##            RENDER FIT RESULTS                ##
##                                              ##
##################################################

  observe({
    input$selectPlot
    input$plotLimY
    input$plotLimX
    isolate({
      if(is.null(dim(vals$xlim)) || is.null(dim(vals$ylim)) || dim(vals$xlim)!=c(vals$nB,2) || dim(vals$ylim)!=c(vals$nB,2))
        vals$xlim <- vals$ylim <- matrix(NA, nrow=vals$nB, ncol=2)
      for(i in 1:vals$nB){      
        if(!is.null(input$selectPlot) && input$selectPlot==paste("fit", i, sep="")){
          vals$xlim[i,] <- input$plotLimX
          vals$ylim[i,] <- input$plotLimY
        }
      }
    })
  })


####################################
##  ==  FIT RESULTS PLOT -- SQ ==
  output$fitResPlot <- renderPlot({
    if( (length(vals$fitRes[[1]]) > 1) ){
      fit.res <- vals$fitRes
      xlim <- vals$xlim
      ylim <- vals$ylim
      N <- vals$nB
      if(N>1)
        return(mPlot.results.banks(fit.res, xlim=xlim, ylim=ylim))
      else
        return(mPlot.results(fit.res[[1]], xlim=xlim, ylim=ylim))
    }  
    else
      return(NA)  
  })
#################
# DOWNLOAD BUTTON
  output$downloadFitResPlotR <- renderUI({
    if( (length(vals$fitRes[[1]]) > 1))    
      return(downloadButton('downloadFitResPlot', 'Download plot'))
    else
      return(NULL)    
   })  

##################
## DOWNLOAD HADLER
  output$downloadFitResPlot <- downloadHandler(
      filename = function() { 'fitPlot.png' }, 
      content = function(file) {
        xlim <- vals$xlim
        ylim <- vals$ylim
        png(file, width=12, height=8, units="in", res=600,  pointsize=12)
        print(if(vals$nB>1) {mPlot.results.banks(vals$fitRes, xlim=xlim, ylim=ylim)} 
               else {mPlot.results(vals$fitRes[[1]], xlim=xlim, ylim=ylim)} )
        dev.off()
     }
  )
  
####################################
##  ==  FIT RESULTS PLOT -- Gr ==
  output$GrPlot <- renderPlot({      
    if( (length(vals$fitRes[[1]]) > 1) && (vals$nB==1) && (length(vals$Gr)!=0)){
      PDF <- vals$Gr
      stdev <- PDF$stdev*2
      gr <- PDF$gr
      r <- PDF$r
      rho.0 <- 0
      if(!is.null(input$rhoInclGr) && is.numeric(input$rhoInclGr))
        rho.0 <- input$rhoInclGr
      if(!is.null(vals$datGr[[1]]$rho.0))
        rho.0 <- vals$datGr[[1]]$rho.0
      
      xlim <- ylim <- NA
      if(!is.null(input$selectPlot) && input$selectPlot==paste("gr")){
        xlim <- input$plotLimX
        ylim <- input$plotLimY
      }
      fplot.Gr(r=r, gr=gr, stdev=stdev, rho.0=rho.0, xlim=xlim, ylim=ylim)
    }  
    else
      return(NA) 
  })            

####################################
##  ==  FIT RESULTS PLOT -- Gr ==
  output$prelimGrPlot <- renderPlot({      
    if(length(vals$estGr)!=0){
      PDF <- vals$estGr
      stdev <- PDF$stdev*2
      gr <- PDF$gr
      r <- PDF$r
      rho.0 <- 0
      if(!is.null(input$rhoInclGr) && is.numeric(input$rhoInclGr))
        rho.0 <- input$rhoInclGr
      if(!is.null(vals$datGr[[1]]$rho.0))
        rho.0 <- vals$datGr[[1]]$rho.0
   
      xlim <- ylim <- NA
      if(!is.null(input$selectPlot) && input$selectPlot==paste("estgr")){
        xlim <- input$plotLimX
        ylim <- input$plotLimY
      }   
      fplot.Gr(r=r, gr=gr, stdev=stdev, rho.0=rho.0, xlim=xlim, ylim=ylim, title="Estimated G(r)")
    }  
    else
      return(NA) 
  })
  
#################
# DOWNLOAD BUTTON
  output$downloadGrPlotR <- renderUI({
    if( (length(vals$fitRes[[1]]) > 1) && (vals$nB==1) && (length(vals$Gr)!=0))    
      return(downloadButton('downloadGrPlot', 'Download plot'))
    else
      return(NULL)    
   })  

##################
## DOWNLOAD HADLER
  output$downloadGrPlot <- downloadHandler(
      filename = function() { 'Gr.png' }, 
      content = function(file) {
        rho.0 <- if(!is.null(vals$datGr[[1]]$rho.0)) {vals$datGr[[1]]$rho.0} else {input$rhoInclGr}
        xlim <- ylim <- NA
        if(!is.null(input$selectPlot) && input$selectPlot==paste("gr")){
          xlim <- input$plotLimX
          ylim <- input$plotLimY
        }        
        png(file, width=12, height=8, units="in", res=600,  pointsize=12)
        print(fplot.Gr(r=vals$Gr$r, gr=vals$Gr$gr, stdev=vals$Gr$stdev*2, 
                       rho.0=rho.0, xlim=xlim, ylim=ylim))
        dev.off()
     }
  )

#  observe({
#      # Initially will be empty
#      if (is.null(input$mainClick)){
#        return(NULL)
#      }
#      if (input$selectRegion==0){
#        return(NULL)
#      }
#      isolate({
#        vals$xlim[vals$selInd] <- input$mainClick$x
#        vals$ylim[vals$selInd] <-input$mainClick$y
#        if(vals$selInd==1)
#          vals$selInd <- 2
#        else 
#          vals$selInd <- 1
#      })    
#  })
    
#  observe({
#  input$resetRegion
#    if (input$resetRegion==0)
#      return(NULL)

#    isolate({
#      if(!is.null(vals$dat[[1]]$x) && !is.null(vals$dat[[1]]$y)){
#        x <- vals$dat[[1]]$x
#        y <- vals$dat[[1]]$y
#        vals$xlim <- c(min(x), max(x))
#        vals$ylim <- c(min(y), max(y))
#      }
#    })    
#  })   
    
#   output$lims <- renderTable({
#    if (length(vals$dat[[1]])==0)
#      return(data.frame()) 
   
#    dat.table <- matrix(c(vals$xlim[1], vals$xlim[2], vals$ylim[1], vals$ylim[2]), nrow=2, ncol=2, byrow=FALSE)
    
#    dat.table  <- data.frame(dat.table)
#    return(dat.table)
#  })

####################################
##  ==  PLOT OPTIONS ==  
  output$selectPlotR <- renderUI({
    if(length(vals$dat[[1]])==0)
      return(NULL)   
     
    choices <- list()

      # BANKS    
    if (vals$nB>1){
      for(i in 1:vals$nB){
        name <- paste("Data bank #", i)
        id <- paste("bank", i, sep="")
        choices[[name]] <- id      
      }
      if(length(vals$fitRes[[1]]) > 1){
        for(i in 1:vals$nB){
          id <- paste("fit", i, sep="")
          name <- paste("Background estimation for bank #", i)
          choices[[name]] <- id 
        }
      }      
    } # SINGLE DATASET
    else{
      choices <- list("Data plot"=paste("bank", 1, sep=""))
      if(length(vals$fitRes[[1]]) > 1)
        choices[["Background estimation"]] <- paste("fit", 1, sep="")
      
      if(length(vals$Gr)!=0)
        choices[["Corrected G(r)"]] = "gr" 
      if(length(vals$estGr)!=0)
        choices[["Estimated G(r)"]] = "estgr"         
    }
 
    return(
      selectInput("selectPlot", 
                  label = strong("Select plot to change"),
                  choices = choices,
                  width="100%")
    )     
  }) 
  
  output$youCanSeePlot <- renderUI({
    if(length(vals$dat[[1]])==0 || is.null(input$selectPlot))
      return(NULL)
    
    selectPlot <- substr(input$selectPlot, 1, 3)
    if(selectPlot=="ban" || selectPlot=="est")
      s1 <- div(span("(you can find it on the"), span(em("'Data Plot'"), style = "color:#0000FF;"), span("inset)"))
    else    
      s1 <- div(span("(you can find it on the"), span(em("'Fit Results Plot'"), style = "color:#0000FF;"), span("inset)"))
      
    return(s1)
  }) 
   
  output$axisLimsTxt <- renderUI({
    if(length(vals$dat[[1]])==0 || is.null(input$selectPlot))
      return(NULL)
        
    return(strong("Set axis limits"))
  }) 
  
  output$plotLimXR <- renderUI({
    if(length(vals$dat[[1]])==0 || is.null(input$selectPlot))
      return(NULL)  
    wis <- whatIsSpecified(vals$dat)
    if (!wis[[1]]$x)
      return(NULL)  
    ps <- input$selectPlot    
    if (vals$nB>1){
      for(i in 1:vals$nB){
        fitN <- paste("fit", i, sep="")      
        bankN <- paste("bank", i, sep="")
        if(ps==fitN || ps==bankN){
          minX <- min(vals$dat[[i]]$x)
          maxX <- max(vals$dat[[i]]$x)
        } 
      }
    }
    else{
      if(ps=="gr"){
        minX <- min(vals$Gr$r)
        maxX <- max(vals$Gr$r)        
      }
      else if(ps=="estgr"){
        minX <- min(vals$estGr$r)
        maxX <- max(vals$estGr$r)        
      }
      else{
        minX <- min(vals$dat[[1]]$x)
        maxX <- max(vals$dat[[1]]$x)      
      }
    
    }
    
    dx=(maxX-minX)/1000
    return(sliderInput("plotLimX", strong("x limits"), 
             min = (minX-0.1*abs(minX)), max = (maxX+0.1*abs(maxX)), 
             step=dx, value = c(minX, maxX))
    )     

  })   
  
  output$plotLimYR <- renderUI({
    if(length(vals$dat[[1]])==0 || is.null(input$selectPlot))
      return(NULL)  
    wis <- whatIsSpecified(vals$dat)
    if (!wis[[1]]$y)
      return(NULL) 

    ps <- input$selectPlot    
    if (vals$nB>1){
      for(i in 1:vals$nB){
        fitN <- paste("fit", i, sep="")      
        bankN <- paste("bank", i, sep="")
        if(ps==fitN || ps==bankN){
          if(wis[[i]]$SB){
            minY <- min(vals$dat[[i]]$y-vals$dat[[i]]$SB)
            maxY <- max(vals$dat[[i]]$y-vals$dat[[i]]$SB)
          }
          else{          
            minY <- min(vals$dat[[i]]$y)
            maxY <- max(vals$dat[[i]]$y)
          }
        }        
      }
    }
    else{
      if(ps=="gr"){
        minY <- min(vals$Gr$gr)
        maxY <- max(vals$Gr$gr)        
      }
      else if(ps=="estgr"){
        minY <- min(vals$estGr$gr)
        maxY <- max(vals$estGr$gr)        
      }
      else{
        if(wis[[1]]$SB){
          minY <- min(vals$dat[[1]]$y-vals$dat[[1]]$SB)
          maxY <- max(vals$dat[[1]]$y-vals$dat[[1]]$SB)
        }
        else{          
          minY <- min(vals$dat[[1]]$y)
          maxY <- max(vals$dat[[1]]$y)
        }    
      }
    
    }
    maxYnew <- vals$yRescale[2]*(maxY-minY) + minY
    minYnew <- vals$yRescale[1]*(maxY-minY) + minY
    maxY <- maxYnew
    minY <- minYnew
    dy=(maxY-minY)/2500
    return(sliderInput("plotLimY", strong("y limits"), 
             min = minY-0.4*abs(minY), max = maxY+0.4*abs(maxY), 
             step=dy, value = c(minY, maxY))
    )

  })   
        
        
  observe({
    if(input$rescaleY==0)
      return(NULL)
    isolate({    ## react on change
      if(length(vals$dat[[1]])==0 || is.null(input$selectPlot))
        return(NULL)  
      wis <- whatIsSpecified(vals$dat)
      if (!wis[[1]]$y)
        return(NULL) 

      ps <- input$selectPlot    
      if (vals$nB>1){
        for(i in 1:vals$nB){
          fitN <- paste("fit", i, sep="")      
          bankN <- paste("bank", i, sep="")
          if(ps==fitN || ps==bankN){
            if(wis[[i]]$SB){
              minY <- min(vals$dat[[i]]$y-vals$dat[[i]]$SB)
              maxY <- max(vals$dat[[i]]$y-vals$dat[[i]]$SB)
            }
            else{          
              minY <- min(vals$dat[[i]]$y)
              maxY <- max(vals$dat[[i]]$y)
            }
          }        
        }
      }
      else{
        if(ps=="gr"){
          minY <- min(vals$Gr$gr)
          maxY <- max(vals$Gr$gr)        
        }
        else if(ps=="estgr"){
          minY <- min(vals$estGr$gr)
          maxY <- max(vals$estGr$gr)        
        }
        else{
          if(wis[[1]]$SB){
            minY <- min(vals$dat[[1]]$y-vals$dat[[1]]$SB)
            maxY <- max(vals$dat[[1]]$y-vals$dat[[1]]$SB)
          }
          else{          
            minY <- min(vals$dat[[1]]$y)
            maxY <- max(vals$dat[[1]]$y)
          }    
        }
      
      }    
      ylim <- input$plotLimY
      vals$yRescale <- (ylim - minY)/(maxY-minY)      
    }) 
  })
  
   observe({
    if(input$resetY==0)
      return(NULL)
    isolate({  
      vals$yRescale <- c(0,1)     
    }) 
  })
  
})


  
  
  