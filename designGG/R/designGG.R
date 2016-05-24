# Function name: designGG 
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


designGG <- function( genotype, nSlides, nTuple, nEnvFactors, nLevels,
                      Level=NULL,  bTwoColorArray=TRUE, initial=NULL, weight=1,
                      region=NULL, optimality="A", method="SA",nIterations=3000,
                      n.search=2, endTemp=1e-10, startTemp=1, maxTempStep=0.9,
                      plotScores=TRUE, directory=NULL  ,fileName=NULL, envFactorNames=NULL,
                      writingProcess=TRUE)
{
    #parameter setting and checking
    nRILs       <- ncol(genotype)
    if( is.null(colnames(genotype) ) )
        colnames(genotype) <- paste( "Strain", 1:ncol(genotype), sep="")
    if( is.null(rownames(genotype) ) )
        rownames(genotype) <- paste( "marker", 1:nrow(genotype), sep="")
    
    if(("A" %in% unique(as.vector(as.matrix(genotype)))) |
      ("B" %in% unique(as.vector(as.matrix(genotype)))))  {
        temp                <- genotype
        temp[genotype=="A"] <- 1
        temp[genotype=="B"] <- 0
        if("H" %in% unique(as.vector(as.matrix(genotype)))) {
           temp[genotype=="H"] <- 0.5
        }
        temp                <- matrix(as.numeric(as.matrix(temp)),
                                          nrow=nrow(genotype))
        dimnames(temp)      <- dimnames(genotype)
        genotype            <- data.frame(temp)
      }
    genotype  <- data.frame(genotype)
        
    if(!is.null(Level)){
         if( any(!(unlist(lapply(Level,length) )==nLevels)))
         stop("The number of levels specified for each environmenal factor (Level) 
         is not consistent to nLevels.\n Please check your settings of nLevels and Level.")    
         if( any(unlist(lapply(Level,duplicated))))
         stop("The levels for certain environmenal factor cannot be duplicated.")
    }
    
    if( is.null(envFactorNames))
    {
        envFactorNames <- paste( "F", seq(1:nEnvFactors), sep="")
    }    
    if( nEnvFactors==1 & all(nLevels==1)){
        if ( is.null(nSlides) )     nSlides <- nTuple/2
        if ( is.null(nTuple) )      nTuple  <- nSlides *2
        if ( bTwoColorArray == FALSE )  stop ("For single channel experiment with no
                                          environmental factors, there is no need
                                          to optimize allocation of samples to
                                          array and different conditions using
                                          this algorithm")
        nEnvFactors <- 0
        nConditions <- 1
    }else{
        if( any(nLevels<=1) )
            stop("The number of levels for certain environmental factor shold be
                    an integer larger than 1.")                   
        
        nConditions <- prod(nLevels)
       
        if (is.null(nSlides) & is.null(nTuple))
           stop( "Either nSlides ( number of slides ) or nTuple( number of strains per condition ) \n is required for this algorithm." )
        
        if (!is.null(nSlides) & !is.null(nTuple)){
          warning("Either nSlides or nTuple is required. When both are specified, nSlides will used only.")
          nTuple  <-  NULL
        }
        if(!is.null(nSlides)) {
          if(as.integer(nSlides) != nSlides){
          nSlides  <- ceiling(nSlides)
          warning(paste("The argument nSlides should be an integer,nSlides of being ",nSlides, " \n has been used by designGG",sep="")) 
          }        
        }
        out     <- nSlidesANDnTuple (nLevels,nSlides,nTuple,bTwoColorArray)
        nSlides <- out[[1]]
        nTuple  <- out[[2]]
    }
    if( length(weight)!=1 )
    {
        varNumber <- variableNumber ( nEnvFactors )
        if( length(weight) != varNumber )
        {
          stop( "The length of weight should be equal to the number of
                  pamameters in the full model. Please check your setting of weight.")
        }
        weight  <- c(0, weight)
    }
    if( !(optimality == "D" | optimality == "A"))
    {
        stop( "This optimality criterion has not been implemented.
              Select `A' or `D'." )
    }
        
    if( is.null(region) )    region <- seq( 1, nrow(genotype) )
    
    genotype <- genotype[region,]

    if( is.null(directory) ) directory  <- getwd()
    if( is.null(fileName) )  fileName   <- "myDesignGG"

    temperature.step <- temperatureStep (startTemp, maxTempStep,endTemp,
                            nIterations)

    out.all <- end.score.all <- NULL
    for (i.search in 1:n.search){
        if ( !is.null(initial) )
        {
            now.array.allocation     <- now.array.allocation0 <- initial[[1]]
            now.condition.allocation <- now.condition.allocation0<- initial[[2]]
        }else{
            out                      <- initialDesign (genotype, nRILs, nSlides,
                                         nConditions, nTuple, bTwoColorArray)
            now.array.allocation     <- now.array.allocation0 <- out[[1]]
            now.condition.allocation <- now.condition.allocation0  <- out[[2]]
        }
        nowScore <- nowScore0 <- designScore( genotype, now.array.allocation,
                                          now.condition.allocation, nEnvFactors,
                                          nLevels, Level, nConditions, weight,
                                          optimality, bTwoColorArray,
                                          envFactorNames)
        accepted    <- 0
        scores      <- nowScore
        cooling     <- NULL

        for (i in 1:nIterations) # S.A. starts
        {
            proposal                 <- updateDesign ( now.array.allocation,
                                        now.condition.allocation, nRILs,
                                        nSlides, nEnvFactors, nTuple,
                                        bTwoColorArray )
            new.array.allocation     <- proposal [[1]]
            new.condition.allocation <- proposal [[2]]
            newScore                 <- designScore ( genotype,
                                        new.array.allocation,
                                        new.condition.allocation, nEnvFactors,
                                        nLevels, Level, nConditions, weight,
                                        optimality, bTwoColorArray,
                                        envFactorNames )
            temperature              <- startTemp * temperature.step^i

            acceptance.probability   <- acceptanceProbability( nowScore,
                                        newScore, method, temperature )
            if( runif(1) < acceptance.probability )
            {
                scores      <- c(scores, newScore)
                if( optimality == "A" )
                {
                    cooling <- c(cooling,
                                    ((1/newScore)-(1/nowScore))/(1/nowScore))
                }else{
                    cooling <- c(cooling, (newScore-nowScore)/nowScore)
                }
                now.array.allocation        <- new.array.allocation
                now.condition.allocation    <- new.condition.allocation
                nowScore                    <- newScore
                accepted                    <- accepted + 1
            }
            if( writingProcess )
            {
                processFile <- file("processing.txt","w")
                cat(format(round(100*(i+nIterations*(i.search-1))/(n.search*
                            nIterations),digits=2),nsmall=2), file = processFile)
                close(processFile)
            }
        }        #S,A ends

        if( scores[1] > scores[length(scores)] )
        {
            end.score                <- nowScore0
            end.array.allocation     <- now.array.allocation0
            end.condition.allocation <- now.condition.allocation0
        }else{
            end.score               <- nowScore
            end.array.allocation    <- now.array.allocation
            end.condition.allocation<- now.condition.allocation
        }

        end.score.all   <- c(end.score.all,end.score)
        out.isearch <- list(optimal.array.allocation = end.array.allocation,
                            optimal.condition.allocation =
                            end.condition.allocation,
                            scores = scores,optimality = optimality,
                            f.accepted = accepted/nIterations,
                            nIterations = nIterations, startTemp = startTemp,
                            temperature.step = temperature.step,
                            final.temp =  temperature)
        out.all     <- c(out.all,list(out.isearch))
    }
    out                       <- out.all[[which.max(end.score.all)[1]]]
    best.array.allocation     <- out[[1]]
    best.condition.allocation <- out[[2]]
    if( optimality == "A" ){
        best.scores               <- 1/out[[3]]
    }else{
        best.scores               <- out[[3]]
    }
    out.design                <- experimentDesignTable (best.array.allocation,
                                 best.condition.allocation, nEnvFactors,
                                 nLevels, Level,fileName, envFactorNames,
                                 directory)
   
    plot.obj                  <- list(scores=best.scores, cooling=cooling, 
                                  startTemp=startTemp, temperature=temperature,
                                  temperature.step=temperature.step, 
                                  nIterations=nIterations,optimality=optimality)
    
    if( all(!is.finite(best.scores)) | all(is.na(best.scores))){
      stop( "The design matrix is always singular. Please either try the 
      algorithm again by increasing the nIterations or change the input parameter e.g. nSlides.")  
    } 
    if( optimality == "A" )  best.scores <- 1/best.scores
    if ( plotScores )   plotAllScores(plot.obj,fileName)
    
    design.obj  <-  list(arrayDesign = out.design$arrayDesign, 
                            conditionDesign = out.design$conditionDesign,
                            plot.obj =plot.obj)
    return( design.obj)
}

