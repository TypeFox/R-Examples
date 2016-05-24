## Disabling recessive and dominant, for the moment
## UPDATE: This code was right all along, the bug was not on my end! YAY!
powerDisabled <- function()
  return( !TRUE )


################################################
## Routines for calculating continuous offset ##
#fy_g <- function( y, x, aa ) {
#  dnorm( y-aa*x );
#}
Xg <- function( gCode, model ) {
  if( model=="additive" )
    return(gCode);
  if( model=="dominant" )
    return( as.integer(gCode==1 || gCode==2) )
  if( model=="recessive" )
    return( as.integer(gCode==2) )
  stop( "Xg, gCode not recognized." )
}
#fy <- function( y, p, aa, model ) {
#  return( y * ( fy_g(y,Xg(0,model),aa) * (1-p) * (1-p)
#                + fy_g(y,Xg(1,model),aa) * 2 * (1-p) * p
#                + fy_g(y,Xg(2,model),aa) * p * p ) )
#}

## this is taken directly from powerSim.cpp, if we ever find any errors.
calc_a <- function( model, p, heritability ) {
  h <- sqrt(heritability)
  SIGMA <- 1;
  varXi <- 0;
  if( model=="additive" ) {
    varXi <- 2*p*(1-p)
  }else if(model=="dominant" ){
    varXi <- 2*p - 5*p*p + 4*p*p*p - p*p*p*p
  }else if(model=="recessive" ){
    varXi <- p*p - p*p*p*p
  }

  aa <- SIGMA * h / sqrt( varXi * ( 1 - h*h ) )
  #cat( "a", aa, "\n" )

  ## newest code x=0,2 for dom/rec so comparable to additive?
  #if( model!="additive" )
  #  aa <- aa / 2

  ##print( aa ) ## debug only

  return( aa )
}

contsOffset <- function( model, p, heritability ) {
  aa <- calc_a( model, p, heritability )
  #return( integrate( fy, -Inf, Inf, p=p, aa=aa, model=model )$value )
  ex <- Xg(0,model)*(1-p)^2 + Xg(1,model)*2*p*(1-p) + Xg(2,model)*p^2;
  return( ex * aa )
}

##########################################################
## new power function, always Monte-Carlo, native style ##
pbat.powerCmd <- function( numOffspring=1, numParents=2, numFamilies=500,
                           additionalOffspringPhenos=TRUE,  ## don't remember messing w/ this in the c++ code
                           ascertainment="affected",
                           modelGen="additive", modelTest=modelGen,
                           afreqMarker=NA,
                           penAA=0.8, penAB=0.5, penBB=0.3,
                           heritability=0.0, contsAscertainmentLower=0.0, contsAscertainmentUpper=1.0,
                           pDiseaseAlleleGivenMarkerAllele=1.0, afreqDSL=0.1,
                           alpha=0.01,
                           offset="default",
                           numSim=1000,
                           ITERATION_KILLER=200 ) {
  ##if( powerDisabled() ) stop("Coming soon!")

  if( powerDisabled() && (modelGen!="additive" || modelTest!="additive" ) )
    stop( "Currently recessive and dominant are disabled." );

  ## Make sure afreqs are set OK
  if( is.na(afreqMarker) || pDiseaseAlleleGivenMarkerAllele==1.0 ) {
    pDiseaseAlleleGivenMarkerAllele <- 1.0;
    afreqMarker <- as.numeric(afreqDSL);
  }
  ##print( afreqMarker )

  ## calculate an offset?
  if( is.character(offset) && offset=="default" ) {
    if( heritability==0.0 ) {
      ## dichotomous trait - popn prevalence
      penAA <- as.numeric(penAA)
      penAB <- as.numeric(penAB)
      penBB <- as.numeric(penBB)
      afreqDSL <- as.numeric(afreqDSL)
      offset <- penAA*afreqDSL*afreqDSL + penAB*2*afreqDSL*(1-afreqDSL) + penBB*(1-afreqDSL)*(1-afreqDSL)
    }else{
      ## it's continuous
      afreqDSL <- as.numeric(afreqDSL)
      heritability <- as.numeric(heritability)
      offset <- contsOffset( modelGen, afreqDSL, heritability )
    }
  }
  ##cat( "offset calculated (R)", offset, "\n" )  ## debug only

  ascertainment <- match(ascertainment,c("affected","unaffected","na")) - 1
  modelGen <- match(modelGen,c("additive","dominant","recessive")) - 1
  modelTest <- match(modelTest,c("additive","dominant","recessive")) - 1
  additionalOffspringPhenos <- additionalOffspringPhenos==TRUE

  MODELS <- c("Additive", "Dominant", "Recessive")
  ASCERTAINMENT <- c("Affected", "Unaffected", "NA")

  if( heritability==0.0 ) {
    cat( "*** P2BAT power simulation (dichotomous trait) ***", "\n" )
  }else{
    cat( "*** P2BAT power simulation (continuous trait) ***", "\n" )
  }
  cat( "# offspring", numOffspring, ", # parents", numParents, ", # fams", numFamilies, "\n" )
  cat( "Phenotype additional offspring ", additionalOffspringPhenos, "\n" )
  cat( "Ascertainment", ASCERTAINMENT[ascertainment+1], "\n" )
  cat( "Model generation", MODELS[modelGen+1], "\n" )
  cat( "Model test", MODELS[modelTest+1], "\n" )
  cat( "Marker allele frequency", afreqMarker, "\n" )
  if( heritability==0.0 ) {
    cat( "Penetrances (AA,AB,BB)", penAA, penAB, penBB, "\n" )
  }else{
    cat( "Heritability", heritability, "\n" )
    cat( " Lower ascertainment", contsAscertainmentLower, "\n" )
    cat( " Upper ascertainment", contsAscertainmentUpper, "\n" )
  }
  cat( "P(Disease A|Marker A)", pDiseaseAlleleGivenMarkerAllele, " Disease allele frequency", afreqDSL, "\n" )
  cat( "Alpha", alpha, "\n" )
  cat( "Offset", offset, "\n" )
  cat( "Number of simulation iterations", numSim, "\n" )

  ## actually, we need to translate the contsAscertainments into 0-1, as they are _quantiles_
  if( any( contsAscertainmentLower>=contsAscertainmentUpper ) ) {
    print( "Continuous ascertainment conditions - the upper must have some spread between them." )
    return( 0 )
  }
  contsAscertainmentLower <- qnorm( as.numeric( contsAscertainmentLower ) )
  contsAscertainmentUpper <- qnorm( as.numeric( contsAscertainmentUpper ) )
  ## but R won't pass infinities very nice (is my guess)
  if( any( contsAscertainmentLower==-Inf ) )
    contsAscertainmentLower[contsAscertainmentLower==-Inf] <- -.Machine$double.xmax
  if( any( contsAscertainmentUpper==Inf ) )
    contsAscertainmentUpper[contsAscertainmentUpper==Inf] <- .Machine$double.xmax

  ##print( contsAscertainmentLower )
  ##print( contsAscertainmentUpper )

  result <- as.double(0.0)
  .C("powerR",
     as.integer(numOffspring), as.integer(numParents), as.integer(numFamilies),
     as.integer(additionalOffspringPhenos),
     as.integer(ascertainment),
     as.integer(modelGen), as.integer(modelTest),
     as.double(afreqMarker),
     as.double(penAA), as.double(penAB), as.double(penBB),
     as.double(heritability), as.double(contsAscertainmentLower), as.double(contsAscertainmentUpper),
     as.double(pDiseaseAlleleGivenMarkerAllele), as.double(afreqDSL),
     as.double(qnorm(1-alpha/2)),
     as.double(offset),
     as.integer(numSim),
     as.integer(ITERATION_KILLER),
     result,
     DUP=FALSE)
  if( result==-1 ) {
    #warning( "Invalid parameters, power is zeroed out." );
    print( "ERROR: Invalid parameters, power is zeroed out." );
    ##result <- 0;
  }else if( result==-2 ) {
    print( "ERROR: Too many iterations. Probability calculation is not converging." );
    ##result <- 0;
  }
  cat( "power", result, "\n\n" )
  return( result ) ## the power!
}

## for debug only
#pbat.powerCmd(numFamilies=500,numSim=1000, penAA=0.2, penAB=0.4, penBB=0.6);
#pbat.powerCmd(numFamilies=500,numSim=1000, penAA=0.6, penAB=0.4, penBB=0.2);
