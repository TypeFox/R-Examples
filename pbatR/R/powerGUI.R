####################################################################
# Thomas Hoffmann                                                  #
# CREATED:     05/29/2007                                          #
# MODIFIED:    06/09/2007                                          #
#                                                                  #
# DESCRIPTION:                                                     #
#  The major pbat-gui interface for power                          #
####################################################################

####################################################################
#                                                                  #
# CONSTANTS                                                        #
#                                                                  #
####################################################################
SLIDER_LENGTH <- 500
POWER_TEXT_LENGTH <- 80

FIRE_UPDATE_PEN <- 1
FIRE_UPDATE_AOR <- 2
FIRE_UPDATE_OR <- 3
FIRE_UPDATE_AF <- 4

FAMILY_MIN <- 1 ## for debug, maybe 20 otherwise?

####################################################################
#                                                                  #
# SET UP THE GLOBAL VARIABLES INTERFACE                            #
#                                                                  #
####################################################################
powerEnv <- new.env();
setPower <- function( x, value )
  assign( x, value, envir=powerEnv );
getPower <- function( x, mode="any", ifnotfound=NA )
  get( x, envir=powerEnv, mode=mode, inherits=FALSE )
getPowerSafe <- function( x, ifnotfound=NA )
  mget( x, envir=powerEnv, ifnotfound=ifnotfound )[[1]]

getPowerV <- function( ... )
  return( tclvalue( getPower( ... ) ) )
getPowerN <- function( ... )
  return( as.numeric(tclvalue( getPower( ... ) ) ) )
getPowerT <- function( ... )
  return( tclvalue( getPower( ... ) ) == TRUE )
getPowerVSafe <- function( ..., default=NA ) {
  res <- getPowerSafe( ... )
  if( is.na(res) )
    return(default)
  return( tclvalue( res ) )
}
####################################################################
#                                                                  #
# Graphical helper functions.                                      #
#                                                                  #
####################################################################

## Creates a slider, etc.
newSlider <- function( text, default, min, max, step=(max-min)/100, grid, update=fireUpdate, state="enabled" )
{
  sliderVal <- tclVar( default )
  ##sliderValLabel <- tklabel( grid, text=as.character(tclvalue(sliderVal)) );
  sliderValLabel <- tkentry( grid, width=10, textvariable=sliderVal )
  slider <- tkscale( grid, from=min, to=max, resolution=step,
                     showvalue=FALSE, variable=sliderVal, orient="horizontal",
                     length=SLIDER_LENGTH )
  tkgrid( slider,
          tklabel(grid, text=text),
          sliderValLabel )
  tkconfigure( sliderValLabel,textvariable=sliderVal )
  if( state=="disabled" ) {
    tkconfigure( slider, state="disabled" )
    tkconfigure( sliderValLabel, state="disbled" )
  }

  ## bind the update firing
  tkbind(sliderValLabel, "<Return>", update)
  tkbind(slider, "<ButtonRelease-1>", update)

  #tkgrid.configure( slider, sticky="news" )

  return( sliderVal )  ## return the tclvariable corresponding to the new object
}

choicesEntry <- function( main, label, choices, defaultChoice=1 ){
  ## Create a frame for everything
  frame <- tkframe( main, relief="groove", borderwidth=2 );
  tkgrid( frame );
  ##tkgrid.configure( frame, sticky="news" );
  tkgrid.configure( frame, sticky="nws" );

  ## The label
  lbl <- tklabel( frame, text=label );

  ## The choices
  rb <- list();
  rb.lbl <- list();
  rb.subframe <- list();

  ## Create the subframe for each choice
  for( i in 1:length(choices) )
    rb.subframe[[i]] <- tkframe( frame, relief='groove', borderwidth=1 );

  ## grid the subframes
  if( length(choices)==1 ){ ## Not really a 'choice' then, but necessary
    tkgrid( lbl, rb.subframe[[1]] );
  }else if( length(choices)==2 ){
    tkgrid( lbl, rb.subframe[[1]], rb.subframe[[2]] );
  }else if( length(choices)==3 ){
    tkgrid( lbl, rb.subframe[[1]], rb.subframe[[2]], rb.subframe[[3]] );
  }else if( length(choices)==4 ){
    tkgrid( lbl, rb.subframe[[1]], rb.subframe[[2]], rb.subframe[[3]], rb.subframe[[4]] );
  }

  ## Create the variable
  rbVal <- tclVar( choices[defaultChoice] );

  ## Create each choice
  for( i in 1:length(choices) ) {
    rb[[i]] <- tkradiobutton( rb.subframe[[i]] );
    tkconfigure( rb[[i]], variable=rbVal, value=choices[i]  );
    rb.lbl[[i]] <- tklabel( rb.subframe[[i]], text=choices[i] );
    tkgrid( rb[[i]], rb.lbl[[i]] );

    ## how about binding the fire event?
    tkbind( rb[[i]], "<ButtonRelease-1>", fireUpdate )
  }

  ## Return the variable
  return( rbVal );
}

## sort of clean-up code if you will - it sets the tclvariable, but doesn't draw it
choicesEntryEmpty <- function( main, label, choices, defaultChoice=1 ) {
  rbVal <- tclVar( choices[defaultChoice] )
  return( rbVal )
}

dblGrid <- function( main ) {
  f <- tkframe(main, borderwidth=0 );
  f1 <- tkframe( f, borderwidth=0 );
  f2 <- tkframe( f, borderwidth=0 );
  fblank <- tkframe( f, borderwidth=0 );
  tkgrid( f );
  tkgrid.configure( f, sticky="nws" );
  tkgrid( f1, fblank, f2 );
  tkgrid( tklabel(fblank,text="    ") );
  return( list(f1=f1,f2=f2) );
}

#newTE <- function( tclVar, option, helps=NULL, width=10, gridframe=form ) {
#  frame <- newframe( borderwidth=3, gridframe=gridframe );
#  but <- NULL;
#  if( !is.null(helps) ) {
#    but <- tkbutton( frame, text="?",
#                      command=function(){
#                        msgText <- paste( option, "-", helps[1] );
#                        msg(msgText)
#                      } );
#  }
#  lab <- tklabel( frame, text=option );
#  entry <- tkentry( frame, width=width, textvariable=tclVar );
#  if( !is.null(but) ) {
#    tkgrid( lab, entry, but );
#  }else{
#    tkgrid( lab, entry );
#  }
#}

newframe <- function( gridframe, grid=TRUE, relief="groove", borderwidth=2, sticky="nws" ) {
  frame <- tkframe( gridframe, relief=relief, borderwidth=borderwidth );
  if( !grid ) return(frame);
  tkgrid( frame );
  tkgrid.configure( frame, sticky=sticky );
  return(frame)
};


####################################################################
#                                                                  #
# Drawing the main form.                                           #
#                                                                  #
####################################################################
createPowerGUI <- function( mode="continuous" )
{
  require( tcltk )
  main <- tktoplevel()
  setPower( "main", main )
  tkwm.title( main, paste("Power calculations (",mode," trait)",sep="") )

  ## newSlider( text, default, min, max, step, main )
  f <- newframe( main )
  setPower( "gNumOffspring", newSlider( "Number of Offspring", 1, 1, 20, 1, f ) )
  setPower( "gNumParents", newSlider( "Number of parents", 2, 0, 2, 1, f ) )
  setPower( "gNumFamilies", newSlider( "Number of families", 500, FAMILY_MIN, 10000, 10, f ) )
  setPower( "gAdditionalOffspringPhenos", choicesEntry( main, "Phenotype non-proband offspring?", c("TRUE","FALSE") ) )

  if( mode=="dichotomous" ) {
    f <- newframe( main )
    f2 <- dblGrid( f )
    setPower( "gAscertainment", choicesEntry( f2$f1, "Ascertainment", c("affected","unaffected","na") ) )
    tkgrid( tkbutton( f2$f2, text="Additional offspring ascertainment", command=createDichotAscertainmentGUI ) )
  }

  f <- newframe( main )
  setPower( "gTrait", choicesEntryEmpty( f, "Trait", c(mode,mode) ) )

  if( mode=="dichotomous" ) {
    #tkgrid( tklabel( f, text="Dichotomous trait options" ) )
    f2 <- newframe( f, relief="sunken", borderwidth=5 )
    tkgrid( tklabel( f2, text="Option 1: Genotype penetrances" ) )
    setPower( "gPenAA", newSlider( "AA Penetrance", 0.8, 0, 1, 0.001, f2, fireUpdatePen ) )
    setPower( "gPenAB", newSlider( "AB Penetrance", 0.5, 0, 1, 0.001, f2, fireUpdatePen ) )
    setPower( "gPenBB", newSlider( "BB Penetrance", 0.3, 0, 1, 0.001, f2, fireUpdatePen ) )

    f3 <- newframe( f, relief="sunken", borderwidth=5 )
    tkgrid( tklabel( f3, text="Option 2: MOI, one of AF/OR_het/AOR, and population prevalence (must select AF/OR_het/AOR first)" ) )
    if( !powerDisabled() ) {
      setPower( "gModelGen", choicesEntry( f3, "Genetic Model for Data Generation", c("additive","dominant","recessive") ) )
      setPower( "gModelTest", choicesEntry( f3, "Genetic Model for Testing", c("additive","dominant","recessive") ) )
    }else{
      setPower( "gModelGen", choicesEntry( f3, "Genetic Model for Data Generation", c("additive") ) )
      setPower( "gModelTest", choicesEntry( f3, "Genetic Model for Testing", c("additive") ) )
    }

    f4 <- newframe( f3 )
    #setPower( "gAF", newSlider( "Allelic Genetic Attributable Fraction", 0.2, 0, 1, 0.0001, f4, fireUpdateAF ) )
    setPower( "gAF", tclVar(0.2) )
    setPower( "gOR1", newSlider( "Heterozygous odds ratio", 0.2, 0, 50, 0.01, f4, fireUpdateOR ) )
    setPower( "gOR2", newSlider( "Homozygous odds raito", 0.2, 0, 50, 0.01, f4, state="disbaled" ) ) ## Want this _disabled_
    setPower( "gAOR", newSlider( "Allelic odds ratio", 0.2, 0, 50, 0.01, f4, fireUpdateAOR ) )

    f5 <- newframe( f3 )
    setPower( "gPopPrev", newSlider( "Population prevalence", 0, 0, 1, 0.001, f5 ) )
  }else{
    #f2 <- newframe( f, borderwidth=0 )
    #f3 <- dblGrid( f2 )
    #tkgrid( tklabel( f3$f1, text="Continuous trait options" ) )
    #tkgrid( tkbutton( f3$f2, text="Additional offspring ascertainment", command=createContsAscertainmentGUI ) )
    tkgrid( tkbutton( f, text="Additional offspring ascertainment criterion", command=createContsAscertainmentGUI ) )
    setPower( "gHeritability", newSlider( "heritability (h^2)", 0.1, 0, 1, 0.001, f ) )
    setPower( "gContsAscertainmentLower", newSlider( "Lower ascertainment criterion", 0, 0, 1, 0.01, f ) )
    setPower( "gContsAscertainmentUpper", newSlider( "Upper ascertainment criterion", 1, 0, 1, 0.01, f ) )

    f <- newframe( main )
    if( !powerDisabled() ) {
      setPower( "gModelGen", choicesEntry( f, "Genetic Model for Data Generation", c("additive","dominant","recessive") ) )
      setPower( "gModelTest", choicesEntry( f, "Genetic Model for Testing", c("additive","dominant","recessive") ) )
    }else{
      setPower( "gModelGen", choicesEntry( f, "Genetic Model for Data Generation", c("additive") ) )
      setPower( "gModelTest", choicesEntry( f, "Genetic Model for Testing", c("additive") ) )
    }
  }

  f <- newframe( main )
  setPower( "gUseOffset", choicesEntry( f, "Offset", c("default (population prevalence)","manual") ) )
  setPower( "gOffset", newSlider( "Offset", 0, 0, 1, 0.001, f ) )

  f <- newframe( main )
  setPower( "gAfreqDSL", newSlider( "DSL Allele Frequency", 0.1, 0, 1, 0.001, f ) )
  setPower( "gIsDSL", choicesEntry( f, "Marker is DSL", c("TRUE","FALSE") ) )
  setPower( "gDiseaseAlleleGivenMarkerAllele", newSlider( "P(Disease Allele|Marker Allele)", 1, 0, 1, 0.01, f ) )
  setPower( "gAfreqMarker", newSlider( "Marker Allele Frequency", 0.1, 0, 1, 0.001, f ) )

  f <- newframe( main )
  setPower( "gAlpha", newSlider( "Significance level (alpha)", 0.01, 0, 0.1, 0.0001, f ) )
  setPower( "gNumSim", newSlider( "Number of Simulations", 1000, 100, 10000, 100, f ) )

  f <- newframe( main )
  #setPower( "gPower", newSlider( "Power (Monte-Carlo simulation)", 0, 0, 1, 0.0001, f) ) ## would be _really_ neat if could alter the sample size
  setPower( "gPower", newSlider( "Power (Monte-Carlo simulation)", 0, 0, 1, 0.0001, f, fireUpdatePower) )

  ## NEW - ITERATIONS before giving up completely
  f <- newframe( main )
  setPower( "gITERATION_KILLER", newSlider( "Max Gen Draw Iter before halt (0 disables)", 500, 0, 10000, 1, f ) )

  ## new - how about a status window?
  f <- newframe( main )
  gStatus=tclVar("")
  statusEntry <- tkentry( f, width=POWER_TEXT_LENGTH, textvariable=gStatus )
  tkgrid( tklabel( f, text="Status" ),
          statusEntry )
  try( tkconfigure( statusEntry, state="readonly" ) )
  setPower( "gStatus", gStatus )
}

createDichotAscertainmentGUI <- function()
{
  numOffspring <- as.integer(tclvalue(getPower("gNumOffspring")))
  if( numOffspring==1 ) {
    tkmessageBox( title="No additional offspring", message="No additional offspring to alter the ascertainment status of.", icon="info", type="ok" )
    return(invisible())
  }

  main <- tktoplevel()
  for( i in 2:numOffspring )
    setPower( paste("gAscertainment",i,sep=""), choicesEntry( main, paste("Ascertainment",i), c("affected","unaffected","na"), defaultChoice=3 ) )

}
createContsAscertainmentGUI <- function()
{
  numOffspring <- as.integer(tclvalue(getPower("gNumOffspring")))
  if( numOffspring==1 ) {
    tkmessageBox( title="No additional offspring", message="No additional offspring to alter the ascertainment status of.", icon="info", type="ok" )
    return(invisible())
  }

  main <- tktoplevel()
  for( i in 2:numOffspring ) {
    setPower( paste("gContsAscertainmentLower",i,sep=""), newSlider( paste("Lower ascertainment criterion",i), 0, 0, 1, 0.1, main ) )
    setPower( paste("gContsAscertainmentUpper",i,sep=""), newSlider( paste("Upper ascertainment criterion",i), 1, 0, 1, 0.1, main ) )
  }
}

####################################################################
#                                                                  #
# Called on the update of any variable to recompute power.         #
#                                                                  #
####################################################################
fireUpdate <- function()
{
  ## New, we need to fix up things according to what was last done!
  if( getPowerV("gTrait") != "continuous" ) {
    updateLast <- getPowerSafe("updateLast", ifnotfound=FIRE_UPDATE_PEN)
    if( updateLast==FIRE_UPDATE_PEN ) fireUpdatePen_recalc()
    if( updateLast==FIRE_UPDATE_AOR ) fireUpdateAOR_recalc()
    if( updateLast==FIRE_UPDATE_OR )  fireUpdateOR_recalc()
    if( updateLast==FIRE_UPDATE_AF )  fireUpdateAF_recalc()
  }
  ## fires an update to the power call
  ##print( "fire in the hole!" ) ## debugging

  ## get all the global values
  numOffspring <- getPowerV( "gNumOffspring" )
  numParents <- getPowerV( "gNumParents" )
  numFamilies <- getPowerV( "gNumFamilies" )
  additionalOffspringPhenos <- getPowerT( "gAdditionalOffspringPhenos" )


  ascertainment <- rep(0,numOffspring)
  ascertainment[1] <- getPowerVSafe( "gAscertainment", default="na" )
  if( as.integer(numOffspring) > 1 ) {
    for( i in 2:as.integer(numOffspring) ) {
      ascertainment[i] <- getPowerVSafe( paste("gAscertainment",i,sep=""), default="na" )
    }
  }
  modelGen <- getPowerV( "gModelGen" )
  modelTest <- getPowerV( "gModelTest" )
  trait <- getPowerV( "gTrait" )

  penAA <- getPowerVSafe( "gPenAA", default=0.0 )
  penAB <- getPowerVSafe( "gPenAB", default=0.0 )
  penBB <- getPowerVSafe( "gPenBB", default=0.0 )

  offset <- "default"
  if( getPowerV("gUseOffset")=="manual" )
    offset <- getPowerV("gOffset")

  heritability <- getPowerVSafe( "gHeritability" )
  contsAscertainmentLower <- rep(0,numOffspring)
  contsAscertainmentUpper <- rep(0,numOffspring)
  contsAscertainmentLower[1] <- getPowerVSafe( "gContsAscertainmentLower", default=0.0 )
  contsAscertainmentUpper[1] <- getPowerVSafe( "gContsAscertainmentUpper", default=1.0 )
  if( as.integer(numOffspring)>1 ) {
    for( i in 2:as.integer(numOffspring) ) {
      contsAscertainmentLower[i] <- getPowerVSafe( paste("gContsAscertainmentLower",i,sep=""), default=0.0 )
      contsAscertainmentUpper[i] <- getPowerVSafe( paste("gContsAscertainmentUpper",i,sep=""), default=1.0 )
    }
  }

  afreqDSL <- getPowerV( "gAfreqDSL" )
  isDSL <- getPowerT( "gIsDSL" )
  diseaseAlleleGivenMarkerAllele <- getPowerV( "gDiseaseAlleleGivenMarkerAllele" )
  afreqMarker <- getPowerV( "gAfreqMarker" )

  alpha <- getPowerN( "gAlpha" )
  numSim <- getPowerV( "gNumSim" )

  ITERATION_KILLER <- getPowerV( "gITERATION_KILLER" )

  ## fix up some things
  if( isDSL==TRUE )
    afreqMarker <- NA;
  if( trait!="continuous" ) {
    heritability <- 0.0
  }else if( heritability == 0.0 ) {
    heritability <- 0.0001 ## It can't be zero, or it won't do continuous...
  }

  ## and update the power
  tkconfigure(getPower("main"),cursor="watch")
  newPower <- pbat.powerCmd( numOffspring, numParents, numFamilies,
                             additionalOffspringPhenos,
                             ascertainment,
                             modelGen, modelTest,
                             afreqMarker,
                             penAA, penAB, penBB,
                             heritability, contsAscertainmentLower, contsAscertainmentUpper,
                             diseaseAlleleGivenMarkerAllele,
                             afreqDSL,
                             alpha,
                             offset,
                             numSim,
                             ITERATION_KILLER )
  tkconfigure(getPower("main"),cursor="arrow")

  gStatus = getPower("gStatus")
  if( newPower<0 ) {
    if( newPower==-1 )
      tclvalue(gStatus) <- "ERROR: Invalid parameters. Generally this occurs when marker allele frequency is not possible."
    if( newPower==-2 )
      tclvalue(gStatus) <- "ERROR: Too many iterations without viable family. Try changing some of the parameters to be less extreme."
    newPower <- 0
  }else{
    tclvalue(gStatus) <- ""
  }

  ggPower <- getPower("gPower")
  tclvalue( ggPower ) <- as.character(newPower)
  #cat( "Power", newPower, "\n" )

}

fireUpdatePen <- function() {
  setPower( "updateLast", FIRE_UPDATE_PEN )
  fireUpdate()
}
fireUpdateAOR <- function() {
  setPower( "updateLast", FIRE_UPDATE_AOR )
  fireUpdate()
}
fireUpdateOR <- function() {
  setPower( "updateLast", FIRE_UPDATE_OR )
  fireUpdate()
}
fireUpdateAF <- function() {
  setPower( "updateLast", FIRE_UPDATE_AF )
  fireUpdate()
}

fireUpdatePen_recalc <- function(do_AOR=TRUE, do_OR=TRUE, do_AF=TRUE, do_pop=TRUE) {
  penAA <- as.numeric( getPowerVSafe( "gPenAA", default=0.0 ) )
  penAB <- as.numeric( getPowerVSafe( "gPenAB", default=0.0 ) )
  penBB <- as.numeric( getPowerVSafe( "gPenBB", default=0.0 ) )

  afreq <- as.numeric( getPowerV( "gAfreqDSL" ) )

  aor <- allelic_OR( penAA, penAB, penBB, afreq )
  or <- calc_or( penAA, penAB, penBB, afreq )
  ar <- calc_ar( penAA, penAB, penBB, afreq )

  gAOR <- getPower("gAOR")
  if( do_AOR ) tclvalue( gAOR ) <- as.character(aor)

  gOR1 <- getPower( "gOR1" )
  if( do_OR ) tclvalue( gOR1 ) <- as.character(or[1])
  gOR2 <- getPower( "gOR2" )
  tclvalue( gOR2 ) <- as.character(or[2])

  gAF <- getPower( "gAF" )
  if( do_AF ) tclvalue( gAF ) <- as.character(ar)

  popPrev <- getPower( "gPopPrev" )
  if( do_pop ) tclvalue( popPrev ) <- calc_popPrev( penAA, penAB, penBB, afreq )
}

getMOI <- function() {
  model <- getPowerVSafe( "gModel", default="additive" )
  if( model=="dominant" )
    return( MOI_DOMINANT )
  if( model=="recessive" )
    return( MOI_RECESSIVE )
  if( model=="additive" )
    return( MOI_ADDITIVE )
  if( model=="multiplicative" )
    return( MOI_MULTIPLICATIVE )
  stop( paste("getMOI(), moi (",model,") misunderstood!",sep="") )
}

fireUpdate_recalc_updatePens <- function( pen ) {
  penAA <- getPower( "gPenAA" )
  penAB <- getPower( "gPenAB" )
  penBB <- getPower( "gPenBB" )

  tclvalue( penAA ) <- as.character( pen[1] )
  tclvalue( penAB ) <- as.character( pen[2] )
  tclvalue( penBB ) <- as.character( pen[3] )
}

fireUpdateAOR_recalc <- function()
{
  moi <- getMOI()
  afreq <- as.numeric( getPowerVSafe( "gAfreqDSL" ) )
  popPrev <- as.numeric( getPowerVSafe( "gPopPrev" ) )
  gAOR <- as.numeric( getPowerVSafe( "gAOR" ) )

  pen <- pen_allelicOR( moi, afreq, popPrev, gAOR )
  fireUpdate_recalc_updatePens( pen )

  fireUpdatePen_recalc(do_AOR=FALSE,do_pop=FALSE)
}

fireUpdateOR_recalc <- function()
{
  moi <- getMOI()
  afreq <- as.numeric( getPowerVSafe( "gAfreqDSL" ) )
  popPrev <- as.numeric( getPowerVSafe( "gPopPrev" ) )
  or1 <- as.numeric( getPowerVSafe( "gOR1" ) )

  pen <- pen_or1( moi, afreq, popPrev, or1 );
  fireUpdate_recalc_updatePens( pen )

  fireUpdatePen_recalc(do_OR=FALSE,do_pop=FALSE)
}

fireUpdateAF_recalc <- function()
{
  moi <- getMOI()
  afreq <- as.numeric( getPowerVSafe( "gAfreqDSL" ) )
  popPrev <- as.numeric( getPowerVSafe( "gPopPrev" ) )
  af <- as.numeric( getPowerVSafe( "gAF" ) )

  pen <- pen_ar( moi, afreq, popPrev, af );
  fireUpdate_recalc_updatePens( pen )

  fireUpdatePen_recalc(do_AF=FALSE,do_pop=FALSE)
}

powerGrid <- function( pseq, power ) {
  powerCalc <- rep(0,length(seq))

  gNumFamilies <- getPower( "gNumFamilies" )

  for( i in 1:length(pseq) ) {
    tclvalue(gNumFamilies) <- pseq[i]
    fireUpdate()
    powerCalc[i] <- as.numeric( getPowerV("gPower") )
    if( powerCalc[i] > power ) {
      if( i==1 ) return( pseq[i] )
      return( mean( pseq[(i-1):i] ) )
    }
    #Sys.sleep(0.01)
    tcl("update")
  }

  return( max(pseq) )
}

fireUpdatePower <- function()
{
  power <- as.numeric( getPowerVSafe( "gPower" ) )
  numSim <- as.numeric( getPowerVSafe( "gNumSim" ) )

  ## okay, first set numSim a little coarser (if possible)
  gNumSim <- getPower("gNumSim")
  tclvalue(gNumSim) <- as.character(numSim/10)

  ## coarse grid search (by 100?)
  pseq <- seq( from=100, to=10000, by=100 )
  grid <- powerGrid( pseq, power )
  tclvalue(gNumSim) <- as.character(numSim/5)
  pseq <- seq( from=grid-50, to=grid+50, by=10 )
  grid <- powerGrid( pseq, power )
  tclvalue(gNumSim) <- as.character(numSim)
  pseq <- seq( from=grid-5, to=grid+5, by=1 )
  grid <- powerGrid( pseq, power )

  ## set the best
  gNumFamilies <- getPower("gNumFamilies")
  tclvalue(gNumFamilies) <- as.character(ceiling(grid))

  ## and make sure didn't fall off the face of the grid
  fireUpdate()
}

###########################
## The exported function ##
###########################
pbat.power <- function(mode="continuous")
{
  ## 04/17/2008 -- bringing in the additive routine
  ##if( powerDisabled() ) stop( "Coming soon!" )

  if( mode!="continuous" && mode!="dichotomous" ) {
    tkmessageBox( title="continuous/binary", message="mode must be 'continuous' or 'binary'", icon="info", type="ok" )
    return(invisible())
  }

  createPowerGUI(mode)
  fireUpdate() ## set some defaults accross the board
}