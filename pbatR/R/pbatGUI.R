####################################################################
# Thomas Hoffmann                                                  #
# CREATED:     06/21/2005                                          #
# MODIFIED:    07/31/2006                                          #
#                                                                  #
# DESCRIPTION:                                                     #
#  The major pbat-interface commands.                              #
####################################################################


##################################################################
## An attempt to allow more than one interface to run at once,  ##
##  without holding the user. This doesn't work.                ##
##################################################################

pbatGUI.isRunning <- function() {
  # get the globals
  globs <- getPbatGUI( "globs" );

  if( is.null(globs) )
    return(FALSE);

  if( is.null(globs$LOCK) )
    return(FALSE);

  if( globs$LOCK )
    return(TRUE);
  return(FALSE);
}

####################################################################
#                                                                  #
# CONSTANTS                                                        #
#                                                                  #
####################################################################
LISTHEIGHT <- 15; # 5 for debugging, 20 in practice?
LISTWIDTH <- 15;
TEXTWIDTH <- 65;

####################################################################
#                                                                  #
# SET UP THE GLOBAL VARIABLES INTERFACE                            #
#                                                                  #
####################################################################
pbatGUIenv <- new.env();
setPbatGUI <- function( x, value )
  assign( x, value, envir=pbatGUIenv );
getPbatGUI <- function( x, mode="any" )
  get( x, envir=pbatGUIenv, mode=mode, inherits=FALSE );

####################################################################
#                                                                  #
# SET UP THE GLOBAL VARIABLES                                      #
#                                                                  #
####################################################################
pbatGUI.setglobs <- function() {
  loadTclTkOrDie()  ## has to be here to pass the check

  globs <- list();
  globs$form <- 0;
  globs$rbVal.pbat <- tclVar("gee");

  globs$tclVar.pbat <- 0;
  globs$te.pbat <- 0;
  globs$tclVar.pbatwine <- 0;
  globs$te.pbatwine <- 0;

  globs$tclVar.ped <- 0;
  globs$te.ped <- 0;
  globs$tclVar.phe <- 0;
  globs$te.phe <- 0;

  globs$tclVar.group <- 0;
  globs$te.group <- 0;

  globs$phe <- 0; # NULL makes it get lost?
  globs$ped <- 0;
  globs$phefile <- "";
  globs$pedfile <- "";
  # We're getting some strange things - this should fix it - yup.
  globs$pheset <- FALSE;
  globs$pedset <- FALSE;

  globs$phenos <- c();
  globs$preds <- c();
  globs$snps <- c();
  globs$group <- c();
  globs$censor <- "";

  ##globs$allPhenosOrder <- c();
  ##globs$allPhenosMI <- c();
  globs$mi <- c();
  globs$order <- c();

  globs$cb.phenos <- list();
  globs$cbValue.phenos <- list();
  globs$cbInitialized.phenos <- 0;

  globs$rbVal.time <- tclVar();
  globs$rbVal.censor <- tclVar();

  globs$cb.preds <- list();
  globs$cbValue.preds <- list();
  globs$cbInitialized.preds <- 0;
  globs$tclVar.predsOrder <- list();
  globs$te.predsOrder <- list();
  globs$cb.predsInter <- list();
  globs$cbValue.predsInter <- list();

  globs$blocks <- c();

  # all the globs for the options
  globs$max.pheno <- tclVar("1");
  globs$min.pheno <- tclVar("1");
  globs$null <- tclVar("no linkage, no association");
  globs$alpha <- tclVar(0.05);
  globs$trans.pheno <- tclVar("none");
  globs$trans.pred <-  tclVar("none");
  globs$trans.inter <- tclVar("none");
  globs$scan.pred <- tclVar("all");
  globs$scan.inter <- tclVar("all");
  globs$scan.genetic <- tclVar("additive");
  globs$offset <- tclVar("gee");
  globs$screening <- tclVar("conditional power");
  globs$distribution <- tclVar("continuous");
  globs$max.gee <- tclVar("1");
  globs$max.ped <- tclVar("14");
  globs$min.info <- tclVar("0");
  globs$incl.ambhaplos <- tclVar("TRUE");
  globs$infer.mis.snp <- tclVar("FALSE");
  globs$sub.haplos <- tclVar("FALSE");
  globs$length.haplos <- tclVar("2");
  globs$adj.snps <- tclVar("FALSE");
  globs$overall.haplo <- tclVar("FALSE");
  ###globs$cutoff.haplo <- tclVar("FALSE");
  globs$cutoff.haplo <- tclVar("0");
  globs$max.mating.types <- tclVar("10000");
  globs$future.expansion <- tclVar("");
  ## newest (12/29/2006)
  globs$monte <- tclVar("0");
  globs$mminsnps <- tclVar("0");
  globs$mmaxsnps <- tclVar("0");
  globs$mminphenos <- tclVar("0");
  globs$mmaxphenos <- tclVar("0");
  globs$env.cor.adjust <- tclVar("FALSE");
  globs$gwa <- tclVar("FALSE");
  ## 01/19/2006
  globs$snppedfile <- tclVar("FALSE");
  ## 01/26/2006
  globs$extended.pedigree.snp.fix <- tclVar("FALSE");
  ## 03/22/2006
  globs$distribution <- tclVar("default");
  ## 04/17/2008
  #globs$new.ped.algo <- tclVar("TRUE");
  globs$new.ped.algo <- tclVar("FALSE");  ## 08/31/2010

  globs$res <- NULL;

  globs$but.process <- NULL;
  globs$but.plot <- NULL;
  globs$but.results <- NULL;
  globs$but.write <- NULL;

  ### CHECK IF ONE IS RUNNING ###
  globs$LOCK <- TRUE;

  ## Newest addition
  globs$rbVal.modes <- tclVar(pbat.getmode()$mode);
  #globs$rbVal.modes <- tclVar('single');
  globs$te.cluster <- NULL;
  globs$tclVar.cluster <- NULL;
  globs$te.refresh <- NULL;
  globs$tclVar.refresh <- NULL;

  globs$tclVar.loadInput <- NULL;

  globs$env.cor.adjust <- tclVar("FALSE");

  ## 04/28/2008
  globs$cnv.intensity <- tclVar(1) ## tclVar(2)
  globs$cnv.intensity.num <- tclVar(1) ## tclVar(3)

  setPbatGUI( "globs", globs );
}
##pbatGUI.setglobs();  ## this shouldn't be necessary here


####################################################################
#                                                                  #
# CONSTANTS                                                        #
#                                                                  #
####################################################################
ENTRYWIDTH <- 40;
##PBATGUIDEBUGVERBOSE <- TRUE;
PBATGUIDEBUGVERBOSE <- FALSE;
BLOCKENTRYWIDTH <- 100;
pbatGUI.debug <- function( msg )
{
  if( PBATGUIDEBUGVERBOSE==TRUE )
    cat( msg );
}

####################################################################
#                                                                  #
#                                                                  #
#  GUI PBAT INTERFACE FUNCTIONS                                    #
#                                                                  #
#                                                                  #
####################################################################

# Start up the GUI when the user calls pbat!
pbat <- function() {
  pbatGUI.setglobs();
  return( pbatGUI.mainForm() );
}
pbatGUI.errorMessage <- function( message ) {
  loadTclTkOrDie()  ## has to be here to pass the check

  tkmessageBox( title="ERROR",
                message=message,
                icon="error", type="ok" );
}

####################################################################
#                                                                  #
# TK HELPER ROUTINES                                               #
#                                                                  #
####################################################################
pbatGUI.populateList <- function( lst, items ) {
  loadTclTkOrDie()  ## has to be here to pass the check

  if( length(items) > 0 ) { ## semi-bug fix 01/20/2006
    for( i in 1:length(items) )
      tkinsert( lst, "end", items[i] );
  }
}
pbatGUI.tkClearText <- function( obj ) {
  loadTclTkOrDie()  ## has to be here to pass the check

  tkdelete( obj, 0, 9999999 );
}
pbatGUI.tkSetText <- function( obj, text, READONLY=TRUE ) {
  loadTclTkOrDie()  ## has to be here to pass the check

  if( READONLY ) try( tkconfigure( obj, state="normal" ) );  ## try
  pbatGUI.tkClearText( obj );
  tkinsert( obj, "end", text );
  if( READONLY ) try( tkconfigure( obj, state="readonly" ) );  ## try
}


####################################################################
#                                                                  #
# PHENOTYPES...                                                    #
#                                                                  #
####################################################################
pbatGUI.phenotypesForm <- function(){
  loadTclTkOrDie()  ## has to be here to pass the check

  ## get the globals
  globs <- getPbatGUI( "globs" );

  ## create a modal dialog
  form <- tktoplevel();
  tkwm.deiconify(form);
  tkgrab.set(form); # make it modal
  tkfocus(form);
  tkwm.title( form, "P2BAT - Phenotypes" );

  ## get the possible/impossible arrays of phenos
  allPhenos <- c();
  if( class(globs$phe) == 'phe' )
    allPhenos <- names( globs$phe[-c(1,2)] );
  posPhenos <- vectorSubtraction( allPhenos, globs$preds );
  posPhenos <- vectorSubtraction( posPhenos, globs$group );
  posPhenos <- c("AffectionStatus",posPhenos);

  if( length(posPhenos) < 1 ) { ## why was it 2?
    tkmessageBox( title="ERROR",
                  message="Not enough phenotypes ",
                  icon="error", type="ok" );
    tkdestroy(form);
    return();
  }

  ## First get the header onto it
  tkgrid( tklabel( form, text="Select Phenotype(s):" ) );

  ## Draw a list box for the phenos stuff ( no order here... )
  scr.phenos <- tkscrollbar( form, repeatinterval=5,
                          command=function(...)tkyview(lst.phenos,...) );
  lst.phenos <- tklistbox( form, height=LISTHEIGHT, width=LISTWIDTH,
                           selectmode="multiple", background="white",
                           yscrollcommand=function(...)tkset(scr.phenos,...) );
  tkgrid( lst.phenos, scr.phenos );
  tkgrid.configure( scr.phenos, sticky="ns" );

  ## Add in only the possible phenotypes
  pbatGUI.populateList( lst.phenos, posPhenos );

  ## And select some of them which might have been selected
  if( length(globs$phenos) > 0 ) {
    ## find the phenotypes that have been selected
    selPhenos <- c();
    for( i in 1:length(globs$phenos) )
      selPhenos <- c( selPhenos, which(globs$phenos[i]==posPhenos) );

    ## and select them!
    for( i in 1:length(selPhenos) )
      tkselection.set( lst.phenos, selPhenos[i]-1 );
  }

  on.exit <- function() {
    ## Need to translate any selected phenotypes to globs$phenos
    phenosIndex <- as.numeric(tkcurselection(lst.phenos));
    if( length(phenosIndex)<1 ) {
      globs$phenos <- c();
    }else{
      globs$phenos <- posPhenos[phenosIndex+1];

      ## 5/17 if AffectionStatus, turn it the offset to 'none'
      if( globs$phenos[1] == "AffectionStatus" )
        globs$offset <- tclVar("none")
    }

    ## set the global variables
    setPbatGUI( "globs", globs );
  }

  ## Bind mouse presses to on.exit(), since we can't capture closing the form!
  tkbind( lst.phenos, "<ButtonRelease>", on.exit );

  ## lastly, an OK button
  cmdOK <- function() {
    phenosIndex <- as.numeric(tkcurselection(lst.phenos));
    ##print( phenosIndex );
    if( phenosIndex==-999 ) print( "Strange." );

    on.exit();
    tkdestroy(form);
  }
  but.ok <- tkbutton( form, text="   OK   ", command=cmdOK );
  tkgrid(but.ok);

  ## Now make the form go modal
  tkfocus( form );
  tkwait.window( form ); # this makes it go modal
  ##on.exit(); # run after going modal --> doesn't work anymore!
}
pbatGUI.logrankForm <- function(){
  loadTclTkOrDie()  ## has to be here to pass the check

  # get the globals
  globs <- getPbatGUI( "globs" );

  # create a modal dialog
  form <- tktoplevel();
  tkwm.deiconify(form);
  tkgrab.set(form); # make it modal
  tkfocus(form);
  tkwm.title( form, "P2BAT - Time / Censor" );

  # get the possible/impossible arrays of phenos
  allPhenos <- c();
  if( class(globs$phe) == 'phe' )
    allPhenos <- names( globs$phe[-c(1,2)] );
  posPhenos <- vectorSubtraction( allPhenos, globs$preds );
  posPhenos <- vectorSubtraction( posPhenos, globs$group );

  # First get the header onto it
  tkgrid( tklabel( form, text="Time:" ),
          tklabel( form, text="" ),
          tklabel( form, text="Censor:" ) );

  # draw out the phenos stuff
  scr.time <- tkscrollbar( form, repeatinterval=5,
                          command=function(...)tkyview(lst.time,...) );
  lst.time <- tklistbox( form, height=LISTHEIGHT, width=LISTWIDTH,
                         selectmode="single", background="white",
                         yscrollcommand=function(...)tkset(scr.time,...) );
  scr.censor <- tkscrollbar( form, repeatinterval=5,
                             command=function(...)tkyview(lst.censor,...) );
  lst.censor <- tklistbox( form, height=LISTHEIGHT, width=LISTWIDTH,
                           selectmode="single", background="white",
                           yscrollcommand=function(...)tkset(scr.censor,...) );

  # populate the lists with possible phenotypes
  pbatGUI.populateList( lst.time, posPhenos );
  pbatGUI.populateList( lst.censor, posPhenos );

  # and select anything if it's been chosen before
  if( length( globs$phenos )==2 ) {
    tkselection.set( lst.time, which(globs$phenos[1]==posPhenos) );
    tkselection.set( lst.censor, which(globs$phenos[2]==posPhenos) );
  }

  # add a couple of text entry stuff so we can see what was selected...
  tclvar.time <- tclVar("");
  tclvar.censor <- tclVar("");
  if( length(globs$phenos)==2 ) {
    tclvalue(tclvar.time) <- globs$phenos[1];
    tclvalue(tclvar.censor) <- globs$phenos[2];
  }
  te.time <- tkentry( form, width=LISTWIDTH, textvariable=tclvar.time );
  try( tkconfigure( te.time, state="readonly" ) );  ## try
  te.censor <- tkentry( form, width=LISTWIDTH, textvariable=tclvar.censor );
  try( tkconfigure( te.censor, state="readonly" ) );

  # and grid everything
  tkgrid( te.time, tklabel(form,text=""), te.censor );
  tkgrid( lst.time, scr.time, lst.censor, scr.censor );
  tkgrid.configure( scr.time, sticky="ns" );
  tkgrid.configure( scr.censor, sticky="ns" );

  # handle pressing selection doesn't persist! stupid!
  on.time <- function() {
    # get the globals
    globs <- getPbatGUI( "globs" );
    # get and set the time
    strTime <- posPhenos[ as.numeric(tkcurselection(lst.time)) + 1 ];
    if( !is.null(strTime) )
      globs$phenos[1] <- strTime;
    # set the global variables
    setPbatGUI( "globs", globs );
    tclvalue(tclvar.time) <- globs$phenos[1];
  }
  on.censor <- function(){
    # get the globals
    globs <- getPbatGUI( "globs" );
    # and set the censor
    strCensor <- posPhenos[ as.numeric(tkcurselection(lst.censor)) + 1 ];
    if( !is.null(strCensor) )
      globs$phenos[2] <- strCensor;
    # set the global variables
    setPbatGUI( "globs", globs );
    tclvalue(tclvar.censor) <- globs$phenos[2];
  }
  on.exit <- function() {
    # get the globals
    globs <- getPbatGUI( "globs" );

    # check consistency
    if( length(globs$phenos) != 2 ) {
      globs$phenos="";
    }else if( globs$phenos[1] == globs$phenos[2] ) {
      tkmessageBox( title="ERROR",
                    message="Phenotype and time must be different.  Nothing has been set.",
                    icon="error", type="ok" );
      globs$phenos="";
    }else if( globs$phenos[1]=="" || globs$phenos[2]=="" ) {
      tkmessageBox( title="ERROR",
                    message="Phenotype and time must be specified.  Nothing has been set.",
                    icon="error", type="ok" );
      globs$phenos="";
    }

    # set the global variables
    setPbatGUI( "globs", globs );
  }


  # Bind mouse presses to functions since selection won't persist
  tkbind( lst.time, "<ButtonRelease>", on.time );
  tkbind( lst.censor, "<ButtonRelease>", on.censor );

  # lastly, an OK button
  cmdOK <- function() {
    tkdestroy(form);
  }
  but.ok <- tkbutton( form, text="OK", command=cmdOK );
  tkgrid( tklabel(form,text=""),tklabel(form,text=""), but.ok );
  tkgrid.configure( but.ok, sticky="we" );
  tkfocus( form );
  tkwait.window( form ); # this makes it go modal
  on.exit(); # run after going modal -- changed, so okay to run
}

####################################################################
#                                                                  #
# PREDICTORS...                                                    #
#                                                                  #
####################################################################
pbatGUI.predictorsForm <- function(){
  loadTclTkOrDie()  ## has to be here to pass the check

  # get the globals
  globs <- getPbatGUI( "globs" );

  # create a modal dialog
  form <- tktoplevel();
  tkwm.deiconify(form);
  tkgrab.set(form); # make it modal
  #tkfocus(form);
  tkwm.title( form, "P2BAT - Predictors" );

  # get the possible/impossible arrays of phenos
  allPhenos <- c();
  if( class(globs$phe)=="phe" )
    allPhenos <- names( globs$phe[-c(1,2)] );
  posPhenos <- vectorSubtraction( allPhenos, globs$phenos );
  posPhenos <- vectorSubtraction( posPhenos, globs$group );

  ## ensure order has been set
  if( length(globs$order) != length(allPhenos) ){
    globs$order <- rep(1,length(allPhenos));
    globs$mi <- rep(FALSE,length(allPhenos));
    setPbatGUI( "globs", globs );
  }

  # First get the header onto it
  tkgrid( tklabel( form, text="Select Phenotypes:" ),
          tklabel( form, text="" ),
          tklabel( form, text="Modify Selected Phenotypes: " ) );

  ## Create the list of available phenotypes
  scr.aphenos <- tkscrollbar( form, repeatinterval=5,
                              command=function(...)tkyview(lst.aphenos,...) );
  lst.aphenos <- tklistbox( form, height=LISTHEIGHT, selectmode="multiple", background="white",
                            yscrollcommand=function(...)tkset(scr.aphenos,...) );
  pbatGUI.populateList( lst.aphenos, posPhenos );

  ## Create the list for selected phenotypes
  scr.phenos <- tkscrollbar( form, repeatinterval=5,
                             command=function(...)tkyview(lst.phenos,...) );
  lst.phenos <- tklistbox( form, height=LISTHEIGHT, selectmode="single", background="white",
                           yscrollcommand=function(...)tkset(scr.phenos,...) );
  ## Populate the list later
  ## Grid everything so far
  tkgrid( lst.aphenos, scr.aphenos, lst.phenos, scr.phenos );
  tkgrid.configure( scr.phenos, sticky="ns" );
  tkgrid.configure( scr.aphenos, sticky="ns" );
  tkgrid.configure( lst.phenos, sticky="news" );
  tkgrid.configure( lst.aphenos, sticky="news" );

  popList <- function(firsttime=FALSE) {
    ## See if anything is selected (reselect later)
    cursel <- tkcurselection(lst.phenos);

    ## clear the list --> need to take into account order and mi...
    for( i in 1:length(allPhenos) )
      tkdelete( lst.phenos, "end" ); # this might error...

    ## Then populate it
    globs <- getPbatGUI( "globs" );
    if( length(globs$preds) >= 1 ) {
      for( i in 1:length(globs$preds) ) {
        ## Get the full index
        idx <- which( globs$preds[i]==allPhenos );
        str <- allPhenos[idx];
        if( globs$order[idx]!=1 )
          str <- paste(str,"^",globs$order[idx] );
        if( globs$mi[idx] )
          str <- paste("mi(",str,")");

        pbatGUI.populateList( lst.phenos, str );
      }
    }

    ## reselect
    ## - this has a tendency to error - can we try() it?
    if( !firsttime ) {
      ## We can't afford to error on the first pass...
      ##  but nothing should be selected then anyway, so we would error!
      try(
          {
            if( length(cursel) > 0 && as.numeric(cursel)+1<=length(globs$preds) )
              tkselection.set( lst.phenos, cursel )
          },
          silent=TRUE );
    }
  }

  ## Create the 'add' button
  cmdAdd <- function() {
    ## get the globals
    globs <- getPbatGUI( "globs" );

    ## set globs$preds
    phenosChosen <- posPhenos[ as.numeric(tkcurselection(lst.aphenos))+1 ];
    if( length(phenosChosen) < 1 ) return; # nothing was selected
    globs$preds <- sort( unique( c(globs$preds, phenosChosen) ) );
    setPbatGUI( "globs", globs );

    ## finally repopulate the list
    popList();
  }

  ## Create the 'remove' ability
  cmdRemove <- function() {
    # get the globals
    globs <- getPbatGUI( "globs" );

    phenosIndex <- as.numeric(tkcurselection(lst.phenos));
    if( length(phenosIndex)<1 ) return;
    tkdelete( lst.phenos, phenosIndex );
    globs$preds <- globs$preds[-phenosIndex];
    setPbatGUI( "globs", globs );
  }

  ## Grid the buttons
  but.add <- tkbutton( form, text="  Add -->  ", command=cmdAdd );
  but.remove <- tkbutton( form, text="  <-- Delete  ", command=cmdRemove );
  tkgrid( but.add, tklabel(form,text=""), but.remove );

  ## Now marker interaction and order buttons...
  cmdMI <- function() {
    # get the globals
    globs <- getPbatGUI( "globs" );

    selPheno <- globs$preds[as.numeric(tkcurselection(lst.phenos))+1];
    if( length(selPheno)!=1 )
      return;
    idx <- which( selPheno==allPhenos );

    ## 01/20/2006 forbidden pbat command fix
    if( tclvalue(globs$rbVal.pbat)=="logrank" ) {
      ## Error message - not possible
      globs$mi[idx] <- FALSE;
      tkmessageBox( title="Error", message="Marker interaction not available for time-to-onset data.", type="ok" );
    }else{
      globs$mi[idx] <- !globs$mi[idx];
    }
    setPbatGUI( "globs", globs );

    ## repopulate list...
    popList();

    tkfocus(but.mi)  ## So we can process key events...
  }
  but.mi <- tkbutton( form, text="Toggle (M)arker Interaction", command=cmdMI );
  ## Order buttons
  cmdOrder <- function( val=-1, plus=FALSE, minus=FALSE ) {
    # get the globals
    globs <- getPbatGUI( "globs" );

    selPheno <- globs$preds[as.numeric(tkcurselection(lst.phenos))+1];
    if( length(selPheno)!=1 )
      return;
    idx <- which( selPheno==allPhenos );

    # finally set it
    if( val > -1 ) {
      globs$order[idx] <- val;
    }else if(plus){
      globs$order[idx] <- globs$order[idx] + 1;
    }else if(minus){
      globs$order[idx] <- globs$order[idx] - 1;
    }
    if( globs$order[idx] < 1 ) globs$order[idx] <- 1;

    # and set globs
    setPbatGUI( "globs", globs );

    ## repopulate list...
    popList();

    tkfocus(but.mi);
  }
  cmdOrderP <- function() cmdOrder(plus=TRUE);
  cmdOrderM <- function() cmdOrder(minus=TRUE);
  cmdOrder0 <- function() cmdOrder(val=0);
  cmdOrder1 <- function() cmdOrder(val=1);
  cmdOrder2 <- function() cmdOrder(val=2);
  cmdOrder3 <- function() cmdOrder(val=3);
  cmdOrder4 <- function() cmdOrder(val=4);
  cmdOrder5 <- function() cmdOrder(val=5);
  cmdOrder6 <- function() cmdOrder(val=6);
  cmdOrder7 <- function() cmdOrder(val=7);
  cmdOrder8 <- function() cmdOrder(val=8);
  cmdOrder9 <- function() cmdOrder(val=9);

  ## can we get a keypress? - have to stick it to a command button
  tkbind( but.mi, "p", cmdOrderP );
  tkbind( but.mi, "s", cmdOrderM );
  tkbind( but.mi, "0", cmdOrder0 );
  tkbind( but.mi, "1", cmdOrder1 );
  tkbind( but.mi, "2", cmdOrder2 );
  tkbind( but.mi, "3", cmdOrder3 );
  tkbind( but.mi, "4", cmdOrder4 );
  tkbind( but.mi, "5", cmdOrder5 );
  tkbind( but.mi, "6", cmdOrder6 );
  tkbind( but.mi, "7", cmdOrder7 );
  tkbind( but.mi, "8", cmdOrder8 );
  tkbind( but.mi, "9", cmdOrder9 );

  tkbind( but.mi, "m", cmdMI );

  but.orderP <- tkbutton( form, text="+ Order (p)", command=cmdOrderP );
  but.orderM <- tkbutton( form, text="- Order (s)", command=cmdOrderM );

  tkbind( lst.phenos, "<ButtonRelease>", function(){tkfocus(but.mi)} );

  ## Don't forget to populate the list for the first time!
  popList(firsttime=TRUE);

  ## lastly, an OK button
  cmdOK <- function() {
    tkdestroy(form);
  }

  but.ok <- tkbutton( form, text="   OK   ", command=cmdOK );

  tkgrid( tklabel(form,text=""),tklabel(form,text=""), but.mi );
  tkgrid( tklabel(form,text=""),tklabel(form,text=""), but.orderP );
  tkgrid( but.ok, tklabel(form,text=""), but.orderM );

  tkfocus( form );

  tkwait.window( form ); # this makes it go modal
}

####################################################################
#                                                                  #
# SNPS / BLOCKS ...                                                #
#                                                                  #
####################################################################
pbatGUI.snpsForm <- function() {
  loadTclTkOrDie()  ## has to be here to pass the check

  ## get the globals
  globs <- getPbatGUI( "globs" );

  # create a modal dialog
  form <- tktoplevel();
  tkwm.deiconify(form);
  tkgrab.set(form); # make it modal
  tkfocus(form);
  tkwm.title( form, "P2BAT - SNPs / Blocks" );

  # get the names of possible snps
  allSnps <- names( as.pedlist( globs$ped ) ); # horribly inefficient
  allSnps <- allSnps[7:length(allSnps)];

  # Grid the headers
  tkgrid( tklabel( form, text="Choose SNPs:" ),
          tklabel( form, text="" ),
          tklabel( form, text="Blocks selected" ) );
  # Create the list of the snps
  scr.snp <- tkscrollbar( form, repeatinterval=5,
                          command=function(...)tkyview(lst.snp,...) );
  lst.snp <- tklistbox( form, height=LISTHEIGHT, selectmode="multiple", background="white",
                        yscrollcommand=function(...)tkset(scr.snp,...) );
  ##lst.snp <- tklistbox( form, height=10, selectmode="multiple", background="white" );
  pbatGUI.populateList( lst.snp, allSnps );
  ####tkselection.set(globs$lst.pbat,0); # indexes from zero!
  # Create the list for the block list
  scr.block <- tkscrollbar( form, repeatinterval=5,
                            command=function(...)tkyview(lst.block,...) );
  lst.block <- tklistbox( form, height=LISTHEIGHT, width=BLOCKENTRYWIDTH, selectmode="single",
                          background="white",
                          yscrollcommand=function(...)tkset(scr.block,...) );
  ##lst.block <- tklistbox( form, height=LISTHEIGHT, width=BLOCKENTRYWIDTH, selectmode="single", background="white" );
  ##print( 'populating list with the following' ); ##print( globs$blocks );
  pbatGUI.populateList( lst.block, globs$blocks );

  # Create the 'add' button
  cmdAddBlock <- function() {
    # get the globals
    globs <- getPbatGUI( "globs" );

    snpsChosen <- allSnps[ as.numeric(tkcurselection(lst.snp)) + 1 ];
    if( nchar(snpsChosen[1]) < 1 ) return; ##############
    #print( snpsChosen );  # DEBUG only
    tkselection.clear( lst.snp, 0, 'end' );

    #newEntry <- pasteVector( snpsChosen, SQUOTE=FALSE, COMMASEP=TRUE );
    newEntry <- pasteVector2( snpsChosen, sep=" + " );
    if( sum(newEntry==globs$blocks) > 0 ) {
      ##print( "Item is already in the list!" );
    }else {
      globs$blocks <- c( globs$blocks, newEntry );
      tkinsert( lst.block, "end", newEntry );

      # set the global variables
      setPbatGUI( "globs", globs );
    }
  }
  cmdAddSnp <- function() {
    # get the globals
    globs <- getPbatGUI( "globs" );

    snpsChosen <- allSnps[ as.numeric(tkcurselection(lst.snp)) + 1 ];
    if( length(snpsChosen<1) || nchar(snpsChosen)<1 ) return;
    tkselection.clear( lst.snp, 0, 'end' );

    for( i in 1:length(snpsChosen) ) {
      newEntry <- snpsChosen[i]; #pasteVector( snpsChosen, SQUOTE=FALSE, COMMASEP=TRUE );
      if( sum(newEntry==globs$blocks) > 0 ) {
        ##print( "Item is already in the list!" );
      }else {
        globs$blocks <- c( globs$blocks, newEntry );
        tkinsert( lst.block, "end", newEntry );

        ;# set the global variables
        setPbatGUI( "globs", globs );
      }
    }
  }
  cmdRemove <- function() {
    # get the globals
    globs <- getPbatGUI( "globs" );

    snpIndex <- as.numeric(tkcurselection(lst.block));
    if( length(snpIndex)<1 ) return;

    tkdelete( lst.block, snpIndex );
    globs$blocks <- globs$blocks[-(snpIndex+1)];

    ;# set the global variables
    setPbatGUI( "globs", globs );
  }
  but.addBlock <- tkbutton( form, text= "  Add Block -->  ", command=cmdAddBlock );
  but.addSnp <- tkbutton( form, text  = "  Add SNPs / CNVs -->  ", command=cmdAddSnp );
  but.remove <- tkbutton( form, text  = "  <-- Delete entry  ", command=cmdRemove );

  if( is.cped(globs$ped) ) tkconfigure( but.addBlock, state="disabled" ); ## who knows what this does...

  # and put everything on the grid
  tkgrid( lst.snp, scr.snp, lst.block, scr.block );
  tkgrid.configure( scr.snp, sticky="ns" );
  tkgrid.configure( scr.block, sticky="ns" );
  tkgrid( but.addBlock, tklabel(form,text=""), but.remove );
#  tkgrid( but.addSnp );


  # lastly, an OK button
  cmdOK <- function() {
    # Store anything on the form that we need!

    # set the global variables
    #setPbatGUI( "globs", globs ); # NO - DON't set them!!!

    # and lastly kill the form
    tkdestroy(form);
  }
  but.ok <- tkbutton( form, text="   OK   ", command=cmdOK );
  tkgrid( but.addSnp, tklabel(form,text=""), but.ok );
  tkfocus( form );
  tkwait.window( form ); # this makes it go modal
}

####################################################################
#                                                                  #
# GROUP...                                                         #
#                                                                  #
####################################################################
pbatGUI.groupForm <- function() {
  loadTclTkOrDie()  ## has to be here to pass the check

  # get the globals
  globs <- getPbatGUI( "globs" );

  # create a modal dialog
  form <- tktoplevel();
  tkwm.deiconify(form);
  tkgrab.set(form); # make it modal
  tkfocus(form);
  tkwm.title( form, "P2BAT - Group" );

  # get the possible/impossible arrays of phenos
  allPhenos <- names( globs$phe[-c(1,2)] );
  posPhenos <- vectorSubtraction( allPhenos, globs$preds );
  posPhenos <- vectorSubtraction( posPhenos, globs$phenos );

  if( length(posPhenos) < 1 ) {
    tkmessageBox( title="ERROR",
                  message="Not enough available phenotypes to select a group.",
                  icon="error", type="ok" );
    tkdestroy(form);
    return();
  }

  # First get the header onto it
  tkgrid( tklabel( form, text="Select Phenotype(s):" ) );

  # Draw a list box for the phenos stuff ( no order here... )
  scr.phenos <- tkscrollbar( form, repeatinterval=5,
                          command=function(...)tkyview(lst.phenos,...) );
  lst.phenos <- tklistbox( form, height=LISTHEIGHT, width=LISTWIDTH,
                           selectmode="single", background="white",
                           yscrollcommand=function(...)tkset(scr.phenos,...) );
  tkgrid( lst.phenos, scr.phenos );
  tkgrid.configure( scr.phenos, sticky="ns" );

  # Add in only the possible phenotypes
  pbatGUI.populateList( lst.phenos, posPhenos );

  # And select some of them which might have been selected
  if( length(globs$group)==1 && globs$group!="" ) {
    selPhenos <- which(globs$group==posPhenos);
    tkselection.set( lst.phenos, selPhenos-1 );
  }

  on.exit <- function() {
    ## Need to translate any selected phenotypes to globs$phenos
    phenosIndex <- as.numeric(tkcurselection(lst.phenos));
    if( length(phenosIndex)<1 ) {
      globs$group <- "";
    }else{
      globs$group <- posPhenos[phenosIndex+1];
    }
    ## set the global variables
    setPbatGUI( "globs", globs );
    ## And set the text on the main form
    pbatGUI.tkSetText( globs$te.group, globs$group );
  }

  # Bind mouse presses to on.exit(), since we can't capture closing the form!
  tkbind( lst.phenos, "<ButtonRelease>", on.exit );

  # draw out the button for NONE
  but.none <- tkbutton( form, text="NONE",
                       command=function(){
                         globs <- getPbatGUI('globs');
                         globs$group <- "";
                         setPbatGUI('globs',globs);
                         pbatGUI.tkSetText( globs$te.group, globs$group );
                         tkdestroy(form);
                       } );
  tkgrid( but.none );
  tkgrid.configure( but.none, sticky="we" );

  # lastly, an OK button
  cmdOK <- function() {
    phenosIndex <- as.numeric(tkcurselection(lst.phenos));
    ##print( phenosIndex );
    if( phenosIndex==-999 ) print( "Strange." );

    on.exit();
    tkdestroy(form);
  }
  but.ok <- tkbutton( form, text="OK", command=cmdOK );
  tkgrid(but.ok);
  tkgrid.configure( but.ok, sticky="we" );

  # Now make the form go modal
  tkfocus( form );
  tkwait.window( form ); # this makes it go modal
  ##on.exit(); # run after going modal --> doesn't work anymore!

}

####################################################################
#                                                                  #
# OPTIONS...                                                       #
#                                                                  #
####################################################################
pbatGUI.optionsForm <- function( whichForm=0 ) {
  ## whichForm - 0 means them all, otherwise now we're cutting it into 2 pieces

  loadTclTkOrDie()  ## has to be here to pass the check

  # get the globals
  globs <- getPbatGUI( "globs" );

  ######isP <- FALSE;
  isG <- FALSE;
  isL <- FALSE;
  ######if( tclvalue(globs$rbVal.pbat) == "pc" ) isP <- TRUE;
  if( tclvalue(globs$rbVal.pbat) == "gee" ) isG <- TRUE;
  if( tclvalue(globs$rbVal.pbat) == "logrank" ) isL <- TRUE;

  # create a modal dialog
  form <- tktoplevel();
  tkwm.deiconify(form);
  tkgrab.set(form); # make it modal
  tkfocus(form);
  tkwm.title( form, "P2BAT - Options" );

  # helper gui functions
  msg <- function( message, title="Information" ) {
    tkmessageBox( title="",
                 message=message,
                 type="ok" );
  }
  newframe <- function( gridframe=form, grid=TRUE, relief="groove", borderwidth=2, sticky="nws" ) {
    frame <- tkframe( gridframe, relief=relief, borderwidth=borderwidth );
    if( !grid ) return(frame);
    tkgrid( frame );
    tkgrid.configure( frame, sticky=sticky );
    return(frame)
  };
  newOpt <- function( tclVar, option, options=c("FALSE","TRUE"), helps=NULL, gridframe=form ) {
    if( length(options)<2 ) return( "ERROR" );
    frame <- newframe( borderwidth=3, gridframe=gridframe );
    but <- NULL;
    if( !is.null(helps) ) {
      but <- tkbutton( frame, text="?",
                       command=function(){
                         msgText <- paste( option, "-", helps[1] );
                         if( length(helps)>1 ) {
                           for( i in 2:length(helps) )
                             msgText <- paste( msgText, "\n\n", options[i-1], ": ", helps[i], sep="" )
                         }
                         msg(msgText)
                       } );
    }

    lab <- tklabel( frame, text=paste(option, ":  ", sep="" ) );

    subframe <- list();
    for( i in 1:length(options) )
      subframe[[i]] <- newframe(gridframe=frame,grid=FALSE)

    if( length(options)==2 ) {
      if( is.null(but) ) {
        tkgrid( lab, subframe[[1]], subframe[[2]] );
      }else{
        tkgrid( lab, subframe[[1]], subframe[[2]], but );
      }
    }else if( length(options)==3 ) {
      if( is.null(but) ) {
        tkgrid( lab, subframe[[1]], subframe[[2]], subframe[[3]] );
      }else{
        tkgrid( lab, subframe[[1]], subframe[[2]], subframe[[3]], but );
      }
    }else if( length(options)==4 ) {
      if( is.null(but) ) {
        tkgrid( lab, subframe[[1]], subframe[[2]], subframe[[3]], subframe[[4]] );
      }else{
        tkgrid( lab, subframe[[1]], subframe[[2]], subframe[[3]], subframe[[4]], but );
      }
    }else if( length(options)==5 ) {
      if( is.null(but) ) {
        tkgrid( lab, subframe[[1]], subframe[[2]], subframe[[3]], subframe[[4]], subframe[[5]] );
      }else{
        tkgrid( lab, subframe[[1]], subframe[[2]], subframe[[3]], subframe[[4]], subframe[[5]], but );
      }
    }else {
      if( is.null(but) ) {
        tkgrid(lab);
      }else{
        tkgrid(lab,but);
      }

      for( i in 1:length(options) )
        tkgrid(subframe[[i]]);
    }

    for( i in 1:length(options) ) {
      rb <- tkradiobutton( subframe[[i]] );
      tkconfigure( rb, variable=tclVar, value=options[i] );
      tkgrid( rb, tklabel(subframe[[i]], text=options[i]) );
    }
  }
  newTE <- function( tclVar, option, helps=NULL, width=10, gridframe=form ) {
    frame <- newframe( borderwidth=3, gridframe=gridframe );
    but <- NULL;
    if( !is.null(helps) ) {
      but <- tkbutton( frame, text="?",
                       command=function(){
                         msgText <- paste( option, "-", helps[1] );
                         msg(msgText)
                       } );
    }
    lab <- tklabel( frame, text=option );
    entry <- tkentry( frame, width=width, textvariable=tclVar );
    if( !is.null(but) ) {
      tkgrid( lab, entry, but );
    }else{
      tkgrid( lab, entry );
    }
  }

  dblGrid <- function() {
    f <- tkframe(form, borderwidth=0 );
    f1 <- tkframe( f, borderwidth=0 );
    f2 <- tkframe( f, borderwidth=0 );
    fblank <- tkframe( f, borderwidth=0 );
    tkgrid( f );
    tkgrid.configure( f, sticky="nws" );
    tkgrid( f1, fblank, f2 );
    tkgrid( tklabel(fblank,text="    ") );
    return( list(f1=f1,f2=f2) );
  }

  #junk <- tclVar();
  #newOpt( junk, "Infer Missing SNPs",
  #        helps=c(
  #          "Handling of missing genotype information in the haplotypes analysis",
  #          "Individuals with missing genotype information are excluded from the analysis.  This analysis is also implemented in the HBAT option of the FBAT program.",
  #          ""
  #         ) );
  #junk2 <- tclVar();
  #newTE( junk2, "Max Matring Types", helps="Yes, wouldn't you like help?" )


  if( whichForm==0 || whichForm==1 ){
    ;##############################
    ;# draw in all of the options #
    ;##############################

    if( !isL ) {
      dg <- dblGrid();

      newTE( globs$max.pheno, "Max Phenotypes",
            helps="The maximum number of phenotypes that will be analyzed in the FBAT-statistic.",
            gridframe=dg$f1 );

      newTE( globs$min.pheno, "Min Phenotypes",
            helps="The minimum number of phenotypes that will be analyzed in the FBAT-statistic.",
            gridframe=dg$f2);
    }

    {
      dg <- dblGrid();

      newOpt( globs$null, "Null Hypothesis",
             options=c("no linkage, no association", "linkage, no association"),
             helps=c("Specification of the null-hypothesis.",
               "Null-hypothesis of no linkage and no association.",
               "Null-hypothesis of linkage, but no association." ),
             gridframe=dg$f1 );

      newTE( globs$alpha, "Alpha",
            helps="Specification of the significance level.",
            gridframe=dg$f2 );
    }
    newOpt( globs$trans.pheno, "Phenotype Transformation",
           options=c("none","ranks","normal score"),
           helps=c("Transformation of the selected phenotypes.",
             "no transformation (default)",
             "transformation to ranks",
             "transformation to normal score (recommended for quantitative phenotypes)" ) );

    newOpt( globs$trans.pred, "Predictor Transformation",
           options=c("none","ranks","normal score"),
           helps=c("Transformation of the selected predictor variables/covariates.",
             "no transformation (default)",
             "transformation to ranks",
             "transformation to normal score (recommended for quantitative phenotypes)" ) );

    newOpt( globs$trans.inter, "Interaction Transformation",
           options=c("none","ranks","normal score"),
           helps=c("Transformation of the selected interaction variables",
             "no transformation (default)",
             "transformation to ranks",
             "transformation to normal score (recommended for quantitative phenotypes)" ) );

    if( !isL ) {
      dg <- dblGrid();

      newOpt( globs$scan.pred, "Covariate Model",
             options=c("all","subsets"),
             helps=c("Computation of all covariate sub-models",
               "The selected FBAT statistic is computed with adjustment for all selected covariates/predictors.",
               "The selected FBAT statistic is computed for all posible subsets of the selected covariates/predictor variables.  The command is particularly useful to examine the dependence of significant results on the selection of a covariate model." ),
             gridframe=dg$f1 );
      newOpt( globs$scan.inter, "Interaction Model",
             options=c("all","subsets"),
             helps=c("Computation of all interaction sub-models",
               "The selected FBAT statistic is computed including all selected interaction variables.",
               "The selected FBAT statistic is computed for all posible subsets of the interaction variables."),
             gridframe=dg$f2);
    }

    newOpt( globs$scan.genetic, "Inheritance Mode",
           options=c("additive","dominant","recessive","heterozygous advantage", "all"),
           helps=c("Specification of the mode of inheritance",
             "Additive model",
             "Dominant model",
             "Recessive model",
             "Heterozygous advantage model",
             "The FBAT-statistics are computed for all 4 genetic models") );

    newOpt( globs$offset, "Covariate Offset",
           options=c("none","max power","gee + marker score","gee"),
           helps=c("Specification of the covariate/predictor variables adjustment",
             "No adjustments for covariates/predictor variables",
             "Offset (=FBAT adjustment for covariates and interaction  variables) that maximizes the power of the FBAT-statistic (computationally slow, efficiency dependent on the correct choice of the mode of inheritance)",
             "FBAT adjustment for covariates and interaction variables) based on standard phenotypic residuals obtained by GEE-estimation including the expected marker score (E(X|H0)), all covariates and interaction variables.",
             "Offset (=FBAT adjustment for covariates and interaction variables) based on standard phenotypic residuals obtained by GEE-estimation including all covariates and interaction variables. The default choice is 'gee' ('no' for dichotomous traits).") );

    {
      ##dg <- dblGrid();

      newOpt( globs$screening, "Screening",
             options=c("conditional power","wald"),
             helps=c("Specification of the screening methods to handle the multiple comparison problem for multiple SNPs/haplotypes and a set of phenotypes.",
               "Screening based on conditional power (parametric approach)",
               "Screening based on Wald-tests (non-parametric approach)") );
      ##       gridframe=dg$f1 );

      #newOpt( globs$distribution, "Phenotype Distribution",
      #       options=c("continuous","categorical"),
      #       helps=c("Specification of the phenotypic distribution",
      #         "Phenotypes are treated as continuous phenotypes in the power calculation",
      #         "Phenotypes are treated as categorical/integer variables. This option is especially recommended for analysis of time-to-onset data and affection status."),
      #       gridframe=dg$f2 );
      newOpt( globs$distribution, "Empirical Pheno Distn",
              options=c("default","jiang","murphy","naive","observed"),
              helps=c(
                "Specification of the phenotype distribution for screening.",
                "Default",
                "Approach by Jiang et. al 2006",
                "Approach by Murphy et. al 2006",
                "Naive allele frequency estimators",
                "Observed allele frequencies") );
    }
  }

  if( whichForm==0 || whichForm==2 ) {
    if( isG )
      newTE( globs$max.gee, "Max GEE Iterations",
            helps="Specification of the maximal number of iterations in the GEE-estimation procedure." );
    {
      dg <- dblGrid();

      newTE( globs$max.ped, "Max Pedigree Iterations",
            helps="Specification of the maximal number of proband in one extended pedigrees.",
            gridframe=dg$f1 );

      newTE( globs$min.info, "Min Families",
            "Specification of the minimum number of informative families required for the computation of the FBAT-statistics.",
            gridframe=dg$f2 );
    }

    {
      dg <- dblGrid();

      newOpt( globs$incl.ambhaplos, "Ambiguous Haplotypes",
             helps=c("This command defines the handling of ambiguous haplotypes in the haplotypes analysis.",
               "Ambiguous haplotypes (phase can not be inferred) are included in the analysis and are weighted according to their estimated frequencies in the probands.",
               "Ambiguous haplotypes are excluded from the analysis."),
             gridframe=dg$f1 );

      newOpt( globs$infer.mis.snp, "Include Missing",
             helps=c("Handling of missing genotype information in the haplotypes analysis.",
               "Individuals with missing genotype information are excluded from the analysis. This is the analysis also implemented in the HBAT option of the FBAT-program.",
               "Individuals with missing genotype information are included in the analysis. The algorithm of Horvath et al (2004) is applied to all individuals, even if they have missing genotype information. This results in more ambiguous haplotypes."),
             gridframe=dg$f2 );
    }

    {
      dg <- dblGrid();

      newOpt( globs$sub.haplos, "Sub-Haplotypes",
             helps=c("",
               "The haplotypes defined by the all SNPs given in the haplotype-block definition are analyzed.",
               "All haplotypes are analyzed that are defined by any subset of SNPs in the haplotypes block definition."),
             gridframe=dg$f1 );

      newTE( globs$length.haplos, "Haplotype length",
            helps="Defines the haplotype length when subhaplos=TRUE",
            gridframe=dg$f2 );
    }

    {
      dg <- dblGrid();

      newOpt( globs$adj.snps, "Adjacent Haplotypes",
             helps=c("Takes effect when subhaplos=TRUE.",
               "All sub-haplotypes are analyzed",
               "Only the sub-haplotypes are analyzed for which the first constituting SNPs are adjacent." ),
             gridframe=dg$f1 );

      newOpt( globs$overall.haplo, "Overall Haplotypes",
             helps=c("Specification of an overall haplotypes test. When this command is included in the batch-file, only one level of the 'groups' variable can be specified.",
               "no overall test",
               "an overall test is computed testing all haplotypes defined by the same set of SNPs simultaneously. This option can not be applied when sub.haplos=TRUE"),
             gridframe=dg$f2);
    }

    {
      dg <- dblGrid();
      newTE( globs$cutoff.haplo, "Min Haplotype Freq.",
            helps="The minimum haplotypes frequency so that a haplotypes is included in the overall test.", gridframe=dg$f1 );

      newTE( globs$max.mating.types, "Max Mating Types",
            helps="Maximal number of mating types in the haplotype analysis.", gridframe=dg$f2 );
    }

    newTE( globs$future.expansion, "(Future Expansion)", width=40,
          helps="(Only included for future expansion of pbat.) Lines to write to the batchfile for pbat." );

    ;##############################
    ;# draw in all of the options #
    ;##############################

    ## new options 12/29/2007
    newTE( globs$monte, "Monte Carlo Iterations",
           helps="When this is nonzero, monte-carlo based methods are used to compute the p-values instead, according to the number of iterations supplied. 1000 iterations is suggested." );
    {
      dg <- dblGrid();
      newTE( globs$mminsnps, "MM SNP min",
             helps="Multi-marker multi-phenotype tests: the minimum number of snps to be tested.",
             gridframe=dg$f1 );
      newTE( globs$mmaxsnps, "MM SNP max",
             helps="Multi-marker multi-phenotype tests: the maximum number of snps to be tested.",
             gridframe=dg$f2 );

    }

    {
      dg <- dblGrid();
      newTE( globs$mminphenos, "MM phenotype min",
             helps="Multi-marker multi-phenotype tests: the minimum number of phenotypes to be tested.",
             gridframe=dg$f1 );
      newTE( globs$mmaxphenos, "MM phenotype max",
             helps="Multi-marker multi-phenotype tests: the maximum number of phenotypes to be tested.",
             gridframe=dg$f2 );

    }

    {
      dg <- dblGrid();
      newOpt( globs$env.cor.adjust, "Environment Correlation Adjust",
              helps=c("Environment Correlation",
                "Do not adjust for",
                "Adjust for"),
              gridframe=dg$f1 );
      newOpt( globs$gwa, "Genome-Wide Accelerated mode",
              helps=c("Whether to use (g)enome (w)ide (a)cceleration mode.  This is faster for genome-wide association tests, and has slightly less output.",
                "Don''t use",
                "Use"),
              gridframe=dg$f2 );
    }

    {
      dg <- dblGrid();
      newOpt( globs$snppedfile, "snppedfile",
              helps=c("snppedfile","The pedigree file does not contain just snps.","The pedigree file does just contain snps. This is advantageous to specify as the storage mode is much more compact, and the program will use much less memory."),
              gridframe=dg$f1 );
      newOpt( globs$extended.pedigree.snp.fix, "Extended pedigree snps fix",
              helps=c("Set to TRUE when you have more extended pedigrees in your dataset, as the pedigree reconstruction will be more accurate. Note this mode is only compatible with 'single' mode, so be sure to set that as well.",
                "Haplotype accelerated mode - faster and good when this is not the case.",
                "Slow, but good."),
              gridframe=dg$f2 );
    }
    {
      newOpt( globs$new.ped.algo, "New pedigree algorithm",
              helps=c("Set to TRUE to use the new 10-100 times faster, more memory efficient algorithm for extended pedigrees.",
                "Do not use the new method.",
                "Use the new method!") );
    }
  }




  # lastly, an OK button
  cmdOK <- function() {
    tkdestroy(form);
  }
  but.ok <- tkbutton( form, text="   Close   ", command=cmdOK );
  tkgrid(but.ok);
  tkfocus( form );
  tkwait.window( form ); # this makes it go modal
  on.exit(); # run after going modal
}

####################################################################
#                                                                  #
# MAIN FORM FUNCTIONS                                              #
#                                                                  #
####################################################################
pbatGUI.pbatset <- function() {
  loadTclTkOrDie()  ## has to be here to pass the check

  # get the globals
  globs <- getPbatGUI( "globs" );
  pbat.set();  # pop's up the GUI selection tool
  pbatGUI.tkSetText( globs$te.pbat, pbat.get() );
}

pbatGUI.pbatwineset <- function() {
  loadTclTkOrDie()
  globs <- getPbatGUI( "globs" );
  pbat.setwine();  # pop's up the GUI selection tool
  pbatGUI.tkSetText( globs$te.pbatwine, pbat.getwine() );
}

pbatGUI.pedFileChoice <- function() {
  loadTclTkOrDie()  ## has to be here to pass the check

  # get the globals
  globs <- getPbatGUI( "globs" );

  tkdelete( globs$te.ped, 0, 999999 );

  # See if we can get the filename; return if cancelled
  tempstr <- tclvalue(tkgetOpenFile(filetypes="{{Un/compressed cped/ped File} {.ped .pped .cped}} {{Pedigree File} {.ped}} {{Compressed Pedigree File} {.pped}} {{Copy Number Variant Pedigree File} {.cped}}",title="Pedigree file - NO SPACES in path"));

  ## 5/17 - make sure there is no space in the filename
  if( spaceInFilename( tempstr ) ) {
    tkmessageBox( title="ERROR - space in path",
                  message=spaceInFilenameError(tempstr) );
    return();
  }

  if( !nchar(tempstr) ) return();
  globs$pedfile <- tempstr;
  pbatGUI.tkSetText( globs$te.ped, tempstr );

  # Load in the data file
  tempstrExtension <- file.extension(tempstr)
  print(tempstrExtension)
  if( tempstrExtension=="ped" ) {
    globs$ped <- read.ped( globs$pedfile );
  }else if( tempstrExtension=="pped" ){
    globs$ped <- read.pped( globs$pedfile );
  }else{
    ## cped! argh!
    globs$ped <- read.cped( globs$pedfile );
    ##print( globs$ped ) ## DEBUG ONLY

    ## New -- can we also set the offset to be zero?
    #globs$offset <- tclVar("none")
    #wrong -- this should only be for AffectionStatus
  }
  globs$pedset <- TRUE;

  # Now, also set the phefile
  phefile <- paste( substring(globs$pedfile,1,nchar(globs$pedfile)-3), "phe", sep="" );
  if( file.exists(phefile) ) {
    # File exists!  Assume this is probably what the user wants...
    globs$phefile <- phefile;
    pbatGUI.tkSetText( globs$te.phe, phefile );

    # Load in the phefile
    globs$phe <- read.phe( globs$phefile );
    globs$pheset <- TRUE;
  }else{
    ## Try again
    phefile <- paste( substring(globs$pedfile,1,nchar(globs$pedfile)-4), "phe", sep="" );
    ## (and just a copy above the code now that we are trying again...)
    if( file.exists(phefile) ) {
      ## File exists!  Assume this is probably what the user wants...
      globs$phefile <- phefile;
      pbatGUI.tkSetText( globs$te.phe, phefile );

      ## Load in the phefile
      globs$phe <- read.phe( globs$phefile );
      globs$pheset <- TRUE;
    }
  }

  # Set the globals
  setPbatGUI( "globs", globs );

  return();
}
pbatGUI.pheFileChoice <- function() {
  loadTclTkOrDie()  ## has to be here to pass the check

  # get the globals
  globs <- getPbatGUI( "globs" );

  # See if we can get the filename, return if cancelled
  tempstr <- tclvalue(tkgetOpenFile(filetypes="{{Phenotype File} {.phe}}",title="Phenotype file - NO SPACES in path"));
  if( !nchar(tempstr) ) return();

  ## 5/17 - make sure there is no space in the filename
  if( spaceInFilename( tempstr ) ) {
    tkmessageBox( title="ERROR - space in path",
                  message=spaceInFilenameError(tempstr) );
    return();
  }

  globs$phefile <- tempstr;
  pbatGUI.tkSetText( globs$te.phe, tempstr );

  # Load in the data file
  globs$phe <- read.phe( globs$phefile );  ## phe, not phefile!!!!
  globs$pheset <- TRUE;

  # Set the globals
  setPbatGUI( "globs", globs );
  return();
}

pbatGUI.ensureDataLoaded <- function() {
  loadTclTkOrDie()  ## has to be here to pass the check

  # get the globals
  globs <- getPbatGUI( "globs" );

  # make sure 'phe' and 'ped' aren't null
  ##if( !globs$pheset | !globs$pedset ) {
  ##  tkmessageBox( title="ERROR",
  ##                message="You need to load in the pedigree and phenotype files first.",
  ##                icon="error", type="ok" );
  ##  return( FALSE );
  ##}
  if( !globs$pedset ) {
    tkmessageBox( title="ERROR",
                  message="You need to at least load in a pedigree file first, and also a phenotype file unless you plan on just using AffectionStatus.",
                  icon="error", type="ok" );
    return( FALSE );
  }
  return( TRUE );
}

pbatGUI.phenotypes <- function() {
  pbatGUI.debug("phenotypes\n");
  if( !pbatGUI.ensureDataLoaded() ) return( FALSE );

  if(  tclvalue( getPbatGUI("globs")$rbVal.pbat )  !=  "logrank"   ) {
    pbatGUI.phenotypesForm();
  }else {
    pbatGUI.logrankForm();
  }
}

pbatGUI.predictors <- function() {
  pbatGUI.debug("predictors\n");
  if( !pbatGUI.ensureDataLoaded() ) return( FALSE );

  pbatGUI.predictorsForm();
}

pbatGUI.snps <- function() {
  pbatGUI.debug("snps / blocks\n" );
  if( !pbatGUI.ensureDataLoaded() ) return( FALSE );

  pbatGUI.snpsForm();
}

pbatGUI.group <- function() {
  pbatGUI.debug("group\n");
  if( !pbatGUI.ensureDataLoaded() ) return( FALSE );

  pbatGUI.groupForm();
}

pbatGUI.options <- function(whichForm=0) {
  #pbatGUI.debug("options\n");
  #if( !pbatGUI.ensureDataLoaded() ) return( FALSE );

  pbatGUI.optionsForm(whichForm=whichForm);
}

pbatGUI.options1 <- function()
  pbatGUI.options(1);
pbatGUI.options2 <- function()
  pbatGUI.options(2);


pbatGUI.write <- function() {
  loadTclTkOrDie()  ## has to be here to pass the check

  # get the globals
  globs <- getPbatGUI( "globs" );

  if( !is.null(globs$res) ) {
    # get a filename to write out
    outfile <- tclvalue(tkgetSaveFile(filetypes="{{Comma Separated Values} {.csv}} {{Text File} {.txt}}"));
    if( !is.null(outfile) && outfile!="" ) {
      if( strfindf(outfile,'.')==-1 )
        outfile <- paste( outfile, ".csv", sep="" );
      write.pbat( globs$res, outfile );
      return(TRUE);
    }
  }
  return(FALSE); # failed to write
}

## addition -- compression
pbatGUI.compress <- function() {
  globs <- getPbatGUI( "globs" );

  if( is.cped(globs$ped) ) {
    print( "CNV files cannot be compressed." )
    return();
  }

  ## copy from onProcess() 02/02/2007 -- rewrite 05/24/06 for modes
  numProcesses <- tclvalue( globs$tclVar.pbatNP );
  mode <- tclvalue( globs$rbVal.modes );
  if( mode=='single' ) numProcesses <- 1;
  if( numProcesses==1 ) mode <- 'single';
  cluster <- tclvalue( globs$tclVar.cluster );
  refresh <- tclvalue( globs$tclVar.refresh );
  pbat.setmode( mode=mode, jobs=numProcesses, clusterCommand=cluster, clusterRefresh=refresh );

  ## and then do the work
  if( is.pped(globs$ped) ){
    print( "File already is compressed." );
    return();
  }
  globs$ped <- as.pped( globs$ped );
  pbatGUI.tkSetText( globs$te.ped, get.sym( globs$ped ) );
  setPbatGUI( "globs", globs );
}

# draw the main form, wait until everything is all done
pbatGUI.mainForm <- function() {
  loadTclTkOrDie()  ## has to be here, period.

  ## get the globals
  globs <- getPbatGUI( "globs" );

  ## Create the window
  globs$form <- tktoplevel();
  tkwm.deiconify( globs$form );
  tkfocus( globs$form );
  tkwm.title( globs$form, "P2BAT" );

  ## Create all of the buttons and objects on the form

  {
    ## Frame 1
    frame.prelim <- tkframe( globs$form, relief="groove", borderwidth=2 );
    tkgrid( frame.prelim );
    tkgrid.configure( frame.prelim, sticky="news" );

    ## Pbat executable
    but.pbat <- tkbutton( frame.prelim, text="Pbat exe...", command=pbatGUI.pbatset );
    globs$tclVar.pbat <- tclVar( pbat.get() );
    globs$te.pbat <- tkentry( frame.prelim, width=ENTRYWIDTH, textvariable=globs$tclVar.pbat );
    try( tkconfigure( globs$te.pbat, state="readonly" ) ); ## try
    tkgrid( but.pbat, globs$te.pbat );
    tkgrid.configure( but.pbat, sticky="ew" );

    if( !isWindows() ) {
      but.pbatwine <- tkbutton( frame.prelim, text="Wine exe...", command=pbatGUI.pbatwineset );
      globs$tclVar.pbatwine <- tclVar( pbat.getwine() );
      globs$te.pbatwine <- tkentry( frame.prelim, width=ENTRYWIDTH, textvariable=globs$tclVar.pbatwine );
      try( tkconfigure( globs$te.pbatwine, state="readonly" ) ); ## try
      tkgrid( but.pbatwine, globs$te.pbatwine );
      tkgrid.configure( but.pbatwine, sticky="ew" );
    }

    ## - pedigree file line
    but.ped <- tkbutton( frame.prelim, text="Pedigree File ...", command=pbatGUI.pedFileChoice );
    globs$tclVar.ped <- tclVar();
    globs$te.ped <- tkentry( frame.prelim, width=ENTRYWIDTH, textvariable=globs$tclVar.ped );
    ## -- addi for compression
    but.compress <- tkbutton( frame.prelim, text="Compress", command=pbatGUI.compress );
    ## -- idda
    ##tkgrid( but.ped, globs$te.ped );
    tkgrid( but.ped, globs$te.ped, but.compress );
    try( tkconfigure( globs$te.ped, state="readonly" ) ); ## try
    tkgrid.configure( but.ped, sticky="ew" );

    ## - phenotype file line
    but.phe <- tkbutton( frame.prelim, text="Phenotype File ...", command=pbatGUI.pheFileChoice );
    globs$tclVar.phe <- tclVar();
    globs$te.phe <- tkentry( frame.prelim, width=ENTRYWIDTH, textvariable=globs$tclVar.phe );
    try( tkconfigure( globs$te.phe, state="readonly" ) ); ## try
    tkgrid( but.phe, globs$te.phe );
  }

  {
    # Frame 1.5 -- cnv.intensity and cnv.intensity.num choices...
    frame.cnv <- tkframe( globs$form, relief="groove", borderwidth=2 );
    tkgrid( frame.cnv )
    tkgrid.configure( frame.cnv, sticky="nws" )
    lbl1 <- tklabel(frame.cnv,text=" CNV only options: intensity number to analyze" )
    lbl2 <- tklabel(frame.cnv,text="# intensities" )
    te.intensity <- tkentry( frame.cnv, width=3, textvariable=globs$cnv.intensity )
    te.intensity.num <- tkentry( frame.cnv, width=5, textvariable=globs$cnv.intensity.num )
    tkgrid( lbl1, te.intensity, lbl2, te.intensity.num )
  }

  {
    # Frame 2 - pbat choice
    ;# - pbat choice
    frame.pbatchoice <- tkframe( globs$form, relief="groove", borderwidth=2 );
    tkgrid( frame.pbatchoice );
    tkgrid.configure( frame.pbatchoice, sticky="nws" );
    pbatchoices <- c("gee","pc","logrank");
    rb.pbat <- list();
    rb.lbl <- list();
    rb.subframe <- list();

    ## Create the subframe for each choice
    for( i in 1:length(pbatchoices) )
      rb.subframe[[i]] <- tkframe( frame.pbatchoice, relief='groove', borderwidth=1 );
    lbl <- tklabel(frame.pbatchoice,text='PBAT:');
    tkgrid( lbl, rb.subframe[[1]], rb.subframe[[2]], rb.subframe[[3]] );

    ## Create each choice
    for( i in 1:length(pbatchoices) ) {
      rb.pbat[[i]] <- tkradiobutton( rb.subframe[[i]] );
      tkconfigure( rb.pbat[[i]], variable=globs$rbVal.pbat, value=pbatchoices[i]  );
      rb.lbl[[i]] <- tklabel( rb.subframe[[i]], text=pbatchoices[i] );
      tkgrid( rb.pbat[[i]], rb.lbl[[i]] );
    }
    ####tkgrid.configure( rb.lbl[[i]], sticky="w" );
    ####tkgrid.configure( rb.pbat[[i]], sticky="w" );
  }

  {
    ## Frame 3 - all other options
    frame.misc <- tkframe( globs$form, relief="groove", borderwidth=2 );
    tkgrid( frame.misc );
    tkgrid.configure( frame.misc, sticky="news" );

    ## - phenotypes
    but.phenotypes <- tkbutton( frame.misc, text="Phenotypes / Censor ... ", command=pbatGUI.phenotypes );
    tkgrid( but.phenotypes );
    tkgrid.configure( frame.misc, sticky="we" );

    ## - predictors
    but.predictors <- tkbutton( frame.misc, text="Predictors ...", command=pbatGUI.predictors );
    tkgrid( but.predictors );
    tkgrid.configure( but.predictors, sticky="we" );

    ## - snps / blocks
    but.snps <- tkbutton( frame.misc, text="SNPs / Blocks / CNVs ...", command=pbatGUI.snps );
    tkgrid( but.snps );
    tkgrid.configure( but.snps, sticky="we" );

    ## - group
    but.group <- tkbutton( frame.misc, text="Group ...", command=pbatGUI.group );
    globs$tclVar.group <- tclVar();
    globs$te.group <- tkentry( frame.misc, width=ENTRYWIDTH, textvariable=globs$tclVar.group );
    try( tkconfigure( globs$te.group, state="readonly" ) ); ## try
    tkgrid( but.group, globs$te.group );
    tkgrid.configure( but.group, sticky="we" );

    ## - options
    ## -- old, all on one form
    #but.options <- tkbutton( frame.misc, text="Options ...", command=pbatGUI.options );
    #tkgrid( but.options );
    #tkgrid.configure( but.options, sticky="we" );

    ## -- two options forms
    but.options1 <- tkbutton( frame.misc, text="Options (1) ...", command=pbatGUI.options1 );
    tkgrid( but.options1 );
    tkgrid.configure( but.options1, sticky="we" );

    but.options2 <- tkbutton( frame.misc, text="Options (2) ...", command=pbatGUI.options2 );
    tkgrid( but.options2 );
    tkgrid.configure( but.options2, sticky="we" );
  }

  {
    ## Frame 3.5a
    frame.mode <- tkframe( globs$form, relief="groove", borderwidth=2 );
    tkgrid( frame.mode );
    tkgrid.configure( frame.mode, sticky='nws' );
    modes <- c('single','multi','cluster');
    rb.modes <- list();
    rb.lbl <- list();
    rb.subframe <- list();

    ## create the subgrids for each choice
    for( i in 1:length(modes) )
      rb.subframe[[i]] <- tkframe( frame.mode, relief='groove', borderwidth=1 );
    lbl <- tklabel(frame.mode,text='Multiprocessor Mode: ');
    tkgrid( lbl, rb.subframe[[1]], rb.subframe[[2]], rb.subframe[[3]] );

    ## create each choice
    for( i in 1:length(modes) ){
      rb.modes[[i]] <- tkradiobutton( rb.subframe[[i]] );
      tkconfigure( rb.modes[[i]], variable=globs$rbVal.modes, value=modes[i] );
      rb.lbl[[i]] <- tklabel( rb.subframe[[i]], text=modes[i] );
      tkgrid( rb.modes[[i]], rb.lbl[[i]] );
    }

  }

  {
    ## Frame 3.5b - number of _jobs_
    frame.np <- tkframe( globs$form, relief="groove", borderwidth=2 );
    tkgrid( frame.np );
    tkgrid.configure( frame.np, sticky="news" );

    globs$tclVar.pbatNP <- tclVar( pbat.getNumProcesses() );
    globs$te.pbatNP <- tkentry( frame.np, width=ENTRYWIDTH, textvariable=globs$tclVar.pbatNP );
    ##tkconfigure( globs$te.pbatNP, state="readonly" );
    lbl.pbatNP <- tklabel( frame.np, text="Number of jobs:" );
    tkgrid( lbl.pbatNP, globs$te.pbatNP );
    tkgrid.configure( lbl.pbatNP, sticky="ew" );

    ## now add the cluster command and refresh
    {
      #frame.np2 <- tkframe( frame.np, relief='groove', borderwidth=2 );
      #tkgrid( frame.np2 );
      #tkgrid.configure( frame.np2, sticky="news" );

      globs$tclVar.cluster <- tclVar( pbat.getmode()$cluster );
      globs$te.cluster <- tkentry( frame.np, width=ENTRYWIDTH, textvariable=globs$tclVar.cluster );
      lbl.cluster <- tklabel( frame.np, text="Cluster Command:" );

      globs$tclVar.refresh <- tclVar( pbat.getmode()$refresh );
      globs$te.refresh <- tkentry( frame.np, width=ENTRYWIDTH, textvariable=globs$tclVar.refresh );
      lbl.refresh <- tklabel( frame.np, text="Cluster Refresh:" );

      tkgrid( lbl.cluster, globs$te.cluster );
      tkgrid( lbl.refresh, globs$te.refresh );
      tkgrid.configure( lbl.cluster, sticky='ew' );
      tkgrid.configure( lbl.refresh, sticky='ew' );
      tkgrid.configure( globs$te.cluster, sticky='w' );
      tkgrid.configure( globs$te.refresh, sticky='w' );
    }
  }

  {
    ## Frame 3.6 - optional loading the input in
    frame.li <- tkframe( globs$form, relief="groove", borderwidth=2 );
    tkgrid( frame.li );
    tkgrid.configure( frame.li, sticky='news' );

    globs$tclVar.loadInput <- tclVar('1');
    checkbut <- tkcheckbutton(frame.li);
    tkconfigure( checkbut, variable=globs$tclVar.loadInput );
    lbl <- tklabel( frame.li, text='Load PBAT results into R object' );
    tkgrid( checkbut, lbl );
  }

  {
    ## Frame 4 - processing
    frame.process <- tkframe( globs$form, relief="groove", borderwidth=3 );
    tkgrid( frame.process );
    tkgrid.configure( frame.process, sticky="news" );

    ;####################################################################
    ;#  onProcess(...)                                                  #
    ;#                                                                  #
    ;####################################################################
    onProcess <- function() {
      ## Clear the results (01/25/2006)
      pbat.last.clear();

      ## Set the number of processes (01/08/2006):
      globs <- getPbatGUI( "globs" );
      ##if( 0 != pbat.setNumProcesses( tclvalue(globs$tclVar.pbatNP) ) ) {
      ##  tkmessageBox( title="ERROR",
      ##                message="Number of processes must be a positive integer.",
      ##                icon="error", type="ok" );
      ##  return(FALSE);
      ##}

      ## -- rewrite 05/24/06 for modes
      numProcesses <- tclvalue( globs$tclVar.pbatNP );
      mode <- tclvalue( globs$rbVal.modes );
      if( mode=='single' ) numProcesses <- 1;
      if( numProcesses==1 ) mode <- 'single';
      cluster <- tclvalue( globs$tclVar.cluster );
      refresh <- tclvalue( globs$tclVar.refresh );
      pbat.setmode( mode=mode, jobs=numProcesses, clusterCommand=cluster, clusterRefresh=refresh );
      ## -- and an additon for repressing loading the output in
      LOAD.OUTPUT <- tclvalue(globs$tclVar.loadInput)==1;
      ##print( LOAD.OUTPUT );
      ##return();  ## DEBUG ONLY

      ##

      # ensure the data was actually loaded in!
      if( !pbatGUI.ensureDataLoaded() ) return( FALSE );

      globs <- getPbatGUI( "globs" );
      allPhenos <- c();
      if( class(globs$phe)=="phe" )
        allPhenos <- names( globs$phe[-c(1,2)] );
      allPhenos <- c( "AffectionStatus", allPhenos );

      ;# First put the R command together...
      ;# We'll be calling pbat.m(...)

      # First fill in the phenotypes...
      formula <- "";
      if(  tclvalue( getPbatGUI("globs")$rbVal.pbat )  !=  "logrank"   ) {
        ;# pbat pc or pbat gee
        if( length(globs$phenos) == 0 ) {
          formula <- "ALL";
        }else{
          formula <- pasteVector2( globs$phenos, sep="+" );
        }
      }else{
        ;# pbat logrank
        if( length(globs$phenos) != 2 ) {
          pbatGUI.errorMessage("Time & Censor must be set for pbat logrank." );
          return();
        }
        formula <- pasteVector2( globs$phenos, sep=" & " );
      }

      # next, we have a tilde '~' in it
      formula <- paste( formula, "~ " );

      # now, put each of the prediction vars into the formula
      if( length(globs$preds)>=1 ) {
        ## MI INTERACTION MI INTERACTION
        #print( "globs$preds" )
        #print( globs$preds )
        #print( "globs$mi" )
        #print( globs$mi )
        mis <- NULL
        print( "allPhenos" )
        print( allPhenos )
        try( {
          if( any( globs$mi ) )
            mis <- allPhenos[c(FALSE,globs$mi)] ## it's AffectionStatus, and then everything else...
        } )
        print( mis )
        for( i in 1:length(globs$preds) ) {
          if( i>1 )
            formula <- paste( formula, " + ", sep="" );

          ## $allPhenosOrder, $allPhenosMI
          phenosIndex <- which(globs$preds[i]==allPhenos);
          ###order <- as.numeric(tclvalue(globs$tclVar.predsOrder[[phenosIndex]]));
          order <- globs$order[phenosIndex];
          ###if( as.numeric(tclvalue(globs$cbValue.predsInter[[phenosIndex]])) ) {

          print( "globs$mi" )
          print( globs$mi )

          #if( globs$mi[phenosIndex] ) {
          if( any(globs$preds[i]==mis) ) { ## is it in the interactions???
            if( order>1 ) {
              formula <- paste( formula, "mi(",globs$preds[i],"^",order,")", sep="" );
            }else{
              formula <- paste( formula, "mi(",globs$preds[i],")", sep="" );
            }
          }else {
            formula <- paste( formula, globs$preds[i], sep="" );
            if( order>1 ) formula <- paste( formula, "^", order, sep="" );
          }

          print( formula )
        }
      }else{
        formula <- paste( formula, "NONE" );
      }

      # and finally follow with the blocks
      if( length(globs$blocks) > 0 ) {
        for( i in 1:length(globs$blocks) )
          formula <- paste( formula, "|", globs$blocks[i] );
      }

      # lastly, if there was a group
      if( length(globs$group)==1 && globs$group!="" )
        formula <- paste( formula, "/", globs$group );

      ##print( paste("Formula:", formula ) );
      # FORMULA HAS BEEN CREATED (now all the other options!)

      ##print( tclvalue(globs$rbVal.pbat) );

      ## NEWEST! - wd set
      cur <- pbat.work(globs$ped)

      #################
      # CALL PBAT!!!! #
      #################
      zero.to.null <- function( x ){ if(x<=0) return(NULL); return(x); };
      globsphe <- NULL; ## fix for AffectionStatus in GUI
      try( {
        globsphe <- globs$phe;
        if( is.na(globsphe) || is.numeric(globsphe) || nchar(globsphe)==0 )
          globsphe <- NULL;
      } )

      ## OK, NEW tryCatch - this would be awesome, awesome, so awesome...
      PBATMERRORED <- FALSE;
      tryCatch( {
      globs$res <- pbat.m(
                          formula,
                          globsphe, globs$ped,
                          fbat=tclvalue(globs$rbVal.pbat),
                          max.pheno=tclvalue(globs$max.pheno),
                          min.pheno=tclvalue(globs$min.pheno),
                          null=tclvalue(globs$null),
                          alpha=tclvalue(globs$alpha),
                          trans.pheno=tclvalue(globs$trans.pheno),
                          trans.pred=tclvalue(globs$trans.pred),
                          trans.inter=tclvalue(globs$trans.inter),
                          scan.pred=tclvalue(globs$scan.pred),
                          scan.inter=tclvalue(globs$scan.inter),
                          scan.genetic=tclvalue(globs$scan.genetic),
                          offset=tclvalue(globs$offset),
                          screening=tclvalue(globs$screening),
                          distribution=tclvalue(globs$distribution),
                          max.gee=tclvalue(globs$max.gee),
                          max.ped=tclvalue(globs$max.ped),
                          min.info=tclvalue(globs$min.info),
                          incl.ambhaplos=tclvalue(globs$incl.ambhaplos),
                          infer.mis.snp=tclvalue(globs$infer.mis.snp),
                          sub.haplos=tclvalue(globs$sub.haplos),
                          length.haplos=tclvalue(globs$length.haplos),
                          adj.snps=tclvalue(globs$adj.snps),
                          overall.haplo=tclvalue(globs$overall.haplo),
                          cutoff.haplo=tclvalue(globs$cutoff.haplo),
                          max.mating.types=tclvalue(globs$max.mating.types),
                          future.expansion=tclvalue(globs$future.expansion),
                          LOAD.OUTPUT=LOAD.OUTPUT,
                          monte=tclvalue(globs$monte),
                          mminsnps=zero.to.null(tclvalue(globs$mminsnps)),
                          mmaxsnps=zero.to.null(tclvalue(globs$mmaxsnps)),
                          mminphenos=zero.to.null(tclvalue(globs$mminphenos)),
                          mmaxphenos=zero.to.null(tclvalue(globs$mmaxphenos)),
                          env.cor.adjust=tclvalue(globs$env.cor.adjust),
                          gwa=tclvalue(globs$gwa),
                          snppedfile=tclvalue(globs$snppedfile),
                          extended.pedigree.snp.fix=tclvalue(globs$extended.pedigree.snp.fix),
                          new.ped.algo=tclvalue(globs$new.ped.algo),
                          cnv.intensity=tclvalue(globs$cnv.intensity),
                          cnv.intensity.num=tclvalue(globs$cnv.intensity.num)
                          );
      }, error=function(e){tkinsert( globs$statusText, "end", paste("GUI ERROR:",e$message,"\n",sep="") ); PBATMERRORED<-TRUE; if(PBATMERRORED){} } ); #,
      ##}, error=function(e){tkinsert( globs$statusText, "end", paste("GUI ERROR:",e$message,"\n",sep="") ); PBATMERRORED<-TRUE } ); #,
      #warning=function(e2){print( paste("GUI Warning:",e2$message,"\n",sep="")) } );
      #warning=function(e){tkinsert( globs$statusText, "end", paste("GUI Warning:",e$message,"\n",sep="") ) } );
      ## Warning above overrides the error! Terrible! Left alone...

      ## Newest! -- unwork
      pbat.unwork(cur)

      if( !PBATMERRORED ) {
        ## -- and go ahead and print the status
        st <- pbat.status(workFirst=TRUE)
        for( i in 1:length(st) )
          tkinsert( globs$statusText, "end", paste(st[i],"\n",sep="") )

        # save any results
        setPbatGUI( "globs", globs );

        # enable the other routines
        tkconfigure( globs$but.results, state="normal" );
        tkconfigure( globs$but.write, state="normal" );
        if(  tclvalue( getPbatGUI("globs")$rbVal.pbat )  ==  "logrank"   )
          tkconfigure( globs$but.plot, state="normal" );

        ## And kill it if refreshing is zero
        if( mode=='cluster' && refresh==0 )
          tkdestroy( globs$form );
      }
    }

    onPlot <- function() {
      res <- pbat.last();
      if( !is.null(res) )
        plot( res );
    }

    onResults <- function() {
      if( !is.null(pbat.last()$results) ) {
        print( pbat.last()$results );
      }else{
        warning( "Results file could not be read in properly; printing output of file. Type 'pbat.last()$results.logfile' to get the name of this file." );
        pbat.last.rawResults(); ## Not ideal, but at least it gives them something...
      }
    }

    onWrite <- function() {
      pbatGUI.write();
    }

    globs$but.process <- tkbutton( frame.process, text="PROCESS", command=onProcess );
    globs$but.plot <- tkbutton( frame.process, text="Plot (logrank)", command=onPlot );
    globs$but.results <- tkbutton( frame.process, text="Results", command=onResults );
    globs$but.write <- tkbutton( frame.process, text="Write", command=onWrite );
    tkgrid( globs$but.process, globs$but.plot, globs$but.results, globs$but.write );
    tkconfigure( globs$but.plot, state="disabled" );
    tkconfigure( globs$but.results, state="disabled" );
    tkconfigure( globs$but.write, state="disabled" );
  }

  ## New: draw a little status window
  {
    stFrame <- tkframe( globs$form, relief="groove", borderwidth=2 )
    tkgrid( stFrame )
    #globs$statusText <- tktext( stFrame, bg="black", fg="white" )
    #tkgrid( globs$statusText )
    #try( tkconfigure( globs$statusText, state="disabled" ) )

    yscr <- tkscrollbar( stFrame, repeatinterval=5,
                         command=function(...)tkyview(globs$statusText,...))
    globs$statusText <- tktext( stFrame, bg="black", fg="white", height=5, width=TEXTWIDTH,
                                yscrollcommand=function(...)tkset(yscr,...) )
    tkgrid(globs$statusText, yscr)
    tkgrid.configure(yscr, sticky="nse")
    tkgrid.configure(globs$statusText, sticky="nse" )

    tkinsert( globs$statusText, "end", "Status of run (last line, enter 'pbat.status(n=0)' from the command prompt for more):\n" )
  }


  ## set the globals before pausing for all the gui operations to finish
  setPbatGUI( "globs", globs );
  ## No longer halt in case you want to do post-processing
  ## Nope - don't know how to catch closing it with the (X)
  tkwait.window( globs$form );

  #return(invisible());
  return( pbat.last() );
}

####################################################################
# pbat.last()                                                      #
# Returns the results of the test run by pbat()                    #
####################################################################
pbat.last <- function() {
  # get the globals
  globs <- getPbatGUI( "globs" );
  # and return the structure
  return( globs$res );
}

####################################################################
## pbat.last.clear()                                               #
## Clears out the pbat.last() result.                              #
####################################################################
pbat.last.clear <- function() {
  globs <- getPbatGUI( "globs" );
  globs$res <- NULL;
  setPbatGUI( "globs", globs );
}

####################################################################
# pbat.last.rawResults()                                           #
# Prints out the direct output of the files... useful especially   #
#  when the program fails to read in Christoph's data file         #
#  properly.                                                       #
####################################################################
pbat.last.rawResults <- function() {
  basefname <- pbat.last()$results.logfile;

  fname <- paste(basefname,"",sep="");
  if( file.exists(fname) )
    printFile(fname);

  fname <- paste(basefname,".hdr",sep="");
  if( file.exists(fname) )
    printFile(fname);

  fname <- paste(basefname,".dat",sep="");
  if( file.exists(fname) )
    printFile(fname);
}

####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
####################################################################

