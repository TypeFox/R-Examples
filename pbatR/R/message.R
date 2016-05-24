####################################################################
# Thomas Hoffmann                                                  #
# CREATED:     03/27/2006                                          #
# MODIFIED:    03/31/2006                                          #
#                                                                  #
# DESCRIPTION:                                                     #
#  Part of removing tcl/tk (making it optional).                   #
####################################################################

#############################################################
# textMessageBox(...)                                       #
# PARAM msg      message to output to user (w/out choices)  #
#       choices  options for the user                       #
# RETURN  the choice selected by user (the actual text)     #
#############################################################
textMessageBox <- function( msg, choices ) {
  msg <- paste(msg,
               " (Please answer ",
               csPasteVector( choices ),
               ")? ",
               sep="");
  answer <- readline( msg )
  counter <- 0;
  while( !isVecContained(answer, choices) ){
    print( "The answer is case sensative. Please answer the question with one of the choices verbatim." );
    answer <- readline( msg )
    counter <- counter + 1;
  }

  return( answer );
}  ## Debugged

#############################################################
# textMessageBox2(...)                                      #
# PARAM msg      message to output to user (w/out choices)  #
#       choices  options for the user                       #
# RETURN  the choice selected by user (the actual text)     #
#############################################################
textMessageBox2 <- function( msg, choices, numericAnswer=FALSE ) {
  printMsg <- function() {
    cat( msg, '\n', sep='' );
    for( i in 1:length(choices) )
      cat( i, ') ', choices[i], '\n', sep='' );
  }

  printMsg();
  answer <- as.numeric( readline( '? ' ) );
  counter <- 0;
  while( is.na(answer) || answer<=0 || answer>length(choices) ){
    print( "Please answer with the number (e.g. '1')." );
    printMsg();
    answer <- as.numeric( readline( '? ' ) );
    counter <- counter + 1;
  }

  if( numericAnswer )
    return( answer );         ## returns index into the choices
  return( choices[answer] );  ## returns the text of the choice
}  ## Not Debugged

#############################################################
# isPackageLoaded(...)                                      #
# PARAM package  name of package                            #
# RETURN  whether that package has been loaded              #
# WARNING:                                                  #
#  isPackageLoaded("tcltk") is fine, but                    #
#   isPackageLoaded(tcltk) is not so fine                   #
#############################################################
isPackageLoaded <- function( package ){
  if( length( grep( package, search() ) ) > 0 )
    return( TRUE );
  return( FALSE );
}  ## Want to test on other OS'es to be sure

#############################################################
## loadTclTkOrDie()                                         #
## Does just that.                                          #
#############################################################
loadTclTkOrDie <- function(){
  if( !require("tcltk") )
    stop('You need to have loaded tcl/tk to use this function {e.g. library(tcltk)}.');
}


################################
## Debug
##choices <- c("yes", "no", "maybe")
##textMessageBox2( "Are you hungry?", choices )  
##textMessageBox2( "Choose a file or directory", dir() )
