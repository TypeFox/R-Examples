## Thomas Hoffmann
## Created: ??
## Last modified: ??
## Functions to pass arguments in to batch mode.

## newest modification allows us to even do arg1 c(1,1) e.g. to pass in vectors...!
parseCommandArgs <- function(evaluate = TRUE) {
  ## start R off in 'batch' with
  ## R --vanilla --args arg1 arg1value arg2 arg2value ... <rfile.R> rfile.Rout
  ## Set any default values before calling this function.

  arglist <- list() ## output parameter...

  args <- commandArgs();
  i <- which( args=="--args" ); ## gets the args...
  ##if( length(i)==0 || i <= length(args) || length(args)<1 )
  if( length(i)==0 || length(args)<1 )
    return(invisible());

  args <- args[(i+1):length(args)];
  for(i in seq(1, length(args), by=2)){
    value <- NA;
    tryCatch(value <- as.double( args[i+1] ), warning=function(e){})
    if(is.na(value)) {
      value <- args[i+1];
      if( substr(value, 1, 2) == "c(" )
        value <- eval(parse(text=args[i+1]));
    }
    if(evaluate)
      assign( args[i], value, inherits=TRUE );
    #print( "Args" ); print( args[i] ); print( value );

    arglist[[length(arglist) + 1]] <- value
    names(arglist)[length(arglist)] <- args[i]
  }

  return(arglist) ## And return the list...
}

# ## Returns a dataframe of the options that are set, for ease of use
# parseCommandArgsDF <- function() {
#   ## start R off in 'batch' with
#   ## R --vanilla --args arg1 arg1value arg2 arg2value ... <rfile.R> rfile.Rout
#   ## Set any default values before calling this function.
#
#   obj <- data.frame();
#
#   args <- commandArgs();
#   i <- which( args=="--args" );
#   if( length(i)==0 || length(args)<1 )
#     return(obj);
#
#   args <- args[(i+1):length(args)];
#
#   obj <- data.frame( t( as.matrix(  args[ seq( from=2, to=length(args), by=2 ) ]  ) ) );
#   obj.names <- args[ seq( from=1, to=length(args), by=2 ) ];
#   names( obj ) <- obj.names;
#
#   return( obj );
# }



globalSet <- function( arg, value, inherits=TRUE ) {
  assign( arg, value, inherits=inherits );
}

##print( commandArgs() );
##parseCommandArgs();
##print( arg1 ) ## debug only

#print( parseCommandArgsDF() );
