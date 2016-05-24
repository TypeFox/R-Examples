#' @rdname opt_get
#' @export

opt_grab <- function( 
  flag,
  n           = 1,
  # required    = FALSE, 
  opts        = commandArgs()
) {  

   flag.str <- Reduce( function(...) paste(..., sep=", " ), flag )
  
  # EXPAND opts
  opts <- opt_expand(opts=opts)
  
    
  # IDENTIFY name/alias FLAG(s)
  wh.alias <- c() 
  for ( alias in flag ) {
    pattern   <- alias 
    wh.alias  <- union(wh.alias, which( opts==pattern ))
  }

  # FLAG OR ALIAS NOT SUPPLIED  
  if( length(wh.alias) == 0 ) # {
    if( n == 0 ) return(FALSE) else return(NA)  
    
#     # The parameter is required and n!=0.
#     if (required) 
#       stop( 
#         call. = FALSE 
#         , "\n\tOption(s): [", flag.str, "] is required, but was not supplied."
#       )
    
    
#   }
  
 
  
  # MULTIPLE MATCHING FLAGS OR ALIAS FOUND
  # allow.multiple (-tk)
  if ( length(wh.alias) > 1 ) {
    stop( 
      call.= FALSE , 
      "\n\tMultiple values supplied for options [", flag.str, "]" 
    )
  }
  
  # SINGLE FLAG OR ALIAS FOUND.
    
  # CHANGE opt.str TO THE PARTICULAR OPTION FOUND    
  op.str <- opts[wh.alias]  

  if ( n == 0 ) vals <- TRUE 
  
  # N is deterministic
  if( n > 0 ) {
    # TEST FOR AVAILABILITY OF ENOUGH ARGUMENTS
    if( wh.alias+n > length(opts) )
      stop(
        .call = FALSE ,
        "\n\tEnd of arugments reached. [", op.str, "] requires ", n, 
        " arguments but ", max(0, length(opts) - wh.alias), " is/are available."
      )

    # The permissable values are from the value above, taking n values.
    val_rng <- (wh.alias+1):(wh.alias+n)
  
    # TEST: enough values available make sure we don't 
    # encounter any other options

    if( any( is.flag( opts[val_rng] ) ) ) {
      wh <- intersect( which.flag(opts), val_rng )   
      flags <- Reduce( paste, opts[wh] )
      stop( 
        "\n\tUnexpected option flag(s) encountered: ", flag, "\n\t", 
        "[", flag, "] requires ", n, " arguments."
      )
    }  
    
    vals <- opts[val_rng]
    vals <- sub("^\\\\", "", vals) # if value was escaped, remove escape character

  }
  
  return(vals) 
    
}
