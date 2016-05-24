#' detect 
#'
#' which elements matching the search pattern 
#'
#' @param str seachable target
#' @param pattern method for searching 
#'
#' @return logical; vector indicating which elements match
#'
#' @seealso 
#'   .matches
#'   
#' @note 
#'   - may not export this function
#' @export 

detect  <- function(str, pattern) { 

  type <- if( ! pattern %>% is('SearchableOrPattern') ) "std" else pattern@type
  
  # ITERATE OVER patterns, ANY MATCH WILL DO  
  m <- pattern@.Data %>%
    sapply( 
      function(x) {
        
        stri_detect_regex( str, x )
  
        switch( type
          , regex    = stri_detect_regex( str, x, opts_regex = stri_opts_regex(pattern@options) )      
          , fixed    = stri_detect_fixed( str, x, opts_fixed = stri_opts_fixed(pattern@options) )    
          , coll     = stri_detect_coll(  str, x, opts_coll  = stri_opts_collator(pattern@options) )
          , std      = stri_detect_std(   str, x, opts_std   = stri_opts_std(pattern@options) )
          , stop( "Unknown search type : ", type )
        )
      
      }
    )
   
  apply(m,1,any)

}
