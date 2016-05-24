xmlblock <- function( tag, xmllines, capitalize = TRUE ) {
# extracts from a loaded XML file the info of a tag
  if( capitalize ) tag <- toupper( tag ) # tags in capital
  tag_ini   <- paste( "<", tag, ">", sep = "" )
  tag_end   <- paste( "</", tag, ">", sep = "" )
  tag_empty <- paste( "<", tag, "/>", sep = "" )

# is empty?
  if( length( grep( tag_empty, xmllines ) ) > 0 ) return( NULL )
# was it a wrong tag?
  if( length( grep( tag_ini, xmllines ) ) == 0 ) return( NA ) # tag does not exist
  
  line_ini <- grep( tag_ini, xmllines )[1]
  line_end <- grep( tag_end, xmllines )[1]
# does it have only one line? if so, return NA
  if( line_ini == line_end ) return( Inf )

# otherwise it is all fine, then return block
  line_ini = line_ini + 1
  line_end = line_end - 1

  return( xmllines[line_ini:line_end] )
}