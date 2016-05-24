xmlitem <- function( tag, xmllines, capitalize = TRUE ) {
# extracts from a loaded XML file the info of a tag
  if( capitalize ) tag <- toupper( tag ) # tags in capital
  tag_ini   <- paste( "<", tag, ">", sep = "" )
  tag_end   <- paste( "</", tag, ">", sep = "" )
  tag_empty <- paste( "<", tag, "/>", sep = "" )

# is empty?
  if( length( grep( tag_empty, xmllines ) ) > 0 ) return( NA )
# get content
  xmlitem   <- sub( tag_end, "", sub( tag_ini, "", xmllines[grep( tag_ini, xmllines )[1]] ) )
# was it a wrong tag?
  if( length( xmlitem ) == 0 ) return( NA ) # tag does not exist

  return( xmlitem )
}