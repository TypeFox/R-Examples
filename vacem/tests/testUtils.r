#*  ----------------------------------------------------------------------------
#*  Copyright (C) 2011-2012 - Justin Lessler, Jessica Metcalf
#*  
#*  This program is free software; you can redistribute it and/or modify
#*  it under the terms of the GNU General Public License as published by
#*  the Free Software Foundation; version 2 of the License.
#*  
#*  This program is distributed in the hope that it will be useful,
#*  but WITHOUT ANY WARRANTY; without even the implied warranty of
#*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#*  GNU General Public License for more details.
#*  
#*  You should have received a copy of the GNU General Public License
#*  along with this program; if not, write to the Free Software
#*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#*  
#*  $Id: testUtils.r 374 2012-01-26 22:28:34Z ken $
#*  ----------------------------------------------------------------------------


# TBD: Add file header
# TBD: Do I need some 'library( pkg )' or 'source( file )' calls here???

# TBD: Do I need/want this 'cat' (or 'print') here?
verbose <- ( getOption( "RUnit" )$verbose && TRUE )


## TBD: Move these constants to a utility file or somewhere in the
##      main vacem code
DAY <- 1
MONTH <- 30.41667  # 31


## TODO: 1/6/2012
## ====
##    * RUnit 'check' functions use an absolute tolerance.  A better
##      design would be to use ULPs (units of least/last precision).
##      For further rationale, see comments in ISLJ (Java) code that
##      I'm writing for TDI/Hunt Lab.
##      Create new wrapper 'check' functions that use ULPs instead of
##      absolute tolerance, ala the ISLJ unit test code
##
##    * 
##


# TBD: Add doc
# TBD: Move to util package
# TBD: Should we use 'cat' instead of 'paste'?
formatVector <- function ( x, prefix = "[", suffix = "]", ... ) {
  #cat( "mode(x): ", mode(x), ", is(x): ", is(x), "summary(x): ", summary(x), "\n" )

  if ( length(x) > 1 ) {
    x <- paste( x, collapse = ", " )
    # x <- paste( "(", x, ")", sep = "" )
    x <- paste( prefix, x, suffix, sep = "" )
  }

  # x <- as.character( x, ... )
  x <- format( x, ... )
  return ( x )
}

# TBD: Add doc
# TBD: Move to util package
# TBD: Should we use 'cat' instead of 'paste'?
formatDataFrame <- function ( x, prefix = "", suffix = "\n", ... ) {
  #cat( "mode(x): ", mode(x), ", is(x): ", is(x), "summary(x): ", summary(x), "\n" )

  # x <- format.data.frame( x )
  #x <- format( x )
  # cat( "mode(x): ", mode(x), ", is(x): ", is(x), "summary(x): ", summary(x), "\n" )
  #cat( "x[0]: ", x[0], "\n" )
  #cat( "mode(x[1]): ", mode(x[1]), ", is(x[1]): ", is(x[1]), "summary(x[1]): ", summary(x[1]), "\n" )
  #cat( "x[1]: ", unlist(x[1]), "\n" )
  #cat( "x[2]: ", unlist(x[2]), "\n" )
  #cat( "x: ", unlist(x), "\n" )

  buf <- textConnection( ".buffer", open = "w" )
  write.table( x, buf, sep = "  ", eol = "\n", quote = FALSE )
  x <- textConnectionValue( buf )
  x <- paste( prefix, x, suffix, collapse = "", sep="" )
  close( buf )

  #cat( "mode(x): ", mode(x), ", is(x): ", is(x), "summary(x): ", summary(x), "\n" )
  return ( x )
}


# TBD: Add doc
# TBD: Move to util package
# TBD: Should we use 'cat' instead of 'paste'?
toStr <- function ( x, ... ) {
  #cat( "mode(x): ", mode(x), ", is(x): ", is(x), "summary(x): ", summary(x), "\n" )

  if ( is(x,"vector") ) {
    return ( formatVector(x,...) )
  }

  if ( is(x,"data.frame") ) {
    return ( formatDataFrame(x,...) )
  }
  
  # x <- as.character( x, ... )
  x <- format( x, ... )
  return ( x )
}


