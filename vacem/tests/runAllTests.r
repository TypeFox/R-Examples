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
#*  $Id: runAllTests.r 387 2012-01-27 17:35:48Z ken $
#*  ----------------------------------------------------------------------------


# TBD: Add header
# TBD: Add doc


# TBD: Add doc
# TBD: Move to a utilities file
# TBD: How 'paste' should work
# TBD: Create unit tests for this function
join <- function ( ..., sep=", ", prefix="", suffix="" ) {
  listing <- paste( ..., collapse=sep )
  listing <- paste( prefix, listing, suffix, sep="" )
  return ( listing )
}


# TBD: Add doc
# TBD: Move to a utilities file
# TBD: Is there a base function that performs this function?
# TBD: Create unit tests for this function
isEmpty <- function ( x ) {

  # cat( sprintf( "Checking if argument x ('%s') is empty, i.e. NA, NULL, \"\", c(), etc\n", x ) )

  try( if ( is.null(x) ) return ( TRUE ) , silent=TRUE )
  try( if ( is.na(x)   ) return ( TRUE ) , silent=TRUE )

  # if ( ! any(x)   ) return ( TRUE )
  try( if ( is.logical(x)     && ! any(x)      ) return ( TRUE ) , silent=TRUE )

  try( if ( is.character(x)   &&  nchar(x) < 1 ) return ( TRUE ) , silent=TRUE )  # empty string
  try( if ( is.character(x)   && length(x) < 1 ) return ( TRUE ) , silent=TRUE )  # null string
  try( if ( is.vector(x)      && length(x) < 1 ) return ( TRUE ) , silent=TRUE )  # empty vector
  try( if ( is.list(x)        && length(x) < 1 ) return ( TRUE ) , silent=TRUE )  # empty list
  try( if ( is.data.frame(x)  && length(x) < 1 ) return ( TRUE ) , silent=TRUE )  # empty data.frame

  return ( FALSE )

}  # end of function isEmpty( x )


# TBD: Add doc
# TBD: Move to a utilities file
# TBD: Create unit tests for this function
# TBD: Test with subs = { NA, NULL, "", c( "dir1", "dir2" ), etc
#
getPkgPaths <- function ( pkg, subs="" ) {
  # TBD: Should we use 'getwd()' instead of '.' and 'getwd()/..' instead of '..'?
  defaultRoots <- c( ".", ".." )

  # cat( sprintf( "Getting path to subdirectories '%s' for package '%s'\n",
  #               join(subs,sep="', '"), pkg ) )

  if ( isEmpty(pkg) ) {
    stop( "No package name specified" )
  }

  root <- system.file( package = pkg )
  # cat( sprintf( "Package '%s' is located in '%s'\n", pkg, root ) )

  if ( isEmpty(root) ) {
    # cat( sprintf( paste( "No root directory found for package '%s'",
    #               "-- Perhaps code is not being run from within a package.",
    #               "\nSearching default root locations ('%s') for the",
    #               "subdirectories: '%s'\n"),
    #               pkg, join(defaultRoots,sep="', '"), join(subs,sep="', '") ) )

    # TBD: We should be testing for _more_ than just _exists_.  We should check
    #      if it (a) exists, (b) is a directory and (c) is accessible
    # TBD: Is there a 'select' function in R that is like 'which' except instead
    #      of returning an index, it returns the matched items?
    #      It would be nice if 'select' (and 'which') supported location flags
    #      like 'first' and last so we could something like:
    #         select( paths, file.exists, first=TRUE )
    paths <- file.path( defaultRoots, subs )
    root <- defaultRoots[ which(file.exists(paths)) ][ 1 ]
    #cat( "   Subs:  '", join(subs,sep="', '"), "'\n",
    #     "   Paths: '", join(paths,sep="', '"), "'\n",
    #     "   Root: '", root, "'\n",
    #     "   file.exists(paths): ", file.exists(paths), "\n",
    #     "   which(file.exists(paths)): ", which(file.exists(paths)), "\n",
    #     "   defaultRoots[ which(file.exists(paths)) ][1]: ", defaultRoots[ which(file.exists(paths)) ][1], "\n"
    #   )

    if ( isEmpty(root) ) {
      cat( sprintf( paste( "Paths not found for package '%s'",
                  " -- Unable to find root for subdirectories: '%s';",
                  " returning NULL\n" ),
                  pkg, join(subs,sep="', '") ) )
      return ( NULL )
    }

  }

  paths <- file.path( root, subs )

  # FIX: Rewrite this code to handle a vector/list of paths
  #      E.g. print warnings about all elements that (a) do not exist
  #      or (b) are not directories or (c) are not accessible.
  #      Return all valid paths.  Remove all "invalid" paths.
  #      As long as one path is valid, just print warning and return
  #      If none are valid then return null or throw an exception?
  if ( ! file.exists(paths) ) {
    cat( sprintf( paste( "Paths not found for package '%s'",
                  " -- Subdirectories '%s' do not exist;",
                  " returning NULL\n"),
                  pkg, join(paths,sep="', '") ) )
    return ( NULL )
  }

  # cat( sprintf( "Paths for package '%s': '%s'\n", pkg, join(paths,sep="', '") ) )
  return ( paths )

}  # end of function getPkgPaths(pkg,subs=NULL)


# TBD: Add doc
# TBD: Move to a utilities file
# TBD: Create unit tests for this function
# TBD: Test with dirs = { NA, NULL, "", c( "dir1", "dir2" ), etc
#
getTestPaths <- function ( pkg, dirs = "tests" ) {
  return ( getPkgPaths( pkg, dirs ) )
}  # end of function getTestPaths(pkg,dirs="tests")



# TBD: Add doc
# TBD: Add parameters with defaults
# TBD: What call this method: 'runAllTests', 'runAllVacemTests',
#      'runVacemTests' (with "all" as the default specifier for which tests to run) ?
# TBD: What are the verbose levels?  min is 0 (I think), what is the max?
#      Note, from RUnit 'options.r' file:
#              ##  integer: == 0: silent, >= 1: add comments to test case log
#
runAllTests <- function ( pkg=NULL, dirs="tests", verbose=1L, showDetails=FALSE ) {
  library( "RUnit" )

  options( warn = 1 )

  # TBD: Create function that gets the package name (does R have a reflection api?)
  #      Obviously, we could assign a default in function signature above, but
  #      if we used some type of reflection, then the code would work no matter
  #      what package we put it in.
  if ( isEmpty(pkg) ) {
    # cat( "No package name specified -- Defaulting to 'vacem'\n" )
    pkg = "vacem"
  }

  paths <- getTestPaths( pkg, dirs )
  if ( is.null(paths) ) {
    msg <- sprintf( paste( "Testing failed for package '%s'",
                    " -- Unable to find testing directory '%s'; Exiting\n" ),
                    pkg, join(dirs,sep="', '") )
    stop( msg )
  }
  
  suite <- defineTestSuite( name = pkg,
                            # dirs = ".",
                            # dirs = "/home/kcline/jhu/epi/dev/measles/coverage/pkgs/vacem/tests",
                            # dirs = "../tests",
                            dirs = paths,
                            testFileRegexp = ".*Test\\.r$",
                            rngKind = "default",
                            rngNormalKind = "default" )

  results <- runTestSuite( suite, verbose=verbose )
  printTextProtocol( results, showDetails=showDetails )

  # HACK HACK HACK
  # Note: This works, just decided not to use it yet. (KC, 1/4/12)
  #    TBD: Add parameter to function to control if HTML is generated
  #    TBD: Add code to generate a file name based on suite name and with a timestamp
  #    TBD: Add code to to test if file exists, add overwrite flag parameter
  #         and add code to generate a unique (temp) file name
  #    TBD: Add parameter to specify output directory and optional file name _pattern_
  #         for the output file
  #
  # printHTMLProtocol( results, fileName = "./vacem-tests.html" )

  # Saw something like this in the Roxygen package's tests directory;
  # still deciding if I need/want to use it...   -- K.Cline, 1/20/2012
  # 
  # errors <- getErrors( results )
  # if ( errors$nFail > 0 || errors$nErr > 0 ) {
  #   stop( 'vacem unit suite failed' )
  # }
  
}  # end of function runAllTests()

# Saw something like this in the Roxygen package's tests directory;
# still deciding if I need/want to use it...   -- K.Cline, 1/20/2012
# 
# if ( require('RUnit') ) {
#   library( 'vacem' )
#   runAllTests( verbose=0L )
# }

runAllTests( verbose=0L )


#
# RUnit print methods: (remove this eventually, just here for reference)
# printTextProtocol( testData, fileName = "", separateFailureList = TRUE, showDetails = TRUE, traceBackCutOff = 9 ) 
# printHTMLProtocol( testData, fileName = "", separateFailureList = TRUE,                     traceBackCutOff = 9,
#                    testFileToLinkMap = function(x) x ) 
#
#
