######################################################################################
#  makeTdmRandomSeed
#
#' Factory method to make a function generating pseudo-random random number seeds.
#'
#' Create first an object of type \code{\link{makeTdmRandomSeed}} and then call the 
#' returned value of that object (a function) as many times as you like. (It is necessary
#' to create the function object first to have in its environment the private storage for 
#' the number of calls to that object.)
#' 
#' @param ID [0] each random seed genarator with a different ID will generate different seeds.
#'    In this way it is possible that parallel jobs ID=0,1,2,... starting in the same second 
#'    are initialized with different seeds and thus produce different results.
#' @return A function object which can be invoked without any arguments and returns
#'    each time a different integer in 0...100001+nCall. This is true even if it is called
#'    many times within the same second (where Sys.time() will return the same integer).
#'    nCall is the number of calls to this function object.
#'
#' @examples
#'
#' tdmRandomSeed = makeTdmRandomSeed();
#' for (i in 1:10) print(c(as.integer(Sys.time()), tdmRandomSeed()));
#'
#' @seealso \code{\link{tdmRandomSeed}}
#' @author Wolfgang Konen, Patrick Koch \email{wolfgang.konen@@fh-koeln.de}
#' @export
#' @keywords internal
makeTdmRandomSeed <- function(ID=0) {
  # this provides private local storage for the function getSeed below, it remains 
  # there even after leaving getSeed. It stores the number of calls to getSeed.
  seedModBuf <- ID;   
  
  getSeed <- function() {
    # with '<<-' assignment to seedModBuf one level above:
    seedModBuf <<- seedModBuf+1;
    #
    # Sys.time() stays the same for 1 sec. By incrementing seedModBuf in each
    # call we ensure that seed will be different in each call (a number from
    # {0,1,...,seedModBuf+100001}), even if called multiple times within a second
    seed <- as.integer(Sys.time()) %% (seedModBuf+100001)
  }
  getSeed;
}

######################################################################################
#  tdmRandomSeed
#
#' Generates pseudo-random random number seeds.
#'
#' To use this mechanism, create first an object \code{tdmRandomSeed} with a call to 
#' \code{\link{makeTdmRandomSeed}}. 
#'
#' @return In each call to this function a different integer 
#'    in 0...100001+nCall is returned. This is true even if it is called
#'    many times within the same second (where Sys.time() will return the same integer).
#'    nCall is the number of calls to this function object.
#'
#' @examples
#'
#' tdmRandomSeed = makeTdmRandomSeed();
#' for (i in 1:10) print(c(as.integer(Sys.time()), tdmRandomSeed()));
#'
#' @seealso \code{\link{makeTdmRandomSeed}}
#' @author Wolfgang Konen, Patrick Koch \email{wolfgang.konen@@fh-koeln.de}
#' @export
#
tdmRandomSeed <- makeTdmRandomSeed();
