#' Get x,y coordinates from a grid reference
#'
#' Given an OSGB or OSNI grid reference string, get the x,y coordinates of the
#' OSGB or OSNI grid for the bottom, left-hand corner of the grid square. The
#' \code{units} parameter controls the units (metres \code{m} or kilometres
#' \code{km}) in which the coordinates should be returned.
#' 
#' @export
#' @param grid an OSGB or OSNI grid reference string
#' @param units metres or kilometres (default)
#' @return a list of class "gridref" with the following contents:
#' \item{grid}{the original grid reference}
#' \item{system}{the grid reference system, either "OSGB" or "OSNI"}
#' \item{x}{the x coordinate (easting) in requested units}
#' \item{y}{the y coordinate (northing)in requested units}
#' \item{units}{"m" or "km"}
#' \item{precision}{the prcision of the original grid reference in metres}
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk}
#' @seealso \code{\link{gridRef}}
#' @examples
#' gridCoords("NZ1265") # returns x=412, y=565 (in km)
#' gridCoords("NZ123654", units="m") # returns x=412300, y=565400 (in m)
#' 
gridCoords <- function(grid=NULL, units=c("km", "m")) {
    
    ##-------------------------------------------------
    ## Bottom left hand corner of a tetrad within its
    ## 10km square in metres
    decodeTetrad <- function(letter) {
        l = as.integer(charToRaw(letter)) - 65
        if (l>13) l <- l-1 ## "O" is not used
        coord <- list()
        coord$x <- (l %/% 5) * 2000
        coord$y <- (l %% 5) * 2000 
        return(coord)
    }

    ## identify output format type
    units <- match.arg(units)
    
    ## clean up input grid ref
    ## make sure there are no spaces and letters are upper-case
    gr <- toupper(gsub(" ", "", grid))
    
    ## regular expression to check that a grid ref consists of
    ## one or two letters, 2-10 numbers and an optional tetrad letter
    ## The positions of matches to the letters, numbers and tetrad letter
    ## are returned in v. The first element of v is zero if no match found
    v <- regexec("^([H,N,O,S,T][A-H,J-Z]|[B-D,F-J,L-O,Q-T,V-X])([0-9]{2,10})([A-N,P-Z]{0,1})$", gr)
    if (v[[1]][[1]] > 0) {
        ## get the letters numbers and tetrad letter
        letters <- unlist(regmatches(gr, v))[2]  
        nums <- unlist(regmatches(gr, v))[3]
        tetrad <- unlist(regmatches(gr, v))[4]
        n <- nchar(nums)
        if ((n %% 2) == 0) {
            n <- n %/% 2
            precision <- 10^(5-n)
            gridref <- list()
            class(gridref) <- "gridref"
            gridref$grid <- grid
            if (nchar(letters) == 2) {
                gridref$system <- "OSGB"
            } else {
                gridref$system <- "OSNI"
            }
            ## now we can extract the coords
            x <- 0
            y <- 0
            l <- 1 ## set up to extract first letter
            if (gridref$system == "OSGB") {
                ## first letter of the OSGB 100km square code
                l = as.integer(charToRaw(substr(letters, l, l))) - 65
                if (l>7) l <- l-1 ## "I" is not used
                x <- ((l %% 5) - 2) * 500000
                y <- (3 - (l %/% 5)) * 500000
                l <- 2 ## set up to extract second letter
            }
            ## second letter of the OSGB 100km square code
            ## or the only letter for OSNI grid refs
            l = as.integer(charToRaw(substr(letters, l, l))) - 65
            if (l>7) l <- l-1 ## "I" is not used
            x <- x + (l %% 5) * 100000
            y <- y + (4 - (l %/% 5)) * 100000
            ## deal with the numbers
            x <- x + as.integer(substr(nums, 1, n)) * precision
            y <- y + as.integer(substr(nums, n+1 ,n+1+n)) * precision
            ## deal with tetrad letter
            if (nchar(tetrad) == 1) {
                c <- decodeTetrad(tetrad)
                x <- x + c$x
                y <- y + c$y
                precision <- 2000
            }
            
            ## prepare list to return
            if (units == "km") {
                x <- x / 1000
                y <- y / 1000
            }
            gridref$x <- x
            gridref$y <- y
            gridref$precision <- precision
            gridref$units <- units
            return(gridref)
        } else {
            stop("must be an even number of digits")
        }   
    } else {
        stop("not a valid grid reference string")
    }
}
