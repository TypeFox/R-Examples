#' Manipulate OSGB or OSNI grid reference string
#' 
#' Extracts grid reference strings at various precisions from the supplied grid
#' reference string - if possible! For example, if you supply a 1km square 
#' reference \code{TL2998}, then you could get the 10km square \code{TL29}, but
#' not a 100m square grid reference.
#' 
#' @export
#' @param grid the grid reference to be manipulated
#' @param format the format you want back. The possibilities are: \code{sq10km},
#'   \code{sq5km}, \code{tetrad}, \code{sq1km}, \code{sq100m}, \code{sq10m}
#' @return a list of class "gridref" with the following contents:
#' \item{grid}{the original grid reference}
#' \item{gfmt}{the grid reference string formatted as requested}
#' \item{system}{the grid reference system, either "OSGB" or "OSNI"}
#' \item{precision}{the precision of the formatted grid reference in metres}
#' @note \bold{Tetrads and 5kms}\cr
#' Tetrads are 2x2km squares and are often used for mapping distributions at a
#' the scale of a county or similar sized local area. They are labelled using 
#' the 10km square followed by a single, upper-case letter (since there are 25
#' tetrads in a 10km square, the letter "O" is not used to avaiod confusion with
#' zero). For example: \code{TL29S}. Letters are used as follows:\cr
#' 
#' \code{0 +---+---+---+---+---+}\cr
#' \code{. | E | J | P | U | Z |}\cr
#' \code{8 +---+---+---+---+---+}\cr
#' \code{. | D | I | N | T | Y |}\cr
#' \code{6 +---+---+---+---+---+}\cr
#' \code{. | C | H | M | S | X |}\cr
#' \code{4 +---+---+---+---+---+}\cr
#' \code{. | B | G | L | R | W |}\cr
#' \code{2 +---+---+---+---+---+}\cr
#' \code{. | A | F | K | Q | V |}\cr
#' \code{0 +---+---+---+---+---+}\cr
#' \code{..0...2...4...6...8...0}\cr
#' 
#' This is named the DINTY system after the letters in the second row of this
#' table.\cr
#' 
#' 5x5km squares (sometimes called "pentads") are used for mapping at a regional
#' scale. They are labelled using the name of the 10km square followed by two
#' upper-case letters as follows:\cr
#' 
#' \code{ 0 +----+----+}\cr
#' \code{ . | NW | NE |}\cr
#' \code{ 5 +----+----+}\cr
#' \code{ . | SW | SE |}\cr
#' \code{ 0 +----+----+}\cr
#' \code{ ..0....5....0}\cr
#' 
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk}
#' @seealso \code{\link{gridCoords}}
#' @examples
#' gridRef("TL2998", "sq10km") # Returns TL29
#' gridRef("TL29", "sq1km") # Returns NULL - you cannot get a 1km from a 10km!
#' 
gridRef <- function(grid=NULL, format=c("sq10km","sq5km","tetrad","sq1km","sq100m","sq10m")) {
    
    ##----------------------------------------------------------
    tetradLetter <- function(nums, n) {
        x2 <- as.integer(substr(nums, 2, 2))
        y2 <- as.integer(substr(nums, n+2, n+2))
        l <- (x2 %/% 2) * 5 + (y2 %/% 2) 
        return(substr("ABCDEFGHIJKLMNPQRSTUVWXYZ", l+1, l+1))
    }

    ##----------------------------------------------------------
    pentadLetter <- function(nums, n) {
        x2 <- as.integer(substr(nums, 2, 2))
        y2 <- as.integer(substr(nums, n+2, n+2))
        l <- (x2 %/% 5) * 2 + (y2 %/% 5) 
        return(c("SW","NW","SE","NE")[l+1])
    }
    
    ## identify output format type
    format <- match.arg(format)
    
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
            gret <- list()
            class(gret) <- "gridref"
            gret$grid <- grid
            ifelse(nchar(letters) == 2, gret$system <- "OSGB", gret$system <- "OSNI")
            switch(format,
                   sq10km={
                       if (n>0) {
                           gret$gfmt <- paste(letters, substr(nums, 1, 1), substr(nums, n+1, n+1), sep="") 
                           gret$precision <- 10000
                       } else {
                           gret <- NULL
                       }
                   },
                   sq5km={
                       if (n>1) {
                           gret$gfmt <- paste(letters, substr(nums, 1, 1), substr(nums, n+1, n+1), 
                                       pentadLetter(nums, n), sep="") 
                           gret$precision <- 5000
                       } else {
                           gret <- NULL
                       }
                   },
                   tetrad={
                       if (n>1) {
                           gret$gfmt <- paste(letters, substr(nums, 1, 1), substr(nums, n+1, n+1), 
                                       tetradLetter(nums, n), sep="") 
                           gret$precision <- 2000
                       } else {
                           gret <- NULL
                       }
                   },
                   sq1km={
                       if (n>1) {
                           gret$gfmt <- paste(letters, substr(nums, 1, 2), substr(nums, n+1, n+2), sep="") 
                           gret$precision <- 1000
                       } else {
                           gret <- NULL
                       }
                   },
                   sq100m={
                       if (n>2) {
                           gret$gfmt <- paste(letters, substr(nums, 1, 3), substr(nums, n+1, n+3), sep="") 
                           gret$precision <- 100
                       } else {
                           gret <- NULL
                       }
                   },
                   sq10m=if (n>3) {
                           gret$gfmt <- paste(letters, substr(nums, 1, 4), substr(nums, n+1, n+4), sep="") 
                           gret$precision <- 10
                   } else {
                           gret <- NULL
                   },
                   stop("requested output format not recognised")
                   )
            return(gret)
        } else {
            stop("must be an even number of digits")
        }   
    } else {
        stop("not a valid grid reference string")
    }
}
