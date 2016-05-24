#' read_px
#'
#' Fixed version of read.px from pxR package
#' 
#' @param filename filename
#' @param encoding encoding
#' @param na.strings na.strings 
#'
#' @return px object
#' @references
#' See citation("statfi") 
#' @author Contact: Leo Lahti \email{louhos@@googlegroups.com}
#' @export
#' @importFrom pxR read.px
#' @examples \dontrun{px <- read_px("http://pxweb2.stat.fi/database/StatFin/vrm/synt/080_synt_tau_203.px")}
#' @keywords utilities

read_px <- function(filename, encoding = "latin1", 
                    na.strings = c('"."', '".."', '"..."', '"...."')) {


  label <- NULL
  attribute <- NULL
  value <- NULL

  # filename <- "http://pxweb2.stat.fi/database/StatFin/asu/ashi/003_ashi_tau_108.px"; encoding = "latin1"; na.strings = c('"."', '".."', '"..."', '"...."')

    ## auxiliary functions ##
    unquote <- function(x){
        gsub('\\"', "", x)
    }

    clean.spaces <- function(x){
        gsub("^[[:space:]]+|[[:space:]]+$", "", x) # elimina blancos por delante|detrás
    }

    get.attributes <- function(x){
        x <- gsub( "([A-Z-]*)\\((.*)\\).*", "\\1;\\2", x ) ## parte etiqueta y atributo con ";"
        x <- strsplit( x, ";" )
        x <- lapply( x, function( x ) c( x, rep( "value", 2 - length( x ) ) ) )
        x <- do.call( rbind, x )
        x[,2] <- unquote( x[,2] )
        clean.spaces( x )
    }

    break.clean <- function( x, sep = '\\"' ) {
        x <- strsplit( x, sep )[[1]]
        if (sep != " ") x <- clean.spaces( x )
        x <- x[ x != "" ]
        x <- x[ x != "," ]
        x
    }

    make.list <- function( dat, my.label ){
        dat <- subset( dat, label == my.label, select = c( attribute, value ) )

        my.list <- as.list( dat$value )
        names( my.list ) <- dat$attribute
        my.list
    }

    ## end: auxiliary functions ##
    a <- scan(filename, what = "character", sep = "\n", quiet = TRUE, fileEncoding = encoding)

    # Recognize field start points
    fb.starts <- c(1, grep(";$", a) + 1)
    fb.starts <- fb.starts[fb.starts <= length(a)]
    fb.ends <- c(fb.starts[-1] - 1, length(a))
    breakpoints <- rbind(fb.starts, fb.ends)

    # add last index separately to avoid its separate handling in
    # the following for loop
    #fb.starts <- c(fb.starts, length(a)) 
    a2 <- lapply(1:ncol(breakpoints), function (i) {
      gsub(";$", "", paste(a[breakpoints[1, i]:breakpoints[2, i]], collapse = " "))
    })

    # ------------------------------------------

    # CRASH:
    # a <- do.call(rbind, strsplit(a, "//=//" ))
    # FIX by antagomir 23.9.2012
    a2 <- sub( "=", "//=//", a2 )
    a2 <- a2[!a2 == " "]
    a2 <- lapply(a2, function (x) {strsplit(x, "//=//")})
    attribs <- get.attributes(sapply(a2, function (x) {x[[1]][[1]]}))
    vals <- sapply(a2, function (x) {x[[1]][[2]]})

    # ------------------------------------------

    a <- data.frame(cbind(attribs, vals))
    colnames(a) <- c("label", "attribute", "value")

    a$label     <- make.names(a$label)
    a$attribute <- make.names(a$attribute)
    a$value     <- as.character(a$value)

    ## build a px object: list with px class attribute ##
    px <- sapply(unique( a$label ), function(label) make.list(a, label), simplify = FALSE)

    # turns data values into an R vector
    px$STUB$value    <- make.names(break.clean(px$STUB$value))
    px$HEADING$value <- make.names(break.clean(px$HEADING$value))
    px$VALUES <- lapply(px$VALUES, break.clean )
    px$CODES  <- lapply(px$CODES,  break.clean )

    tmp <- gsub('"-"', 0, px$DATA$value) # 0 can be encoded as "-"
    dat <- textConnection(tmp) # much faster than with cleanDat (strsplit)

    # This fix was needed to interpret for instance:
    # "http://pxweb2.stat.fi/database/StatFin/hin/pthi/004_pthi_tau_004_fi.px"
    # px$DATA$value <- scan(dat, na.strings = na.strings, quiet = TRUE)
    px$DATA$value <- as.numeric(scan(dat, na.strings = na.strings, quiet = TRUE, what = "character"))

    close(dat)
    
    class(px) <- "px"

    px

}


