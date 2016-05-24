#' Constructor of ecdq class
#' 
#' Construct an ecdq class by providing the required parameters.
#'
#' @param ecd An object of ecd class
#' @param verbose logical, display timing information, for debugging purpose.
#'
#' @return An object of ecdq class 
#'
#' @keywords constructor
#'
#' @author Stephen H. Lihn
#'
#' @export ecdq
#'
#' @importFrom stats na.omit lm
#'
#' 
#' @examples
#' \dontrun{
#'     d <- ecd()
#'     dq <- ecdq(d)
#' }
### <======================================================================>
"ecdq" <- function(ecd, verbose=FALSE)
{
    call <- match.call()

    if(class(ecd) != "ecd"){
      stop("Parameter 'ecd' must be an ecd object!\n")
    }
    if (length(ecd@stats)==0) {
        stop("stats is not computed in ecd object (ecd)")
    }

    d <- ecd
    mean <- d@stats$mean
    dv <- d@stats$stdev
    kurt <- d@stats$kurtosis
    xe <- ellipticity(d)
    
    
    N_dv <- 4 # number of stdev to cover on each side of mean
    S_dv <- if(kurt > 12.258) 4 else 2 # divisions within each stdev
    N_fit <- 8 # number of points to fit in each division
        
    x_dv <- mean + dv/S_dv*seq(-N_dv*S_dv, N_dv*S_dv)
    points <- ecd.mp2f(c(x_dv, xe$xe1, xe$xe2)) # x doesn't need mpfr

    if (verbose) print(paste(Sys.time(), "ecdq: x points=", length(points)))
    
    xseg0 <- sort(na.omit(unique(points)))
    xseg.to <- ecd.lag(xseg0, -1, na.omit=TRUE)
    xseg.from <- xseg0[1:length(xseg.to)]
    
    cseg0 <- ecd.cdf(d, xseg0, piece.wise=TRUE)
    cseg.to <- ecd.lag(cseg0, -1, na.omit=TRUE)
    N_seg <- length(cseg.to)
    cseg.from <- cseg0[1:N_seg]
    
    cseg.min <- min(cseg.from)
    cseg.max <- max(cseg.to)
    
    conf <- list(mean=mean, stdev=dv, xe1=xe$xe1, xe2=xe$xe2, 
                 N_fit=N_fit, N_stdev=N_dv, S_stdev=S_dv)
    
    if (min(abs(xseg.from-xseg.to)) < .Machine$double.eps) {
        stop("ERROR: segment(s) too fragmented!")
    }
    gen_fit <- function(idx) {
        x <- ecd.mp2f(xseg.from[idx] + seq(0,N_fit)*(xseg.to[idx] - xseg.from[idx])/N_fit)
        cdf <- ecd.mp2f(ecd.cdf(d, x, piece.wise=TRUE))
        lm(x ~ poly(cdf,6))
    }
    cdf.fit <- parallel::mclapply(seq(1,N_seg), gen_fit, mc.allow.recursive=FALSE)
    if (verbose) print(paste(Sys.time(), "ecdq: cdf.fit generated"))
    
    # tail region, direction: left: -1, right: +1
    get_x_tail <- function(direction) {
        x_tail <- NULL
        if (direction == -1) {
            x_tail <- min(xseg.from)
        }
        if (direction == 1) {
            x_tail <- max(xseg.to)
        }
        x_tail
    }
    gen_tail_fit <- function(direction, N, S) {
        x_tail <- get_x_tail(direction)
        dx <- abs(x_tail - mean)
        x <- ecd.mp2f(seq(0,N*S)*dx/S * direction + x_tail)
        if (direction == -1) {
            # left tail
            lcdf <- ecd.mp2f(log(ecd.cdf(d,x, piece.wise=FALSE)))
            # note, in tail region, do not use piece.wise=T
            lm(x ~ poly(lcdf,6))
        } else {
            # right tail
            lccdf <- ecd.mp2f(log(ecd.ccdf(d,x, piece.wise=FALSE)))
            # note, in tail region, do not use piece.wise=T
            lm(x ~ poly(lccdf,6))
        }
    }
    N <- 10
    S <- 4

    conf$N_tail <- N
    conf$S_tail <- S
    
    fit.tail <- parallel::mclapply(c(-1,1), function(dir) gen_tail_fit(dir, N, S))
    if (verbose) print(paste(Sys.time(), "ecdq: fit.tail generated"))
    
    dq <- new("ecdq", call = call,
              xseg.from = xseg.from,
              xseg.to = xseg.to,
              cseg.from = cseg.from,
              cseg.to  = cseg.to,
              cseg.min = cseg.min,
              cseg.max = cseg.max,
              N_seg = N_seg,
              cdf.fit = cdf.fit,
              x_left_tail = get_x_tail(-1),
              x_right_tail = get_x_tail(1),
              fit.left = fit.tail[[1]],
              fit.right = fit.tail[[2]],
              conf = conf)
    
    invisible(dq)
}
### <---------------------------------------------------------------------->
