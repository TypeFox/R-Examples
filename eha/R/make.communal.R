make.communal <-
function (dat, com.dat, communal = TRUE, start, period = 1, lag = 0, 
    surv = c("enter", "exit", "event", "birthdate"), tol = 1e-04, 
    fortran = TRUE) 
{
    if (!is.data.frame(dat)) 
        stop("dat must be a data frame")
    if (!is.data.frame(com.dat)) 
        stop("com.dat must be a data frame")
    if (length(surv) != 4) 
        stop("surv must have length 4")
    fixed.names <- names(dat)
    surv.indices <- match(surv, fixed.names)
    if (length(which(is.na(surv.indices)))) {
        x <- which(is.na(surv.indices))
        stop(paste(surv[x], " is not a name in the fixed data frame."))
    }
    com.names <- names(com.dat)
    if (length(x <- which(!is.na(match(com.names, c(surv, fixed.names)))))) 
        stop(paste(com.names[x], "are names in fixed data frame."))
    nn <- nrow(dat)
    n.years <- nrow(com.dat)
    n.com <- NCOL(com.dat)
    if (n.com > 1) 
        stop("Only one communal covariate at a time!")
    cuts <- start + c(0, (1:n.years) * period) - lag
    beg.per <- cuts[1]
    end.per <- cuts[n.years + 1]
    ##cat("beg.per = ", beg.per, "end.per = ", end.per, "\n")
    iv.length <- period
    if (communal){ ## Added 26 March 2004 ##
        spell.tot <- sum(dat[, surv.indices[2]] - dat[, surv.indices[1]])
        dat <- cal.window(dat, c(beg.per, end.per), surv)
        if (sum(dat[, surv.indices[2]] - dat[, surv.indices[1]]) < 
            spell.tot) 
            warning("Spells are cut")
    }#######################################
    nn <- nrow(dat)
    if (!communal) {
        get.per <- function(dates) pmin(pmax(1, ceiling((dates - 
            beg.per)/iv.length)), n.years)

        dates <- dat[, surv.indices[4]]
        
        ppp <- get.per(dates)

        yy <- matrix(0, ncol = n.com, nrow = nn)
        for (i in 1:n.com) {
            yy[, i] <- com.dat[ppp, i]
        }
        ## Added 26 March 2004: ##
        yy <- ifelse((dat[, surv.indices[4]] <= beg.per) |
                     (dat[, surv.indices[4]] > end.per),
                     NA, yy)
        ##########################       
        yy <- as.data.frame(yy)
        names(yy) <- com.names
        yy <- cbind(dat, yy)
    }
    else {
        get.iv <- function(dates) cbind(pmin(pmax(1, floor((dates[, 
            1, drop = FALSE] - beg.per)/iv.length) + 1), n.years), 
            pmin(pmax(1, ceiling((dates[, 2, drop = FALSE] - 
                beg.per)/iv.length)), n.years))

        event.boolean <- is.logical(dat[, surv.indices[3]])
        xx <- cbind(dat[, surv.indices, drop = FALSE], 1:nn)
        xx[, 3] <- as.numeric(xx[, 3])
        xx <- as.matrix(xx)
        if (!is.numeric(xx)) 
            stop("Internal error in [make.communal]: xx not numeric")
        ind.date <- cbind(xx[, 1, drop = FALSE] + xx[, 4, drop = FALSE], 
            xx[, 2, drop = FALSE] + xx[, 4, drop = FALSE])
        cases <- ((ind.date[, 1] < end.per) & (ind.date[, 2] > 
            beg.per))
        xx <- xx[cases, , drop = FALSE]
        ind.date <- ind.date[cases, , drop = FALSE]
        ind.iv <- get.iv(ind.date)
        
        nn <- nrow(xx)
        nn.out <- sum(ind.iv[, 2] - ind.iv[, 1] + 1)
        yy <- matrix(0, nrow = nn.out, ncol = ncol(xx) + 1)
        nn.out <- ind.iv[, 2] - ind.iv[, 1] + 1
        cur.row <- 0
        com.dat <- as.matrix(com.dat)
        split <- function(i) {
            n.rows <- nn.out[i]
            if (n.rows == 1) {
                return(list(c(xx[i, ], ind.iv[i, 1])))
            }
            else {
                x.i <- xx[i, ]
                out <- matrix(0, nrow = n.rows, ncol = ncol(yy))
                out[n.rows, 3] <- x.i[3]
                out[, 4] <- x.i[4]
                out[, 5] <- x.i[5]
                out[1, 1] <- x.i[1]
                out[1, 2] <- cuts[ind.iv[i, 1] + 1] - x.i[4]
                out[1, 6] <- ind.iv[i, 1]
                out[n.rows, 1] <- cuts[ind.iv[i, 2]] - x.i[4]
                out[n.rows, 2] <- x.i[2]
                out[n.rows, 6] <- ind.iv[i, 2]
                if (n.rows > 2) {
                  for (j in 2:(n.rows - 1)) {
                    out[j, 1] <- out[j - 1, 2]
                    out[j, 2] <- out[j - 1, 2] + iv.length
                    out[j, 6] <- ind.iv[i, 1] + j - 1
                  }
                }
                return(out)
            }
        }
        beg.row <- end.row <- 0
        if (!fortran) {
            for (j in 1:nn) {
                beg.row <- end.row + 1
                end.row <- end.row + nn.out[j]
                yy[beg.row:end.row, ] <- split(j)
                if (j%/%100 * 100 == j) 
                  cat("j = ", j, "\n")
                NULL
            }
        }
        if (fortran) {
            yy <- .Fortran("split",
                           as.double(xx),
                           as.integer(nn), 
                           as.integer(ncol(xx)),
                           yy = as.double(yy),
                           as.integer(nrow(yy)), 
                           as.integer(ncol(yy)),
                           as.integer(nn.out),
                           as.integer(ind.iv), 
                           as.double(cuts),
                           as.integer(n.years),
                           ## DUP = FALSE,
                           PACKAGE = "eha")$yy
            yy <- matrix(yy, ncol = ncol(xx) + 1)
        }
        yy <- cbind(yy[, 1:4, drop = FALSE], dat[yy[, 5], -surv.indices, 
            drop = FALSE], com.dat[yy[, 6], , drop = FALSE])
        names(yy)[1:4] <- surv
        all.names <- c(surv, fixed.names[-surv.indices], com.names)
        row.names(yy) <- as.character(1:nrow(yy))
        yy <- as.data.frame(yy)
        names(yy) <- all.names
        if (event.boolean) 
            yy[, 3] <- as.logical(yy[, 3])
    }
    yy
}
