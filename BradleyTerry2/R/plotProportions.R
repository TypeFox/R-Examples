## P(win|not tie) in terms of expit(lambda_i - lambda_j)
GenDavidsonTie <- function(p){
    scale <- match("tie.scale", substring(names(coef), 1, 9), 0)
    if (scale != 0) scale <- exp(coef[scale])
    else scale <- 1
    tie.mode <- match("tie.mode", substring(names(coef), 1, 8), 0)
    if (tie.mode != 0) tie.mode <- coef["tie.mode"]
    delta <- coef[match("tie.max", substring(names(coef), 1, 7))]
    ## first player is at home
    weight1 <- plogis(tie.mode)
    weight2 <- 1 - weight1
    ## plogis = expit
    plogis(delta - scale * (weight1 * log(weight1) + weight2 * log(weight2)) +
           scale * (weight1 * log(p) + weight2 * log(1-p)))
}

#tmp <- eval(substitute(player1), data, parent.frame())
plotProportions <- function(win, tie = NULL, loss,
                            player1,
                            player2,
                            abilities = NULL,
                            home.adv = NULL,
                            tie.max = NULL,
                            tie.scale = NULL,
                            tie.mode = NULL,
                            at.home1 = NULL,
                            at.home2 = NULL,
                            data = NULL,
                            subset = NULL,
                            bin.size = 20,
                            xlab = "P(player1 wins | not a tie)",
                            ylab = "Proportion",
                            legend = NULL,
                            col = 1:2,
                            ...){
    call <- as.list(match.call())
    var <- intersect(names(call), c("win", "tie", "loss",
                                    "player1", "player2",
                                    "at.home1", "at.home2"))
    var <- var[!sapply(call[var], is.null)]

    dat <- with(data, do.call("data.frame", call[var]))
    if (!missing(subset)){
        subset <- eval(substitute(subset), data, parent.frame())
        dat <- subset(dat, subset)
    }
    if (!missing(tie) && sum(dat$tie) == 0) dat$tie <- NULL
    if (!is.null(home.adv) && (missing(at.home1) || missing(at.home2)))
        stop("at.home1 and at.home2 must be specified")
    if (!is.null(home.adv)){
        ## exclude neutral contests, make sure home player is first
        dat <- subset(dat, at.home1 | at.home2)
        swap <- which(as.logical(dat$at.home2))
        if (length(swap)) {
            dat$win[swap] <- dat$loss[swap]
            if (is.null(dat$tie)) dat$loss[swap] <- !dat$win[swap]
            else dat$loss[swap] <- !(dat$win[swap] | dat$tie[swap])
            tmp <- dat$player1[swap]
            dat$player1[swap] <- dat$player2[swap]
            dat$player2[swap] <- tmp
            dat$at.home1[swap] <- TRUE
            dat$at.home2[swap] <- FALSE
        }
    } else home.adv <- 0
    ### get proportions
    p <- with(dat, plogis(home.adv + abilities[as.character(player1)] -
                          abilities[as.character(player2)]))
    ## Depending on the distribution of p_ij (across all matches),
    ## divide the range of probabilities p_ij into discrete "bins", each
    ## of which has at least (say) 20 matches in it
    getBins <- function(p, bin.size) {
        ## alternatively estimate bins to same size intervals
        ## at least bin.size - distribute extra evenly over range
        min.size <- bin.size
        n <- length(p)
        r <- n %% min.size
        size <- rep(min.size, n %/% min.size)
        if (r > 0) {
            step <- length(size)/r
            extra <- round(seq(from = step/2 + 0.01,
                               to = step/2 + 0.01 + (r - 1)*step, by = step))
            size[extra] <- min.size + 1
        }
        bin <- factor(rep(seq(length(size)), size))[match(p, sort(p))]
        low <- sort(p)[cumsum(c(1, size[-length(size)]))] #first
        high <- sort(p)[cumsum(size)] #last
        mid <- (high - low)/2 + low
        list(bin = bin, mid = mid)
    }
    winBin <- getBins(p, bin.size)
    ## Within each bin b, calculate
    ## d_b = proportion of matches in that bin that were drawn
    if (!is.null(dat$tie)) {
        tieBin <- winBin
        tri <- with(dat, win - (!win & !tie))
        d_b <- tapply(tri, tieBin$bin, function(x) sum(x == 0)/length(x))
        ## recompute bins omitting ties
        winBin <- getBins(p[!dat$tie], bin.size)
    }

    ## h_b = proportion of *non-drawn* matches in that bin that were won
    ## by the home team
    if (!is.null(dat$tie)) {
        h_b <- tapply(tri[!dat$tie], winBin$bin,
                      function(x) sum(x == 1)/length(x))
    }
    else h_b <- tapply(dat$win, winBin$bin, function(x) sum(x == 1)/length(x))

    ## Plot d_b and h_b against the bin midpoints, in a plot with
    ## axis limits both (0,1)
    plot(h_b ~ winBin$mid, xlim = c(0, 1), ylim = c(0, 1),
         xlab = xlab, ylab = ylab, ...)
    if (missing(legend)) {
        if (is.null(dat$tie)) legend <- "Matches won"
        else legend <- c("Non-tied matches won", "Matches tied")
    }
    legend("topleft", legend, col = col[c(1, 2[!missing(tie)])],
           pch = 1)
    if (!is.null(dat$tie)) points(d_b ~ tieBin$mid, col = col[2])

    ## Add to the plot the lines/curves
    ## y = x
    ## y = expit(log(nu * sqrt(p_ij * (1 - p_ij))))
    ## The d_b should lie around the latter curve, and the h_b should
    ## lie around the former line.  Any clear patterns of departure are
    ## of interest.
    curve(I, 0, 1, add = TRUE)

    env <- new.env()
    environment(GenDavidsonTie) <- env
    coef <- na.omit(c(home.adv = unname(home.adv),
                      tie.max = unname(tie.max),
                      tie.scale = unname(tie.scale),
                      tie.mode = unname(tie.mode)))
    assign("coef", coef, envir=env)
    curve(GenDavidsonTie, 0, 1, col = col[2], add = TRUE)
    out <- list(win = data.frame(prop.win = h_b, bin.win = winBin$mid))
    if (!is.null(dat$tie))
        out <- c(out, tie = data.frame(prop.tie = d_b, bin.tie = tieBin$mid))
    invisible(out)
}
