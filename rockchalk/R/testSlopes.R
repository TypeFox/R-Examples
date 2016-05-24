##' Hypothesis tests for Simple Slopes Objects
##'
##' Conducts t-test of the hypothesis that the "simple slope" line for
##' one predictor is statistically significantly different from zero
##' for each value of a moderator variable. The user must first run
##' \code{plotSlopes()}, and then give the output object to
##' \code{plotSlopes()}. A plot method has been implemented for
##' testSlopes objects. It will create an interesting display, but
##' only when the moderator is a numeric variable.
##'
##' This function scans the input object to detect the focal values of
##' the moderator variable (the variable declared as \code{modx} in
##' \code{plotSlopes}). Consider a regression with interactions
##'
##' y <- b0 + b1*x1 + b2*x2 + b3*(x1*x2) + b4*x3 + ... + error
##'
##' If \code{plotSlopes} has been run with the argument plotx="x1" and
##' the argument modx="x2", then there will be several plotted lines,
##' one for each of the chosen values of x2.  The slope of each of
##' these lines depends on x1's effect, b1, as well as the interactive
##' part, b3*x2.
##'
##' This function performs a test of the null hypothesis of the slope
##' of each fitted line in a \code{plotSlopes} object is statistically
##' significant from zero. A simple t-test for each line is offered.
##' No correction for the conduct of multiple hypothesis tests (no
##' Bonferroni correction).
##'
##' When \code{modx} is a numeric variable, it is possible to conduct
##' further analysis. We ask "for which values of \code{modx} would
##' the effect of \code{plotx} be statistically significant?"  This is
##' called a Johnson-Neyman (Johnson-Neyman, 1936) approach in
##' Preacher, Curran, and Bauer (2006). The interval is calculated
##' here.  A plot method is provided to illustrate the result.
##'
##' @param object Output from the plotSlopes function
##' @return A list including 1) the hypothesis test table, 2) a copy of
##' the plotSlopes object, and, for numeric
##' modx variables, 3) the Johnson-Neyman (J-N) interval boundaries.
##' @export
##' @import car
##' @seealso plotSlopes
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @references
##' Preacher, Kristopher J, Curran, Patrick J.,and Bauer, Daniel J. (2006).
##' Computational Tools for Probing Interactions in Multiple Linear
##' Regression, Multilevel Modeling, and Latent Curve Analysis.
##' Journal of Educational and Behavioral Statistics. 31,4, 437-448.
##'
##' Johnson, P.O. and Neyman, J. (1936). "Tests of certain linear
##' hypotheses and their applications to some educational problems.
##' Statistical Research Memoirs, 1, 57-93.
##' @example inst/examples/testSlopes-ex.R
testSlopes <-
    function(object)
{
    if (!inherits(object, "plotSlopes"))
        stop("use only with \"plotSlopes\" objects")
    model <-  eval(parse(text=object$call$model))
    plotx <- object$call$plotx
    modx <- object$call$modx

    ivs <- attr(terms(model), "term.labels")
    relevantInteractions <- c(paste(plotx, ":", modx, sep = ""),
                              paste(modx, ":", plotx, sep = ""))

    interactionsIn <- relevantInteractions[which(relevantInteractions %in% ivs)]

    if (!any(relevantInteractions %in% ivs)) {
        cat("There were no interactions in the plotSlopes object, so testSlopes can't offer any advice.\n")
        return(NULL)
    }

    ## Whew. We did not return, so let's get down to business.
    modxVar <- object$newdata[ , modx]
    modxVals <- object$modxVals
    bs <- coef(model)
    V <- vcov(model)

    ## If modx is a factor, only need to test a particular set of lines.
    if (is.factor(modxVar)) {
        mcoef <- coef(model)
        modxContrastNames <- c(grep(plotx,  grep(modx, names(mcoef), value = TRUE),  value=TRUE))
        slope <- mcoef[modxContrastNames] + coef(model)[plotx]

        ## 2013-02-19 Sunthud Pornprasertmanit spots bug here:
        ## Problem: diag doesn't work when argument is a single real number.
        ## Fix by inserting drop=FALSE (wrote a blog post about the "drop gotcha"
        seslope <- sqrt(V[plotx, plotx, drop = FALSE] +  diag(V[modxContrastNames, modxContrastNames, drop = FALSE]) + 2* V[modxContrastNames, plotx, drop = FALSE])

        modxContrastNames <- c(plotx, modxContrastNames)
        slope <- c(mcoef[plotx], slope)
        seslope <- c(sqrt(V[plotx, plotx]) , seslope)
        tbsimple <- slope/seslope

        testSlopes <- data.frame(modx = modxContrastNames, b = slope,
                                 se = seslope, t = tbsimple, p = 2 *
                                 pt(abs(tbsimple), df =
                                 model$df.residual, lower.tail =
                                 FALSE))

        row.names(testSlopes) <-  levels(model.frame(model)[ , modx])

        colnames(testSlopes) <- c(deparse(modx), "slope",
                                  "Std. Error", "t value", "Pr(>|t|)")

        cat(paste("These are the straight-line \"simple slopes\" of the variable", plotx, " \n for the selected moderator values. \n"))
        print(testSlopes)
        res <- list("hypotests" = testSlopes, pso = object)
        class(res) <- "testSlopes"
        return(invisible(res))
    }

    ## Aha! We did not return. So this input is numeric. Continue.

    bmodx <- NULL
    bplotx <- bs[plotx]

    bmodx <- bs[interactionsIn]
    bsimple <- bplotx + bmodx * modxVals
    covbsimple <- cbind(1, modxVals^2, 2 * modxVals) %*%
        c(V[plotx, plotx], V[names(bmodx), names(bmodx)], V[plotx, names(bmodx)])
    tbsimple <- bsimple/sqrt(covbsimple)

    testSlopes <- data.frame( modx = modxVals, b = bsimple,
                             se = sqrt(covbsimple), t = tbsimple,
                             p = 2 * pt(abs(tbsimple),
                             df = model$df.residual, lower.tail = FALSE))

    colnames(testSlopes) <- c(deparse(modx), "slope", "Std. Error", "t value", "Pr(>|t|)")

    ## Just for numeric variables, the J-N calculation
    roots <- NULL

    mm <- model.matrix(model)
    b2 <- bmodx <- NULL
    b1 <- bplotx <- bs[plotx]

    b2 <- bs[interactionsIn]

    if(is.null(b2)) stop("b2 is null, there's no interation! Logic error")

    Tcrit <- qt(0.975, model$df)

    ## Quadratic formula. Solve
    ## Tcrit < T = (b1 + b2*modx)/s.e.(b1 + b2*modx)
    ## Tcrit < T = (b1 + b2*modx)/sqrt(Var(b1) + modx^2*Var(b2) + 2 Var[b1,b2])

    jn <- list()
    jn$a <- b2^2 - (Tcrit^2) * V[interactionsIn, interactionsIn]
    jn$b <- 2*b1*b2 - (Tcrit^2) * 2 * V[plotx, interactionsIn]
    jn$c <- b1^2 - (Tcrit^2) * V[plotx, plotx]
    jn$inroot <- (jn$b^2) - 4 * jn$a * jn$c

    ##complex root check, if yes, exit with message.
    if (jn$inroot <= 0) {
        print(paste("There are no real roots to the quadratic equation that represents regions of statistical significance."))
            if (jn$a > 0) {
                print("That means the slope (b1 + b2modx)plotx is statistically significant for all values of modx")
            } else {
                print("In this case, that means the  slope (b1 + b2modx)plotx is never statistically significant")
            }
        res <- list("hypotests" = testSlopes, "jn" = jn, pso = object)
        class(res) <- "testSlopes"
        return(invisible(res))
    }


    ## else (!inroot <= 0)
    ## Whew. Roots not complex, otherwise would have returned
    jn$roots <- c( (-jn$b - sqrt(jn$inroot))/(2*jn$a),
                  (-jn$b + sqrt(jn$inroot))/(2*jn$a) )
    jn$roots <- sort(jn$roots)
    names(jn$roots) <- c("lo","hi")
    if (jn$a > 0) {
        cat(paste("Values of", modx, "OUTSIDE this interval:\n"))
        print(jn$roots)
        cat(paste("cause the slope of (b1 + b2*", modx,")", plotx, " to be statistically significant\n", sep = ""))
    } else {
        cat(paste("Values of", modx, "INSIDE this interval:\n"))
        print(jn$roots)
        cat(paste("cause the slope of (b1 + b2*", modx,")", plotx, " to be statistically significant\n", sep = ""))
    }


    ## print(paste("b1 = b2x at", -b1/b2))
    ## print(paste("quadratic minimum/maximum point at", -jn$b/(2*jn$a)))
    ## if(jn$a > 0) print("that is a minimum") else print("that is a maximum")
    res <- list("hypotests" = testSlopes, "jn" = jn, pso = object)
    class(res) <- "testSlopes"
    return(invisible(res))
}

NULL

##' Plot testSlopes objects
##'
##' plot.testSlopes is a method for the
##' generic function plot. It has been revised so that it creates a plot
##' illustrating the marginal effect, using the
##' Johnson-Neyman interval calculations to highlight the
##' "statistically significantly different from zero" slopes.
##' @return \code{NULL}
##' @author <pauljohn@@ku.edu>
##' @method plot testSlopes
##' @export
##' @param x   A testSlopes object.
##' @param ... Additional arguments that are ignored currently.
##' @param shade Optional. Create colored polygon for significant regions.
##' @param col Optional. Color of the shaded area. Default transparent pink.
plot.testSlopes <-
    function(x, ..., shade = TRUE, col = rgb(1, 0, 0, 0.10))
{
    ## Following should be unnecessary.
    if (!inherits(x, "testSlopes"))
         stop("use only with \"testSlopes\" objects")
    tso <- x
    model <-  eval(parse(text = tso$pso$call$model))
    modx <- tso$pso$call$modx
    modxVar <- model$model[, modx]
    modxRange <- magRange(modxVar, 1.1)

    if (x$jn$inroot <= 0) {
        print("There are no real roots. There is nothing worth plotting")
        return(NULL)
    }

    if ((x$jn$roots["hi"] < modxRange[1]) || (x$jn$roots["lo"] > modxRange[2])) {
            ## both roots outside observed x
    }

    ## tedious exercise to figure out which type of interval of significant values we have
    if (x$jn$roots["lo"] > modxRange[1] && x$jn$roots["hi"] < modxRange[2]){
        ##2 roots inside
        idxStart <- c(modxRange[1], x$jn$roots)
        idxEnd <- c(x$jn$roots, modxRange[2])
        intervals <- if(x$jn$a > 0) c(1,3) else c(2)
        markAt <- x$jn$roots
    } else if (modxRange[1] < x$jn$roots["lo"]){
        ## hi root inside
        idxStart <- c(modxRange[1], x$jn$roots["lo"])
        idxEnd <- c(x$jn$roots["lo"], modxRange[2])
        intervals <- if(x$jn$a > 0) c(1) else c(2)
        markAt <- x$jn$roots["lo"]
    } else if (x$jn$roots["hi"] < modxRange[2]){
        ## lo root inside
        idxStart <- c(modxRange[1], x$jn$roots["hi"])
        idxEnd <- c(x$jn$roots["hi"], modxRange[2])
        intervals <- if(x$jn$a > 0) c(2) else c(1)
        markAt <- x$jn$roots["hi"]
    } else {
        stop("Your observed moderator data does not overlap with
              the region on which the slope would be significant.
              So we are stopping now.")
    }


    if (!is.numeric(modxVar)){
        print("Sorry, but I can't see how it makes sense to use a non-numeric moderator in these plots")
        warning("testSlopes: non-numeric moderator")
        return(NULL)
    }

    plotx <- tso$pso$call$plotx

    modxSeq <- plotSeq(modxRange, length.out=100)

    depVarVals <- model.response(model.frame(model))
    modyRange <- magRange(depVarVals, 1.1)
    Tcrit <- qt(0.975, model$df)

    ivs <- attr(terms(model), "term.labels")
    relevantInteractions <- c(paste(plotx, ":", modx, sep = ""),
                              paste(modx, ":", plotx, sep = ""))
    interactionsIn <- relevantInteractions[which(relevantInteractions %in% ivs)]
    if (length(interactionsIn) > 1) stop("sorry, I haven't figured out what do to do if interactionsIn > 1")

    bs <- coef(model)
    V <- vcov(model)

    bsslope <-  bs[plotx] + bs[interactionsIn]*modxSeq
    bsse <- sqrt(V[plotx,plotx] + modxSeq^2 * V[interactionsIn,interactionsIn] + 2 * modxSeq * V[plotx, interactionsIn])

    ## MM marginal matrix, similar to return of the predict() function
    ## This makes reasonable looking confidence intervals
    MM <- cbind(fit = bsslope,
                lwr = bsslope - Tcrit * bsse,
                upr = bsslope + Tcrit * bsse,
                modxSeq = modxSeq, se = bsse,
                p = 2*pt(abs(bsslope/bsse),
                lower.tail = FALSE, df = model$df))


    ylab <- substitute("Marginal Effect of" ~~ AA: ~~ "(" * hat(b)[AA] + hat(b)[BB:AA]*BB[i] *")", list(AA=plotx, BB=modx))
    ylim <- magRange(range(c(MM[ , "lwr"],MM[, "upr"])), c(1.15, 1))

    op1 <- par(no.readonly = TRUE)
    par(mar = par("mar") + c(0,0.2,0,0))

    plot(fit ~ modxSeq, data = MM, type="l", xlab = paste("The Moderator: ", modx), ylim = ylim, ylab = ylab)
    par <- op1

    abline(h = 0, col = gray(.80))
    lines(upr ~ modxSeq, data = MM, lty = 2, col = gray(.50))
    lines(lwr ~ modxSeq, data = MM, lty = 2, col = gray(.50))


    for (i in intervals){
        modxSeq <- seq(idxStart[i], idxEnd[i],
                   length.out = as.integer(40*(idxEnd[i] - idxStart[i])/diff(modxRange)))
        bsslope <-  bs[plotx] + bs[interactionsIn]*modxSeq
        bsse <- sqrt(V[plotx,plotx] + modxSeq^2 * V[interactionsIn,interactionsIn] + 2 * modxSeq * V[plotx, interactionsIn])

        ## MM marginal matrix, similar to return of the predict() function
        ## This makes reasonable looking confidence intervals
        MMsm <- cbind(fit = bsslope,
                       lwr = bsslope - Tcrit * bsse,
                       upr = bsslope + Tcrit * bsse,
                       modxSeq = modxSeq, se = bsse,
                       p = 2*pt(abs(bsslope/bsse),
                       lower.tail = FALSE, df = model$df))

        arrows(x0 = MMsm[ ,"modxSeq"], y0 = MMsm[ ,"lwr"],
               x1 = MMsm[ ,"modxSeq"], y1 = MMsm[ ,"upr"],
               angle = 90, length = 0.05, code = 3, col = gray(.70))

        if (shade == TRUE){
            polygon(x = c(MMsm[ ,"modxSeq"], rev(MMsm[ ,"modxSeq"])),
                    y = c(MMsm[ ,"upr"], rev(MMsm[ , "lwr"]) ),
                    col = col, border = gray(.80))
        }
    }

    ## draw mtext markers for internal roots
    for (i in markAt){
        lines(x = rep(i, 2),
              y = c(0, ylim[1]), lty=4, col = gray(.70))
        mtext(text = round(i, 2),
              at = i,  side = 1, line = 2, col = rgb(1, 0, 0, 0.70))
    }


    legend("topleft", legend = c("Marginal Effect", "95% Conf. Int."),
           lty = c(1, 2), col = c(1, gray(.50)), bg = "white")

    if (shade == TRUE){
        legend("bottomright", title = "Shaded Region: Null Hypothesis",
               legend = substitute(b[AA] + b[BB:AA]*BB[i] == 0 ~~ "rejected", list(AA=plotx, BB=modx)),
               fill = c(rgb(1,0,0, 0.10)), bg = "white")
    }






    if (0){
        ## Now draw the quadratic equation, to show the solution of slope/se - Tcrit = 0

        ##print("A plot of the quadratic equation should appear now, just for fun")
        if (tso$jn$a > 0) {
            xps <- plotSeq(magRange(tso$jn$roots, 1.3), length.out=100)
            ## cat(paste("Values of modx OUTSIDE this interval:\n"))
        } else {
            xps <- plotSeq(magRange(modxVar, 1.5), length.out=100)
        }

        y <- tso$jn$a * xps^2 + tso$jn$b*xps + tso$jn$c
        plot(xps, y, type="l", xlab=paste("The Moderator:", modx), ylab="slope/se - Tcrit")
        if( !is.null(tso$jn$roots) ){
            if(tso$jn$a < 0 ){
                arrows( tso$jn$roots[1], 0, tso$jn$roots[2], 0, col = "red", angle = 90, lwd = 3, code = 3, length = 0.1)
                text( mean(range(xps)), range(y)[1], pos = 3, label = expression(paste((b[plotx] + b[modx:plotx]*modx)*plotx, " is significant in the red zone")))
            } else {
                arrows(min(xps), 0, tso$jn$roots[1], 0, col = "red", angle = 90, lwd = 3, code = 2, length = 0.1)
                arrows(tso$jn$roots[2], 0, max(xps), 0, col = "red", angle = 90, lwd = 3, code = 1, length = 0.1)
                text( mean(range(xps)), range(y)[2], pos = 1,
                     label = expression(paste((b[plotx] + b[modx:plotx]*modx)*plotx,
                         " is significant in the red zone")))
            }
        }
        abline(h = 0, col = "gray80")
    }

    invisible(MM)
}


