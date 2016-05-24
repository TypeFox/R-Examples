##' Model formula construction
##' 
##' Interactive prompt for constructing model formulas.
##'
##' All of the staistical models in the \pkg{games} package require the
##' specification of multiple model formulas, as each player's utility is a
##' function of potentially different regressors.  The number of equations to
##' specify ranges from two in \code{\link{ultimatum}} to eight in
##' \code{\link{egame123}}.  \code{makeFormulas} is an interactive tool to
##' simplify the specification process.
##'
##' To use \code{makeFormulas}, specify the model you want to fit (\code{model})
##' and descriptions of the outcomes of the game (\code{outcomes}).  The order
##' of the descriptions in \code{outcomes} should match the numbering in the
##' game tree in the help page for \code{model}.  For example, with
##' \code{\link{egame122}}, the order is:
##' \enumerate{
##' \item Player 1 moves Left, Player 2 moves Left
##' \item Player 1 moves Left, Player 2 moves Right
##' \item Player 1 moves Right, Player 2 moves Left
##' \item Player 1 moves Right, Player 2 moves Right}
##' If the dependent variable in the dataset (\code{dat}) is a factor (\code{y})
##' whose levels contain the descriptions, then either \code{outcomes = dat$y}
##' or \code{outcomes = levels(dat$y)} will work.
##'
##' As an example, consider the following use of \code{\link{egame122}}.  Player
##' 1 is the legislature, which can propose budget cuts (left on the game tree)
##' or increases (right).  Player 2 is the president, who can sign or veto the
##' legislature's budget proposal (left and right respectively).  The variables
##' of interest are the president's party (\code{presparty}), the legislature's
##' party (\code{legparty}), and the year's percentage GDP growth
##' (\code{growth}).  To construct the formulas for this case, run
##' \code{makeFormulas(egame122, outcomes = c("budget cuts passed",
##' "budget cuts vetoed", "budget increase passed", "budget increase vetoed"))}.
##' The first set of options that appears is
##' \preformatted{Equation for player 1's utility from budget cuts passed: 
##' 
##' 1: fix to 0
##' 2: intercept only
##' 3: regressors, no intercept
##' 4: regressors with intercept}
##' To specify this utility as a function of a constant, the legislature's
##' party, and GDP growth, select option \code{4} and enter \code{legparty growth}
##' at the prompt.  \code{makeFormulas} will then move on to ask about Player
##' 1's utility for the other three outcomes, followed by Player 2's utility for
##' the outcomes for which her utility is not fixed to 0 (see
##' \code{\link{egame122}}).  See "Examples" below for a full example of the
##' input and constructed formula in this case.
##'
##' It is \strong{not} necessary to use \code{makeFormulas} to specify model
##' formulas.  See the help file of each model for examples of "manually" making
##' the formula.
##' @param model name of the model (must be from the \pkg{games} package) for
##' which to make a formula.
##' @param outcomes character vector with descriptions of the possible outcomes
##' of the game (see "Details" for a more precise explanation)
##' @return An object of class \code{"Formula"}, typically to be used as the
##' \code{formulas} argument in a statistical model from the \pkg{games}
##' package.
##' @seealso \code{\link{Formula}} (and the \pkg{Formula} package generally) for
##' the details of how \pkg{games} deals with multiple model formulas.
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
##' @example inst/examples/makeFormulas.r
makeFormulas <- function(model, outcomes)
{
    ## convert 'model' to character if supplied without quotation marks
    ischar <- tryCatch(is.character(model) && length(model) == TRUE, error =
                       identity)
    if (inherits(ischar, "error"))
        ischar <- FALSE
    if (!ischar)
        model <- deparse(substitute(model))

    ## extract outcome names from 'outcomes' if not directly supplied
    if (missing(outcomes) || !is.character(outcomes))
        outcomes <- levels(as.factor(outcomes))

    ## which players can have utilities in which equations
    eqByPlayer <- switch(model,
                         egame12 = list(rep(TRUE, 3), c(FALSE, FALSE, TRUE)),
                         egame122 = list(rep(TRUE, 4), c(FALSE, TRUE, FALSE,
                         TRUE)),
                         egame123 = list(rep(TRUE, 4), c(FALSE, rep(TRUE, 3)),
                         c(rep(FALSE, 3), TRUE)),
                         ultimatum = list(c(TRUE, FALSE), c(FALSE, TRUE)))
    if (is.null(eqByPlayer))
        stop("'", model, "' is not a model in the 'games' package; see 'help(package=\"games\")'")

    ## set of equations relevant to identification for a player's utilities
    ## (i.e. no equations for actions before their move)
    idByPlayer <- switch(model,
                         egame12 = list(1:3, 2:3),
                         egame122 = list(1:4, 1:4),
                         egame123 = list(1:4, 2:4, 3:4),
                         ultimatum = list(1:2, 1:2))  # not actually relevant
                                                      # for ultimatum

    nplayers <- length(eqByPlayer)
    neqs <- length(eqByPlayer[[1]])
    ans <- vector("list", nplayers)

    choices <- c("fix to 0", "intercept only", "regressors, no intercept",
                 "regressors with intercept")

    for (i in seq_len(nplayers)) {
        ## storing selections so options with intercepts can be removed when
        ## necessary for identification.  the last is set to 2 to avoid making
        ## 'allowCond' (below) even more convoluted
        sel <- c(rep(1, neqs - 1), 2)
        noconstant <- rep(FALSE, neqs)
        ans[[i]] <- vector("list", neqs)
        
        for (j in seq_len(neqs)) {
            ## skip if player i doesn't have utility for outcome j (either fixed
            ## to 0 or occurs before i's move)
            if (!eqByPlayer[[i]][j])
                next

            if (model == "ultimatum") {
                title <- paste("\n---\nPlayer ", i, "'s reservation value:", sep = "")
            } else {
                title <- paste("\n---\nEquation for player ", i, "'s utility from ",
                               outcomes[j], ":", sep = "")
            }

            ## disallow options with intercepts if this is the last of the
            ## player's identification-relevant equations and all previous
            ## equations contain an intercept
            allowCond <- j < neqs || !all(sel[idByPlayer[[i]]] %in% c(2, 4))
            allowed <- if (allowCond) 1:4 else c(1,3)
            eqtype <- menu(choices[allowed], title = title)

            if (eqtype == 0L) {
                stop("stopped by user")
            } else if (eqtype == 1L) {
                ans[[i]][[j]] <- "0"
            } else if (eqtype == 2L && allowCond) {
                ans[[i]][[j]] <- "1"
            } else {
                if (eqtype == 3L || !allowCond)
                    noconstant[j] <- TRUE

                ## keep asking for variable names until a valid list is supplied
                ## (nothing that would cause identification problems)
                repeat {
                    varnames <-
                        readline("\nEnter variable names (separated by spaces):\n")
                    varnames <- str_split(varnames, " ")[[1]]
                    varnames <- varnames[varnames != ""]
                    varnames <- varnames[!duplicated(varnames)]
                    ans[[i]][[j]] <- varnames
                    bad <- do.call(intersectAll, ans[[i]][idByPlayer[[i]]])
                    if (j == neqs && length(bad)) {
                        cat("The following variables cannot be used due to identification problems: ",
                            paste(bad, collapse = ", "), ". Try again.\n", sep =
                            "")
                    } else {
                        break
                    }
                }
            }

            sel[j] <- eqtype
        }

        ## need to close the old loop and start a new one before combining
        ## regressors into equation form, else the identification check won't
        ## work
        for (j in seq_len(neqs)) {
            if (length(ans[[i]][[j]]) >= 1) {
                ans[[i]][[j]] <- paste(ans[[i]][[j]], collapse = " + ")
                if (noconstant[j])
                    ans[[i]][[j]] <- paste(ans[[i]][[j]], "- 1")
            }
        }
    }

    ## get the dependent variable name
    if (model == "ultimatum") {
        yname <-
            readline("\n---\nWhat is the variable name for the level of offer made?\n")
        aname <-
            readline("\nWhat is the variable name for the offer acceptance indicator? (Leave blank if none will be used.)\n")
        if (str_trim(aname) != "")
            yname <- paste(yname, "+", aname)
    } else {
        yname <-
            readline("\n---\nWhat is the name of the dependent variable in the dataset? (If stored as action indicators/dummies, separate their names with spaces.)\n")
        yname <- str_split(yname, " ")[[1]]
        yname <- yname[yname != ""]
        yname <- paste(yname, collapse = " + ")
    }

    ## combine the strings obtained above and convert them to a Formula object
    ans <- unlist(ans)
    ans <- paste(ans, collapse = " | ")
    ans <- paste(yname, "~", ans)
    ans <- as.Formula(ans)

    return(ans)
}
