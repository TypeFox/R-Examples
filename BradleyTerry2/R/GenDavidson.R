GenDavidson <- function(win, # TRUE/FALSE
                        tie, # TRUE/FALSE
                        loss, # TRUE/FALSE
                        player1, # player1 in each contest
                        player2, # ditto player2
                        home.adv = NULL,
                        tie.max = ~1,
                        tie.mode = NULL,
                        tie.scale = NULL,
                        at.home1 = NULL,
                        at.home2 = NULL){
    call <- as.expression(sys.call()[c(1,5:6)])
    extra <- NULL
    if (is.null(tie.max)) stop("a formula must be specified for tie.max")
    if (!is.null(home.adv) & is.null(at.home1))
        stop("at.home1 and at.home2 must be specified")
    has.home.adv <- !is.null(home.adv)
    has.tie.mode <- !is.null(tie.mode)
    has.tie.scale <- !is.null(tie.scale)
    if (has.home.adv) extra <- c(extra, list(home.adv = home.adv))
    if (has.tie.mode) extra <- c(extra, list(tie.mode = tie.mode))
    if (has.tie.scale) extra <- c(extra, list(tie.scale = tie.scale))
    i <- has.home.adv + has.tie.mode + has.tie.scale
    a <- match("home.adv", names(extra), 1)
    b <- match("tie.mode", names(extra), 1)
    c <- match("tie.scale", names(extra), 1)
    adv <- has.home.adv | has.tie.mode
    list(predictors = {c(extra,
                         list(tie.max = tie.max,
                              substitute(player1), # player1 & 2 are homogeneous
                              substitute(player2)))},
         ## substitutes "result" for "outcome", but also substitutes all of code vector
         variables = {c(list(loss = substitute(loss),
                             tie = substitute(tie),
                             win = substitute(win)),
                        list(at.home1 = substitute(at.home1),
                             at.home2 = substitute(at.home2))[adv])},
         common =  c(1[has.home.adv], 2[has.tie.mode], 3[has.tie.scale], 4, 5, 5),
         term = function(predLabels, varLabels){
             if (has.home.adv) {
                 ability1 <- paste("(", predLabels[a], ") * ", varLabels[4],
                                   " + ", predLabels[i + 2], sep = "")
                 ability2 <- paste("(", predLabels[a], ") * ", varLabels[5],
                                   " + ", predLabels[i + 3], sep = "")
             }
             else {
                 ability1 <- predLabels[i + 2]
                 ability2 <- predLabels[i + 3]
             }
             tie.scale <- ifelse(has.tie.scale, predLabels[c], 0)
             scale <- paste("exp(", tie.scale, ")", sep = "")
             if (has.tie.mode) {
                 psi1 <- paste("exp((", predLabels[b], ") * ",  varLabels[4],
                               ")", sep = "")
                 psi2 <- paste("exp((", predLabels[b], ") * ",  varLabels[5],
                               ")", sep = "")
                 weight1 <- paste(psi1, "/(", psi1, " + ", psi2, ")", sep = "")
                 weight2 <- paste(psi2, "/(", psi1, " + ", psi2, ")", sep = "")
             }
             else {
                 weight1 <- weight2 <- "0.5"
             }
             nu <- paste(predLabels[i + 1], " - ", scale, " * (",
                         weight1, " * log(", weight1, ") + ",
                         weight2, " * log(", weight2, "))", sep = "")
             paste(varLabels[1], " * (", ability2, ") + ",
                   varLabels[2], " * (", nu, " + ",
                   scale, " * ", weight1, " * (", ability1, ") + ",
                   scale, " * ", weight2, " * (", ability2, ") + ",
                   "(1 - ", scale, ") * ",
                   "log(exp(", ability1, ") + exp(", ability2, "))) + ",
                   varLabels[3], " * (", ability1, ")", sep = "")
         },
         start = function(theta) {
             init <- runif(length(theta)) - 0.5
             init[c] <- 0.5
         }
         )
}
class(GenDavidson) <- "nonlin"
