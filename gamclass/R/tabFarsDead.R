tabFarsDead <-
function (restrict = "age>=16&age<998&inimpact%in%c(11,12,1)",
              fatal = 4, statistics = c("airbagAvail", "airbagDeploy",
                         "Restraint"))
{
    data('FARS', package='gamclass', envir=environment())
    FARS <- get("FARS", envir=environment())
    vars2 <- paste("D_", statistics, sep = "")
    yrs <- sort(unique(FARS[,"year"]))
    tabAa <- tabAd <- tabRe <- array(0, c(length(yrs), 2, 4),
                                     dimnames = list(years = yrs, D_airbagAvail = levels(FARS[,'D_airbagAvail'])[1:2],
                                     injury = c("P_injury", "D_injury", "tot", "prop")))
    names(tabAd)[2] <- vars2[2]
    names(tabRe)[2] <- vars2[3]
    i <- 0
    for (yr in yrs) {
        i <- i + 1
        yrdat <- subset(FARS, FARS[,"year"]==yr)
        yrdat <- yrdat[eval(parse(text = restrict), yrdat), ]
        subdat <- subset(yrdat[, c(vars2[1], "injury", "D_injury")],
                         yrdat[, statistics[1]] == "no")
        if (yr %in% c(2009, 2010))
            tabAa[i, , 1:2] <- NA
        else tabAa[i, , 1:2] <- as.matrix(aggregate(subdat[,
                                                           c("injury", "D_injury")], by = list(subdat[, vars2[1]]),
                                                    FUN = function(x) sum(x %in% fatal))[-3, -1])
        tabAa[i, , 3] <- tabAa[i, , 1] + tabAa[i, , 2]
        tabAa[i, , 4] <- tabAa[i, , 2]/tabAa[i, , 1]
        subdat <- subset(yrdat[, c(vars2[2], "injury", "D_injury")],
                         yrdat[, statistics[2]] == "no")
        tabAd[i, , 1:2] <- as.matrix(aggregate(subdat[, c("injury",
                                                          "D_injury")], by = list(subdat[, vars2[2]]), FUN = function(x) sum(x %in%
                                                                                                       fatal))[-3, -1])
        tabAd[i, , 3] <- tabAd[i, , 1] + tabAd[i, , 2]
        tabAd[i, , 4] <- tabAd[i, , 2]/tabAd[i, , 1]
        subdat <- subset(yrdat[, c(vars2[3], "injury", "D_injury")],
                         yrdat[, statistics[3]] == "no")
        tabRe[i, , 1:2] <- as.matrix(aggregate(subdat[, c("injury",
                                                          "D_injury")], by = list(subdat[, vars2[3]]), FUN = function(x) sum(x %in%
                                                                                                       fatal))[-3, -1])
        tabRe[i, , 3] <- tabRe[i, , 1] + tabRe[i, , 2]
        tabRe[i, , 4] <- tabRe[i, , 2]/tabRe[i, , 1]
    }
    df <- data.frame(years = yrs, airbagAvail = tabAa[, 2, 4]/tabAa[,
                                  1, 4], airbagDeploy = tabAd[, 2, 4]/tabAd[, 1, 4], Restraint = tabRe[,
                                                                                     2, 4]/tabRe[, 1, 4])
    invisible(list(airbagAvail = tabAa, airbagDeploy = tabAd,
                   restraint = tabRe))
}
