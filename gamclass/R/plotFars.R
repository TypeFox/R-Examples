plotFars <-
function (restrict = "age>=16&age<998&inimpact%in%c(11,12,1)",
              fatal = 4,
              statistics = c("airbagAvail", "airbagDeploy", "Restraint"))
{
tabDeaths <- tabFarsDead(restrict = restrict,
              fatal = 4,
              statistics = statistics)
tabAa <- tabDeaths[['airbagAvail']]
tabAd <- tabDeaths[['airbagDeploy']]
tabRe <- tabDeaths[['restraint']]
yrs <- as.numeric(dimnames(tabAa)[[1]])
    df <- data.frame(years = yrs, airbagAvail = tabAa[, 2, 4]/tabAa[,
                                  1, 4], airbagDeploy = tabAd[, 2, 4]/tabAd[, 1, 4], Restraint = tabRe[,
                                                                                     2, 4]/tabRe[, 1, 4])
        form <- formula(paste(paste(statistics, collapse="+"), "~ years"))
        gph <- xyplot(form, data = df,
                      par.settings=simpleTheme(pch = 16),
                      auto.key=list(columns=length(statistics)))
gph
}
