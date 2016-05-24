## setup sleepstudy data
require(robustbase)

data(sleepstudy, package="lme4")
sleepstudy$intercept <- 1
sleepstudy$slope <- 1
for (Subj in levels(sleepstudy$Subject)) {
    idx <- sleepstudy$Subject == Subj
    fit <- lmrob(Reaction ~ Days, sleepstudy, subset=idx, setting="KS2011")
    sleepstudy[idx, "intercept"] <- coef(fit)[1]
    sleepstudy[idx, "slope"] <- coef(fit)[2]
}
sleepstudy <- within(sleepstudy, {
    ## center
    ## DaysC <- Days - 4.5
    ## reorder by increasing mean value
    ## Subject <- reorder(Subject, Reaction)
    ## reorder Subject by intercept
    Subject <- reorder(Subject, intercept)
    attr(Subject, "scores") <- NULL
    intercept <- NULL
    slope <- NULL
})

## fit using heavy
## require(heavy)
## hfm <- heavyLme(Reaction ~ Days, random =~ Days, groups =~ Subject, sleepstudy)

## compare ranef:
## plot(ranef(robust)[[1]][,1], hfm$ranef[,1])
## plot(ranef(robust)[[1]][,1], ranef(classical)[[1]][,1])
## plot(ranef(robust)[[1]][,2], hfm$ranef[,2])
## plot(ranef(robust)[[1]][,2], ranef(classical)[[1]][,2])
