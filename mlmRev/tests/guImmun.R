library(mlmRev)
options(digits=6, useFancyQuotes = FALSE)# signif.stars for once..
fm <- glmer(immun ~ kid2p + mom25p + ord + ethn + momEd +
            husEd + momWork + rural + pcInd81 + (1|mom) + (1|comm),
            data = guImmun, family = binomial)
print(fm, symbolic.cor = TRUE)

fm.h <- update(fm, ~ . - husEd)
print(fm.h, corr = FALSE)
fm.ho <- update(fm.h, ~ . - ord)
## FIXME: shows 53 outer iterations (+ probably IRLS ones) --
##        but no such info is kept stored
print(fm.ho, corr = FALSE)

anova(fm, fm.h, fm.ho)

(fm.hoe <- update(fm.ho, ~ . - ethn))

(fm.hoem <- update(fm.hoe, ~ . - mom25p))

(AN <- anova(fm, fm.h, fm.ho, fm.hoe, fm.hoem))

AN[, "logLik"] + 1362                   # an inversion in the first two models
## FIXME: AN doesn't have a deviance column!
## AN[, "deviance"] - 2711                 # deviance scale shows this more clearly
stopifnot(AN[,"Df"] == c(9,10,12,15,18),
#          all.equal(AN[,"logLik"] + 1362,
#                    c(0.6072186497422, 0.6289103306312, 0.8541186984307,
#                      2.725550814599, 6.299084917162), tol = 1e-6),
#          all.equal(fixef(fm.hoem)[-1],
#                    c("kid2pY" = 1.2662536,  "momEdP"= 0.35116180,
#                      "momEdS"= 0.3487824136, "momWorkY"=0.2672759992340,
#                      "ruralY"=-0.678846606719, "pcInd81"=-0.9612710104134),
#                    tol = 1e-4),
          TRUE
          )


cat('Time elapsed: ', proc.time(),'\n') # "stats"
