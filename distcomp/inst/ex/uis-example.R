library(distcomp)

## First we do the conventional thing.

uis <- readRDS("uis.RDS")
coxOrig <- coxph(formula=Surv(time, censor) ~ age + becktota + ndrugfp1 + ndrugfp2 +
                     ivhx3 + race + treat + strata(site), data=uis)
summary(coxOrig)

## We define the computation

coxDef <- data.frame(compType = names(availableComputations())[1],
                     formula = "Surv(time, censor) ~ age + becktota + ndrugfp1 + ndrugfp2 + ivhx3 + race + treat",
                     id = "UIS",
                     stringsAsFactors=FALSE)

library(opencpu)
## We split the data by site
siteData <- with(uis, split(x=uis, f=site))
nSites <- length(siteData)
sites <- lapply(seq.int(nSites),
                function(x) list(name = paste0("site", x),
                                 url = opencpu$url()))

ok <- Map(uploadNewComputation, sites,
          lapply(seq.int(nSites), function(i) coxDef),
          siteData)

stopifnot(all(as.logical(ok)))

master <- CoxMaster$new(defnId = coxDef$id, formula=coxDef$formula)

for (site in sites) {
  master$addSite(name = site$name, url = site$url)
}

result <- master$run()

master$summary()

print(master$summary(), digits=5)

sessionInfo()
