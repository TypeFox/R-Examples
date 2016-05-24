simulateLatinsquare.fnc = function (dat, with = TRUE, trial = 0, 
                                    nruns = 100, nsub = NA, nitem = NA, ...) 
{
    stop("this function is not working due to changes in lme4.
          an update using lmerTest is in progress")
    require("lme4", quietly = TRUE)
    if (with) {
        model.lmer = lmer(RT ~ SOA + (1 | Subject) + (1 | Word), 
            data = dat)
    }
    else {
        model.lmer = lmer(RT ~ 1 + (1 | Subject) + (1 | Word), 
            data = dat)
    }
    res = matrix(0, nruns, 3)
    if (is.na(nsub)) {
        nsub = nlevels(dat$Subject)
    }
    else {
        if (nsub%%3 != 0) {
            cat("number of subjects should be multiple of 3\n")
            return(-1)
        }
    }
    if (is.na(nitem)) {
        nitem = nlevels(dat$Word)
        if (nitem%%3 != 0) {
            cat("number of items should be multiple of 3\n")
            return(-1)
        }
    }
    for (run in 1:nruns) {
        subjects = data.frame(Subject = paste("S", 1:nsub, sep = ""), 
            rSubIntercept = rnorm(nsub, 0, sd(unlist(ranef(model.lmer)$Subject))))
        items = data.frame(Word = paste("W", 1:nitem, sep = ""), 
            rItemIntercept = rnorm(nitem, 0, sd(unlist(ranef(model.lmer)$Word))))
        simdat = expand.grid(Subject = paste("S", 1:nsub, sep = ""), 
            Word = paste("W", 1:nitem, sep = ""), SOA = c("short", 
                "medium", "long"))
        simdat$SOA = relevel(simdat$SOA, ref = "long")
        itemsvec = as.character(simdat$Word)
        simdat$ItemCnt = as.numeric(substr(itemsvec, 2, nchar(itemsvec)))
        subjectsvec = as.character(simdat$Subject)
        simdat$SubjectCnt = as.numeric(substr(subjectsvec, 2, 
            nchar(subjectsvec)))
        g1 = simdat[simdat$SOA == "short" & simdat$ItemCnt %in% 
            1:(nitem/3) & simdat$SubjectCnt %in% 1:(nsub/3), 
            ]
        g1$List = rep("L1", nrow(g1))
        g2 = simdat[simdat$SOA == "medium" & simdat$ItemCnt %in% 
            ((nitem/3) + 1):(2 * nitem/3) & simdat$SubjectCnt %in% 
            1:(nsub/3), ]
        g2$List = rep("L2", nrow(g1))
        g3 = simdat[simdat$SOA == "long" & simdat$ItemCnt %in% 
            ((2 * nitem/3) + 1):(3 * nitem/3) & simdat$SubjectCnt %in% 
            1:(nsub/3), ]
        g3$List = rep("L3", nrow(g3))
        g4 = simdat[simdat$SOA == "medium" & simdat$ItemCnt %in% 
            1:(nitem/3) & simdat$SubjectCnt %in% ((nsub/3) + 
            1):(2 * nsub/3), ]
        g4$List = rep("L1", nrow(g4))
        g5 = simdat[simdat$SOA == "long" & simdat$ItemCnt %in% 
            ((nitem/3) + 1):(2 * nitem/3) & simdat$SubjectCnt %in% 
            ((nsub/3) + 1):(2 * nsub/3), ]
        g5$List = rep("L2", nrow(g5))
        g6 = simdat[simdat$SOA == "short" & simdat$ItemCnt %in% 
            ((2 * nitem/3) + 1):(3 * nitem/3) & simdat$SubjectCnt %in% 
            ((nsub/3) + 1):(2 * nsub/3), ]
        g6$List = rep("L3", nrow(g6))
        g7 = simdat[simdat$SOA == "long" & simdat$ItemCnt %in% 
            1:(nitem/3) & simdat$SubjectCnt %in% ((2 * nsub/3) + 
            1):(3 * nsub/3), ]
        g7$List = rep("L1", nrow(g1))
        g8 = simdat[simdat$SOA == "short" & simdat$ItemCnt %in% 
            ((nitem/3) + 1):(2 * nitem/3) & simdat$SubjectCnt %in% 
            (2 * (nsub/3) + 1):(3 * nsub/3), ]
        g8$List = rep("L2", nrow(g8))
        g9 = simdat[simdat$SOA == "medium" & simdat$ItemCnt %in% 
            ((2 * nitem/3) + 1):(3 * nitem/3) & simdat$SubjectCnt %in% 
            ((2 * nsub/3) + 1):(3 * nsub/3), ]
        g9$List = rep("L3", nrow(g9))
        simdat = rbind(g1, g2, g3, g4, g5, g6, g7, g8, g9)
        simdat$List = as.factor(simdat$List)
        rownames(simdat) = 1:nrow(simdat)
        simdat$Group = as.factor(paste(rep("G", nrow(simdat)), 
            rep(1:3, rep(nrow(simdat)/3, 3)), sep = ""))
        simdat = simdat[order(simdat$Subject), ]
        n = nrow(simdat)/nsub
        trials = sample(1:n)
        for (i in 2:nsub) {
            trials = c(trials, sample(1:n))
        }
        simdat$Trial = trials
        simdat$population = fixef(model.lmer)["(Intercept)"]
        if (with) {
            simdat[simdat$SOA == "medium", ]$population = simdat[simdat$SOA == 
                "medium", ]$population + fixef(model.lmer)["SOAmedium"]
            simdat[simdat$SOA == "short", ]$population = simdat[simdat$SOA == 
                "short", ]$population + fixef(model.lmer)["SOAshort"]
        }
        simdat = merge(simdat, subjects, by.x = "Subject", by.y = "Subject")
        simdat = merge(simdat, items, by.x = "Word", by.y = "Word")
        simdat$error = rnorm(nrow(simdat), 0, sd(resid(model.lmer)))
        simdat$RTsim = apply(simdat[, 9:12], 1, sum)
        if (trial > 0) {
            simdat$RTsim = simdat$RTsim + trial * simdat$Trial
            sim.lmer = lmer(RTsim ~ Trial + SOA + (1 | Subject) + 
                (1 | Word), data = simdat)
        }
        else {
            sim.lmer = lmer(RTsim ~ SOA + (1 | Subject) + (1 | 
                Word), data = simdat)
        }
        # mcmc = mcmcsamp(sim.lmer, n = 10000)
        # aov = aovlmer.fnc(sim.lmer, as.data.frame(as.matrix(mcmc)), 
        #     c("SOAmedium", "SOAshort"))
        res[run, 2] = NA #aov$p
        res[run, 3] = subjects.latinsquare.fnc(simdat)$p
        cat(".")
    }
    cat("\n")
    res = data.frame(res)
    res = res[, -1]
    colnames(res) = c("MCMC", "F1")
    alpha05 = apply(res < 0.05, 2, sum)/nrow(res)
    alpha01 = apply(res < 0.01, 2, sum)/nrow(res)
    return(list(alpha05 = alpha05, alpha01 = alpha01, res = res, 
        with = with))
}
