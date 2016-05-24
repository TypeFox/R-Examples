# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.
xp.boot.par.est <- function (bootgam.obj = NULL,
                             sd.norm = TRUE,
                             by.cov.type = FALSE,
                             abs.values = FALSE,
                             show.data = TRUE,
                             show.means = TRUE,
                             show.bias = TRUE,
                             dotpch = c(1,19),
                             labels = NULL,
                             pch.mean = "|",
                             xlab = NULL,
                             ylab = NULL,
                             col = c(rgb(.8, .5, .5), rgb(.2, .2, .7), rgb(.2,.2,.7), rgb(.6,.6,.6)),
                             ...) {
    boot.type <- "bootscm"
    bootgam.obj <- get.boot.obj(bootgam.obj, boot.type)
    if (is.null(bootgam.obj)) {
        return()
    }
    if (bootgam.obj$group.by.cov == TRUE) {
        cat ("This plot cannot be created when imported bootscm results are grouped by covariate.\nPlease re-import the bootscm results.")
        return()
    }
    if (!("par.est.first" %in% names(bootgam.obj))) {
        cat ("The required data is not available. Please check that all necessary PsN data was imported.\n")
        cat ("Note: If you've used the bootscm.import function, please set 'skip.par.est.import' to FALSE.\n\n")
        return(NULL)
    }
    if (is.null(xlab)) {
        xlab <- "Relative parameter estimate (from 1st scm-step)"
    }
    if (is.null(ylab)) {
        ylab <- "Covariate"
    }

    if (sd.norm == TRUE) {
        pl.dat <- bootgam.obj$par.est.long.norm
        bias.dat <- bootgam.obj$bias.dat.norm
    } else {
        pl.dat <- bootgam.obj$par.est.long
        bias.dat <- bootgam.obj$bias.dat
    }

    ## order by inclusion frequency
    rem <- seq(along = bootgam.obj$results.tab[,1])[bootgam.obj$failed == 1]
    cleaned.data <- bootgam.obj$results.tab
    if (length(rem)>0) {
        cleaned.data <- cleaned.data[-rem,]
    }
    incl.freq <- apply (cleaned.data, 2, sum)
    lev.ord <- names(incl.freq)[order(incl.freq)]

    abs.fun <- function (dat) {return(dat)}
    if (abs.values == TRUE) {
        abs.fun <- abs
        xlab <- paste("Absolute", xlab)
    }
    if (by.cov.type == TRUE) {
        formula <- factor(cov, levels=lev.ord) ~ abs.fun(value) | cov.type
    } else {
        formula <- factor(cov, levels=lev.ord) ~ abs.fun(value)
    }
    if (!is.null(labels)) {
        labels <- rev(labels)
        if (length(labels)==length(lev.ord)) {
            idx1 <- match(bias.dat$cov, lev.ord)
            idx2 <- match(names(incl.freq), lev.ord)
            idx3 <- match(pl.dat$cov, lev.ord)
            bias.dat$cov <- labels[idx1]
            names(incl.freq) <- labels[idx2]
            pl.dat$cov <- labels[idx3]
            lev.ord <- names(incl.freq)[order(incl.freq)]
        } else {
            cat ("Length of specified labels-vector not equal to number of covariate-parameter relationships. Returning to default.")
        }
    }
    legend <- list(text = list("Selected", cex=.75),
                   points = list(pch=dotpch[2], col=col[3], cex=1),
                   text = list("Not selected", cex=.75),
                   points = list(pch=dotpch[1], col=col[1], cex=1) )
    if (show.means == TRUE) {
        legend <- list(text = list("Selected", cex=.75),
                       points = list(pch=dotpch[2], col=col[3], cex=1),
                       text = list("Not selected", cex=.75),
                       points = list(pch=dotpch[1], col=col[1], cex=1),
                       text = list("mean (selected)", cex=.75),
                       lines = list(lwd=1.5, span=0.1, col=col[3]),
                       text = list("mean (all)", cex=.75),
                       lines = list(lwd=1.5, span=0.1, col=col[4])
                   )
    }
    p <- stripplot (formula,
                    data = pl.dat,
                    ylab = ylab,
                    xlab = xlab,
                    groups = factor(eval(as.name("incl")), levels = c("Not included", "Included")),
                    par.settings = simpleTheme (col=col, pch=dotpch),
                    key = legend,
                    levels = lev.ord,
                    panel = function (...) {
                        panel.abline (v=0, lty=3)
                        if (show.data == TRUE) {
                            panel.stripplot (jitter.data=TRUE, ...)
                        }
                        if (show.means == TRUE) {
                            panel.xyplot (y = factor(bias.dat[bias.dat$incl == "Included",]$cov, levels=lev.ord),
                                          x = abs.fun(as.num(bias.dat[bias.dat$incl == "Included",]$mean)),
                                          bias.data=bias.dat, pch = pch.mean, cex=2.5, col=col[3]
                                          )
                            panel.xyplot (y = factor(bias.dat[bias.dat$incl == "Included",]$cov, levels=lev.ord),
                                          x = abs.fun(as.num(bias.dat[bias.dat$incl == "Included",]$All)),
                                          bias.data=bias.dat, pch = pch.mean, cex=2.5, col=col[4]
                                          )
                        }
                        if (show.bias == TRUE) {
                            panel.text (y = factor(bias.dat[bias.dat$incl=="Included",]$cov, levels=lev.ord),
                                        x = abs.fun(max(pl.dat[!is.na(pl.dat$value),]$value)*0.94),
                                        labels = paste (round(bias.dat[bias.dat$incl=="Included",]$bias,0), "%", sep=""), cex=0.8)
                        }
                    }, ...)
    return(p)
}

ask.covs.plot <- function (bootgam.obj = NULL) {
    if (!is.null(bootgam.obj)) {
        cat ("Covariates in database: ")
        covs <- colnames(bootgam.obj$covariate$sd.all)
        cat (covs)
        cat ("\n\nPlot for which covariates (separate by space, return for all): ")
        ans <- readline()
        if (ans == "") {
            return()
        }
        ans.cov <- strsplit(ans, " ")[[1]]

        if (length(ans.cov) < 2) {
            cat("Please choose at least 2 covariatess from the list!\n\n")
            Recall(bootgam.obj)
        } else {
            if (sum((ans.cov %in% covs)*1) == length(ans.cov)) {
                return (ans.cov)
            } else {
                cat("Please choose covariates from the list only!\n\n")
                Recall(bootgam.obj)
            }
        }
    }
}

xp.boot.par.est.corr <- function (bootgam.obj = NULL,
                                  sd.norm = TRUE,
                                  by.cov.type = FALSE,
                                  cov.plot = NULL, # covariates to plot if not all are wanted
                                  ask.covs = FALSE,
                                  dotpch = 19,
                                  col = rgb(.2, .2, .9, .75),
                                  ...) {
    boot.type <- "bootscm"
    bootgam.obj <- get.boot.obj(bootgam.obj, boot.type)
    if (is.null(bootgam.obj)) {
        return()
    }
    if (bootgam.obj$group.by.cov == TRUE) {
        cat ("This plot cannot be created when imported bootscm results are grouped by covariate.\nPlease re-import the bootscm results.")
        return()
    }
    if (!("par.est.first" %in% names(bootgam.obj))) {
        cat ("The required data is not available. Please check that all necessary PsN data was imported.\n")
        cat ("Note: If you've used the bootscm.import function, please set 'skip.par.est.import' to FALSE.\n\n")
        return(NULL)
    }
    tmp <- bootgam.obj$par.est.first
    if (sd.norm == TRUE) { # for non-dichotomous covariates, do correction
        tmp <- bootgam.obj$par.est.first.corr
        xlab <- "Parameter estimate (from 1st scm-step), SD-normalized"
    }
    pl.dat <- tmp
    pl.dat.incl <- (!is.na(bootgam.obj$par.est.final))*1

    ## filter out the desired covariates
    if (is.null(cov.plot)) {
        if (ask.covs==TRUE) {
            cov.plot <- ask.covs.plot (bootgam.obj)
        }
    }
    if ((!is.null(cov.plot))&&(sum(cov.plot %in% colnames(tmp))>0)) {
        pl.dat <- tmp[,cov.plot]
        pl.dat.incl <- pl.dat.incl[,cov.plot]
    }

    p <- splom (pl.dat, pch = dotpch, col=col)
    return(p)
}

bootgam.print <- function(bootgam.obj = NULL) {
    bootgam.obj <- get.boot.obj(bootgam.obj, NULL)
    if (is.null(bootgam.obj)) {
        return()
    }
    boot.type <- get.boot.type (bootgam.obj)
    cat("\n********************************************************************\n")
    if (boot.type == "bootgam") {
        cat("************************* BootGAM results **************************\n")
    } else {
        cat("************************* BootSCM results **************************\n")
    }
    cat("Run number:", bootgam.obj$runno, "\n")
    failed <- NULL
    if (boot.type == "bootgam") {
        if(is.null(startm <- bootgam.obj$start.mod)) {
            cat("No start model specified.\n")
        } else {
            cat("Start model set to:", startm,"\n")
        }
        cat("Seed number:", bootgam.obj$seed,"\n")
        if(length(bootgam.obj$excluded.ids)>0) {
            cat("Excluded individuals:",bootgam.obj$excluded.ids,"\n")
        } else {
            cat("No individuals were excluded.\n")
        }
        cat("\nConvergence algorithm:", bootgam.obj$algo,"\n")
        if(bootgam.obj$algo == "fluct.ratio") {
            cat("Lowest important inclusion frequency:")
            cat("\n  Convergence criterium:", format(bootgam.obj$fluct.ratio.last, digits = 5), "(target=", bootgam.obj$conv.value, ")\n")
        } else {
            cat("Lowest absolute joint inclusion frequency:")
            cat("\n  Convergence criterium:", format(bootgam.obj$ljif.last, digits = 5), "(target=", bootgam.obj$ljif, ")\n")
        }
        failed <- seq(along=eval(as.name("current.bootgam"))$failed)[eval(as.name("current.bootgam"))$failed==1]
        cat ("Failed BootGAM replicates: ", failed, "\n")
    }
    cat("\nTotal number of iterations:", length(bootgam.obj$results.tab[,1]), "\n")
    cat("\nModel size: ")
    res <- bootgam.obj$results.tab
    if (!is.null(failed)) {
        if (length(failed)>0) {
            res <- bootgam.obj$results.tab[-failed,]
        }
    }
    print (summary(apply(res, 1, sum)))
    cat("\nInclusion probabilities:\n")
    tot.prob <- tail(bootgam.obj$incl.freq,1)
    ord <- rev(order(tot.prob))
    print(t(as.list(round(tot.prob[ord],3))))
    cat("********************************************************************\n\n")
}

check.bootgamobj <- function () {
    getit <- function() {
        cat("\nYou have to specify the parameter name and the run number",
            "of the bootgam objects you want to plot. The following",
            "bootgam objects are available:\n", fill = 60)
        if (.Platform$OS == "windows") {
            cat(objects(pattern = "bootgam.xpose*", pos = 1), fill = 60)
        }
        else {
            cat(objects(pattern = "^bootgam.xpose", pos = 1), fill = 60)
        }
        cat("\nParameter (0 to exit): ")
        ans <- readline()
        if (ans == 0) {
            return(ans <- NULL)
        }
        cat("Run number (0 to exit):")
        ans1 <- readline()
        if (ans1 == 0) {
            return(ans1 <- NULL)
        }
        gobjname <- paste("bootgam.xpose.", ans, ".", ans1, sep = "")
        if (!exists(gobjname, where = 1)) {
            cat("\n*There are no objects that matches", gobjname,
                "\n")
            gobjname <- Recall()
        }
        return(gobjname)
    }
    
    if (exists("current.bootgam", where = 1)) {
      cur.boot <- eval(as.name("current.bootgam"))
      cat("\nThe current bootgam object is for", cur.boot$parnam,
          "in run", cur.boot$runno, ".\n")
      cat("\nDo you want to proceed with this bootgam object? y(n) ")
      ans <- readline()
      if (ans != "y" && ans != "") {
        gobjname <- getit()
        if (!is.null(gobjname)) {
          c1 <- call("assign",pos = 1, "current.bootgam", eval(as.name(gobjname)),
                 immediate = T)
          eval(c1)
        }
      } else {
        gobjname <- T
      }
    } else {
        gobjname <- getit()
        if (!is.null(gobjname)) {
            c2 <- call("assign",pos = 1, "current.bootgam", eval(as.name(gobjname)),
                   immediate = T)
            eval(c2)
        }
    }
    return(gobjname)
}

ask.bootgam.bootscm.type <- function () {
    cat ("Both a bootgam and a bootscm object are available, which one\nwould you like to summarize?\n")
    cat ("  1) the current bootgam object\n")
    cat ("  2) the current bootscm object\n")
    ans <- readline()
    if (ans == "") {
        Recall()
    } else {
        if ((ans == 1)|(ans == 2)) {
            if (ans == 1) {return ("bootgam")}
            if (ans == 2) {return ("bootscm")}
        } else {
            cat("Please choose either 1 or 2!\n\n")
            Recall()
        }
    }
}

get.boot.obj <- function (bootgam.obj = NULL,
                          boot.type = NULL
                          ) {
                                        # switch between supplied object or global object, and bootscm/bootgam
    if ((is.null(boot.type))&(is.null(bootgam.obj))) {
        if (("current.bootgam" %in% ls(.GlobalEnv))&(!"current.bootscm" %in% ls(.GlobalEnv))) {
            boot.type <- "bootgam"
        }
        if (("current.bootscm" %in% ls(.GlobalEnv))&(!"current.bootgam" %in% ls(.GlobalEnv))) {
            boot.type <- "bootscm"
        }
        if (("current.bootscm" %in% ls(.GlobalEnv))&("current.bootgam" %in% ls(.GlobalEnv))) {
            boot.type <- ask.bootgam.bootscm.type()
            cat ("\n")
        }
        if (is.null(boot.type)) {
            cat ("No bootgam or bootscm object found!\n")
            return()
        }
    }
    if (is.null(boot.type)) {
        boot.type <- get.boot.type (bootgam.obj)
    }
    if (boot.type == "bootscm") {
        if (is.null(bootgam.obj)) {
            if ("current.bootscm" %in% objects(pos=1)) {
                if (!is.null(eval(as.name("current.bootscm")))) {
                    bootgam.obj <- eval(as.name("current.bootscm"))
                } else {
                    cat ("Data not available. Did you import the bootSCM data?\n")
                }
            } else {
                cat (paste(objects()))
                cat ("Data not available. Did you import the bootSCM data?\n")
            }
        } else {
            c3 <- call("assign",pos = 1, "current.bootscm", bootgam.obj, immediate = T)
            eval(c3)
        }
    } else { # load bootgam object
        if (is.null(bootgam.obj)) {
            if ("current.bootgam" %in% objects()) {
                if (!is.null(bootgam.obj)) {
                    bootgam.obj <- eval(as.name("current.bootgam"))
                }
            } else {
                if (check.bootgamobj()) {
                    bootgam.obj <- eval(as.name("current.bootgam"))
                } else {
                    cat ("Data not available. Did you run the bootGAM data?\n")
                }
            }
        } else {
            c4 <- call("assign",pos = 1, "current.bootgam", bootgam.obj, immediate = T)
            eval(c4)
        }
    }
    return(bootgam.obj)
}

xp.distr.mod.size <- function (bootgam.obj = NULL,
                               boot.type = NULL,
                               main = NULL,
                               bw = 0.5,
                               xlb = NULL,
                               ... ) {
    bootgam.obj <- get.boot.obj(bootgam.obj, boot.type)
    if (is.null(bootgam.obj)) {
        return()
    }
    boot.type <- get.boot.type (bootgam.obj)

    ## Sort out the titles
    if(is.null(main)) {
        main <- paste("Distribution of covariate model sizes", bootgam.obj$runno)
    }
    if (is.null(xlb)) {
        if (boot.type == "bootgam") {
            xlb <- paste ("Covariate model size (on", bootgam.obj$parnam, ")", sep = "")
        } else {
            xlb <- paste ("Covariate model size (on any parameter)")
        }
    }
                                        # Plot
    res <- bootgam.obj$results.tab
    if (!is.null(bootgam.obj$failed)) {
        res <- res[bootgam.obj$failed == 0,]
    }
    sizes <- apply (res, 1, sum)
    pl <- densityplot (sizes,
                       bw = bw,
                       main = main,
                       ... )
    return(pl)
}

ask.incl.range <- function (bootgam.obj = NULL) {
    text <- paste("The plots that show correlations between covariate inclusion\n",
                  "frequencies (inclusion index) are not informative when the inclusion\n",
                  "frequency for a covariate is either very high or very low. Therefore\n",
                  "it is advised to show these plots only for intermediately strong \n",
                  "covariates. The default range is 20% to 80%.\n\n", sep="")
    cat (text)
    cat ("Specify range (e.g.: 20 80): ")
    ans <- readline()
    if (ans == "") {
        range <- c(20,80)
    } else {
        range <- as.numeric(strsplit (ans, " ")[[1]])
    }
    if (length(range) == 2) {
        return (range)
    } else {
        cat("Please specify two numbers, separated by a space!\n\n")
        Recall(bootgam.obj)
    }
}

xp.incl.index.cov <- function (bootgam.obj = NULL,
                               boot.type = NULL,
                               main = NULL,
                               xlb = "Index",
                               ylb = "Covariate",
                               add.ci = FALSE,
                               incl.range = NULL,
                               ... ) {

    bootgam.obj <- get.boot.obj(bootgam.obj, boot.type)
    if (is.null(bootgam.obj)) {
        return()
    }
    boot.type <- get.boot.type (bootgam.obj)
    as.num <- function (dat) { return (as.numeric(as.character(dat))) }

    ## Sort out the titles
    if(is.null(main)) {
        main <- paste("Inclusion index for", bootgam.obj$runno)
    }

    se_idx <- function (p, q, n) {
        A <- (p/n)*(1-(p/n))/n
        B <- (q/n)*(1-(q/n))/n
        rho <- 1  # cor (A,B)
        se <- sqrt ( A + B + 2 * sqrt (A) * sqrt(B) * rho )
        return (se)
    }

    inc_obs <- tail(bootgam.obj$incl.freq,1)
    res <- bootgam.obj$results.tab

    ## filter out only covariates within specified inclusion freq range
    if (is.null(incl.range)) {
        incl.range <- ask.incl.range()
    }
    if (length(incl.range) == 2) {
        filter <- inc_obs > incl.range[1]/100 & inc_obs < incl.range[2]/100
        res <- res[,filter]
        inc_obs <- inc_obs[,filter] # remove dummy columns
    }
    n_cov <- length(inc_obs)
    nam <- names(inc_obs)

    ## filter out failed scms
    if (!is.null(bootgam.obj$failed)) {
        res <- res[bootgam.obj$failed == 0,]
    }
    if (boot.type == "bootscm") {
        cols.dum <- grep("^X.", colnames(res))
        if (length(cols.dum)>0) {
            res <- res[,-cols.dum]
        }
    }

    cov_idx <- c()
    for (i in 1:n_cov) {
        sub <- res[res[,i]==1,]
        obs <- apply (sub, 2, sum)
        n <- length(sub[,1])
        expect <- inc_obs * n
        idx <- as.num((obs/n) / as.num (inc_obs[i]) * as.num(inc_obs)) - 1
        se <- 0
        if (add.ci == TRUE) {
            se <- unlist (se_idx(p = obs, q = expect, n = length(res[,1])))
        }
        cov_idx <- data.frame (rbind (cov_idx, cbind ("COV1" = nam[i], "COV2" = nam, idx, se, "lbnd" = (idx-(1.96*se)), "ubnd"=(idx+(1.96*se)))))
    }

    p <- dotplot (as.factor(COV1) ~ as.num(idx) | as.factor(COV2),
                  data=cov_idx,
                  plot.zero=TRUE,
                  main = main,
                  xlab = xlb,
                  ylab = ylb,
                  lx = as.num(cov_idx$lbnd), ux = as.num(cov_idx$ubnd),
                  prepanel = prepanel.ci,
                  panel = panel.ci,
                  ... )
    return(p)
}

ask.cov.name <- function (bootgam.obj = NULL) {
  if (!is.null(bootgam.obj)) {
    cat ("Covariates in database: ")
    cat (paste (bootgam.obj$covnams))
    cat ("\n\nPlot for which covariate (return to exit): ")
    ans <- readline()
    if (ans == "") {
      return()
    }
    if (ans %in% (bootgam.obj$covnams)) {
      return (ans)
    } else {
      cat("Please choose a covariate from the list!\n\n")
      Recall(bootgam.obj)
    }
  }
}


xp.incl.index.cov.ind <- function (bootgam.obj = NULL,
                                   boot.type = NULL,
                                   cov.name = NULL,
                                   main = NULL,
                                   ylb = "ID",
                                   xlb = "Individual inclusion index",
                                   ... ) {
    bootgam.obj <- get.boot.obj(bootgam.obj, boot.type)
    if (is.null(bootgam.obj)) {
        return()
    }
    boot.type <- get.boot.type (bootgam.obj)

    as.num <- function (dat) { return (as.numeric(as.character(dat))) }

    if (is.null(cov.name)) {
        cov.name <- ask.cov.name(bootgam.obj)
    }
    if (is.null(cov.name)) { return() }

    if(is.null(main)) {
        main <- paste ("Individual inclusion index (", cov.name, " on ", bootgam.obj$parnam, ") for ", bootgam.obj$runno, sep="")
    }

    ids <- colnames(bootgam.obj$oid)
    oid.cnt <- apply (bootgam.obj$oid, 2, sum)
    res <- bootgam.obj$results.tab

    if (!is.null(bootgam.obj$failed)) {
        res <- res[bootgam.obj$failed == 0,]
    }
    oid.rel <- oid.cnt / length(res[,1])
    nam <- names(res)

    cov_idx <- c()
    sub <- bootgam.obj$oid[res[,cov.name == nam]==1,]
    obs <- apply (sub, 2, sum)
    n <- length(sub[,1])
    idx <- (as.num(obs) / (n * as.num(oid.rel))) - 1
    ord <- order(idx)
    ids <- as.num(gsub("X","", ids))

    cov_idx <- data.frame(cbind ("idn" = ids[ord], "idx" = as.num(idx[ord])))
    scales <- list(y = list (labels = rev(cov_idx$idn)), cex=c(0.7,1))
    p <- xyplot (factor(idn, levels=rev(idn)) ~ as.num(idx),
                 data =cov_idx,
                 main = main,
                 xlab = xlb,
                 ylab = ylb,
                 scales = scales,
                 lx = 0, ux = 0, plot.zero=TRUE,
                 prepanel = prepanel.ci,
                 panel = panel.ci,
                 ... )
    return (p)
}

xp.incl.index.cov.comp <- function (bootgam.obj = NULL,
                                    boot.type = NULL,
                                    main = NULL,
                                    xlb = "Individual inclusion index",
                                    ylb = "ID",
                                    ... ) {
    bootgam.obj <- get.boot.obj(bootgam.obj, boot.type)
    if (is.null(bootgam.obj)) {
        return()
    }
    as.num <- function (dat) { return (as.numeric(as.character(dat))) }
    boot.type <- get.boot.type (bootgam.obj)

    if(is.null(main)) {
        main <- paste ("Individual inclusion indices for", bootgam.obj$runno)
    }

    ids <- colnames(bootgam.obj$oid)
    oid.cnt <- apply (bootgam.obj$oid, 2, sum)
    res <- bootgam.obj$results.tab
    if (!is.null(bootgam.obj$failed)) {
        res <- res[bootgam.obj$failed == 0,]
    }
    oid.rel <- oid.cnt / length(res[,1])
    nam <- names(res)

    cov_idx <- c()
    for (i in seq(along=nam)) {
        sub <- bootgam.obj$oid[res[,nam[i] == nam]==1,]
        obs <- apply (sub, 2, sum)
        n <- length(sub[,1])
        idx <- (as.num(obs) / (n * as.num(oid.rel))) - 1
        cov_idx <- data.frame(rbind (cov_idx, cbind ("cov" = nam[i], "idx" = as.num(idx))))
    }

    p <- xyplot (factor(cov) ~ as.num(idx),
                 data=cov_idx,
                 xlab = xlb,
                 ylab = ylb,
                 main = main,
                 lx = 0, ux = 0,
                 plot.zero = TRUE,
                 prepanel = prepanel.ci,
                 panel = panel.ci,
                 ... )
    return (p)
}

xp.inc.prob <- function (bootgam.obj = NULL,
                         boot.type = NULL,
                         main = NULL,
                         col = "#6495ED",
                         xlb = NULL,
                         ylb = "Covariate",
                         ... ) {
    bootgam.obj <- get.boot.obj(bootgam.obj, boot.type)
    if (is.null(bootgam.obj)) {
        return()
    }
    boot.type <- get.boot.type (bootgam.obj)

    ## Sort out the titles
    if(is.null(main)) {
        main <- paste("Total frequency of covariates for", bootgam.obj$runno)
    }

    rem <- seq(along = bootgam.obj$results.tab[,1])[bootgam.obj$failed == 1]
    cleaned.data <- bootgam.obj$results.tab
    if (length(rem)>0) {
        cleaned.data <- bootgam.obj$results.tab[-rem,]
    }
    frac <- function (data) { sum (data) / length(data) }
    se <- function (data) {
        p <-  sum (data) / length(data)
        se <- p * (1-p) / length(data)
        return (se)
    }
    as.num <- function (data) { return (as.numeric(as.character(data)))}

    cov.prob <- apply (cleaned.data, 2, frac)
    cov.prob <- cov.prob[order(cov.prob)]
    cov.se <- apply (cleaned.data, 2, se)
    cov.se <- cov.se[order(cov.prob)]
    cov.ci <- cbind ("ubnd" = cov.prob + 1.96*cov.se, "lbnd" = cov.prob - 1.96*cov.se)
    cov.comb <- data.frame ( cbind ( "cov" = names(cov.prob), "prob" = cov.prob, cov.ci) )
    cov.comb <- cov.comb[order(cov.comb$prob),]

    if (is.null(xlb)) {
        xlb <- paste("Inclusion frequency (%) on ", bootgam.obj$parnam, sep="")
        if (boot.type == "bootscm") {
            xlb <- "Inclusion frequency (%)"
        }
    }

    pl <- xyplot (factor(cov, levels=cov) ~ 100*as.num(prob),
                  lx = as.num(cov.comb$lbnd), ux = as.num(cov.comb$ubnd),
                  data = cov.comb,
                  prepanel = prepanel.ci,
                  panel = panel.ci,
                  main = main,
                  xlim = c(0,100),
                  xlab = xlb,
                  ylab = ylb,
                  ... )
    return(pl)
}

xp.inc.prob.comb.2 <- function (bootgam.obj = NULL,
                                boot.type = NULL,
                                main = NULL,
                                col = "#6495ED",
                                xlb = NULL,
                                ylb = "Covariate combination",
                                ... ) {
    bootgam.obj <- get.boot.obj(bootgam.obj, boot.type)
    if (is.null(bootgam.obj)) {
        return()
    }
    if(is.null(main)) {
        main <- paste("Most common 2-covariate combinations for", bootgam.obj$runno)
    }
    boot.type <- get.boot.type (bootgam.obj)

    rem <- seq(along = bootgam.obj$results.tab[,1])[bootgam.obj$failed == 1]
    cleaned.data <- bootgam.obj$results.tab
    if (length(rem) > 0) {
        cleaned.data <- cleaned.data[-rem,]
    }

    frac <- function (data) { sum (data) / length(data) }
    se <- function (data) {
        p <-  sum (data) / length(data)
        se <- p * (1-p) / length(data)
        return (se)
    }
    as.num <- function (data) { return (as.numeric(as.character(data)))}

    covs <- colnames(cleaned.data)
    cov_all <- c()
    for (i in seq(along=covs)) {
        tmp <- cleaned.data[cleaned.data[,i] == 1, -i]
        cov.prob <- apply (tmp, 2, frac)
        cov_all <- data.frame (rbind (cov_all, cbind ("cov1"=covs[i], "cov2" = names(cov.prob), "idx" = as.num(cov.prob))))
    }
    cov_all$idx <- as.num(cov_all$idx)
    cov_all_10 <- head(cov_all[order(cov_all$idx, decreasing=TRUE),], 10)
    cov_all_10$label <- paste(cov_all_10$cov1, "+", cov_all_10$cov2)

    if (is.null(xlb)) {
        xlb <- paste("Inclusion frequency (%) on ", bootgam.obj$parnam, sep="")
        if (boot.type == "bootscm") {
            xlb <- "Inclusion frequency (%) on any parameter)"
        }
    }
    pl <- dotplot (factor(label, levels=rev(cov_all_10$label)) ~ 100*as.num(idx),
                   lx = 0, ux=0,
                   data = cov_all_10,
                   prepanel = prepanel.ci,
                   panel = panel.ci,
                   xlim = c(0,100),
                   main = main,
                   xlab = xlb,
                   ylab = ylb,
                   ...)
    return(pl)
}

prepanel.ci <- function(x, y, lx, ux, subscripts, ...) {
    x <- as.numeric(x)
    lx <- as.numeric(lx[subscripts])
    ux <- as.numeric(ux[subscripts])
    list(xlim = range(x, ux, lx, finite = TRUE))
}

panel.ci <- function(x, y, lx, ux, subscripts, pch = 16, plot.zero = FALSE, ...) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    lx <- as.numeric(lx[subscripts])
    ux <- as.numeric(ux[subscripts])
    if (plot.zero == TRUE) {
        panel.abline (v = 0, lty = 3, lwd = 1, col="#999999")
    }
    panel.abline(h = unique(y), col = "grey", lwd = 1)
                                        # show SE of estimate. Disabled.
                                        #    panel.arrows(lx, y, ux, y, col = 'black', lwd = 2,
                                        #                 length = 0, unit = "native",
                                        #                 angle = 90, code = 3)
    panel.xyplot(x, y, pch = pch, ...)
}

xp.inc.stab.cov <- function (bootgam.obj = NULL,
                             boot.type = NULL,
                             main = NULL,
                             xlb = "Bootstrap replicate number",
                             ylb = "Inclusion frequency",
                             ...) {
    ## Create a plot of inclusion frequency (y) versus bootstrap replicate number (x)
    bootgam.obj <- get.boot.obj(bootgam.obj, boot.type)
    if (is.null(bootgam.obj)) {
        return()
    }
    boot.type <- get.boot.type (bootgam.obj)

    if(is.null(main)) {
        main <- paste("Inclusion stability for", bootgam.obj$runno)
    }
    freq <- bootgam.obj$incl.freq
    if (!is.null(bootgam.obj$failed)) {
        freq <- freq[bootgam.obj$failed==0,]
    }
    freq <- data.frame (cbind (row = seq(along = freq[,1]), freq ))
    freq.long <- reshape (freq,
                          ids=row.names(freq), varying = names(freq)[-1],
                          idvar = "row", timevar = "var", v.names = "value",
                          times = names(freq)[-1], direction="long")
    pl <- xyplot (value ~ row | var,
                  data = freq.long,
                  main = main,
                  xlab = xlb,
                  ylab = ylb,
                  type = "l",
                  ...)
    return (pl)
}

xp.dofv.plot <- function (bootscm.obj = NULL,
                          main = NULL,
                          xlb = "Difference in OFV",
                          ylb = "Density",
                          ... ) {
    bootscm.obj <- get.boot.obj(bootscm.obj, boot.type = "bootscm")
    if (is.null(bootscm.obj)) {
        return()
    }

    ## Sort out the titles
    if(is.null(main)) {
        main <- paste("Distribution of dOFV for", bootscm.obj$runno)
    }

                                        # Plot
    dofv <- bootscm.obj$dofv[!is.na(bootscm.obj$dofv$dOFV),]$dOFV
    dofv <- dofv[-1]
    pl <- densityplot (dofv, lwd=3,
                       main=main,
                       xlab = xlb,
                       ylab = ylb,
                       panel = function () {
                           panel.abline (v=0, lty=3, col="#888888")
                           panel.densityplot (dofv)
                       },
                       ... )
    return (pl)
}

get.boot.type <- function (bootscm.obj) {
    boot.type <- "bootgam"
    if ("dofv" %in% names(bootscm.obj)) {
        boot.type <- "bootscm"
    }
    return(boot.type)
}
