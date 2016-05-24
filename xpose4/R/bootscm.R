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
bootscm.import <- function (scm.folder = NULL,
                            silent=FALSE,
                            n.bs = NULL,
                            cov.recoding=NULL,
                            group.by.cov=NULL,
                            skip.par.est.import=FALSE,
                            dofv.forward = 3.84,
                            dofv.backward = 6.64
                            )  {
    bootscm.obj <- list()
    cat.s <- function (txt) { if (silent == FALSE) { cat (txt) } }
    if (is.null(scm.folder)) {
        scm.folder <- ask.folder()
    }
    if (is.null(scm.folder)) {
        return()
    }
    if(is.null(group.by.cov)) {
        group.by.cov.num <- ask.group.by.cov()
        if (group.by.cov.num == "1") { group.by.cov <- FALSE} else { group.by.cov <- TRUE}
    }
    if (is.null(NULL)) {
        cat (paste("Please note that if you are manually recoding the covariates in the \n",
                   "scm config file, you will have to import the results manually using\n",
                   "the bootscm.import() function and the 'cov.recoding=...' argument.\n",
                   "See manual for more details.\n\n", sep=""))
    }
    cat.s ("Importing bootstrap...\n")
    if (is.null(n.bs)) {
        scm_dirs <- dir (scm.folder, "scm_dir\\d")
        n.bs <- length(scm_dirs)
        cat.s (paste("* Found ", n.bs, " scm folders in the bootscm folder.\n", sep=""))
    } else {
        cat.s (paste("* Trying to read data from ", n.bs, " scm's in the bootscm folder.\n", sep=""))
    }
    if (n.bs == 0) {
        cat.s ("* No scm folders were found in the bootscm folder. Please check correct execution of PsN.\n\n")
        return()
    }

    cat.s ("* Trying to import covariate inclusion data...")
    if (file.exists(paste(scm.folder,"/covariate_inclusion.csv", sep=""))) {
        cov.incl <- read.csv(paste(scm.folder,"/covariate_inclusion.csv", sep=""))
        cat.s ("Converting...")
        bootscm.obj <- convert.bootscm.cov.inclusion (cov.incl, group.by.cov = group.by.cov)
        bootscm.obj$group.by.cov <- group.by.cov
        bootscm.obj$covnams <- colnames(bootscm.obj$results.tab)
        if (0 %in% cov.incl$bs_n) {
            bootscm.obj$reestimate_final <- TRUE
            bootscm.obj$results.tab.orig <- bootscm.obj$results.tab[1,]
            bootscm.obj$results.tab <- bootscm.obj$results.tab[-1,]
        } else {
            bootscm.obj$reestimate_final <- FALSE
        }
        ## Create table of inclusion frequency
        incl.freq <- bootscm.obj$results.tab
        for (i in 1:length(bootscm.obj$results.tab[,1])) {
            incl.freq[i,] <-  apply (bootscm.obj$results.tab[1:i,], 2, function (data) { sum (as.num(data)) / length(data) })
        }
        bootscm.obj$incl.freq <- incl.freq
        cat.s ("OK.\n")
    } else {
        cat.s ("Data not found.\n")
    }

    cat.s ("* Trying to import individual inclusion data...")
    if (file.exists(paste(scm.folder,"/bs_ids.csv", sep=""))) {
        bs_ids <- read.csv(paste(scm.folder,"/bs_ids.csv", sep=""))
        cat.s ("Converting...")
        bootscm.obj$oid <- convert.bootscm.bs.ids (bs_ids)
        cat.s ("OK.\n")
    } else {
        cat.s ("Data not found.\n")
    }

    if (skip.par.est.import == FALSE) {
        cat.s ("* Trying to import covariate parameter estimates...")
        if (file.exists(paste(scm.folder, "/scm_dir1/raw_results_bsmod_1.csv", sep=""))) {
            cat.s ("Importing...")
            tmp <- read.bootscm.par.est(folder=paste(scm.folder,sep=""),
                                        n.bs = n.bs,
                                        cov.recoding=cov.recoding,
                                        dofv.forward=dofv.forward,
                                        dofv.backward=dofv.backward)
            bootscm.obj <- c(bootscm.obj, tmp)
            bootscm.obj$par.est.first.corr <- tmp$par.est.first * tmp$covariate$sd
            cat.s ("OK.\n")
        } else {
            cat.s ("Data not found.\n")
        }
    }

    cat.s ("* Trying to import final objective function values...")
    if (bootscm.obj$reestimate_final) {
        if (file.exists(paste(scm.folder,"/ofv_final.csv", sep=""))) {
            ofv_final <- read.csv(paste(scm.folder,"/ofv_final.csv", sep=""))
            ofv_final$dOFV <- 0
            ofv_final[2:length(ofv_final$OFV),]$dOFV <- ofv_final[2:length(ofv_final$OFV),]$OFV - ofv_final[1,]$OFV
            ofv_final <- ofv_final[!(is.na(ofv_final[,2]) | is.na(ofv_final[,4])),]
            ofv_final <- ofv_final[ofv_final$OFV!=0,]
            bootscm.obj$dofv <- ofv_final
            cat.s ("OK.\n")
        } else {
            cat.s ("Data not found.\n")
        }
    } else {
        cat.s ("Re-estimation of final models not performed.\n")
    }
    if (silent==FALSE) {
        cat.s ("\nWhat was the run number of the base model for this bootSCM? ")
        bootscm.obj$runno <- readline()
    }
    cat.s ("\n")
    c1<-call("assign", pos = 1, "current.bootscm", bootscm.obj, immediate=T)
    eval(c1)
    return ()
}

read.scm.covariate.sd <- function (file, cov.recoding = NULL) {
    ## calculate standard deviation from PsN's covariate_statistics.txt
    ## and also return most common covariates and number of levels
    tmp <- readLines (file)
    cov.lines <- grep ("=> \\{", tmp)
    cov.lines <- cov.lines[!(cov.lines %in% grep ("[fractions|factors]", tmp))]
    covs <- tmp[cov.lines]
    covs <- gsub("[\\'|\\{|=>| ]", "", covs)
    factor.lines <- grep ("factors", tmp)
    close.lines <- grep ("\\}", tmp)
    stats <- list()
    stats$covs.sd <- list()
    stats$most.common <- list()
    stats$n.levels <- list()
    for (i in seq(along=covs)) {
        tmp.sub <- tmp[((factor.lines[i])+1) : (min(close.lines[close.lines>factor.lines[i]])-1)]
        tmp.stats <- get.cov.stats(tmp.sub, covs[i], cov.recoding)
        stats$covs.sd[[covs[i]]] <- tmp.stats$sd
        stats$most.common[[covs[i]]] <- tmp.stats$most.common
        stats$min.val[[covs[i]]] <- tmp.stats$min.val
        stats$max.val[[covs[i]]] <- tmp.stats$max.val
        stats$n.levels[[covs[i]]] <- tmp.stats$n.levels
    }
    stats$dichot <- names(stats$n.levels)[stats$n.levels == 2]
    return(stats)
}

get.cov.stats <- function (dat, cov.name, cov.recoding) {
    dat <- gsub ("[\\'|\\,| ]","", dat)
    covs <- strsplit (dat, "=>")
    covs.num <- data.frame (t(matrix (nrow=2, as.num(unlist(covs)))))
    covs.num$X3 <- 0
    if (cov.name %in% names(cov.recoding)) {
        recode <- cov.recoding[[cov.name]]
        if (length(covs.num[,1]) == length(recode[,1])) {
            covs.num[,3] <- recode[match(recode[,1], covs.num[,1]),2]
        }
        covs.num$X1 <- covs.num$X3
    }
    covs.list <- rep (covs.num$X1, covs.num$X2)
    most.common <- covs.num[covs.num$X2==max(covs.num$X2),]$X1[1]
    min.val <- min(covs.num$X1)
    max.val <- max(covs.num$X1)
    n.levels <- length(unique(covs.num[,1]))
    stats <- list ("sd" = sd(covs.list), "most.common" = most.common, "min.val" = min.val, "max.val" = max.val, "n.levels" = n.levels)
    return (stats)
}

as.num <- function (dat) {
    return (as.numeric (as.character(dat)))
}

convert.most.common <- function (th) {
    ## For dichotomous covariates, PsN always uses:
    ##   CL = TVCL                     # for most common value
    ##   CL = TVCL * (1 + THETA(1))    # for other covariate value
    ## this script recalculates the value if the other implementation is desired
    return (-(1-(1/(1+th))))
}

get.pars.from.relations.file <- function (folder) {
    par.file <- readLines (paste(folder,"/scm_dir1/relations.txt", sep=""))
    grp <- par.file[grep("->\\{", par.file)]
    p.grp <- gregexpr ("\\'" , grp)
    pars <- c()
    for (i in seq(along=p.grp)) {
        pars <- c(pars, substr(grp[i], p.grp[[i]][1]+1, p.grp[[i]][2]-1 ))
    }
    return(unique(pars))
}

reshape.simple <- function (dat) { # simple reshape function
    comb <- c()
    for (i in seq(along=colnames(dat))) {
        comb <- rbind (comb, cbind ("cov" = colnames(dat)[i], "value" = dat[,i]))
    }
    comb <- data.frame(comb)
    comb$value <- as.numeric(as.character(comb$value))
    return (comb)
}

read.bootscm.par.est <- function (folder, n.bs = 100, cov.recoding = NULL, verbose = TRUE, 
                                  dofv.forward = 3.84, dofv.backward = 6.64){
  bs_final <- c()
  bs_first <- c()
  covariate <- list()
  covariate$sd <- c()
  covariate$most.common <- c()
  covariate$n.levels <- c()
  pars <- get.pars.from.relations.file(folder)
  first.non.na <- function(dat) {
    return(dat[!is.na(dat)][1])
  }
  add.to.table <- function(tab, row, nams) {
    if (is.null(tab)) {
      tab <- data.frame(t(c(row)))
      colnames(tab) <- nams
    }
    else {
      tab <- rbind(tab, NA)
      j <- length(tab[, 1])
      for (k in seq(along = nams)) {
        if (nams[k] %in% colnames(tab)) {
          tab[j, ][[nams[k]]] <- as.numeric(row[k])
        }
        else {
          tab[[nams[k]]] <- NA
          tab[j, ][[nams[k]]] <- as.numeric(row)
        }
      }
    }
    return(tab)
  }
  for (j in 1:n.bs) {
    if (file.exists(paste(folder, "/scm_dir", j, sep = ""))) {
      tmp_full <- read.csv(file = paste(folder, "/scm_dir", 
                                        j, "/raw_results_bsmod_", j, ".csv", sep = ""))
      cov_cols <- c((grep("ofv", colnames(tmp_full))[1] + 
                       1):(grep("^OM", colnames(tmp_full))[1] - 
                             1))
      covs <- tmp_full[, cov_cols]
      if (length(grep("th[[:digit:]]", colnames(covs))) > 
            0) {
        covs <- covs[, -(grep("th[[:digit:]]", colnames(covs)))]
      }
      nams <- cbind(colnames(covs))
      for (i in seq(nams)) {
        splt <- strsplit(nams[i], "\\.")[[1]]
        nams[i] <- paste(splt[1], splt[2], sep = ".")
      }
      for (i in seq(pars)) {
        sel <- grep(pars[i], substr(nams, 1, nchar(pars[i])))
        nams[sel] <- paste(pars[i], substr(nams[sel], 
                                           nchar(pars[i]) + 1, nchar(nams[sel])), sep = ".")
      }
      forward_steps <- unique(tmp_full[tmp_full$action == 
                                         "added", ]$step.num)
      tmp_full$row <- 1:length(tmp_full[, 1])
      ofv.base <- tmp_full[tmp_full$step.number == 0, ]$ofv
      best.model.row <- 0
      best.model.ofv <- ofv.base
      for (i in forward_steps) {
        step.tmp <- tmp_full[tmp_full$step.number == 
                               i & tmp_full$action == "added", ]
        min.tmp <- min(step.tmp$ofv)
        if (!is.na(min.tmp)) {
          if (min.tmp < (best.model.ofv - dofv.forward)) {
            best.tmp <- step.tmp[step.tmp$ofv == min.tmp, 
                                 ]$row
            best.model.row <- best.tmp
            best.model.ofv <- min.tmp
          }
        }
      }
      backward_steps <- unique(tmp_full[tmp_full$action == 
                                          "removed", ]$step.num)
      for (i in backward_steps) {
        step.tmp <- tmp_full[tmp_full$step.number == 
                               i & tmp_full$action == "removed", ]
        min.tmp <- min(step.tmp$ofv)
        if (!is.na(min.tmp)) {
          if (min.tmp < (best.model.ofv + dofv.backward)) {
            best.tmp <- step.tmp[step.tmp$ofv == min.tmp, 
                                 ]$row
            best.model.row <- best.tmp
            best.model.ofv <- min.tmp
          }
        }
      }
      est_final <- tmp_full[best.model.row, cov_cols]
      if (length(grep("th[[:digit:]]", colnames(est_final))) > 
            0) {
        est_final <- est_final[, -(grep("th[[:digit:]]", 
                                        colnames(est_final)))]
      }
      est_first <- apply(covs, 2, first.non.na)
      names(est_first) <- nams
      bs_first <- add.to.table(bs_first, est_first, names(est_first))
      names(est_final) <- nams
      bs_final <- add.to.table(bs_final, as.numeric(est_final), 
                               names(est_final))
      stats <- read.scm.covariate.sd(paste(folder, "/scm_dir", 
                                           j, "/covariate_statistics.txt", sep = ""), cov.recoding)
      std <- stats$covs.sd
      std[!is.na(std)] <- as.num(std[!is.na(std)])
      if (j == 1) {
        nams.sd <- names(std)
        covariate$dichot <- stats$dichot
      }
      covariate$sd <- rbind(covariate$sd, as.num(t(std)))
      covariate$most.common <- rbind(covariate$most.common, 
                                     unlist(stats$most.common))
      covariate$min.val <- rbind(covariate$min.val, unlist(stats$min.val))
      covariate$max.val <- rbind(covariate$max.val, unlist(stats$max.val))
      covariate$n.levels <- rbind(covariate$n.levels, as.num(t(stats$n.levels)))
    }
  }
  if (verbose == TRUE) {
    cat("\n    Importing step done. Processing imported data...")
  }
  colnames(covariate$sd) <- nams.sd
  colnames(covariate$n.levels) <- nams.sd
  tmp <- covariate$sd
  tmp.min <- covariate$min.val
  covariate$sd.all <- c()
  covs.dichot <- c()
  for (i in seq(along = nams)) {
    splt <- strsplit(nams[i], "\\.")[[1]]
    covariate$sd.all <- cbind(covariate$sd.all, covariate$sd[, 
                                                             match(splt[2], colnames(tmp))])
    covariate$min.val.all <- cbind(covariate$min.val.all, 
                                   covariate$min.val[, match(splt[2], colnames(tmp))])
    colnames(covariate$sd.all)[i] <- nams[i]
    colnames(covariate$min.val.all)[i] <- nams[i]
    sel.dichot <- match(splt[2], covariate$dichot)
    if (!is.na(sel.dichot)) {
      covs.dichot <- c(covs.dichot, nams[i])
    }
  }
  for (i in seq(along = covs.dichot)) {
    if (!is.na(match(covs.dichot[i], colnames(bs_first)))) {
      min.val <- min(data.frame(covariate$min.val.all)[[covs.dichot[i]]])
      sel <- data.frame(covariate$most.common)[[covs.dichot[i]]] != 
        min.val
      bs_first[[covs.dichot[i]]][sel] <- convert.most.common(bs_first[[covs.dichot[i]]][sel])
      bs_final[[covs.dichot[i]]][sel] <- convert.most.common(as.numeric(bs_final[[covs.dichot[i]]][sel]))
    }
  }
  bs_first_norm <- bs_first * covariate$sd.all
  bs_final_norm <- bs_final * covariate$sd.all
  mean.na <- function(dat) {
    dat <- as.numeric(dat)
    return(mean(dat[!is.na(dat)]))
  }
  sd.na <- function(dat) {
    dat <- as.numeric(dat)
    return(sd(dat[!is.na(dat)]))
  }
  rse.na <- function(dat) {
    dat <- as.numeric(dat)
    return(sd(dat[!is.na(dat)])/mean(dat[!is.na(dat)]))
  }
  first.step.stats <- rbind(apply(bs_first, 2, mean.na), apply(bs_first, 
                                                               2, sd.na), apply(bs_first, 2, rse.na))
  rownames(first.step.stats) <- c("mean", "sd", "rse")
  final.step.stats <- rbind(apply(bs_final, 2, mean.na), apply(bs_final, 
                                                               2, sd.na), apply(bs_final, 2, rse.na))
  rownames(final.step.stats) <- c("mean", "sd", "rse")
  create.bias.table <- function(first, final) {
    tmp <- reshape.simple(first)
    tmp$cov.type <- "Continuous"
    if (length(covariate$dichot) > 0) {
      tmp[tmp$cov %in% covs.dichot, ]$cov.type <- "Dichotomous"
    }
    tmp$incl <- c((!is.na(final)) * 1)
    tmp <- tmp[order(tmp$incl), ]
    tmp[tmp$incl == 0, ]$incl <- "Not included"
    tmp[tmp$incl == 1, ]$incl <- "Included"
    mn.tmp1 <- aggregate(tmp[tmp$incl == "Included", ]$value, 
                         by = list(tmp[tmp$incl == "Included", ]$cov), mean.na)
    never.incl <- unique(tmp$cov)[!(unique(tmp$cov) %in% 
                                      mn.tmp1$Group.1)]
    if (length(never.incl) > 0) {
      mn.tmp1 <- rbind(mn.tmp1, aggregate(tmp[tmp$cov %in% 
                                                never.incl, ]$value, by = list(tmp[tmp$cov %in% 
                                                                                     never.incl, ]$cov), mean.na))
    }
    mn.tmp1 <- mn.tmp1[match(unique(tmp$cov), mn.tmp1$Group.1), 
                       ]
    b.stats.incl <- data.frame(cbind(cov = as.character(unique(tmp$cov)), 
                                     mean = mn.tmp1$x, All = as.num(aggregate(tmp$value, 
                                                                              by = list(tmp$cov), mean.na)$x), incl = "Included"))
    mn.tmp2 <- aggregate(tmp[tmp$incl == "Not included", 
                             ]$value, by = list(tmp[tmp$incl == "Not included", 
                                                    ]$cov), mean.na)
    always.incl <- unique(tmp$cov)[!(unique(tmp$cov) %in% 
                                       mn.tmp2$Group.1)]
    if (length(always.incl) > 0) {
      mn.tmp2.add <- aggregate(tmp[tmp$cov %in% always.incl, 
                                   ]$value, by = list(tmp[tmp$cov %in% always.incl, 
                                                          ]$cov), mean.na)
      mn.tmp2.add$x <- NA
      mn.tmp2 <- rbind(mn.tmp2, mn.tmp2.add)
    }
    mn.tmp2 <- mn.tmp2[match(unique(tmp$cov), mn.tmp2$Group.1), 
                       ]
    b.stats.nincl <- data.frame(cbind(cov = as.character(unique(tmp$cov)), 
                                      mean = mn.tmp2$x, All = as.num(aggregate(tmp$value, 
                                                                               by = list(tmp$cov), mean.na)$x), incl = "Not Included"))
    b.stats <- data.frame(rbind(b.stats.incl, b.stats.nincl))
    add.nincl <- b.stats.incl[!b.stats.incl$cov %in% b.stats.nincl$cov, 
                              ]$cov
    add.incl <- b.stats.nincl[!b.stats.nincl$cov %in% b.stats.incl$cov, 
                              ]$cov
    for (i in seq(along = add.nincl)) {
      b.stats <- rbind(b.stats, c(as.character(add.nincl[i]), 
                                  "Not included", NA))
    }
    for (i in seq(along = add.incl)) {
      b.stats <- rbind(b.stats, c(as.character(add.incl[i]), 
                                  "Included", NA))
    }
    bias.dat <- data.frame(b.stats)
    bias.dat$bias <- as.num(100 * (as.num(bias.dat$mean) - 
                                     as.num(bias.dat$All))/as.num(bias.dat$All))
    bias.dat$cov.type <- "Continuous"
    if (length(covariate$dichot) > 0) {
      bias.dat[bias.dat$cov %in% covs.dichot, ]$cov.type <- "Dichotomous"
    }
    return(list(table.long = tmp, bias.dat = bias.dat))
  }
  tab <- create.bias.table(bs_first, bs_final)
  tab.norm <- create.bias.table(bs_first_norm, bs_final_norm)
  return(list(par.est.first = bs_first, par.est.first.norm = bs_first_norm, 
              par.est.first.stats = data.frame(first.step.stats), par.est.final = bs_final, 
              par.est.final.norm = bs_final_norm, par.est.final.stats = data.frame(final.step.stats), 
              covariate = covariate, bias.dat = tab$bias.dat, bias.dat.norm = tab.norm$bias.dat, 
              par.est.long = tab$table.long, par.est.long.norm = tab.norm$table.long, 
              pars = pars))
}


convert.bootscm.cov.inclusion <- function (cov.incl, group.by.cov = FALSE) {
    ## convert table: discard info about parameter (CL / V / etc.) and relation type (linear / nonlin / etc.) if
    n <- max(cov.incl$bs_n)
    cov.incl <- cov.incl[,-1]
    covs <- colnames(cov.incl)
    if (group.by.cov == TRUE) {
        covs.names <- c()
        for (i in seq(along=covs)) {
            tmp_cov <- convert.cov.name (covs[i])
            if (!is.na(tmp_cov)) {
                covs.names[i] <- tmp_cov
            }
        }
        covs.unq <- unique(covs.names)
        results.tab <- data.frame (matrix(0, nrow = length(cov.incl[,1]), ncol=length(covs.unq)))
        colnames(results.tab) <- covs.unq
        for (i in seq(along=cov.incl[,1])) {
            for (j in seq(along=cov.incl[i,])) {
                if (cov.incl[i,j] == 1) {
                    results.tab[i,][[covs.names[j]]] <- 1
                }
            }
        }
    } else {
        results.tab <- cov.incl[,-(length(cov.incl[1,]))]
        covs.names <- colnames(results.tab)
    }

    ## separate off Dummy covariates
    res <- list ("n" = n, "results.tab" = results.tab)
    if (group.by.cov == FALSE) {
        patt <- "\\.X"
    } else {
        patt <- "^X."
    }
    cols.dum <- grep(patt, colnames(results.tab))
    if (length(cols.dum)>0) {
        results.tab.dum <- cbind(results.tab[,cols.dum])
        colnames(results.tab.dum) <- gsub(patt, "", colnames(results.tab)[cols.dum])
        res$results.tab <- results.tab[,-cols.dum]
        res$results.tab.dum <- results.tab.dum
    }
    return (res)
}

convert.cov.name <- function (cov) {
    return (strsplit (cov, "\\.")[[1]][2])
}

convert.bootscm.bs.ids <- function (bs_ids) {
    ids <- unique (unlist(bs_ids[,-1]))
    n <- length(bs_ids[,1])
    oid <- data.frame (matrix(0, nrow = n, ncol=length(ids)))
    colnames(oid) <- paste("X", ids, sep="")
    tmp <- bs_ids[,-1]
    for (i in 1:length(tmp[,1])) {
        cols <- match(paste("X", tmp[i,], sep=""), colnames(oid))
        for (k in seq(along=cols)) {
            oid[i,cols[k]] <- oid[i,cols[k]] + 1                   # save original ID numbers selected in the bootstrap
        }
    }
    return(oid)
}

ask.folder <- function () {
    d <- gsub ("\\./", "", list.dirs (path = ".", pattern="scm", full.names = TRUE))
    if (length(d) == 0) { # maybe folder was not named "scm"
        d <- gsub ("\\./", "", list.dirs (path = ".", full.names = TRUE))
    }
    cat ("Import from subfolder (filtered on 'scm'):\n  ")
    cat (paste (d,"\n"), sep = "  ")
    cat ("\nFolder with bootSCM data (Enter to abort): ")
    ans <- readline()
    if (ans == "") {
        return()
    }
    if (!is.na(file.info(ans)$isdir)) {
        if (file.info(ans)$isdir) {
            return (ans)
        } else {
            cat("Please choose a valid folder!\n\n")
            Recall()
        }
    } else {
        cat("Please choose a valid folder!\n\n")
        Recall()
    }
}

ask.group.by.cov <- function () {
    cat ("\nPlease choose how you want to import covariate inclusion frequencies:\n")
    cat ("  1: Not grouped, inclusion frequencies per parameter-covariate relationship\n")
    cat ("  2: Grouped, inclusion frequencies per covariate\n")
    ans <- readline()
    if (ans == "") {
        return()
    }
    if (ans==1|ans==2) {
        return (ans)
    } else {
        cat("Please choose a valid option!\n\n")
        Recall()
    }
}

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
    all <- list.files(path, pattern, all.dirs, full.names, recursive=FALSE, ignore.case)
    return(all[file.info(paste(path, "/", all, sep=""))$isdir])
}
