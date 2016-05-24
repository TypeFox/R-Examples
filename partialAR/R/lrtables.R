# lrtables.R -- functions for building up tables of likelihood ratios
# Copyright (C) 2015 Matthew Clegg

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


# This file builds and loads various pre-computed tables that have been
# generated to calibrate the likelihood ratio tests.

# if(getRversion() >= "2.15.1")  utils::globalVariables(c("PAR.SAMPLES",
#  "PAR.JOINT.CRITICAL.VALUES.DT",
#  "PAR.JOINT.CRITICAL.VALUES.KPSS.DT",
#  "PAR.SAMPLES.DT"))

# if (!exists("PAR.SAMPLES")) PAR.SAMPLES <- NULL
# if (!exists("PAR.SAMPLES.DT")) PAR.SAMPLES.DT <- NULL
if (!exists("PAR.POWER.SAMPLES")) PAR.POWER.SAMPLES <- NULL
if (!exists("PAR.POWER.PVMR.SAMPLES")) PAR.POWER.PVMR.SAMPLES <- NULL
# if (!exists("PAR.JOINT.CRITICAL.VALUES.DT")) PAR.JOINT.CRITICAL.VALUES.DT <- NULL
# if (!exists("PAR.JOINT.CRITICAL.VALUES.KPSS.DT")) PAR.JOINT.CRITICAL.VALUES.KPSS.DT <- NULL
  
build_par_tables <- function (dir="tables", debug=FALSE, nrep=10000, rebuild_samples=TRUE) { 
    # Rebuilds all of the tables that are contained in this file.
    # This function takes several days to a week to run to completion.
    # The tables are written to files in the specified directory,
    # where they can then be loaded back into R and used to update
    # this file.
    
    if (debug) {
        nr <- 1
    } else {
        nr <- nrep
    }

    if (rebuild_samples) {
        cat("Rebuilding PAR likelihood ratio samples ...\n")
        par.generate.likelihood_ratio.samples(nrep=nrep)
    }
    SAMPLES <- par.load.likelihood_ratio.samples()
    
    dir.create(dir, recursive=TRUE)

    par.rwnull.lrqt <- quantile.table.from.samples("rw_lrt", 
        SAMPLES[SAMPLES$rho==1.0 & SAMPLES$robust == 0,])
    dput (par.rwnull.lrqt, sprintf("%s/%s", dir, "PAR.RWNULL.LRQT"))

    par.mrnull.lrqt <- quantile.table.from.samples("mr_lrt", 
        SAMPLES[SAMPLES$rho==0.9 & SAMPLES$robust == 0,])
    dput (par.mrnull.lrqt, sprintf("%s/%s", dir, "PAR.MRNULL.LRQT"))

    par.mrnull.kpss <- quantile.table.from.samples("kpss_nstat", 
        SAMPLES[SAMPLES$rho==0.9 & SAMPLES$robust == 0,])
    dput (par.mrnull.kpss, sprintf("%s/%s", dir, "PAR.MRNULL.KPSS"))

    par.rwnull.rob.lrqt <- quantile.table.from.samples("rw_lrt", 
        SAMPLES[SAMPLES$rho==1.0 & SAMPLES$robust == 1,])
    dput (par.rwnull.rob.lrqt, sprintf("%s/%s", dir, "PAR.RWNULL.ROB.LRQT"))

    par.mrnull.rob.lrqt <- quantile.table.from.samples("mr_lrt", 
        SAMPLES[SAMPLES$rho==0.9 & SAMPLES$robust == 1,])
    dput (par.mrnull.rob.lrqt, sprintf("%s/%s", dir, "PAR.MRNULL.ROB.LRQT"))

    par.mrnull.rob.kpss <- quantile.table.from.samples("kpss_nstat", 
        SAMPLES[SAMPLES$rho==0.9 & SAMPLES$robust == 1,])
    dput (par.mrnull.rob.kpss, sprintf("%s/%s", dir, "PAR.MRNULL.ROB.KPSS"))

    par.joint.critical.values <- par.findall.joint.critical.values()
    dput (par.joint.critical.values, sprintf("%s/%s", dir, "PAR.JOINT.CRITICAL.VALUES"))

    par.joint.critical.values <- par.findall.joint.critical.values(ar1test="kpss")
    dput (par.joint.critical.values, sprintf("%s/%s", dir, "PAR.JOINT.CRITICAL.VALUES.KPSS"))
    
    printf ("%s done\n", Sys.time())  
}

load_table <- function (..., dir="tables") {
    # Loads a table and stores it in a global variable
    for (table_name in list(...)) {
        printf("Loading %s\n", table_name)
        tab <- dget(sprintf("%s/%s", dir, table_name))
        if (exists(table_name, envir=asNamespace("partialAR"))) {
            unlockBinding(table_name, asNamespace("partialAR"))
        }
        assign(table_name, tab, envir = asNamespace("partialAR"))
    }
}

load.lrtables <- function () {
    # PAR: Likelihood ratio tests vs. single hypothesis of AR(1) or Random Walk
    load_table("PAR.RWNULL.LRQT", "PAR.MRNULL.LRQT")

    # PAR: Robust likelihood ratio tests vs. single hypothesis
    load_table("PAR.RWNULL.ROB.LRQT", "PAR.MRNULL.ROB.LRQT")

    # PAR: KPSS tests for AR(1) vs. PAR hypothesis
    load_table("PAR.MRNULL.KPSS", "PAR.MRNULL.ROB.KPSS")
        
    load_table("PAR.JOINT.CRITICAL.VALUES")
#    if (exists("PAR.JOINT.CRITICAL.VALUES.DT", envir=asNamespace("partialAR"))) {
#        unlockBinding("PAR.JOINT.CRITICAL.VALUES.DT", asNamespace("partialAR"))
#    }
#    PAR.JOINT.CRITICAL.VALUES.DT <<- as.data.table(PAR.JOINT.CRITICAL.VALUES)
    
    load_table("PAR.JOINT.CRITICAL.VALUES.KPSS")
#    if (exists("PAR.JOINT.CRITICAL.VALUES.KPSS.DT", envir=asNamespace("partialAR"))) {
#        unlockBinding("PAR.JOINT.CRITICAL.VALUES.KPSS.DT", asNamespace("partialAR"))
#    }
#    PAR.JOINT.CRITICAL.VALUES.KPSS.DT <<- as.data.table(PAR.JOINT.CRITICAL.VALUES.KPSS)
    
}

dump.lrtables <- function(filename="lrdata.R") {
    tables_list <- c("PAR.RWNULL.LRQT", "PAR.MRNULL.LRQT",
        "PAR.RWNULL.ROB.LRQT", "PAR.MRNULL.ROB.LRQT",
        "PAR.MRNULL.KPSS", "PAR.MRNULL.ROB.KPSS",
        "PAR.JOINT.CRITICAL.VALUES",
        "PAR.JOINT.CRITICAL.VALUES.KPSS")
        
    dump(tables_list, filename)

    cat("\n\n", file=filename, append=TRUE)
#    cat("PAR.SAMPLES.DT <- data.table(PAR.SAMPLES)\n", 
#        file = filename, append=TRUE)
    cat("# PAR.JOINT.CRITICAL.VALUES.DT <- data.table(PAR.JOINT.CRITICAL.VALUES)\n", 
        file = filename, append=TRUE)
    cat("# PAR.JOINT.CRITICAL.VALUES.KPSS.DT <- data.table(PAR.JOINT.CRITICAL.VALUES.KPSS)\n", 
        file = filename, append=TRUE)

}

