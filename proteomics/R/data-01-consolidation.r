# Copyright (C) 2012-2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
#
#     This program is free software: you can redistribute it and/or
#     modify it under the terms of the GNU General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see
#     <http://www.gnu.org/licenses/>.

###########################################
## Add Ion and Retention Time Statistics ##
###########################################

#' Summary statistics -- Ion intensities per spectra
#'
#' @param dwide iTRAQ data in wide format including columns corresponding to
#' iTRAQ channels containing their intensities.
#'
#' @export
addIonSatistics <- function(dwide) {
    ch <- attr(dwide, "channels")
    dwide$ions.sum    <- apply(dwide[,ch], 1, sum, na.rm=TRUE)
    dwide$ions.min    <- apply(dwide[,ch], 1, min, na.rm=TRUE)
    dwide$ions.max    <- apply(dwide[,ch], 1, max, na.rm=TRUE)
    dwide$ions.median <- apply(dwide[,ch], 1, median, na.rm=TRUE)
    dwide$ions.mean   <- apply(dwide[,ch], 1, mean, na.rm=TRUE)
    return(dwide)
}

#' Summary statistics -- Generic to calculate summary statistics
#'
#' Calculates generic summary statistics based on a given formula.
#'
#' @param dwide iTRAQ data in wide format
#' @param frm for example:
#'    frm <- value ~ protein + variable
#'    frm <- value ~ peptide + variable
#' @export
responseStatisics <- function (dwide, frm) {
    statistics <- c("median", "mean", "sd", "length", "min", "max")

    val <- all.vars(frm)[1]

    func <- function(stat) {
        dstat <- aggregate(frm, data=dwide, stat)
        names(dstat) <- mapvalues(  names(dstat)
                                  , val
                                  , paste(val, stat, sep="."))
        return(dstat)
    }

    stat=NULL
    statList <- foreach(stat=statistics) %dopar% func(stat)
    return(Reduce(merge, statList))
}

#' Summary statistics -- Calculates retention time statistics at apex
#'
#' Calculates different summary retention time statistics for each peptide (a
#' subsequence of a protein including post translational modifications).  The
#' idea is that each peptide is supposed to have roughly the same retention
#' time.
#'
#' @param dwide iTRAQ data in wide format
#' @param ... Additional arguments passed for ddply
#' @export
addRetentionAtApex <- function(dwide, ...) {
    atApex <-  function(d) {
        return(data.frame(retention.atApex=with(d, retention[which.max(intensity)])))
    }

    dwideAtApex <- ddply(dwide, .(peptide, experiment, injection), atApex, ...)

    rtstats <- merge(dwide, dwideAtApex)
    rtstats$outlier <- FALSE
    attr(rtstats, "channelnames") <- attr(dwide, "channelnames")
    attr(rtstats, "channels") <- match(attr(rtstats, "channelnames"), names(rtstats))
    attr(rtstats, "referencename") <- attr(dwide, "referencename")
    attr(rtstats, "reference") <- match(attr(rtstats, "reference"), names(rtstats))
    return(rtstats)
}

#' Summary statistics -- Calculates index retention time statistics
#'
#' @param dwide iTRAQ data in wide format
#' @param ... Additional arguments passed for ddply
#' @export
addRetentionIndexTimeStatistics <- function (dwide, ...) {
    retention_frm <- retention.atApex ~ peptide + charge
    dwideStats <- responseStatisics(dwide, retention_frm)
    rtstats <- merge(dwide, dwideStats)
    attr(rtstats, "channelnames") <- attr(dwide, "channelnames")
    attr(rtstats, "channels") <- match(attr(rtstats, "channelnames"), names(rtstats))
    attr(rtstats, "referencename") <- attr(dwide, "referencename")
    attr(rtstats, "reference") <- match(attr(rtstats, "reference"), names(rtstats))
    return(rtstats)
}

####################################
## Plot Retention Time Statistics ##
####################################

#' Plot Retention Time Statistics
#'
#' Plot retention times with possible outliers
#'
#' @param rwide iTRAQ data in wide format with retention time information
#' @examples
#' \dontrun{
#' iglobal <- addIonSatistics(pglobal)
#' rglobal <- addRetentionTimeStatistics(iglobal, .parallel=TRUE)
#' rglob$outlier <- with(rglob, abs(retention.atApex - retention) > 4)
#' p <- pRetention(rglobal)
#'
#' p + geom_point(aes(retention.atApex, retention))
#' p + geom_point(aes(retention.atApex, retention-retention.atApex))
#' p + geom_point(aes(ppm, retention-retention.atApex))
#' p + geom_density(aes(x=ppm), alpha=.242)
#' }
#' @export
pRetention <- function(rwide) {
    return(  ggplot(rwide, aes(colour=outlier, fill=outlier))
           + scale_colour_manual(name="Possible\noutlier", values=c("black", "red"))
           + scale_fill_manual(name="Possible\noutlier", values=c("black", "red"))
           + ylab("Density")
           + coord_cartesian()
           + facet_wrap(~injection)
           )
}

##################################################################################
## Tranform data from intestity scales to density histrograms for data analysis ##
##################################################################################

#' Transformation -- From intensity scales to density histrograms
#'
#' @param dwide iTRAQ data in wide format
#'
#' @export
toProportions <- function(dwide) {
    ch <- attr(dwide, "channels")
    cs <- apply(dwide[,ch],1,sum)
    dwide[,ch] <- dwide[,ch] / cs
    return(dwide)
}

##############################################################################
## Measuring stability by evaluating angle of loading vector from identity ###
##############################################################################

alpha <- function (p) {
    return( acos(sum(p) / sqrt( sum(p^2) * length(p))))
}

#' Measuring stability -- angle of loading vector
#'
#' Measuring stability by evaluating angle of loading vector from identity
#'
#' @param dwide iTRAQ data in wide format
#' @export
toAlpha <- function (dwide) {
    ch <- attr(dwide, "channels")
    dwide$alpha <- apply(dwide[,ch],1,alpha)
    return(dwide)
}

####################################
## Plot Retention Time Statistics ##
####################################

#' Plot Retention Time Statistics in violine form
#'
#' @param dat iTRAQ in log format
#' @param target of the norming
#' @export
pVioline <- function (dat, target) {
    return(ggplot(dat, aes(variable, value, fill=variable))
           + geom_violin()
           + geom_hline(yintercept=target)
           + scale_fill_discrete(expression("Channel"))
           + xlab(expression("Channel"))
           + ylab(expression("Standardised abundance"))
           + theme(axis.text = element_text(colour = "black"), legend.position="none")
           )
}

###################################################################
## Adjust for confounding due to differences in channel loadings ##
###################################################################

avrgLoadingCalculation <- function(dprop, chan) {
    return(as.data.frame(t(apply(dprop[,chan], 2, median, na.rm=T))))
}

#' Adjust for confounding -- calculates the average loading
#'
#' @param dwide iTRAQ data in wide format
#' @export
avrgLoading <- function(dwide) {
    chan <- attr(dwide, "channels")
    avrg <- ddply(dwide, .(experiment), avrgLoadingCalculation, chan=chan)
    attr(avrg, "channelnames") <- attr(dwide, "channelnames")
    attr(avrg, "channels") <- match(attr(avrg, "channelnames"), names(avrg))
    attr(avrg, "referencename") <- attr(dwide, "referencename")
    attr(avrg, "reference") <- match(attr(avrg, "reference"), names(avrg))
    return(avrg)
}

#' Adjust for confounding -- add an appropiate  target
#'
#' @param dwide iTRAQ data in wide format
#' @param byRef schould the average be calculated from the loading of the
#' reference channel. Default is FALSE and this is recommended.
#' @export
addLoadings <- function (dwide, byRef=F) {
    avrg <- avrgLoading(dwide)
    trgt <- ifelse(  is.na(attr(avrg, "reference") | !byRef)
                   , 0
                   , mean(as.vector(log(avrg[,attr(avrg, "reference")]))))
    attr(dwide, "loadings") <- avrg
    attr(dwide, "logtargetloading") <- trgt
    return(dwide)
}

#' Adjust for confounding -- copy loadings from one experiment to another
#'
#' This is important when analysing enriched samples. Here, use the loading
#' averages from the corresponig non-eriched sample.
#'
#' @param fromWide iTRAQ data in wide format
#' @param toWide iTRAQ data in wide format
#' @export
copyLoadings <- function (fromWide, toWide) {
    attr(toWide, "loadings") <- attr(fromWide, "loadings")
    attr(toWide, "logtargetloading") <- attr(fromWide, "logtargetloading")
    return(toWide)
}

#' Adjust for confounding -- State of the art adjustments for confounding
#'
#' Compensate for possible confounding due the transformed values for the
#' statistical analysis.
#'
#' @param dwide iTRAQ data in wide format.
#' @export
adjusting <- function (dwide) {
    e = NULL # otherwise foreach throws a global variable
    avrg <- attr(dwide, "loadings")
    trgt <- attr(dwide, "logtargetloading")

    adjustV <- function (v, logloading, trgt) {
        exp(log(v) - logloading + trgt)
    }

    adjustExperiment <- function (exper, logloading, trgt) {
        unadjchannels <- as.matrix(exper[,attr(exper, "channels")])

        adjchan <- apply(unadjchannels, 1, adjustV, logloading, trgt)
        rownames(adjchan) <- colnames(unadjchannels)
        adjchan <- as.data.frame(t(adjchan))

        attr(adjchan, "channelnames") <- attr(dwide, "channelnames")
        attr(adjchan, "channels") <- match(attr(adjchan, "channelnames"), names(adjchan))
        exper[,attr(exper, "channels")] <- adjchan[,attr(adjchan, "channels")]
        return(exper)
    }

    lwide <- foreach(e=levels(dwide$experiment)) %dopar% {
        exper <- dwide[dwide$experiment == e,]
        logloading <- log(as.matrix(avrg[avrg$experiment == e,attr(avrg, "channels")]))
        return(adjustExperiment(exper, logloading, trgt))
    }
    awide <- do.call(rbind,lwide)
    attr(awide, "loadings") <- avrg
    attr(awide, "loadingtarget") <- exp(trgt)
    return(awide)
}

#' Adjust for confounding -- In one single experiment only
#'
#' Simple code when only one iTRAQ-experiment has been performed. (Code not used anymore.)
#'
#' @param dwide iTRAQ data in wide format.
#' @export
adjustOne <- function(dwide) {
    chan <- attr(dwide, "channels")
    refe <- attr(dwide, "reference")
    chan <- c(refe, setdiff(chan,refe))
    medi <- apply(dwide[,chan], 2, median, na.rm=T)
    effe <- exp(log(medi) - log(medi)[1])

    dwide[,chan] <- t(t(dwide[,chan])/effe)
    attr(dwide, "loadingsadjustments") <- effe
    return(dwide)
}

#' Adjust for confounding -- Generic function for centring data
#'
#' This function calculates from given adjusting factors that compensate for
#' possible confounding due the transformed values for the statistical analysis.
#'
#' Can be used to perfome custom adjustments.  (Code not used anymore.)
#'
#' @param dwide iTRAQ data in wide format.
#' @param effect estimated effects which may yield to confounding.
#' @param ch names of the channel columns.
#' @export
adjustBy <- function(dwide, effect, ch) {
    dwide[ch] <- t(t(dwide[ch]) / effect)
    return(dwide)
}

######################################################################
## From spectrum to protein level: Response variable calculation
######################################################################

#' Response calculation
#'
#' Calculates needed sample size accumulation from iTRAQ data which is given on
#' spectrum level.
#'
#' @param dwide iTRAQ data in wide format including columns corresponding to
#' iTRAQ channels containing their intensities.
#' @export
accum <- function (dwide) {
    accumOne <- function (subwide) {
        subwide <- droplevels(subwide)
        return(data.frame(peptide=nlevels(subwide$peptide), id=nlevels(subwide$id)))
    }

    return(ddply(dwide, .(protein), accumOne, .parallel=TRUE))
}

#' Response calculation
#'
#' From spectrum to protein level -- Response variable calculation
#'
#' @param dwide iTRAQ data in wide format including columns corresponding to
#' iTRAQ channels containing their intensities.
#' @param acc result of an accumulation of sample sizes
#' @export
channelResponses <- function (dwide, acc) {
    n <- attr(dwide, "channels")

    estimate <- function (subwide) {
        d <- log2(subwide[,n])
        prt <- unique(subwide$protein)
        return(data.frame(variable=names(d),
                          mean=apply(d, 2, mean),
                          var=apply(d, 2, var),
                          accum.peptide=acc$peptide[acc$protein == prt],
                          accum.id=acc$id          [acc$protein == prt]
                          ))
    }

    dw <- ddply(dwide, .(protein), estimate, .parallel=TRUE)

    attr(dw, "referencename") <- attr(dwide, "referencename")
    return(dw)
}

#' Response calculation
#'
#' Norming the responses of a single iTRAQ to a given reference channel.
#'
#' @param dlong iTRAQ data in long format.
#' @export
norm2Reference <- function (dlong) {
    ref <- attr(dlong, "referencename")
    estimate <- function (sublong) {
        sublong$mean <- sublong$mean - sublong$mean[sublong$variable == ref]
        sublong$var  <- sublong$var  + sublong$var [sublong$variable == ref]
        return(sublong)
    }

    dl <- ddply(dlong, .(protein), estimate, .parallel=TRUE)
    return(dl[!(dl$variable == ref),])
}

######################################################################
## Setting the stage (sampling design) for the statistical analysis ##
######################################################################

#' Sample design -- Generating multiple factor designs from one-dimensional factor
#'
#' Making a multiple-factor ANOVA from the single channel variable of an iTRAQ
#' experiment.
#'
#' This function uses a matrix convmat to convert the single channel into a
#' full fledged multiple factor ANOVA.
#'
#' @param dwide iTRAQ data in wide format including columns corresponding to
#' iTRAQ channels containing their intensities.
#'
#' @param cvmat a matrix that hold the information on which channel is mapped
#' to which factor.
#'
#' @examples
#' channels <- c("X113", "X114", "X115", "X116", "X117", "X118", "X119") #, "X121")
#' typus     <- c(rep(c("A", "B", "C"), each=2), "reference")
#' treatment <- c(rep(c("I", "II"), 3), "mixed")
#' convmat   <- data.frame(channels=channels, typus=typus, treatment=treatment)
#' print(convmat)
#' \dontrun{factoring(dwide, cvmat=convmat)}
#' @export
factoring <- function(dwide, cvmat) {
    if (!all(paste(cvmat$channels) == attr(dwide, "channelnames"))) {
        error("Channel names do not match!")
    }

    dat <- melt(dwide, measure.vars=attr(dwide, "channels"))

    chs <- cvmat$channels

    addFactor <- function(name, dat) {
        dat[[name]] <- with(dat, mapvalues(variable, paste(chs), paste(cvmat[[name]])))
        return(dat)
    }

    n <- setdiff(names(cvmat), "channels")

    # this for-loop may be replaced by a complicated
    # 'reduce' function.  Would make the whole thing
    # unnecessarily complicated, though.
    for (i in n) {
        dat <- addFactor(i, dat)
    }
    return(dat)
}

#####################################################
## Plot interaction plots of peptides and proteins ##
#####################################################

#' Plot interaction plots of peptides
#'
#' @param datP subframe of peptide data
#' @export
plotMePeptide <- function (datP) {
    return(ggplot(datP, aes(x=as.numeric(variable), y=value))
           + facet_wrap(~peptide)
           + geom_line(aes(colour=id, linetype=outlier))
           #+ geom_line(aes(size=abundance.wt, colour=id, linetype=outlier))
           )
}

#' Plot interaction plots of proteins
#'
#' @param datP subframe of protein data
#' @export
plotMeProtein <- function (datP) {
    return(ggplot(datP, aes(x=as.numeric(variable), y=value))
           + facet_wrap(~protein)
           + geom_line(aes(colour=id, linetype=outlier))
           #+ geom_line(aes(size=abundance.wt, colour=id, linetype=outlier))
           )
}
