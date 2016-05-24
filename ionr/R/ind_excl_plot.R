#' Plot indicator exlusion results with and without excluded indicators
#' @description Provides an overview of the indicator exclusion results. Marked(x) indicators are excluded in the indicator exclusion procedure. See \code{\link{ind_excl}} for details.
#' \describe{
#'   \item{left}{correlations between single indicator and outcome}
#'   \item{right}{correlations between sum-score and outcome with and without the marked indicators}
#' }
#'
#' @inheritParams ind_excl
#' @param tagged items to be marked as excluded by the indicator exclusion procedure
#' @param tagged2 same as 'tagged' for second scale (e.g., informant report)
#'
#' @encoding utf-8


#' @export
#' @return See \code{\link{ind_excl}}

# for debugging
#  a<-scale_sim(n=2500, to_n=2, tn_n=6); indicators=a[[1]] ;outcome=a[[2]] ;indicators2 = vector(); scalename =
#  'scale'; outcomename = 'outcome'; indicatornames = 1:ncol(indicators); tagged = c('8', '7'); tagged2 = vector();
#  location1 = 'topleft'; location2 = 'topright'; pcrit = 0.05; multi = 1; coruse = 'everything'; ci=100

ind_excl_plot <- function(indicators, indicators2 = vector(), outcome, scalename = "scale",
                          outcomename = "outcome", indicatornames = 1:ncol(indicators),
                          tagged = vector(), tagged2 = vector(), location1 = "topleft",
                          location2 = "topright", pcrit = 0.05, multi = 1, coruse = "everything",
                          ci="estimate")
    {
    old <- options(stringsAsFactors = FALSE)

    # plot sumscore with and without excluded indicators

    #PREPARE a function to extract correlation and CI-s
    scalecor <- function(indicators, outcome, coruse) {
        kor <- cor(outcome, rowMeans(indicators, na.rm = T), use = coruse)
        return(kor)
    }

    scalecor_ci <- function(indicators, outcome,coruse="everything",ci="estimate"){
        temp=cbind(outcome, rowMeans(indicators, na.rm = T))
        # if we want CI-s based on normal theory
        if (ci=="estimate") {
            kor <- psych::corr.test(temp,adjust="none",use=coruse)
            korout=as.matrix(t(kor$ci))  # t so it would be similar to the boostrapped version
            return(korout)
        # if we want CI-s to be boostrapped
        } else if (is.numeric(ci)) {
            cis <- psych::cor.ci(temp,plot=FALSE, n.iter=ci)
            kor=scalecor(indicators, outcome, coruse)
            korout=as.matrix(unlist(cis$ci))
            korout=rbind(kor,korout)
            korout=korout[c(2,1,4,6),]
            names(korout)[2]="r"
            return(korout)
        } else {
            print("ci can either equal to 'estimate' or to a number of boostraps in integers")
        }
    }

    # function to add absolute value
    awayzero <- function(x, space) {
        if (x > 0)
            x <- x + space else x <- x - space
        return(x)
    }

    # PREPARE having 2 plots
    par(mfrow = c(1, 2))

    # PREPARE string collections
    exclid <- c("All \nindicators", "Marked (x) \n excluded")
    indicatorsid <- c("Self-report", "Informant")
    y <- paste("Correlation with", outcomename)
    x <- "Indicator number"
    main1 <- "Single indicator-\noutcome correlation"
    main2 <- "Scale-outcome correlation"
    main2.cex <- multi

    star <- paste("x: SONE\n <", pcrit)


    # make them double rowed main1=paste0(scalename,'\n',main1) main2=paste0('\n',main2)

    # barplot colours
    cols <- c("gray", "white")

    # indicators to exclude
    excl_indicators <- grep(paste0(tagged, collapse = "|"), indicatornames)
    excl_indicators2 <- grep(paste0(tagged2, collapse = "|"), indicatornames)


    # PREPARE for single indicator - outcome plot

    # n_indicators <- ncol(indicators)
    cor_indicators <- cor(indicators, outcome, use = coruse)
    if (length(indicators2) != 0) {
        cor_indicators2 <- cor(indicators2, outcome, use = coruse)
        cor_indicators <- cbind(cor_indicators, cor_indicators2)
    }
    rownames(cor_indicators) <- indicatornames
    # PREPARE for sumscore - outcome plot

    # preallocate whole array. correlations only

    cors =ci_l =ci_u <- array(NA, dim = c(2, 2), dimnames = list(indicatorsid, exclid))

    # fill in self-reported original correlation
    cors[1, 1] <- scalecor(indicators, outcome, coruse)

    # if exclusions are present
    if (length(tagged) != 0) {
        cors[1, 2] <- scalecor(indicators[, -excl_indicators], outcome, coruse)
    }

    if (length(indicators2) != 0) {
        cors[2, 1] <- scalecor(indicators2, outcome, coruse)
        if (length(tagged2) != 0) {
            cors[2, 2] <- scalecor(indicators2[, -excl_indicators2], outcome, coruse)

        }
    }
    # calculate cors with CI-s if any are required
    if (ci=="estimate"|is.numeric(ci)) {

        # old version, if the cors_ci dimensions should be dynamic. currently not needed.
        ## self-reported original CI-s
        #cors_ci1 <- scalecor_ci(indicators=indicators, outcome=outcome, coruse=coruse, ci=ci)
        ## preallocate whole matrix. Determines matrix dimensions partly on cors_ci1
        #cors_ci <- array (NA,dim = c(4,length(cors_ci1)),dimnames=list(paste(rep(indicatorsid,each=2),c("All", "Excluded")),c("lower", "r", "upper", "p")))
        #cors_ci[1,]=cors_ci1

        # preallocate matrix

        cors_ci <- array (NA,dim = c(4,4),dimnames=list(paste(rep(indicatorsid,each=2),c("All", "Excluded")),c("lower", "r", "upper", "p")))

        ## self-reported original CI-s
        cors_ci[1,] <- scalecor_ci(indicators=indicators, outcome=outcome, coruse=coruse, ci=ci)

        # if exclusions are present
        if (length(tagged) != 0) {
            cors_ci[2,] <- scalecor_ci(indicators=indicators[, -excl_indicators], outcome=outcome, coruse=coruse,ci=ci)
        }

        if (length(indicators2) != 0) {
            cors_ci[3, ] <- scalecor_ci(indicators=indicators2, outcome=outcome, coruse=coruse, ci=ci)
            if (length(tagged2) != 0) {
                cors_ci[4, ] <- scalecor_ci(indicators=indicators2[, -excl_indicators2], outcome=outcome, coruse=coruse, ci=ci)

            }
        }

        # PREPARE ci-s for plotting

        ci_l[1,] =cors_ci[1:2,1]
        ci_l[2,] =cors_ci[3:4,1]

        ci_u[1,] =cors_ci[1:2,3]
        ci_u[2,] =cors_ci[3:4,3]
    }



    # PREPARE limits

    limits <- range(cors, cor_indicators, na.rm = T)

    ## add CI-s to the limits

    if (ci=="estimate"|is.numeric(ci)) limits <- range(limits, ci_u, ci_l, na.rm = T)


    # reverse limits if the effect is negative

    if (cors[1, 1] < 0)
        limits <- rev(limits)

    # do the limits include zero? if no, include zero!

    haszero <- sum((limits[1] > 0 & limits[2] < 0) | (limits[1] < 0 & limits[2] > 0))
    if (haszero == 0)
        limits[1] <- 0

    # add extra space around limits

    spacing <- 0.02

    if (cors[1, 1] < 0) {
        limits[1] <- limits[1]
        limits[2] <- limits[2] - spacing
    } else {
        limits[1] <- limits[1]
        limits[2] <- limits[2] + spacing

    }

    ### are the scale-outcome correlation bars statistically significantly different?

    # sample sizes
    n=nrow(indicators)
    n2=nrow(indicators2)

    # data frame version.
    difout=data.frame(cor_excl_all=numeric(),t=numeric(), p=numeric())

    if (length(tagged) != 0) {
        cor_excl_all=cor(rowMeans(indicators[, -excl_indicators],na.rm = T),rowMeans(indicators,na.rm=T),use = coruse)
        dif= psych::r.test(n=n, r12=cors[1,1], r13=cors[1,2], r23=cor_excl_all)

        # commented out is the old version with string output
        #if (dif$p<0.0001) p.out=paste0(", p<.0001") else  p.out=paste0(", p=",dif$p)
        #difout=paste0("Test for difference in correlations in self-report: t=",round(dif$t,2), p.out)
        #difreturn=difout

        difout[1,]=c(cor_excl_all,dif$t,dif$p)
        rownames(difout)[1]="Self-report: All vs Excluded"
    }

    if (length(tagged2) != 0) {

        cor_excl_all2=cor(rowMeans(indicators2[, -excl_indicators2],na.rm = T),rowMeans(indicators2,na.rm=T),use = coruse)
        dif2= psych::r.test(n=n2, r12=cors[2,1], r13=cors[2,2], r23=cor_excl_all2)
        #if (dif2$p<0.0001) p.out=paste0(", p<.0001") else  p.out=paste0(", p=",dif2$p)
        #difout2=paste0("Test for difference in correlations in other-report: t=",round(dif2$t,2), p.out)
        #difreturn=list(difout,difout2)
        difout[2,]=c(cor_excl_all2,dif2$t,dif2$p)
        rownames(difout)[2]="Informant: All vs Excluded"
    }

#     } else {
#         difreturn=NULL

    if (length(indicators2) != 0) {
        cors[2, 1] <- scalecor(indicators2, outcome, coruse)
        if (length(tagged2) != 0) {
            excl_indicators <- grep(paste0(tagged2, collapse = "|"), indicatornames)
            cors[2, 2] <- scalecor(indicators2[, -excl_indicators], outcome, coruse)

        }
    }







    # PLOTTING PLOT single indicator-outcome


    par(mar = c(5, 5, 4, 1) + 0.1)
    cor_indicators <- t(cor_indicators)

    # add CI-s?
    ci_add=FALSE
    if (ci=="estimate"|is.numeric(ci)) ci_add=TRUE

    if (length(indicators2) == 0) {
        # if indicators2 are missing, we have a different plottic scheme

        bp <- gplots::barplot2(cor_indicators, ylim = limits, col = cols[1], main = main1, cex.main = main2.cex, axes = F, beside = T, cex.names = multi)
        axis(2, cex.axis = multi)
        if (length(tagged) != 0 | length(tagged2) != 0) {
            # add extra legend if there is something to exclude

            legend(x = location1, star, bty = "n", cex = multi)
        }
    } else {

        # regular plotting
        bp <- gplots::barplot2(cor_indicators, ylim = limits, col = cols, main = main1, cex.main = main2.cex, axes = F, beside = T, cex.names = multi)
        axis(2, cex.axis = multi)
        # add extra legend if there is something to exclude

        if (length(tagged) != 0 | length(tagged2) != 0) {

            legend(x = location1, star, bty = "n", cex = multi)
        }

    }

    # prepare starring locations
    indicatorloc <- cor_indicators
    indicatorloc[] <- vapply(indicatorloc, awayzero, numeric(1), 0.007)  #http://stackoverflow.com/questions/8579257/r-applying-function-over-matrix-and-keeping-matrix-dimensions

    # tagged or not
    taggedm <- matrix("", nrow = nrow(cor_indicators), ncol = ncol(cor_indicators))
    if (length(tagged) != 0) {
        tagged_nr <- grep(paste0(tagged, collapse = "|"), indicatornames)
        taggedm[1, tagged_nr] <- "x"
    }
    if (length(tagged2) != 0) {
        tagged2_nr <- grep(paste0(tagged2, collapse = "|"), indicatornames)
        taggedm[2, tagged2_nr] <- "x"
    }

    text(bp, indicatorloc, taggedm, cex = multi)

    # labels
    title(main = main1, ylab = y, xlab = x, cex.main = multi, cex.lab = multi)


    # PLOT sumscore

    # add space to the right
    par(mar = c(5, 1, 4, 5) + 0.1)

    if (length(indicators2) == 0) {
        # if indicators2 are missing, we have a different plottic scheme

        # cors=cors[,,-2]
        sumplot=gplots::barplot2(cors[1, ], ylim = limits, col = cols[1], main = main2, cex.main = main2.cex, axes = F, beside = T,cex.names = multi, plot.ci=ci_add, ci.l=ci_l[1,], ci.u=ci_u[1,])
        # legend(x=location2,exclid[1],fill=cols )

        axis(4, cex.axis = multi)
        mtext(4, text = y, line = 3, cex = multi)

        # add star if correlation differencei is below 0.05
        if(length(tagged)<0 & difout[1,3]<0.05) text(sumplot[2,1], ci_u[1,2]+0.01, "*", cex = multi)



        # # add extra SONE legend if there is something to exclude if (length(tagged) != 0 | length(tagged2) != 0) {
        # legend(x = location1, star, bty = 'n', cex = multi) }
    } else {
        # regular plotting
        sumplot2=gplots::barplot2(cors, ylim = limits, col = cols, main = main2, cex.main = main2.cex, axes = F, beside = T,
            cex.names = multi, plot.ci=ci_add, ci.l=ci_l, ci.u=ci_u)
        legend(x = location2, indicatorsid, fill = cols, bty = "n", cex = multi)

        # if (length(tagged) != 0 | length(tagged2) != 0) { # add extra legend if there is something to exclude legend(x =
        # location1, c('', '', star), bty = 'n', cex = multi) }

        axis(4, cex.axis = multi)
        mtext(4, text = y, line = 3, cex = multi)

        print(sumplot2)
        print(ci_u)

        # add star if correlation difference is below 0.05
        if(length(tagged)<0 & difout[1,3]<0.05) text(sumplot2[1,2], ci_u[1,2]+0.01, "*", cex = multi)
        if(length(tagged2)<0 & difout[2,3]<0.05) text(sumplot2[2,2], ci_u[2,2]+0.01, "*", cex = multi)



    }
    on.exit(options(old), add = TRUE)

    # prepare output



    # cut indicators2 from output if not needed
    if (length(indicators2) == 0) {
        cors=cors[1,]
        if (ci=="estimate"|is.numeric(ci)) cors_ci=cors_ci[1:2,]
    }
    # return either cors_ci or cors. add diff result, if needed

    if (ci=="estimate"|is.numeric(ci)) {
        if (length(tagged) > 0| length(tagged2) > 0){
            out=list(cors_ci,difout)
            return(out)
        } else {
            return(cors_ci)

        }
    } else {
        if (length(tagged) > 0| length(tagged2) > 0){
            out=list(cors,difout)
            return(out)
        } else {
            return(cors)
        }
    }

}
