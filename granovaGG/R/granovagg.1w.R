#' Elemental Graphic Display for One-Way ANOVA
#'
#' Graphic to display data for a one-way analysis of
#' variance -- that is for unstructured groups. Also to help understand
#' how data play out in the context of the basic one-way model,
#' how the F statistic is generated for the data at hand,
#' etc. The graphic may be called 'elemental' or 'natural'
#' because it is built upon the central question that drives
#' one-way ANOVA (see details below).
#'
#' The one-way ANOVA graphic shows how the comparison of unstructured
#' groups, viz. their means, entails a particular linear combination (L.C.)
#' of the group means. In particular, we use the fact that the numerator of
#' the one-way F statistic, the mean square between (MS.B), is a linear combination
#' of the group means; each weight -- one for each group -- in the L.C. is (principally)
#' a function of the difference between the group's mean and the grand
#' mean, viz., (M_{j} - M..) where M_{j} denotes the jth group's mean, and M.. denotes
#' the grand mean. The L.C. can be written as a sum of products of the form
#' MS.B = Sum((1/df.B)(n_j (M_j - M..) M_j)) for j = 1...J.
#' The denominator of the F-statistic, MS.W
#' (mean square within), can be described as a 'scaling factor'. It is just the (weighted)
#' average of the variances of the J groups (j = 1 ... J). (n_{j}'s are group sizes.)
#' The differences (M_{j} - M..) are themselves the 'effects' in the analysis.
#' When the effects are plotted against the group means (the horizontal and
#' vertical axes) a straight line necessarily ensues. Group means are plotted as
#' triangles along this line. Once the means have been plotted, the data points
#' (jittered) for the groups are displayed (vertical axis) with respect to the
#' respective contrasts. Since the group means are just the fitted values in
#' one-way ANOVA, and the deviations of the scores within groups are the residuals
#' (subsetted by groups), the graphic can be seen as showing fitted vs. residual
#' values for the line that shows the locus of ordered group means -- from the smallest
#' on the left) the the largest (on the right). If desired, the aggregate of all
#' such residuals can be plotted (as a rug plot) on the right margin of the
#' graphic centered on the grand mean (large green dot in 'middle'). The use of
#' effects to locate groups this way yields what we term an 'elemental'
#' graphic because it is based on the central question that drives one-way
#' ANOVA.
#'
#' Note that groups need not have the same size, nor do data need to
#' reflect any particular distributional characteristics. Finally, the gray bars
#' (one for each group) at the bottom of the graphic show the relative sizes of
#' the group standard deviations with referene to the 'average' group s.d. (more precisely,
#' the square root of the MS.W). This 'average' corresponds to the thin white
#' line that runs horizontally across these bars.
#'
#' @param data Dataframe or vector. If a dataframe, the two or more columns
#'   are taken to be groups of equal size (whence \code{group} is NULL).  If
#'   \code{data} is a vector, \code{group} must be a vector, perhaps a factor,
#'   that indicates groups (unequal group sizes allowed with this option).
#' @param group Group indicator, generally a factor in case \code{data} is a
#'   vector.
#' @param h.rng Numeric; controls the horizontal spread of groups, default =
#'   1
#' @param v.rng Numeric; controls the vertical spread of points, default =
#'   1.
#' @param jj Numeric; sets horiz. jittering level of points. \code{jj} gets passed as the
#'   \code{amount} parameter to \code{\link{jitter}}.
#'   When \code{jj = NULL} (the default behavior), the degree of jitter will take on a sensible value.
#'   In addition, if pairs of ordered means are close to one another and \code{jj = NULL},
#'   the degree of jitter will default to the smallest difference between two adjacent contrasts.
#' @param dg Numeric; sets number of decimal points in output display, default = 2
#' @param resid Logical; displays marginal distribution of residuals (as a
#'   'rug') on right side (wrt grand mean), default = FALSE.
#' @param print.squares Logical; displays graphical squares for visualizing the F-statistic as a ratio
#'   of MS-between to MS-within
#' @param xlab Character; horizontal axis label, can be supplied by user, default = \code{"default_x_label"},
#'   which leads to a generic x-axis label ("Contrast coefficients based on group means").
#' @param ylab Character; vertical axis label, can be supplied by user, default = \code{"default_y_label"},
#'   which leads to a generic y-axis label ("Dependent variable (response)").
#' @param main Character; main label, top of graphic; can be supplied by user,
#'   default = \code{"default_granova_title"}, which will print a generic title for graphic.
#' @param plot.theme argument indicating a ggplot2 theme to apply to the
#'   graphic; defaults to a customized theme created for the one-way graphic
#' @param ... Optional arguments to/from other functions
#' @return Returns a plot object of class \code{ggplot}. The function also provides printed output including by-group
#'   statistical summaries and information about groups that might be overplotted (if applicable):
#'      \item{group}{group names}
#'      \item{group means}{means for each group}
#'      \item{trimmed.mean}{20\% trimmed group means}
#'      \item{contrast}{Contrasts (group main effects)}
#'      \item{variance}{variances}
#'      \item{standard.deviation}{standard deviations}
#'      \item{group.size}{group sizes}
#'      \item{overplotting information}{Information about groups that, due to their close means, may be overplotted}
#' @seealso \code{\link{granovagg.contr}},
#'   \code{\link{granovagg.ds}}, \code{\link{granovaGG}}
#'
#' @author Brian A. Danielak \email{brian@@briandk.com}\cr
#'   Robert M. Pruzek \email{RMPruzek@@yahoo.com}
#'
#' with contributions by:\cr
#'   William E. J. Doane \email{wil@@drdoane.com}\cr
#'   James E. Helmreich \email{James.Helmreich@@Marist.edu}\cr
#'   Jason Bryer \email{jason@@bryer.org}
#'
#' @references Fundamentals of Exploratory Analysis of Variance, Hoaglin D.,
#'   Mosteller F. and Tukey J. eds., Wiley, 1991.
#' @include granovagg.1w-helpers.R
#' @include shared-functions.R
#' @include theme-defaults.R
#' @keywords hplot htest
#' @example /demo/granovagg.1w.R
#' @references Wickham, H. (2009). Ggplot2: Elegant Graphics for Data Analysis. New York: Springer.
#' @references Wilkinson, L. (1999). The Grammar of Graphics. Statistics and computing. New York: Springer.
#' @import plyr
#' @import RColorBrewer
#' @import stats
#' @import utils
#' @export
granovagg.1w <- function(data,
                         group         = NULL,
                         h.rng         = 1,
                         v.rng         = 1,
                         jj            = NULL,
                         dg            = 2,
                         resid         = FALSE,
                         print.squares = TRUE,
                         xlab          = "default_x_label",
                         ylab          = "default_y_label",
                         main          = "default_granova_title",
                         plot.theme    = "theme_granova_1w",
                         ...
                )

{

  yy <- data

  CoerceHigherDimensionalDataToMatrix <- function(data) {
    if (is.data.frame(data)) {
      if (length(names(data)) > 1) {
        data <- CoerceToMatrix(data)
      }
    }

    return(data)
  }

  CoerceToMatrix <- function(x) {
    error.message <- "It looks like you've tried to pass in higher-dimension data AND a separate group indicator.
    If your data contains columns of equal numbers of observations, try re-calling granovagg.1w
    on your data while setting group = NULL"

    if (!is.null(group)) {
      message(error.message)
    }
    return(as.matrix(x))
  }

  yy <- CoerceHigherDimensionalDataToMatrix(yy)

  #Testing input data type
  mtx <- is.matrix(yy)
  if (!mtx) {
    # if this executes, yy is already a vector as handed in via data
    yr <- yy
    groupf<-factor(group)
  }

  #If yy is matrix sans colnames, need to give them some
  if(mtx & is.null(colnames(yy))) { #Note changes here;did not work before, and I thought LETTERS looked better (your thoughts?)
       dmyy2<-dim(yy)[2]
       cnms<-LETTERS[1:dmyy2]     #Note that numbers replaced by LETTERS if initial matrix yy does not have col. names
       colnames(yy)<-c(paste(" G",cnms))
  }  #1:dim(yy)[2]))

  # If yy is a matrix, the data represents a balanced case. The code below creates a yr (outcomes) vector and a groupf (factored groups) vector by repeating each of the column numbers (group) of the source matrix for each of the outcomes in that column.
  if(mtx) {
    group <- rep(1:ncol(yy), each = nrow(yy))
    groupf<-factor(group,labels=colnames(yy))
    yr <- as.vector(yy)
  }

  ngroups<-length(unique(groupf))

  # By this point, all data have been transformed to yr and groupf - two vectors of equal length. Think of them as the "score" and "group" columns of an n x 2 dataframe, where n = total number of observations.



  #
  # all we now care about are yr, a vector, and groupf, a vector of group names/labels
  #



  #Basic stats by group; 'stats' is matrix with rows corresponding to groups, columns for effect size contrasts, means, variances & sd's
  mvs <- function(x){c(mean(x),var(x),sd(x))}
  stats <- matrix(unlist(tapply(yr,groupf,mvs)),byrow=T,ncol=3)
  # stats will have as many rows as we have groups, one row per group

  #      [,1]   [,2]  [,3]
  # [1,] 3.75 12.917 3.594
  # [2,] 2.50  1.667 1.291
  # [3,] 6.50  1.667 1.291
  # [4,] 3.00 16.667 4.082


  groupn <- as.numeric(groupf)
  # [1] 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4

  #yrm is a vector of same length as 'groupn' with the appropriate group mean in each cell
  yrm <- stats[,1][groupn]

  # The second bracket notation is a way of repeatedly selecting from the stats matrix
  # > x <- c("a", "b", "c")
  # > y <- c(1,1,1,2,2,2,3,3,3)
  # > x[y]
  # [1] "a" "a" "a" "b" "b" "b" "c" "c" "c"


  #  [1] 3.75 3.75 3.75 3.75 2.50 2.50 2.50 2.50 6.50 6.50 6.50 6.50 3.00 3.00 3.00 3.00


  tabc <- table(groupf)
  #mn.n is mean groupsize
  mn.n <- mean(tabc)
  tabc.dm <- tabc/mn.n
  grandmean <- mean(yr)

  #Stats now has 6 cols, first is group size, second group mean minus grandmean,
  #third is weighted (by group size) mean, then mean, var, sd.

  stats <- cbind(tabc,stats[,1]-grandmean, tabc.dm * stats[,1],stats)

  #Creating x, y ranges for graph
  #Parameters h.rng, v.rng, jj for horizontal, vertical and jitter enabled.

  # generate the contrasts
  stats.vc <- yrm - grandmean

  # We now have four vectors, each of length n = number of observations.
  # yr contains the raw observations
  # groupf is the vector of factors
  # yrm is the vector of group means for each raw score
  # stats.vc is the vector of contrasts

  rng.v <- range(yr)
  rng.h <- h.rng*range(stats.vc)
  rng.vv <- c(rng.v[1]-v.rng*diff(rng.v), rng.v[2] + v.rng * diff(rng.v))
  rng.rt <- diff(rng.vv)/diff(rng.h)
  rng.sts <- range(stats[, 2])


  ammt <- (1/200) * diff(rng.sts)
  stats.vcj<-jitter(stats.vc, amount = ammt)

  #Reordering the stats matrix by the mean of each group
  statsro<-stats[order(stats[,4]),]

  #Calculation of squares etc.
  r.xy.sqd<-cor(yr,stats.vc)^2
  SS.tot<-(length(yr)-1)*var(yr)
  SS.bet<-r.xy.sqd*SS.tot
  df.b<-ngroups-1
  df.w<-length(yr)-1-df.b
  SS.w<-SS.tot-SS.bet
  MS.w<-SS.w/df.w
  MS.b<-SS.bet/df.b
  residuals<-round(yr-stats.vc,3)

  #sdw is standard deviation within, ie of residuals.
  sdw<-sd(residuals)*sqrt((length(yr)-1)/df.w)

  #This interval based on pooled standard error within.
  grandmean.pm.sdw<-c(grandmean-sdw,grandmean+sdw)
  grandmean.pm.sewR<-round(grandmean.pm.sdw,1-1)

    F.stat <- MS.b/MS.w

  p.F <- 1 - pf(F.stat, df.b,df.w)
  sqrF<-sqrt(F.stat)
  sqrs<-2*sqrt(MS.w)/rng.rt

  #Trimmed means marked and outputted if 'trmean = TRUE'
  trmd.mn<-tapply(yr,list(groupf), mean, tr=.2)

  gsum<-array(c(grandmean, df.b, df.w, MS.b, MS.w, F.stat, p.F, round(r.xy.sqd,3)))
  dimnames(gsum)<-list(c('Grandmean', 'df.bet', 'df.with', 'MS.bet', 'MS.with', 'F.stat', 'F.prob', 'SS.bet/SS.tot'))
  stats.out<-cbind(statsro[,1:4],round(as.matrix(trmd.mn[order(stats[,4])]),2),statsro[,5:6])

  dimnames(stats.out)[2]<-list(c('Size','Contrast Coef',
  "Wt'd Mean",'Mean', "Trim'd Mean" , 'Var.','St. Dev.'))
  out<-list(grandsum = round(gsum, 1), stats = round(stats.out, 1))

  if(FALSE){
           if(is.null(pt.lab) & !mtx & !is.null(rownames(yy))){pt.lab<-rownames(yy)}
           if(is.null(pt.lab) & !mtx & is.null(rownames(yy))){pt.lab<-1:length(yy)}
           if(is.null(pt.lab) & mtx){pt.lab<-paste(rep(1:dim(yy)[1],dim(yy)[2]),",", rep(1:dim(yy)[2],ea = dim(yy)[1]), sep="")}
           FALSEify(stats.vcj,yr,labels = pt.lab, ...)
           }


  AdaptVariablesFromGranovaComputations <- function() {
    result  <- list(data = data.frame(score             = yr,
                                      group             = groupf,
                                      group.mean        = yrm,
                                      contrast          = stats.vc
                           )
    )
    result$stats <- list(F.statistic               = F.stat,
                         SS.between                = SS.bet,
                         SS.within                 = SS.w,
                         df.between                = df.b,
                         df.within                 = df.w,
                         grand.mean                = grandmean,
                         square.side.length        = sqrs,
                         sd.within                 = sdw
    )
    result$residuals <- data.frame(within.group.residuals                   = residuals,
                                   within.1.sd.of.the.mean.of.all.residuals =
                                     ConvertBooleanValuesToResidualLabels(abs(residuals - grandmean) < sdw)
    )

    result$range.expansion <- list(horizontal.range.expansion = 1.25 * h.rng,
                                   vertical.range.expansion   = 0.2 * v.rng
                              )
    return(result)
  }

  ConvertBooleanValuesToResidualLabels <- function(boolean.vector) {
    label.vector                          <- as.character(boolean.vector)
    label.vector[label.vector == "TRUE"]  <- "Within 1 SDpooled"
    label.vector[label.vector == "FALSE"] <- "Outside 1 SDpooled"

    return(label.vector)
  }

  GetSummary <- function(owp) {
    # To appease R CMD Check
    score <- NULL
    contrast <- NULL

    return(
      ddply(owp$data, .(group), summarise,
        group              = unique(group),
        group.mean         = mean(score),
        trimmed.mean       = mean(score, trim = 0.2),
        contrast           = unique(contrast),
        variance           = var(score),
        standard.deviation = sd(score),
        maximum.score      = max(score),
        group.size         = length(score)
      )
    )
  }

  PrintGroupSummary <- function(data, digits.to.round) {
    # To appease R CMD Check
    maximum.score <- NULL

    groups <- subset(data, select = group)
    stats  <- subset(data, select = c(-group, -maximum.score))
    rounded.stats <- round(stats, digits = digits.to.round)
    output <- ReorderDataByColumn(cbind(groups, rounded.stats), "group.mean")
    message("\nBy-group summary statistics for your input data (ordered by group means)")

    return(
       print(output)
    )
  }

  GetGroupMeanLine <- function(owp) {
    return(
      data.frame(
        x     = min(owp$data$contrast),
        y     = min(owp$data$group.mean),
        xend  = max(owp$data$contrast),
        yend  = max(owp$data$group.mean),
        color = "blue"
      )
    )
  }

  GetGraphicalParameters <- function(owp) {
    .expanded.contrast.range <- range(owp$summary$contrast) * owp$range.expansion$horizontal.range.expansion
    .score.range             <- range(owp$data$score)
    .expanded.score.range    <- c(  min(.score.range) - (owp$range$vertical.range.expansion/2 * diff(.score.range)),
                                    max(.score.range) + (owp$range$vertical.range.expansion/2 * diff(.score.range))
                                )
    .score.range.distance    <- max(.expanded.score.range) - min(.expanded.score.range)
    .contrast.range.distance <- max(.expanded.contrast.range) - min(.expanded.contrast.range)
    .aggregate.y.breaks      <- c(owp$summary$group.mean, range(owp$data$score))
    .aggregate.x.breaks      <- c(owp$summary$contrast, 0)
    .vertical.percent        <- .score.range.distance / 100
    .horizontal.percent      <- .contrast.range.distance / 100
    .y.range                 <- c(
                                  min(.expanded.score.range) - (10 * .vertical.percent),
                                  max(.expanded.score.range) + (10 * .vertical.percent)
                                )
    .x.range                 <- c(min(.expanded.contrast.range) - (10 * .horizontal.percent),
                                  max(.expanded.contrast.range) + (10 * .horizontal.percent))
    .aspect.ratio            <- .contrast.range.distance / .score.range.distance

    return(list(
             expanded.contrast.range = .expanded.contrast.range,
             score.range.distance    = .score.range.distance,
             aggregate.x.breaks      = .aggregate.x.breaks,
             aggregate.y.breaks      = .aggregate.y.breaks,
             y.range                 = .y.range,
             x.range                 = .x.range,
             vertical.percent        = .vertical.percent,
             horizontal.percent      = .horizontal.percent,
             aspect.ratio            = .aspect.ratio,
             contrast.range.distance = .contrast.range.distance
           )
    )
  }

  GetSquareParameters <- function(owp) {
    return(
      list(
        x.center = max(owp$params$x.range) - (5.5 * (owp$params$horizontal.percent)),
        y.center = min(owp$params$y.range) + (5.5 * (owp$params$vertical.percent)),
        height   = 10 * owp$params$vertical.percent,
        width    = 10 * owp$params$horizontal.percent
      )
    )
  }

  GetMSbetweenColor <- function(owp) {
    if ((IsFSignificant(owp$model.summary)) == TRUE) {
      return(brewer.pal(n = 8, name = "Paired")[5])
    }

    else {
      return(brewer.pal(n = 8, name = "Paired")[2])
    }
  }

  GetMSwithinColor <- function(owp) {
    if ((IsFSignificant(owp$model.summary)) == TRUE) {
      return(brewer.pal(n = 8, name = "Paired")[6])
    }

    else {
      return(brewer.pal(n = 8, name = "Paired")[1])
    }
  }


  GetColors <- function() {
    strokeColors <- c(
      "Grand Mean"         = brewer.pal(n = 8, name = "Set1")[3],
      "Group Mean Line"    = brewer.pal(n = 8, name = "Paired")[2],
      "Within 1 SDpooled"  = "darkblue",
      "Outside 1 SDpooled" = "darkorange"
    )

    fillColors <- c(
      "MS-between"         = GetMSbetweenColor(owp),
      "MS-within"          = GetMSwithinColor(owp),
      "Group Means"        = brewer.pal(n = 8, name = "Paired")[8]
    )
    return(list(stroke = strokeColors, fill = fillColors))
  }

  GetAsTheOuterSquare <- function(owp, name.of.square) {
    return(
      data.frame(
        xmin  = owp$squares$x.center - (owp$squares$width / 2),
        xmax  = owp$squares$x.center + (owp$squares$width / 2),
        ymin  = owp$squares$y.center - (owp$squares$height / 2),
        ymax  = owp$squares$y.center + (owp$squares$height / 2),
        fill  = factor(paste(name.of.square))
      )
    )
  }

  GetMSbetweenAsTheInnerSquare <- function(owp) {
    return(
      data.frame(
        xmin = owp$squares$x.center - (owp$squares$width  / (2 * sqrt(1 / owp$stats$F.statistic))),
        xmax = owp$squares$x.center + (owp$squares$width  / (2 * sqrt(1 / owp$stats$F.statistic))),
        ymin = owp$squares$y.center - (owp$squares$height / (2 * sqrt(1 / owp$stats$F.statistic))),
        ymax = owp$squares$y.center + (owp$squares$height / (2 * sqrt(1 / owp$stats$F.statistic))),
        fill = factor(paste("MS-between"))
      )
    )
  }

  GetMSwithinAsTheInnerSquare <- function(owp) {
    return(
      data.frame(
        xmin = owp$squares$x.center - (owp$squares$width  / (2 * sqrt(owp$stats$F.statistic))),
        xmax = owp$squares$x.center + (owp$squares$width  / (2 * sqrt(owp$stats$F.statistic))),
        ymin = owp$squares$y.center - (owp$squares$height / (2 * sqrt(owp$stats$F.statistic))),
        ymax = owp$squares$y.center + (owp$squares$height / (2 * sqrt(owp$stats$F.statistic))),
        fill = factor(paste("MS-within"))
      )
    )
  }

  GetOuterSquare <- function(owp) {
    if (owp$stats$F.statistic > 1) {
      return(GetAsTheOuterSquare(owp, "MS-between"))
    }

    else {
      return(GetAsTheOuterSquare(owp, "MS-within"))
    }
  }

  GetInnerSquare <- function(owp) {
    if (owp$stats$F.statistic > 1) {
      return(GetMSwithinAsTheInnerSquare(owp))
    }

    else {
      return(GetMSbetweenAsTheInnerSquare(owp))
    }

  }

  GetModelSummary <- function(owp) {
    model <- lm(score ~ group, data = owp$data)

    return(summary(model))
  }

  GetSquaresData <- function(owp) {
    if (print.squares == TRUE) {
      vertical.position <- owp$outer.square$ymax + (2.0 * owp$params$vertical.percent)
    } else {
      vertical.position <- owp$outer.square$ymin
    }

    GetSquaresText(owp, vertical.position)
  }

  GetSquaresText <- function(owp, position) {
    test.statistic         <- owp$model.summary$fstatistic["value"]
    test.statistic.rounded <- round(test.statistic, digits = 2)

    if (length(levels(owp$data$group)) == 2) { # 2-group t-test case
      test.statistic.rounded = round(sqrt(test.statistic), digits = 2)
      text <- paste("t = ", test.statistic.rounded, sep = "")
    } else {
      text <- paste("F = ", test.statistic.rounded, sep = "")
    }

    return(
      data.frame(label     = text,
                 x         = owp$squares$x.center,
                 y         = position,
                 text.size = GetSquaresTextSize(test.statistic.rounded)
      )
    )
  }

  GetSquaresTextSize <- function(number) {
    if (number < 10) {
      return(2.5)
    }

    else {
      return(2.25)
    }

  }

  GetWithinGroupVariation <- function(owp) {
    lower.bound           <- min(owp$params$y.range)
    contrast              <- owp$summary$contrast
    standard.deviation    <- owp$summary$standard.deviation
    root.mean.square.variation <- sqrt(mean(owp$summary$variance))
    return(
      data.frame(
        x                  = contrast,
        ymin               = lower.bound,
        ymax               = lower.bound + RescaleVariationData(standard.deviation),
        baseline.variation = lower.bound + RescaleVariationData(root.mean.square.variation)
      )
    )
  }

  RescaleVariationData <- function(data) {
    scale.factor <- 1/2
    return(scale.factor * data)
  }

  GetBackgroundForGroupSizesAndLabels <- function(owp) {
    return(data.frame(ymin = max(owp$params$y.range) - 15 * owp$params$vertical.percent,
                      ymax = max(owp$params$y.range),
                      xmin = min(owp$params$x.range),
                      xmax = max(owp$params$x.range)
           )
    )
  }

  GetGroupSizes  <- function(owp) {
    return(data.frame(
             y           = max(owp$params$y.range) - (1 * owp$params$vertical.percent),
             x           = owp$overplot$contrast,
             label       = owp$overplot$group.size,
             overplotted = owp$overplot$overplotted,
             angle       = 90
           )
    )

  }

  GetGroupLabels <- function(owp) {
    return(data.frame(
             y           = max(owp$params$y.range) - (10 * owp$params$vertical.percent),
             x           = owp$overplot$contrast,
             label       = owp$overplot$group,
             overplotted = owp$overplot$overplotted,
             angle       = 90
           )
    )
  }

  AddOverplotInformation <- function(data, variable, tolerance) {
    more.than.two.groups <- length(data[[variable]]) > 2
    ordered.data <- ReorderDataByColumn(data, variable)

    if (more.than.two.groups) {
      overplotted  <- OverlapWarning(ordered.data[[variable]], tolerance)
    }
    else {
      overplotted <- rep(FALSE, times = length(data[[variable]]))
    }

    return(data.frame(ordered.data, overplotted))
  }

  ######## Plot Functions Below

  GroupMeanLine <- function(owp) {
    return(
      geom_segment(
        aes_string(
          x      = "x",
          y      = "y",
          xend   = "xend",
          yend   = "yend",
          color  = "factor('Group Mean Line')"
        ),
        alpha = I(1/2),
        data  = owp$group.mean.line
      )
    )
  }

  GroupMeansByContrast <- function(owp) {
    return(
      geom_point(
        aes_string(
          x     = "contrast",
          y     = "group.mean",
          fill  = "factor('Group Means')"
        ),
        size  = I(3),
        shape = 24,
        color = "black",
        alpha = 0.50,
        data  = owp$summary
      )
    )
  }

  Residuals <- function(owp, resid) {
    if (resid == TRUE) {
      return(
        geom_rug(
          aes_string(
            x     = "NULL",
            y     = "within.group.residuals",
            color = "factor(within.1.sd.of.the.mean.of.all.residuals)"
          ),
          alpha = I(1),
          data  = owp$residuals,
          sides = "l"
        )
      )
    }
  }

  SquaresForFstatistic <- function() {
    output <- NULL

    if (print.squares == TRUE) {
      output <- list(OuterSquare(), InnerSquare())
    }

    return(output)
  }

  OuterSquare <- function() {
    return(
      geom_rect(
        aes_string(
          xmin   = "xmin",
          xmax   = "xmax",
          ymin   = "ymin",
          ymax   = "ymax",
          fill   = "fill",
          color  = "NULL"
        ),
        data  = owp$outer.square
      )
    )
  }

  InnerSquare <- function() {
    return(
      geom_rect(
        aes_string(
          xmin   = "xmin",
          xmax   = "xmax",
          ymin   = "ymin",
          ymax   = "ymax",
          fill   = "fill",
          color  = "NULL"
        ),
        data  = owp$inner.square,
      )
    )
  }

  SquaresText <- function(owp) {
    return(
      geom_text(
        aes_string(
          x     = "x",
          y     = "y",
          label = "label"
        ),
        color = "grey20",
        size  = owp$squares.text$text.size,
        data  = owp$squares.text,
        vjust = ifelse(print.squares == TRUE, 0.5, -1)
      )
    )
  }

  WithinGroupVariation <- function(owp) {
    return(
      geom_linerange(
        aes_string(
          x      = "x",
          ymin   = "ymin",
          ymax   = "ymax"
        ),
        color = "grey30",
        size  = GetWithinGroupVariationSize(),
        data  = owp$variation
      )
    )
  }

  MaxWithinGroupVariation <- function(owp) {
    return(
      geom_linerange(
        aes_string(
          x      = "x",
          ymin   = "ymin",
          ymax   = "max(ymax)"
        ),
        color = "grey",
        size  = GetWithinGroupVariationSize(),
        data  = owp$variation
      )
    )
  }

  GetWithinGroupVariationSize <- function () {
    return(1.5)
  }

  BaselineWithinGroupVariation <- function(owp) {
    return(
      geom_hline(
        aes_string(
          yintercept = "baseline.variation"
        ),
        color = "white",
        size  = I(1/4),
        data  = owp$variation
      )
    )
  }

  ColorScale <- function(owp) {
    output <- scale_color_manual(
      values = owp$colors$stroke,
      name = ""
    )
    if(exists("guides")) {
      output <- scale_color_manual(
        values = owp$colors$stroke,
        name = "",
        guide = "legend"
      )
    }
    return(output)
  }

  FillScale <- function() {
    return(
      scale_fill_manual(
        values = owp$colors$fill,
        name = ""
      )
    )
  }

  XLabel <- function(xlab) {
    label.to.output <- xlab

    if ((!is.null(xlab)) && (xlab == "default_x_label")) {
      label.to.output <- "Contrast coefficients based on group means"
    }

    return(xlab(label.to.output))
  }

  YLabel <- function(ylab) {
    label.to.output <- ylab

    if ((!is.null(ylab)) && (ylab == "default_y_label")) {
      label.to.output <- "Dependent variable (response)"
    }

    return(ylab(label.to.output))
  }

  BackgroundForGroupSizesAndLabels <- function(owp) {
    return(
      geom_rect(
        aes_string(
          ymin  = "ymin",
          ymax  = "ymax",
          xmin  = "xmin",
          xmax  = "xmax"
        ),
        fill  = "white",
        data  = owp$label.background
      )
    )
  }

  GroupSizes  <- function(owp) {
    return(
      geom_text(
        aes_string(
          x     = "x",
          y     = "y",
          label = "label",
          angle = "angle"
        ),
        size  = 2.5,
        color = "grey10",
        hjust = 1,
        vjust = 0.5,
        data  = owp$group.sizes
      )
    )
  }

  NonOverplottedGroupLabels <- function(owp) {
    overplotted <- NULL # to appease R CMD check
    if (FALSE %in% owp$group.labels$overplotted) {
      return(
        geom_text(
          aes_string(
            x     = "x",
            y     = "y",
            label = "label",
            angle = "angle"
          ),
          size  = GetGroupLabelSize(),
          color = "grey50",
          hjust = 0.5,
          vjust = 0.5,
          data  = subset(
            owp$group.labels,
            overplotted == FALSE
          )
        )
      )
    }
  }

  OverplottedGroupLabels <- function(owp) {
    overplotted <- NULL # to appease R CMD check
    if (TRUE %in% owp$group.labels$overplotted) {
      return(
        geom_text(
          aes_string(
            x     = "x",
            y     = "y",
            label = "label",
            angle = "angle"
          ),
          size  = GetGroupLabelSize(),
          color = brewer.pal(n = 8, name = "Paired")[6],
          hjust = 0.5,
          vjust = 0.5,
          data  = subset(
            owp$group.labels,
            overplotted == TRUE)
        )
      )
    }
  }

  GetGroupLabelSize <- function() {
    return(3)
  }

  RotateXTicks <- function() {
    return(
      theme(
        axis.text.x = element_text(angle = 90)
      )
    )
  }

  ForceCoordinateAxesToBeEqual <- function(owp) {
    return(
      coord_fixed(
        ratio = owp$params$aspect.ratio
      )
    )
  }

  GetClassicTitle <- function () {
    return(
      paste("One-way ANOVA displaying",ngroups,"groups")
    )
  }

  PlotTitle <- function (main) {
    title.to.output <- main

    if ( !is.null(title.to.output) && (title.to.output == "default_granova_title")) {
      title.to.output <- GetClassicTitle()
    }

    return(
      ggtitle(title.to.output)
    )
  }

  RemoveSizeElementFromLegend <- function() {
    # return(scale_size_continuous(legend = FALSE))
    return(NULL)
  }

  ### Warning Function Below

  PrintOverplotWarning <- function(owp, digits.to.round) {
    # To appease R CMD Check
    overplotted <- NULL
    group.mean <- NULL
    contrast <- NULL

    if (TRUE %in% owp$overplot$overplotted) {
      overplotted.groups <- subset(
        owp$overplot,
        overplotted == TRUE,
        select = c(
          "group",
          "group.mean",
          "contrast"
        )
      )
      overplotted.groups <- transform(
        overplotted.groups,
        group.mean = round(
          group.mean,
          digits = digits.to.round
        ),
        contrast = round(
          contrast,
          digits = digits.to.round
        )
      )
      message("\nThe following groups are likely to be overplotted")
      print(overplotted.groups)
    }
  }

  PrintLinearModelSummary <- function(owp) {

    if (length(levels(owp$data$group)) == 2) {
      PrintTtest(owp$data[, c("score", "group")])
    } else {
      message("\nBelow is a linear model summary of your input data")
      print(owp$model.summary)
    }
  }

  PrintTtest <- function(data) {
    unstacked.data <- unstack(data, score ~ group)
    message("\nBelow is a t-test summary of your input data")
    print(
      t.test(unstacked.data[, 1],
             unstacked.data[, 2],
             var.equal = TRUE
      )
    )
  }

  # Pepare OWP object
  owp                       <- AdaptVariablesFromGranovaComputations()
  owp$summary               <- GetSummary(owp)
  owp$model.summary         <- GetModelSummary(owp)
  owp$group.mean.line       <- GetGroupMeanLine(owp)
  owp$params                <- GetGraphicalParameters(owp)
  owp$overplot              <- AddOverplotInformation(owp$summary, "contrast", 2*owp$params$horizontal.percent)
  owp$squares               <- GetSquareParameters(owp)
  owp$colors                <- GetColors()
  owp$outer.square          <- GetOuterSquare(owp)
  owp$inner.square          <- GetInnerSquare(owp)
  owp$squares.text          <- GetSquaresData(owp)
  owp$variation    <- GetWithinGroupVariation(owp)
  owp$label.background      <- GetBackgroundForGroupSizesAndLabels(owp)
  owp$group.labels          <- GetGroupLabels(owp)
  owp$group.sizes           <- GetGroupSizes(owp)
  PrintGroupSummary(owp$summary, dg)


  #Plot OWP object
  p <- InitializeGgplot_1w()
  p <- p + GrandMeanLine(owp)
  p <- p + GrandMeanPoint(owp)
  p <- p + ScaleX_1w(owp)
  p <- p + ScaleY_1w(owp)
  p <- p + JitteredScoresByGroupContrast(owp, jj)
  p <- p + GroupMeanLine(owp)
  p <- p + GroupMeansByContrast(owp)
  p <- p + Residuals(owp, resid)
  p <- p + MaxWithinGroupVariation(owp)
  p <- p + WithinGroupVariation(owp)
  p <- p + BaselineWithinGroupVariation(owp)
  p <- p + SquaresForFstatistic()
  p <- p + SquaresText(owp)
  p <- p + ColorScale(owp)
  p <- p + FillScale()
  p <- p + XLabel(xlab)
  p <- p + YLabel(ylab)
  p <- p + BackgroundForGroupSizesAndLabels(owp)
  p <- p + GroupSizes(owp)
  p <- p + NonOverplottedGroupLabels(owp)
  p <- p + OverplottedGroupLabels(owp)
  p <- p + RotateXTicks()
  p <- p + Theme(plot.theme)
  p <- p + ForceCoordinateAxesToBeEqual(owp)
  p <- p + PlotTitle(main)
  p <- p + RemoveSizeElementFromLegend()
  PrintOverplotWarning(owp, dg)
  PrintLinearModelSummary(owp)

  return(p)
}
