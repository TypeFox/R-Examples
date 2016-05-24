#' @title Generic function for Enhanced Tipping Point Displays
#'
#' @description Generic function for Enhanced Tipping Point Displays


#' @param \dots Additional arguments, see \code{\link{TippingPoint.default}}, \code{\link{TippingPoint.formula}} for more details.
#' @export
#' @seealso \code{\link{TippingPoint.default}}, \code{\link{TippingPoint.formula}}.
#' @references 1. Liublinska, V. and Rubin, D.B. Enhanced Tipping-Point Displays. In JSM Proceedings, Section on Survey Research Methods, San Diego, CA: American Statistical Association. 3861-3686 (2012).
#' @references 2. Liublinska, V. & Rubin, D.B. Sensitivity analysis for a partially missing binary outcome in a two-arm randomized clinical trial. Stat Med 33, 4170-85 (2014).
#' @references 3. Liublinska, V. (May, 2013) Sensitivity Analyses in Empirical Studies Plagued with Missing Data. PhD Dissertation, Harvard University, \url{https://dash.harvard.edu/handle/1/11124841}
#' @references 4. \url{https://sites.google.com/site/vliublinska/research}

#' @examples
#' TippingPoint(outcome=tippingdata$binary, treat= tippingdata$treat,
#'   group.infor=TRUE, plot.type = "estimate",ind.values = TRUE,
#'   impValuesT  = NA,  impValuesC = NA,
#'   summary.type = "density", alpha = 0.95, S=1.5, n.grid = 100,
#'   HistMeanT = c(0.38,0.4), HistMeanC =  c(0.2,0.55))

TippingPoint <- function(...) UseMethod("TippingPoint")




#' @title Default method for TippingPoint
#'
#' @description The default method for enhanced tipping point displays.

#' @param outcome A  numeric vector of the outcomes, a binary or continuous outcome.
#' @param treat A (non-NA) numeric vector of treatment group.
#' @param group.infor A logical, whether to display the group information.
#' @param plot.type A character, one of "estimate", "p.value" or "both" indicating which one should be represented by a heat-map layer.
#' @param summary.type A character, how to summarize the joint posterior distribution of imputed outcomes for  treated and controls, one of "density", "credible.region" or "convex.hull".  see \code{\link[ggplot2]{geom_density2d}},  \code{\link[stats]{mahalanobis}},  \code{\link[ggplot2]{geom_polygon}}, \code{\link[bayesSurv]{credible.region}} for more details.
#' @param alpha A numeric between 0-1, with alpha  of points in Convex hull, 1-alpha removed by Machalanobis distance. It also specifies the probabilities for credible regions used in \code{\link[bayesSurv]{credible.region}}, in this case, alpha should be above 0.5 and below 1. The default value is 0.95.
#' @param HistMeanT A numeric vector or NULL, historical values or proportions for the treatment group.
#' @param HistMeanC A numeric vector or NULL, historical values or proportions for the control group.
#' @param ind.values A logical, whether or not to display values in heat-map layer.
#' @param impValuesT NA or imputed values for the treatment group, see \code{\link{imputedata}} for more details.
#' @param impValuesC NA or imputed values for the control group, see \code{\link{imputedata}} for more details.
#' @param impValuesColor NA or imputed colors correspond to the columns in impValuesT or impValuesC. The default colors are from \emph{Set1} in \strong{RColorBrewer} allowing up to 9. Specify explicitly if need more colors.  See \code{\link[RColorBrewer]{display.brewer.all}} for more colors.
#' @param show.points A logical, whether to show the points for imputed values.
#' @param point.size  Size of points for imputed values.
#' @param point.shape Shape of points for imputed values.
#' @param S A integer indicating range of plotting, the default value is 3.
#' @param n.grid A integer, number of points in the grid, only for continuous case, the default is 150.
#' @param \dots Additional arguments
#' @import ggplot2 RColorBrewer
#' @importFrom reshape2 melt
#' @importFrom bayesSurv credible.region
#' @importFrom stats sd prop.test cov mahalanobis pt quantile var as.formula get_all_vars terms
#' @importFrom grDevices chull
#' @importFrom utils head
#' @seealso \code{\link{TippingPoint}}, \code{\link{TippingPoint.formula}}.
#'
#'
#' @export
#' @examples
#' #  See more details in vignette using:
#' #  vignette("TippingPoint")
#' TippingPoint(outcome=tippingdata$binary,treat= tippingdata$treat,
#'  plot.type = "p.value",ind.values = TRUE,
#'  impValuesT  = imputedata[,c("MAR_T2","MCAR_T2")],
#'  impValuesC = imputedata[,c("MAR_C2","MCAR_C2")],
#'  impValuesColor = RColorBrewer::brewer.pal(8,"Accent")[c(4,6)],
#'  summary.type = "credible.region", alpha = 0.95,
#'  S=1.5, n.grid = 100, HistMeanT = c(0.38,0.4), HistMeanC =  c(0.2,0.55))


#' @method TippingPoint default
TippingPoint.default <- function(outcome, treat, group.infor=FALSE,
                         plot.type = c("estimate", "p.value", "both"),
                         summary.type = c("density", "credible.region", "convex.hull"),
                         alpha = 0.95, HistMeanT = NULL, HistMeanC = NULL,
                         ind.values = FALSE,
                         impValuesT  = NA,  impValuesC = NA,impValuesColor=NA,
                         show.points = TRUE, point.size=1 , point.shape=19, S=3, n.grid =150, ...) {

  if (is.numeric(outcome)) {
    Yobs<-outcome
  } else stop("Outcome should be a numeric vector!", call. = FALSE)

  if (all(!(treat %in% c(0,1))))
      stop("treat should be a numeric vector with value of 0 or 1 !", call. = FALSE)

  stopifnot(is.logical(group.infor),is.logical(ind.values),is.logical(show.points))

  if (alpha>1||alpha<0) stop ("alpha must be a single number between 0 and 1", call. = FALSE)

  plot.type<-match.arg(plot.type)
  summary.type<-match.arg(summary.type)

  if (length(unique(Yobs[!is.na(Yobs)])) > 2) TPtype <- "continuous" else TPtype <- "binary"


  ## TP display for continuous outcome
  if (TPtype == "continuous") {

    X_t_mis = seq.int(mean(Yobs[treat==1], na.rm = TRUE) - S*sd(Yobs[treat==1], na.rm = TRUE), mean(Yobs[treat==1], na.rm = TRUE) + S*sd(Yobs[treat==1], na.rm = TRUE),length.out = n.grid)
    X_c_mis = seq.int(mean(Yobs[treat==0], na.rm = TRUE) - S*sd(Yobs[treat==0], na.rm = TRUE), mean(Yobs[treat==0], na.rm = TRUE) + S*sd(Yobs[treat==0], na.rm = TRUE),length.out = n.grid)

    # t-test
    Welch.t.test = function(meanYt_mis,meanYc_mis, treat, Yobs){
      Nt = sum(treat==1)
      Nc = sum(treat==0)
      N = length(treat)
      L = length(meanYt_mis)
      Ntmis = sum(is.na(Yobs)&treat==1)
      Ntobs = Nt - Ntmis
      meanYobsT = mean(Yobs[treat==1], na.rm = TRUE)
      Ncmis = sum(is.na(Yobs)&treat==0)
      Ncobs = Nc - Ncmis
      meanYobsC = mean(Yobs[treat==0], na.rm = TRUE)
      #p.val = vector(mode = "numeric", length =L)
      p.val = array(NA, L)
      for (i in 1:L){
        d = (meanYobsT*Ntobs + meanYt_mis[i]*Ntmis)/Nt - (meanYobsC*Ncobs + meanYc_mis[i]*Ncmis)/Nc
        Vart = (var(Yobs[treat==1], na.rm = TRUE)*(Ntobs-1) + Ntobs*Ntmis/Nt*(meanYobsT - meanYt_mis[i])^2 )/(Ntobs)
        Varc = (var(Yobs[treat==0], na.rm = TRUE)*(Ncobs-1) + Ncobs*Ncmis/Nc*(meanYobsC - meanYc_mis[i])^2 )/(Ncobs)
        #CHECK THE DEGREES OF FREEDOM!!!
        Welch.DF = (Vart/Nt + Varc/Nc)^2/((Vart/Nt)^2/(Ntobs)+(Varc/Nc)^2/(Ncobs))
        t.stat = d/sqrt(Vart/Nt + Varc/Nc)
        p.val[i] = 2*pt(-abs(t.stat), df=Welch.DF)
      }
      p.val
    }

    effect.size.mean = function(meanYt_mis,meanYc_mis,treat,Yobs){
      (meanYt_mis*sum(is.na(Yobs[treat==1])) + sum(Yobs[treat==1], na.rm = TRUE))/sum(treat==1) -
        (meanYc_mis*sum(is.na(Yobs[treat==0])) + sum(Yobs[treat==0], na.rm = TRUE))/sum(treat==0)}

    p.values <- outer(X_t_mis, X_c_mis, FUN = Welch.t.test, treat=treat, Yobs=Yobs)
    theta <- outer(X_t_mis, X_c_mis, FUN = effect.size.mean, treat, Yobs)

    colnames(theta) = X_c_mis
    rownames(theta) = X_t_mis

    df <- melt(theta)
    names(df) <- c("MisYt", "MisYc", "value")
    df = data.frame(df)
    df$p.value = melt(p.values)[,3]

    #### Basic plot
    p <- ggplot(df, aes(x=MisYt, y=MisYc))

    if (!is.null(HistMeanT)) p<-p+ geom_vline(xintercept = HistMeanT, colour = "purple",  lwd = 1.01, alpha = 0.7, lty = 1)

    if (!is.null(HistMeanC)) p<-p+geom_hline(yintercept = HistMeanC,  colour = "purple", lwd = 1.01, alpha = 0.7, lty = 1)

    if (plot.type == "estimate") {
      p <- p + geom_tile(aes(fill = value)) +  labs(title = "Treatment effect")
      if(min(df$value)>0)
        p<-p+scale_fill_gradient("estimate",low="white", high="orange", space="Lab")
      if(max(df$value)<0)
        p<-p+scale_fill_gradient("estimate",low="steelblue", high="white", space="Lab")
      if(max(df$value)>0 & min(df$value)<0)
        p<-p+ scale_fill_gradient2("estimate",low="steelblue", high="orange", space="Lab")
    }

    if (plot.type == "p.value"){
      p <- p + geom_tile(aes(fill = p.value)) + labs(title = "P-value from a hypothesis test")
      p <- p + scale_fill_gradient2("p.value",low="white", high="darkolivegreen", space="Lab", limits = c(0,1))
      if (sum(df$p.value <= 0.05) > 0) {
        p <- p  + geom_tile(data = subset(df, p.value <= 0.05), fill = "darkred", alpha = 0.1) + #color = "white", size = 0.1
          stat_contour(aes(z = p.value), breaks = c(0.05), color = "darkred")
      }
    }

    if (plot.type == "both"){
      p <- p + geom_tile(aes(fill = value)) + labs(title = "Treatment effect and significant p-value in red grid if any")
      if(min(df$value)>0)
        p <- p + scale_fill_gradient("estimate",low="white", high="orange", space="Lab")
      if(max(df$value)<0)
        p <- p + scale_fill_gradient("estimate",low="steelblue", high="white", space="Lab")
      if(max(df$value)>0 & min(df$value)<0)
        p <- p +  scale_fill_gradient2("estimate",low="steelblue", high="orange", space="Lab")
      if (sum(df$p.value <= 0.05) > 0) {
        p <- p  + geom_tile(data = subset(df, p.value <= 0.05), fill = "darkred", alpha = 0.1) + #color = "white", size = 0.1,
          stat_contour(aes(z = p.value), breaks = c(0.05), color = "darkred")
      }
    }


    p <- p + labs(x = "Average outcome for nonrespondents in treatment group", y = "Average outcome for nonrespondents in control group")

    # Observed averages
    p <- p + geom_hline(yintercept = mean(Yobs[treat==0], na.rm = TRUE), colour = "darkblue", lty = 2, alpha = 0.7) +
      geom_vline(xintercept = mean(Yobs[treat==1], na.rm = TRUE), colour = "darkblue", lty = 2, alpha = 0.7)

    ## Adding minimum and maximum on the plot (only for continuous outcome)
    ContMinMax = check.range(range(Yobs[treat==0], na.rm = TRUE),X_c_mis)
    if (length(ContMinMax))
      p <- p + geom_hline(yintercept = range(Yobs[treat==0], na.rm = TRUE), colour = "dark blue", alpha = 0.5)

    TreatMinMax = check.range(range(Yobs[treat==1], na.rm = TRUE),X_t_mis)
    if (length(TreatMinMax))
      p <- p + geom_vline(xintercept = range(Yobs[treat==1], na.rm = TRUE), colour = "dark blue", alpha = 0.5)


    if (group.infor) {
      cat("\nGroup Information:\n\n")
      cat("Groups                    \tTreatment \tControl\n")
      cat("Size                     \t",format_cat(sum(treat)),"\t", format_cat(sum(treat==0)),"\n")
      cat("Number of nonrespondents \t",format_cat(sum(is.na(Yobs[treat==1]))),"\t", format_cat(sum(is.na(Yobs[treat==0]))),"\n")
      cat("% of nonrespondents      \t",format_cat(sum(is.na(Yobs[treat==1]))/sum(treat)),"\t", format_cat(sum(is.na(Yobs[treat==0]))/sum(treat==0)),"\n")
      cat("Observed average response\t",format_cat(mean(Yobs[treat==1], na.rm = TRUE)),"\t", format_cat(mean(Yobs[treat==0], na.rm = TRUE)),"\n")
      cat("Observed min response    \t",format_cat(min(Yobs[treat==1], na.rm = TRUE)),"\t", format_cat(min(Yobs[treat==0], na.rm = TRUE)),"\n")
      cat("Observed max response    \t",format_cat(max(Yobs[treat==1], na.rm = TRUE)),"\t", format_cat(max(Yobs[treat==1], na.rm = TRUE)),"\n")
    }

  } # continuous


  ## TP display for binary outcome
  if (TPtype == "binary"){

    X_t_mis = seq(0,sum(is.na(Yobs)&treat==1),1)
    X_c_mis = seq(0,sum(is.na(Yobs)&treat==0),1)

    proportion.test = function(sumYt_mis,sumYc_mis, treat, Yobs){
      L = length(sumYt_mis)
      p.val = array(NA,L)
      for (i in 1:L)
        p.val[i] = prop.test(c(sumYt_mis[i] + sum(Yobs*treat, na.rm = TRUE),sumYc_mis[i] + sum(Yobs*(1-treat), na.rm = TRUE)),
                             c(sum(treat==1),sum(treat==0)))$p.value
      p.val
    }

    effect.size.sum = function(sumYt_mis,sumYc_mis,treat,Yobs){
      (sumYt_mis + sum(Yobs[treat==1], na.rm = TRUE))/sum(treat==1) -
        (sumYc_mis + sum(Yobs[treat==0], na.rm = TRUE))/sum(treat==0)}

    p.values <- outer(X_t_mis, X_c_mis, FUN = proportion.test, treat, Yobs)
    theta <- outer(X_t_mis, X_c_mis, FUN = effect.size.sum, treat, Yobs)

    colnames(theta) = X_c_mis
    rownames(theta) = X_t_mis

    df <- melt(theta)
    names(df) <- c("MisYt", "MisYc", "value")
    df = data.frame(df)

    df$p.value = melt(p.values)[,3]

    # Basic plot
    p <- ggplot(df, aes(x=MisYt, y=MisYc))

    if (!is.null(HistMeanT)) {
      HistValueT = HistMeanT*sum(treat==1) - sum(Yobs[treat==1], na.rm = TRUE)
      HistValueT = check.range(HistValueT,X_t_mis)
      if (!is.null(HistValueT))
        p <- p +geom_vline(xintercept = HistValueT, colour = "purple",  lwd = 1.01, alpha = 0.7, lty = 1)
    }

    if (!is.null(HistMeanC)){
      HistValueC = HistMeanC*sum(treat==0) - sum(Yobs[treat==0], na.rm = TRUE)
      HistValueC = check.range(HistValueC,X_c_mis)
      if (!is.null(HistValueC))
        p <- p +geom_hline(yintercept = HistValueC, colour = "purple",  lwd = 1.01, alpha = 0.7, lty = 1)
    }

    if (plot.type == "estimate"){
      ## If interested in an estimate
      p <- p + geom_tile(aes(fill = value)) + labs(title = "Treatment effect")
      if(min(df$value)>0)
        p <- p + scale_fill_gradient("estimate",low="white", high="orange", space="Lab") #color gradietn
      if(max(df$value)<0)
        p <- p + scale_fill_gradient("estimate",low="steelblue", high="white", space="Lab") #color gradietn
      if(max(df$value)>0 & min(df$value)<0)
        p <- p +  scale_fill_gradient2("estimate",low="steelblue", high="orange", space="Lab") #color gradietn
    }

    if (plot.type == "p.value") {
      ## If interested in an p-value
      p <- p + geom_tile(aes(fill = p.value)) + labs(title = "P-value from a hypothesis test")
      p <- p + scale_fill_gradient2("p.value", low="white", high="darkolivegreen", space="Lab", limits = c(0,1)) #color gradietn

      if (sum(df$p.value > 0.05) > 0) p <- p +
        geom_tile(data = subset(df, p.value > 0.05), alpha = 0, colour = "white")
      if (sum(df$p.value <= 0.05) > 0) p <- p +
        geom_tile(data = subset(df, p.value <= 0.05), color = "darkred", alpha = 0)

    }

    if (plot.type == "both"){
      ## If interested in an estimate and p-value together
      p <- p + geom_tile(aes(fill = value)) + labs(title = "Treatment effect and significant p-value in red grid if any")
      if(min(df$value)>0)
        p <- p + scale_fill_gradient("estimate",low="white", high="orange", space="Lab") #color gradietn
      if(max(df$value)<0)
        p <- p + scale_fill_gradient("estimate",low="steelblue", high="white", space="Lab") #color gradietn
      if(max(df$value)>0 & min(df$value)<0)
        p <- p  + scale_fill_gradient2("estimate",low="steelblue", high="orange", space="Lab") #color gradietn

      if (sum(df$p.value > 0.05) > 0) p <- p +
          geom_tile(data = subset(df, p.value > 0.05), alpha = 0, colour = "white")

      if (sum(df$p.value <= 0.05) > 0)  p <- p +
          geom_tile(data = subset(df, p.value <= 0.05), colour = "darkred", alpha = 0)
    }

    ## Displaying individual values in the grid
    if (ind.values == TRUE&&(plot.type != "p.value")) p <- p + geom_text(aes(x = MisYt, y = MisYc, label = round(value,2)), color = "black", alpha = 0.6, size = 3)

    if (ind.values == TRUE&&(plot.type == "p.value")) p <- p + geom_text(aes(x = MisYt, y = MisYc, label = round(p.value,2)), color = "black", alpha = 0.6, size = 3)


    p <- p + labs(x = "Number of successes among nonrespondents in treatment group", y = "Number of successes among nonrespondents in control group")

    # Adding observed proportions
    p <- p + geom_hline(yintercept = mean(Yobs[treat==0], na.rm = TRUE)*sum(is.na(Yobs[treat==0])), colour = "dark blue", lty = 2, alpha = 0.7) +
      geom_vline(xintercept = mean(Yobs[treat==1]*sum(is.na(Yobs[treat==1])), na.rm = TRUE), colour = "dark blue", lty = 2, alpha = 0.7)

    if (group.infor) {
      cat("\nGroup Information:\n\n")

      cat("Groups                    \tTreatment \tControl\n")
      cat("Size                     \t",format_cat(sum(treat)),"\t", format_cat(sum(treat==0)),"\n")
      cat("Number of nonrespondents \t",format_cat(sum(is.na(Yobs[treat==1]))),"\t", format_cat(sum(is.na(Yobs[treat==0]))),"\n")
      cat("% of nonrespondents      \t",format_cat(sum(is.na(Yobs[treat==1]))/sum(treat)),"\t", format_cat(sum(is.na(Yobs[treat==0]))/sum(treat==0)),"\n")
      cat("Observed proportion      \t",format_cat(mean(Yobs[treat==1], na.rm = TRUE)),"\t", format_cat(mean(Yobs[treat==0], na.rm = TRUE)),"\n")
    }


} # binary

  ########### Adding imputed values
  # Options: Rectangles representing credibility intervals, alpha level
  #          Convex hull around alpha % of points, with (1-alpha)% removed by Machalanobis distance
  #          Fitted density using bivariate normal kernel


  if (!all(is.na(impValuesT))) {
    if (all(is.na(impValuesC)))
        stop("No simulations for control group are specified.")

    n.methods = dim(impValuesT)[2]
    MeanTtrx <- vector("list", n.methods)


    if (!all(is.na(impValuesColor))) Col.pal<-impValuesColor else {
      Col.pal = brewer.pal(9,"Set1")
      if (n.methods<3) {
        Col.pal<-Col.pal[c(5,6)]
      } else  Col.pal<-Col.pal[1:n.methods]

    }


    for (i in 1:n.methods){

      if (summary.type == "density"){
        MeanTtrx[[i]] = cbind(impValuesT[,i],impValuesC[,i])

        p = p +
          geom_density2d(data = data.frame(MeanTtrx[[i]], row.names = NULL), aes(X1,  X2),
                         col = Col.pal[i], alpha = 1, bins = 6, size = 0.7)#, fill = Col.pal[i]), size = 1.05
      }

      if (summary.type == "convex.hull"){

        RemoveOutliers <- function(Mtrx){
          Sx <- cov(Mtrx)
          MahDist <- mahalanobis(Mtrx, colMeans(Mtrx), Sx)
          Mtrx = Mtrx[MahDist < quantile(MahDist,alpha),]
          Mtrx
        }

        MeanTtrx[[i]] = RemoveOutliers(cbind(impValuesT[,i],impValuesC[,i]))

        ######## Plotting convex hull of the set of points
        hpts <- chull(MeanTtrx[[i]])
        hpts <- c(hpts, hpts[1])
        p = p +
          geom_polygon(data = data.frame(MeanTtrx[[i]][hpts,], row.names = NULL),
                       aes(X1,  X2), col = Col.pal[i], alpha = 0.3, fill = Col.pal[i])
      }

      if (summary.type == "credible.region"){

        MeanTtrx[[i]] = cbind(impValuesT[,i],impValuesC[,i])
        Rect = bayesSurv::credible.region(MeanTtrx[[i]], probs = alpha)[[1]]

        p = p +
          geom_rect(aes_string(
            xmin = Rect[1,1],
            ymin = Rect[1,2],
            xmax = Rect[2,1],
            ymax = Rect[2,2]
          ), fill = Col.pal[i], alpha = 0, color = Col.pal[i], size = 1)

      }

      if (show.points) {
        if (TPtype == "continuous") p = p +geom_point(size = point.size, shape= point.shape, data = data.frame(X1=impValuesT[,i], X2=impValuesC[,i]), aes(X1,  X2), col = Col.pal[i])
        if (TPtype == "binary") p = p +geom_jitter(size = point.size, shape= point.shape, data = data.frame(X1=impValuesT[,i], X2=impValuesC[,i]), aes(X1,  X2), col = Col.pal[i])
      }
    }
  }
  p
}






#' @title TippingPoint.formula
#'
#' @description The formula method for enhanced tipping point displays.

#' @param formula A formula of the form outcome ~ treat.
#' @param data A data.frame containing the variables in the formula.
#' @param \dots Additional arguments, see details in \code{\link{TippingPoint.default}}.
#' @seealso \code{\link{TippingPoint}}, \code{\link{TippingPoint.default}}.
#' @export
#' @method TippingPoint formula
#' @examples
#' #  See more details in vignette using:
#' #  vignette("TippingPoint")
#' TippingPoint(binary~treat, data=tippingdata,
#'   plot.type = "both", ind.values = TRUE,
#'   impValuesT  = imputedata[,c("MAR_T2","MCAR_T2")],
#'   impValuesC = imputedata[,c("MAR_C2","MCAR_C2")],
#'   impValuesColor =c("red","blue"),
#'   point.size=0.8,point.shape = 15,
#'   summary.type = "convex.hull", alpha = 0.95, S=1.5, n.grid = 100,
#'   HistMeanT = c(0.38,0.4), HistMeanC =  c(0.2,0.55))



TippingPoint.formula <- function(formula, data, ...) {

  if (missing(formula) || (length(formula) != 3L) ||
      (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect", call. = FALSE )

  formula<-as.formula(formula)
  data <- as.data.frame(data)

  outcome<-unlist(get_all_vars(formula,data)[1])
  treat<-unlist(get_all_vars(formula,data)[2])
  TippingPoint.default(outcome,treat,...)

}





## define global variables for cran submission.
utils::globalVariables(c("MisYt", "MisYc", "value", "p.value", "X1", "X2"))
