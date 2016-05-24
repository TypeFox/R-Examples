#' Produce forest plots for categorical covariates.
#' @return Forestplots by ggplot2.
#' @param x A cdtafit object from \link{fit}.
#' @param graph An optional numeric value indicating which forest to plot(s) to graph. Valid values are:0 - for no graph, 1 - yielding a forest plot of the
#' sensitivity and specificity with a 95 percent exact confidence intervals, 2 - yielding a forest plot of the posterior study-specific sensitivity and specificity
#' and the marginal mean sensitivy and specificity and their corresponding 95 percent credible intervals, 3 - yielding a combination of 1 and 2 in one plot, and NULL(default) - yielding plots of
#' 1, 2 and 3.
#' @param title.1 An optional string indicating the title of graph 1.
#' @param title.2 An optional string indicating the title of graph 2.
#' @param title.3 An optional string indicating the title of graph 3.
#' @param width An optional numeric value to adjust the dogding position. The default is 0.2.
#' @param shape.1 An optional numeric value(0-255) indicating the symbol to plot in graph 1. The defualt is 19 which is a solid circle. See \link[graphics]{points} for more details.
#' @param size.1 An optional positive numeric value indicating the size of symbols in graph 1. The defualt is 2.5.
#' @param shape.2 An optional numeric value(0-255) indicating the symbol to plot in graph 2. The defualt is 8 which is a star. See \link[graphics]{points} for more details.
#' @param size.2 An optional positive numeric value indicating the size of symbols in graph 2. The defualt is 2.5.
#' @param shape.O An optional numeric value(0-255) indicating the symbol representing the posterior marginal mean in graph 2. The defualt is 19 which is a solid circle. See \link[graphics]{points} for more details.
#' @param size.O An optional numeric value indicating the size of symbols representing the posterior marginal means in graph 2.
#' @param cols.1 An optional string vector specifying colours of shapes in graph 1.
#' @param cols.2 An optional string vector specifying colours of shapes in graph 2.
#' @param digits An optional positive value to control the number of digits to print when printing numeric values. The default is 3.
#' @param ... other \link[rstan]{stan} options.
#' @examples
#' \dontrun{
#' fit1 <- fit(data=telomerase,
#'              SID = "ID",
#'              copula="fgm",
#'              iter = 400,
#'              warmup = 100,
#'              seed=1,
#'              cores=1)
#'
#' plot(fit1)
#' }
#' @references {Watanabe S (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular
#' Learning Theory. Journal of Machine Learning Research, 11, 3571-3594.}
#' @references {Vehtari A, Gelman A (2014). WAIC and Cross-validation in Stan. Unpublished, pp. 1-14.}
#' @export
#' @author Victoria N Nyaga <victoria.nyaga@outlook.com>
forestplot.cdtafit <- function(x,
                           title.1=NULL,
                           title.2=NULL,
                           title.3=NULL,
                           graph=NULL,
                           width=0.2,
                           shape.1=19,
                           size.1=2.5,
                           shape.2=8,
                           size.2=2.5,
                           shape.O=9,
                           size.O=3.5,
                           cols.1=NULL,
                           cols.2=NULL,
                           digits=3,
                         ...){
 #==================================================================================#
    ID <- NULL
    SID <- NULL
	COV <- NULL
	Lower <- NULL
	Upper <- NULL
	mean.p <- NULL
	Lower.p <- NULL
	Upper.p <- NULL

 #==================================================================================#
    df <- prep.data(data=x@data,
                    SID = x@SID,
                    formula.se=x@modelargs$formula.se,
                    formula.sp=x@modelargs$formula.sp,
                    formula.omega=x@modelargs$formula.omega)$data


    df$Dis <- df$TP + df$FN
    df$NDis <- df$TN + df$FP
    df$ID <- 1:nrow(df)

    #Long format data
    df <- df[,names(df) %in% c("FN", "FP")==FALSE]
    event <- reshape2::melt(df, id.vars = names(df)[names(df) %in% c("TP", "TN")==FALSE])
    total <- reshape2::melt(df, id.vars = names(df)[names(df) %in% c("Dis", "NDis")==FALSE])

    names(event)[grepl(pattern="variable", x=names(event))] <- "Parameter"
    names(event)[grepl(pattern="value", x=names(event))] <- "Event"
    event$Parameter <- factor(event$Parameter, labels=c("Sensitivity", "Specificity"))

    names(total)[grepl(pattern="value", x=names(total))] <- "Total"
    names(total)[grepl(pattern="variable", x=names(total))] <- "Parameter"
    total$Parameter <- factor(total$Parameter, labels=c("Sensitivity", "Specificity"))

    df<- merge(event, total, by=intersect(names(event), names(total)))

    df$p <- df$Event/df$Total

    formula <- x@modelargs$formula.se

    df1 <- stats::get_all_vars(formula, df)

    if (ncol(df1) > 1){
        df$COV <- df1[,2]
        covname <- names(df1[2])
    } else{
        df$COV <- 1
        covname <- "Intercept"
    }
#====================Exact CI========================================================#
    df$Lower <- NA
    df$Upper <- NA

    for (r in 1:nrow(df)){
        bt <- stats::binom.test(df$Event[r], df$Total[r], conf.level = 0.95)$conf.int
        attr(bt, 'conf.level') <- NULL
        df[r,c("Lower", "Upper")] <- bt
    }
#=======================Extract Model Parameters ===================================#
    sm <- rstan::summary(x@fit,...)

    p <- data.frame(sm$summary[grepl('p_i', rownames(sm$summary)), c("mean", "2.5%", "97.5%")])
    names(p) <- c("mean.p", "Lower.p", "Upper.p")
    p$ID <- rep(1:(nrow(p)/2), each=2)

    df1 <- stats::get_all_vars(formula, data=prep.data(data=x@data,
                                                  SID = x@SID,
                                                  formula.se=x@modelargs$formula.se,
                                                  formula.sp=x@modelargs$formula.sp,
                                                  formula.omega=x@modelargs$formula.omega)$data)

    p$SID <- factor(rep(as.character(df1[,1]), each=2))

    if (ncol(df1) > 1){
        p$COV <- factor(rep(as.character(df1[,2]), each=2))
    } else{
        p$COV <- 1
    }

    p$Parameter <- rep(c("Sensitivity", "Specificity"), length.out=nrow(p))

    MU <- sm$summary[grepl('MU', rownames(sm$summary)), c("mean", "2.5%", "97.5%")]
    colnames(MU) <- c("mean.p", "Lower.p", "Upper.p")
    rownames(MU) <- NULL
    MU <- data.frame(MU)
    MU$Parameter <- rep(c("Sensitivity", "Specificity"), each=nrow(MU)/2)
    MU$ID <- c((nrow(p)/2 + 1):((nrow(p)+nrow(MU))/2))
    MU$SID <- "Overall"

    if (ncol(df1) > 1){
        cov <- c(as.character(levels(df1[,2])))
        MU$COV <- rep(cov, times=2)
    } else {
        MU$COV <- 1
    }

    P <- plyr::rbind.fill(p, MU)
#==============================Combine Data and Model parameters===========================#
    df <- merge(df, p, by=intersect(names(df), names(p)))
    df <- df[order(df$ID),]
#=======================================PLOTTING =============================================#
    if (ncol(df1) > 1) {
        legend <- "bottom"
    }else{
        legend <- "none"
    }

    if (ncol(df1) < 2){
        if (is.null(cols.1)){
            cols.1 <- "magenta"
            }
        if (is.null(cols.2)){
            cols.2 <- "blue"
            }
    }
#=====================================        DATA  ============================================#
    if (is.null(title.1)) title.1 <- paste("Plot of study-specific sensitivity and specificity by\n",
                                        x@SID,
                                       " and ",
                                       covname,
                                       ": mean and 95% exact CI",sep='')

    dodge <- ggplot2::position_dodge(width)
    if (ncol(df1) > 1){
        g1 <-  ggplot2::ggplot(data=df, ggplot2::aes(x = stats::reorder(SID, -ID), y = p, ymax= max(p)*1.05, group = COV, colour = COV)) +
            ggplot2::coord_flip() +
            ggplot2::theme_bw() +
            ggplot2::facet_grid( ~ Parameter) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower, ymax=Upper),size=0.75,
                          width=0,
                          colour="black",
                          position=dodge) +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=10),
                  axis.text.y=ggplot2::element_text(size=10),
                  axis.title.x=ggplot2::element_text(size=10),
                  legend.title=ggplot2::element_text(size=10),
                  legend.position = legend,
                  legend.direction = 'horizontal',
                  legend.text=ggplot2::element_text(size=10)) +
            ggplot2::geom_point(size=size.1, shape=shape.1, position=dodge) +
            ggplot2::scale_x_discrete(name=x@SID) +
            ggplot2::scale_colour_discrete(name=covname) +
            ggplot2::scale_y_continuous(name="", limits=c(0,1)) +
            ggplot2::ggtitle(title.1)
    }else{
        g1 <-  ggplot2::ggplot(data=df, ggplot2::aes(x = stats::reorder(SID, -ID), y = p, ymax = max(p)*1.05)) +
            ggplot2::coord_flip() +
            ggplot2::theme_bw() +
            ggplot2::facet_grid( ~ Parameter) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower, ymax=Upper),size=0.75,
                          width=0,
                          colour="black",
                          position=dodge) +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=10),
                  axis.text.y=ggplot2::element_text(size=10),
                  axis.title.x=ggplot2::element_text(size=10),
                  legend.title=ggplot2::element_text(size=10),
                  legend.position = legend,
                  legend.direction = 'horizontal',
                  legend.text=ggplot2::element_text(size=10)) +
            ggplot2::geom_point(size=size.1, shape=shape.1, position=dodge, colour=cols.1) +
            ggplot2::scale_x_discrete(name=x@SID) +
            ggplot2::scale_y_continuous(name="", limits=c(0,1)) +
            ggplot2::ggtitle(title.1)
    }

#=====================================   Posterior                   ========================#
    if (is.null(title.2)) title.2 <- paste("Plot of study-specific posterior sensitivity and specificity by\n",
                                        x@SID,
                                       " and ",
                                       covname,
                                       ": marginal mean and 95% CI",sep='')
    if (ncol(df1) > 1){
        g2 <-  ggplot2::ggplot(data=P, ggplot2::aes(x = stats::reorder(SID, -ID), y = mean.p, ymax = max(mean.p)*1.05, group = COV, colour = COV)) +
            ggplot2::coord_flip() +
            ggplot2::theme_bw() +
            ggplot2::facet_grid( ~ Parameter) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower.p, ymax=Upper.p),
                          size=0.75, width=0,
                          colour="black",
                          position=dodge) +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=10),
                  axis.text.y=ggplot2::element_text(size=10),
                  axis.title.x=ggplot2::element_text(size=10),
                  legend.title=ggplot2::element_text(size=10),
                  legend.position = legend,
                  legend.direction = 'horizontal',
                  legend.text=ggplot2::element_text(size=10)) +
           ggplot2:: geom_point(data=P[P$ID %in% df$ID,,],ggplot2::aes(x = SID, y = mean.p),
                       size=size.1,
                       shape=shape.1,
                       position=dodge) +
            ggplot2::geom_point(data=P[P$ID %in% df$ID==FALSE,], ggplot2::aes(x = SID, y = mean.p),
                       size=size.O,
                       shape=shape.O,
                       position=dodge) +
            ggplot2::scale_x_discrete(name=x@SID) +
            ggplot2::scale_colour_discrete(name=covname) +
            ggplot2::scale_y_continuous(name="", limits=c(0,1)) +
            ggplot2::ggtitle(title.2)
    } else{
        g2 <-  ggplot2::ggplot(data=P, ggplot2::aes(x = stats::reorder(SID, -ID), y = mean.p, ymax=max(mean.p)*1.05)) +
            ggplot2::coord_flip() +
            ggplot2::theme_bw() +
            ggplot2::facet_grid( ~ Parameter) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower.p, ymax=Upper.p),
                          size=0.75, width=0,
                          colour="black",
                          position=dodge) +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=10),
                  axis.text.y=ggplot2::element_text(size=10),
                  axis.title.x=ggplot2::element_text(size=10),
                  legend.title=ggplot2::element_text(size=10),
                  legend.position = legend,
                  legend.direction = 'horizontal',
                  legend.text=ggplot2::element_text(size=10)) +
            ggplot2::geom_point(data=P[P$ID %in% df$ID,],ggplot2::aes(x = SID, y = mean.p),
                       size=size.1,
                       shape=shape.1,
                       position=dodge,
                       colour=cols.1) +
            ggplot2::geom_point(data=P[P$ID %in% df$ID==FALSE,], ggplot2::aes(x = SID, y = mean.p),
                       size=size.O,
                       shape=shape.O,
                       position=dodge,
                       colour=cols.1) +
            ggplot2::scale_x_discrete(name=x@SID) +
            ggplot2::scale_y_continuous(name="", limits=c(0,1)) +
            ggplot2::ggtitle(title.2)
    }

#=====================================COMBINED Data + Posterior ========================#
    df <- plyr::rbind.fill(df, MU)
    dodge <- ggplot2::position_dodge(width + 0.2)

    if (is.null(title.3)) title.3 <- paste("Plot of study-specific sensitivity and specificity  by\n",
                                           x@SID,
                                       " and ",
                                       covname,
                                       ": marginal mean and 95% CI",sep='')
    if (ncol(df1) > 1){
        g3 <-  ggplot2::ggplot(data=df, ggplot2::aes(x = stats::reorder(SID, -ID), y = p, ymax= max(p)*1.05, group = COV, colour = COV)) +
            ggplot2::coord_flip() +
            ggplot2::theme_bw() +
            ggplot2::facet_wrap( ~ Parameter) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower, ymax=Upper),
                          size=2.5,
                          width=0,
                          colour="gray", position=dodge) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower.p, ymax=Upper.p),
                          size=1,
                          width=0,colour="black",
                          position=dodge) +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=10),
                  axis.text.y=ggplot2::element_text(size=10),
                  axis.title.x=ggplot2::element_text(size=10),
                  legend.title=ggplot2::element_text(size=10),
                  legend.position = legend,
                  legend.direction = 'horizontal',
                  legend.text=ggplot2::element_text(size=10)) +
            ggplot2::scale_x_discrete(name=x@SID) +
            ggplot2::scale_colour_discrete(name=covname) +
            ggplot2::scale_y_continuous(name="", limits=c(0,1)) +
            ggplot2::geom_point(data=df[df$ID<=nrow(x@data),],
                       ggplot2::aes(x = SID, y = p),
                       size=size.1,
                       shape=shape.1,
                       position=dodge) +
            ggplot2::geom_point(data=df[df$ID<=nrow(x@data),],
                       ggplot2::aes(x = SID, y = mean.p),
                       size=size.2,
                       shape=shape.2,
                       position=dodge) +
            ggplot2::geom_point(data=df[df$ID>nrow(x@data),],
                       ggplot2::aes(x = SID, y = mean.p),
                       size=size.O, shape=shape.O,
                       position=dodge) +
            ggplot2::ggtitle(title.3)
    }else{
        g3 <-  ggplot2::ggplot(data=df, ggplot2::aes(x = stats::reorder(SID, -ID), y = p, ymax= max(p)*1.05)) +
            ggplot2::coord_flip() +
            ggplot2::theme_bw() +
            ggplot2::facet_wrap( ~ Parameter) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower, ymax=Upper),
                          size=2.5,
                          width=0,
                          colour="gray", position=dodge) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower.p, ymax=Upper.p),
                          size=1,
                          width=0,
                          position=dodge,
                          colour="black") +
            ggplot2::theme(axis.text.x=ggplot2::element_text(size=10),
                  axis.text.y=ggplot2::element_text(size=10),
                  axis.title.x=ggplot2::element_text(size=10),
                  legend.title=ggplot2::element_text(size=10),
                  legend.position = legend,
                  legend.direction = 'horizontal',
                  legend.text=ggplot2::element_text(size=10)) +
            ggplot2::scale_x_discrete(name=x@SID) +
            ggplot2::scale_colour_discrete(name=covname) +
            ggplot2::scale_y_continuous(name="", limits=c(0,1)) +
            ggplot2::geom_point(data=df[df$ID<=nrow(x@data),],
                       ggplot2::aes(x = SID, y = p),
                       size=size.1,
                       shape=shape.1,
                       position=dodge,
                       colour=cols.1) +
            ggplot2::geom_point(data=df[df$ID<=nrow(x@data),],
                       ggplot2::aes(x = SID, y = mean.p),
                       size=size.2,
                       shape=shape.2,
                       position=dodge,
                       colour=cols.2) +
            ggplot2::geom_point(data=df[df$ID>nrow(x@data),],
                       ggplot2::aes(x = SID, y = mean.p),
                       size=size.O, shape=shape.O,
                       position=dodge,
                       colour=cols.2) +
            ggplot2::ggtitle(title.3)
    }

#========================================= PLOTTING ================================#

    if (is.null(graph)){
        if (grDevices::dev.interactive()) {
            grDevices::dev.new()
            print(g1)
            }
        if (grDevices::dev.interactive()) {
            grDevices::dev.new()
            print(g2)
            }
        if (grDevices::dev.interactive()){
            grDevices::dev.new()
            print(g3)
            }
    } else if(graph==1){
        if (grDevices::dev.interactive()) {
            grDevices::dev.new()
            print(g1)
            }
    } else if(graph==2){
        if (grDevices::dev.interactive()) {
            grDevices::dev.new()
            print(g2)
        }
    } else if(graph==3){
        if (grDevices::dev.interactive()) {
            grDevices::dev.new()
            print(g3)
        }
    }

out <- list(G1 =g1, G2=g2, G3=g3)

return(out)
}

