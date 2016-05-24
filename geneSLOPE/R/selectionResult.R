#' selectionResult class
#'
#' A result of applying SLOPE to matrix of SNPs obtained by
#' clumping produced. Result of function \code{\link{select_snps}}
#'
#' @details Always a named list of eighteen elements
#' \enumerate{
#' \item \code{X} numeric matrix, consists of one snp representative for each clump
#' selected by SLOPE
#' \item \code{effects} numeric vector, coefficients in linear model build on
#' snps selected by SLOPE
#' \item \code{R2} numeric, value of R-squared in linear model build on
#' snps selected by SLOPE
#' \item \code{selectedSNPs} which columns in matrix \code{X_all}
#' are related to snps selected by SLOPE
#' \item \code{y} selectedClumps list of numeric vectors, which columns in SNP matrix
#' \code{X_all} are related to clump members selected by SLOPE
#' \item \code{lambda} numeric vector, lambda values used by SLOPE procedure
#' \item \code{y} numeric vector, phenotype
#' \item \code{clumpRepresentatives} numeric vector, which columns in SNP matrix \code{X_all}
#' are related to clumps representatives
#' \item \code{clumps} list of numeric vectors, which columns in SNP matrix
#' \code{X_all} are related to clump members
#' \item \code{X_info} data.frame, mapping information about SNPs from .map file.
#' Copied from the result of clumping procedure
#' \item \code{X_clumps} numeric matrix, consists of one snp representative for each clump
#' \item \code{X_all} numeric matrix, all the snps that passed screening procedure
#' \item \code{selectedSnpsNumbers} numeric vector, which rows of \code{X_info}
#' data.frame are related to snps that were selected by SLOPE
#' \item \code{clumpingRepresentativesNumbers} numeric vector, which rows of \code{X_info}
#' data.frame are related to snps that are clump represenatives
#' \item \code{screenedSNPsNumbers} numeric vector, which rows of \code{X_info}
#' data.frame are related to snps that passed screening
#' \item \code{numberOfSnps} numeric, total number of SNPs before screening procedure
#' \item \code{pValMax} numeric, p-value used in screening procedure
#' \item \code{fdr} numeric, false discovery rate used by \code{\link{SLOPE}}
#' }
#' @seealso \code{\link{screeningResult}} \code{\link{clumpingResult}}
#' \code{\link{select_snps}} \code{\link[SLOPE]{SLOPE}}
#' @name selectionResult
NULL


#' Print selectionResult class object
#'
#' @param x selectionResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @return Nothing.
#' @export
#'
#' @method print selectionResult
print.selectionResult <- function(x, ...){
  cat("Object of class selectionResult\n")
  cat("$X: numeric matrix\n")
  cat("\t", nrow(x$X), " rows\n")
  cat("\t", ncol(x$X), " columns\n")
  cat("$effects: numeric vector of length ", length(x$effects), "\n")
  cat("$R2: ", x$R2, "\n")
  cat("$selectedSNPs: numeric vector of length",
      length(x$selectedSNPs), "\n")
  cat("$selectedClumps: list of vectors of length",
      length(x$selectedClumps), "\n")
  cat("$lambda: numeric vector of length",
      length(x$lambda), "\n")
  cat("$y: numeric vector\n")
  cat("$X_clump: Matrix after clumping\n")
  cat("\t", nrow(x$X_clump), " rows\n")
  cat("\t", ncol(x$X_clump), " columns\n")
  cat("$X_all: Matrix before clumping\n")
  cat("\t", nrow(x$X_all), " rows\n")
  cat("\t", ncol(x$X_all), " columns\n")
  cat("$X_info: Information about snps\n")
  cat("\t", nrow(x$X_info), " rows\n")
  cat("\t", ncol(x$X_info), " columns\n")
  cat("$clumpRepresentatives: numeric vector of length",
      length(x$clumpRepresentatives), "\n")
  cat("$clumps: list of numeric vectors of length",
      length(x$clumps), "\n")
  cat("$selectedSnpsNumbers: numeric vector of length",
      length(x$selectedSnpsNumbers), "\n")
  cat("$clumpingRepresentativesNumbers: numeric vector of length",
      length(x$clumpingRepresentativesNumbers), "\n")
  cat("$screenedSNPsNumbers: numeric vector of length",
      length(x$screenedSNPsNumbers), "\n")
  cat("$numberOfSnps: number of SNPs before screening:", x$numberOfSnps, "\n")
  cat("$pValMax: p-value threshold: ", x$pValMax, "\n")
  cat("$fdr: false discovery rate: ", x$fdr, "\n")
}

#' Summary selectionResult class object
#'
#' @param object selectionResult class object
#' @param clumpNumber number of clump to be summarized
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method summary selectionResult
summary.selectionResult <- function(object, clumpNumber = NULL, ...){
    if(is.null(clumpNumber)) {
    lambda_diffs <- diff(object$lambda)
    if(any(lambda_diffs==0)){
      kink <- which.min(lambda_diffs==0)
    } else {
      kink <- length(object$lambda)
    }
    cat("Object of class selectionResult\n")
    cat(length(object$selectedSNPs), "snps selected out of",
        length(object$clumpRepresentatives), "clump representatives\n")
    cat("Effect size for selected snps (absolute values)\n")
    cat("\tMin: ", min(abs(object$effects)), "\n")
    cat("\tMean: ", mean(abs(object$effects)), "\n")
    cat("\tMax: ", max(abs(object$effects)), "\n")
    cat("R square of the final model: ", object$R2, "\n")
    cat("Kink value: ", kink, "\n")
    } else {
      if(length(object$selectedClumps)<clumpNumber){
        stop("Number of selected clumps is smaller than ", clumpNumber)
      }
      cat("Summary of", clumpNumber, "selected clump\n")
      print(object$X_info[object$screenedSNPsNumbers[object$selectedClumps[[clumpNumber]]],])
    }
}


#' Plot selectionResult class object
#'
#' @param x selectionResult class object
#' @param chromosomeNumber optional parameter, only selected chromosome will be plotted
#' @param clumpNumber optional parameter, only SNPs from selected clump will be plotted
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
plot.selectionResult <- function(x, chromosomeNumber=NULL, clumpNumber=NULL, ...){
  chromosome <- snp <- val <- clump <- representatives <- NULL #to remove CRAN's NOTE
  if(length(x$selectedSNPs)==0){
    message("No SNPs selected by SLOPE")
    return(NULL)
  }
  if(!is.null(x$X_info)){
    plot.data <- create_slope_plot_data(x)
    granice <- aggregate(x$X_info[,3], list(x$X_info[,1]), max)
    granice_max <- cumsum(granice$x)
    granice$x <- c(0,head(cumsum(granice$x),-1))

    if(!is.null(chromosomeNumber)) {
      plot.data <- subset(plot.data, plot.data$chromosome%in%chromosomeNumber)
      if(nrow(plot.data)==0) {
        message("No SNPs selected in chromosme ", chromosomeNumber)
        return(NULL)
      }
      plot.data$clump <- as.factor(plot.data$clump)
      if(nrow(plot.data[plot.data$representatives,])==0) {
        p <- ggplot(plot.data)
      } else {
        p <- ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = clump, size = 6),
                                            plot.data[plot.data$representatives,])
      }
      p +
        geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=representatives,
                         color=clump)) +
        ggtitle(expression(atop(bold("SLOPE selection result"),
                                atop(italic("Dots indicate clump representatives. Colors indicate different clumps"), "")))) +
        ylab("% of variance explained") + scale_y_continuous() +
        xlab("Genome") +
        scale_x_continuous(limits=c(min(granice$x[chromosomeNumber]),
                                    max(granice_max[chromosomeNumber])),
                           breaks=rowMeans(cbind(granice$x, granice_max)),
                           labels=granice$Group.1,
                           minor_breaks=c(granice$x, max(granice_max))) +
        scale_alpha_manual(guide=FALSE, values = c(0.5, 1)) +
        scale_color_discrete(guide=FALSE, "Clump") +
        scale_size_area(guide=FALSE, max_size = 4) +
        slope_result_theme
    } else if(!is.null(clumpNumber)) {
      plot.data <- subset(plot.data, plot.data$clump%in%clumpNumber)
      if(nrow(plot.data)==0 | nrow(plot.data[plot.data$representatives,])==0) {
        message("No SNPs selected in clump ", clumpNumber)
        return(NULL)
      }
      plot.data$clump <- as.factor(plot.data$clump)
      ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = clump, size = 6),
                                     plot.data[plot.data$representatives,]) +
        geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=representatives,
                         color=clump)) +
        ggtitle(expression(atop(bold("SLOPE selection result"),
                                atop(italic("Dots indicate clump representatives"), "")))) +
        ylab("% of variance explained") + scale_y_continuous() +
        xlab("Genome") +
        scale_x_continuous(limits=c(min(granice$x[plot.data$chromosome]),
                                    max(granice_max[plot.data$chromosome])),
                           breaks=rowMeans(cbind(granice$x, granice_max)),
                           labels=granice$Group.1,
                           minor_breaks=c(granice$x, max(granice_max))) +
        scale_alpha_manual(guide=FALSE, values = c(0.5, 1)) +
        scale_color_discrete(guide=FALSE, "Clump") +
        scale_size_area(guide=FALSE, max_size = 4) +
        slope_result_theme
    } else {
      plot.data$clump <- as.factor(plot.data$clump)
      ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = clump, size = 6),
                                     plot.data[plot.data$representatives,]) +
        geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=representatives,
                         color=clump)) +
        ggtitle(expression(atop(bold("SLOPE selection result"),
                                atop(italic("Dots indicate clump representatives. Colors indicate different clumps"), "")))) +
        ylab("% of variance explained") +
        xlab("Genome") +
        scale_x_continuous(expand = c(0,0),
                           limits=c(0, max(granice_max)+1),
                           breaks=rowMeans(cbind(granice$x, granice_max)),
                           labels=granice$Group.1,
                           minor_breaks=c(granice$x, max(granice_max))) +
        scale_y_continuous(expand = c(0,0),limits=c(0, 1.1*max(plot.data$val))) +
        scale_alpha_manual(guide=FALSE, values = c(0.5, 1)) +
        scale_color_discrete(guide=FALSE, "Clump") +
        scale_size_area(guide=FALSE, max_size = 4) +
        slope_result_theme
    }


    } else {
      plot.data <- NULL
      for(i in 1L:length(x$selectedClumps)){
        plot.data <- rbind(plot.data,
                           cbind(x$selectedSnpsClumpingNumbers[unlist(x$selectedClumps[[i]])],
                                 i, abs(x$effects[i])/2))
      }
      plot.data <- data.frame(plot.data)
      colnames(plot.data) <- c("snp", "clump", "val")
      rownames(plot.data) <- NULL
      ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = "red", size = 6),
                                     plot.data[representatives,]) +
        geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=val/4)) +
        ylab("") +
        xlab("SNP number") +
        scale_alpha_continuous(guide=FALSE) +
        scale_color_discrete(guide=FALSE) +
        scale_size_area(guide=FALSE) +
        theme(panel.background=element_blank(),
              panel.grid.major.y=element_line(colour = "grey80"),
              panel.grid.minor.y=element_line(colour = "grey90"),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank())
    }
}

slope_result_theme <- theme(panel.background=element_blank(),
                        panel.grid.major.y=element_line(colour = "grey80"),
                        panel.grid.minor.y=element_line(colour = "grey90"),
                        panel.grid.major.x=element_blank(),
                        panel.grid.minor.x=element_line(colour = "grey70", linetype = "dotted", size=0.5),
                        axis.ticks.x=element_blank(),
                        legend.text = element_text(size=15),
                        legend.position="bottom",
                        legend.key =element_rect(fill="white"))

