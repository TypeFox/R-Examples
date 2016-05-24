#' clumpingResult class
#'
#' A result of procedure for snp clumping produced by \code{\link{clump_snps}}
#'
#' @details Always a named list of eleven elements
#' \enumerate{
#' \item \code{X} numeric matrix, consists of one snp representative for each clump
#' \item \code{y} numeric vector, phenotype
#' \item \code{SNPnumber} numeric vector, which columns in SNP matrix \code{X_all}
#' are related to clumps representatives
#' \item \code{SNPclumps} list of numeric vectors, which columns in SNP matrix
#' \code{X_all} are related to clump members
#' \item \code{X_info} data.frame, mapping information about SNPs from .map file.
#' Copied from the result of screening procedure.
#' \item \code{selectedSnpsNumbers} numeric vector, which rows of \code{X_info}
#' matrix are related to selected clump representatives
#' \item \code{X_all} numeric matrix, all the snps that passed screening procedure
#' \item \code{numberOfSnps} numeric, total number of SNPs before screening procedure
#' \item \code{selectedSnpsNumbersScreening} numeric vector, which rows of \code{X_info}
#' data.frame are related to snps that passed screening
#' \item \code{pVals} numeric vector, p-values from marginal tests for each snp
#' \item \code{pValMax} numeric, p-value used in screening procedure
#' }
#' @seealso \code{\link{screeningResult}} \code{\link{clump_snps}}
#' @name clumpingResult
NULL

#' Print clumpingResult class object
#'
#' @param x clumpingResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method print clumpingResult
print.clumpingResult <- function(x, ...){
  cat("Object of class clumpingResult\n")
  cat("$X: numeric matrix\n")
  cat("\t", nrow(x$X), " rows\n")
  cat("\t", ncol(x$X), " columns\n")
  cat("$y: numeric phenotype vector of length", length(x$y), "\n")
  cat("$SNPnumber: list with snp representatives for clumps \n")
  cat("\t[", paste(head(x$SNPnumber), collapse=","), "..., ]\n")
  cat("$SNPclumps: list of length", length(x$SNPclumps), " containing numeric vectors\n")
  cat("$X_info: data.frame\n")
  cat("\t", nrow(x$X_info), " rows\n")
  cat("\t", ncol(x$X_info), " columns\n")
  cat("$selectedSnpsNumbers: numeric vector of length" ,
      length(x$selectedSnpsNumbers), "\n")
  cat("$X_all: numeric matrix\n")
  cat("\t", nrow(x$X_all), " rows\n")
  cat("\t", ncol(x$X_all), " columns\n")
  cat("$numberOfSnps: ", x$numberOfSnps, "\n")
  cat("$selectedSnpsNumbersScreening: numeric vector of length" ,
      length(x$selectedSnpsNumbersScreening), "\n")
  cat("$pVals: numeric vector of length ", length(x$pVals), "\n")
  cat("$pValMax: ", x$pValMax, "\n")
}

#' Summary clumpingResult class object
#'
#' @param object clumpingResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method summary clumpingResult
summary.clumpingResult <- function(object, ...){
  cat("Object of class clumpingResult\n")
  cat(ncol(object$X_all), " SNPs grouped in ", length(object$SNPclumps), " clumps\n")
  cat("Mean clump size ", mean(unlist(lapply(object$SNPclumps, length))), "\n")
  cat("Min clump size ", min(unlist(lapply(object$SNPclumps, length))), "\n")
  cat("Max clump size ", max(unlist(lapply(object$SNPclumps, length))), "\n")
}


#' Plot clumpingResult class object
#'
#' @param x clumpingResult class object
#' @param chromosomeNumber optional parameter, only selected chromosome will be plotted
#' @param clumpNumber optional parameter, only SNPs from selected clump will be plotted
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
plot.clumpingResult <- function(x, chromosomeNumber=NULL, clumpNumber=NULL, ...){
  if(!is.null(x$X_info)){
    chromosome <- snp <- val <- clump <- representatives <- NULL #to remove CRAN's NOTE
    plot.data <- create_clumping_plot_data(x)
    granice <- aggregate(x$X_info[,3], list(x$X_info[,1]), max)
    granice_max <- cumsum(granice$x)
    granice$x <- c(0,head(cumsum(granice$x),-1))

    if(!is.null(chromosomeNumber)){
      plot.data <- subset(plot.data, plot.data$chromosome%in%chromosomeNumber)
      if(nrow(plot.data)==0) {
        message("No SNPs selected in chromosme ", chromosomeNumber)
        return(NULL)
      }
      ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = "red", size = 1),
                                     data=plot.data[plot.data$representatives,]) +
        geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=representatives)) +
        ylab("") + scale_y_continuous("Marginal test p-value", breaks=-log(0.1^(1:20)),
                                      labels=0.1^(1:20)) +
        xlab("Genome") +
        scale_x_continuous(limits=c(min(granice$x[chromosomeNumber]),
                                    max(granice_max[chromosomeNumber])),
                           breaks=rowMeans(cbind(granice$x, granice_max)),
                           labels=granice$Group.1,
                           minor_breaks=c(granice$x, max(granice_max))) +
        scale_alpha_manual(guide=FALSE, values = c(0.5, 1)) +
        scale_color_manual("", values = "red", labels="Clump representative") +
        scale_size_area(guide=FALSE, max_size = 4) +
        clumping_theme
    } else if(!is.null(clumpNumber)){
      plot.data <- subset(plot.data, clump%in%clumpNumber)
      if(nrow(plot.data)==0 | nrow(plot.data[plot.data$representatives,])==0) {
        message("No SNPs selected in clump ", clumpNumber)
        return(NULL)
      }
      ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = "red", size = 1),
                                     data=plot.data[plot.data$representatives,]) +
        geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=representatives)) +
        ylab("") + scale_y_continuous("Marginal test p-value", breaks=-log(0.1^(1:20)),
                                      labels=0.1^(1:20)) +
        xlab("Genome") +
        scale_x_continuous(limits=c(min(granice$x[plot.data$chromosome]),
                                    max(granice_max[plot.data$chromosome])),
                           breaks=rowMeans(cbind(granice$x, granice_max)),
                           labels=granice$Group.1,
                           minor_breaks=c(granice$x, max(granice_max))) +
        scale_alpha_manual(guide=FALSE, values = c(0.5, 1)) +
        scale_color_manual("", values = "red", labels="Clump representative") +
        scale_size_area(guide=FALSE, max_size = 4) +
        clumping_theme
    } else{
      ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = "red", size = 1),
                                     data=plot.data[plot.data$representatives,]) +
        geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=representatives)) +
        ylab("") +
        xlab("Genome") +
        scale_x_continuous(expand = c(0,0),
                           limits=c(0, max(granice_max)+1),
                           breaks=rowMeans(cbind(granice$x, granice_max)),
                           labels=granice$Group.1,
                           minor_breaks=c(granice$x, max(granice_max))) +
        scale_y_continuous("Marginal test p-value", expand = c(0,0),
                           limits=c(0, 1.1*max(plot.data$val)),
                           breaks=-log(0.1^(1:20)),
                           labels=0.1^(1:20)) +
        scale_alpha_manual(guide=FALSE, values = c(0.5, 1)) +
        scale_color_manual("", values = "red", labels="Clump representative") +
        scale_size_area(guide=FALSE, max_size = 4) +
        clumping_theme
    }
  } else {
    clumpingResult_no_info_print()
  }

}


clumpingResult_no_info_print <- function(x, ...){
  snp <- val <- clump  <- NULL #to remove CRAN's NOTE
  plot.data <- data.frame(cbind(snp=x$selectedSnpsNumbersScreening[unlist(x$SNPclumps)],
                                val=-log(x$pVals[x$selectedSnpsNumbersScreening[unlist(x$SNPclumps)]])))
  representatives = which(x$selectedSnpsNumbersScreening %in% x$selectedSnpsNumbers)
  ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = "red", size = 6),
                                 plot.data[representatives,]) +
    geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=val/4)) +
    ylab("") + scale_y_continuous("Marginal test p-value", breaks=-log(0.1^(1:5)),
                                  labels=0.1^(1:5)) +
    xlab("SNP number") +
    scale_alpha_continuous(guide=FALSE) +
    scale_color_discrete(guide=FALSE) +
    scale_size_area(guide=FALSE) +
    theme(panel.background=element_blank(),
          panel.grid.major.y=element_line(colour = "grey80"),
          panel.grid.minor.y=element_line(colour = "grey90"),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank())
}


clumping_theme <- theme(panel.background=element_blank(),
                        panel.grid.major.y=element_line(colour = "grey80"),
                        panel.grid.minor.y=element_line(colour = "grey90"),
                        panel.grid.major.x=element_blank(),
                        panel.grid.minor.x=element_line(colour = "grey70", linetype = "dotted", size=0.5),
                        axis.ticks.x=element_blank(),
                        legend.text = element_text(size=15),
                        legend.position="top",
                        legend.key =element_rect(fill="white"))

