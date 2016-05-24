#' @title Decouples highly imbalanced labels
#' @description This function implements the REMEDIAL algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to decouple frequent and rare classes appearing in the same instance. For doing so, it aggregates new instances to the dataset
#' and edit the labels present in them.
#' @source F. Charte, A. J. Rivera, M. J. del Jesus, F. Herrera. "Resampling Multilabel Datasets by Decoupling Highly Imbalanced Labels". Proc. 2015 International Conference on Hybrid Artificial Intelligent Systems (HAIS 2015), pp. 489-501, Bilbao, Spain, 2015
#' @param mld \code{mldr} object with the multilabel dataset to preprocess
#' @return An mldr object containing the preprocessed multilabel dataset
#' @seealso \code{\link{concurrenceReport}}, \code{\link{labelInteractions}}
#' @examples
#'
#' library(mldr)
#'\dontrun{
#' summary(birds)
#' summary(remedial(birds))
#'}
#' @export
remedial <- function(mld) decoupleImbalancedLabels(mld, mld$measures$scumble)

decoupleImbalancedLabels <- function(mld, atkLevel) {
  mldbase <- mld[mld$dataset$.SCUMBLE <= atkLevel]
  mldhigh <- mld[mld$dataset$.SCUMBLE > atkLevel]  # Samples with coocurrence of highly imbalanced labels

  # Indexes of minority and majority labels
  minIndexes <- mld$labels[mld$labels$IRLbl > mld$measures$meanIR, "index"]
  majIndexes <- mld$labels[mld$labels$IRLbl <= mld$measures$meanIR, "index"]

  # Duplicate rows affected by coocurrence of highly imbalanced labels
  ninstances <- mldhigh$measures$num.instances
  mldhigh$dataset[(ninstances+1):(ninstances*2), ] <- mldhigh$dataset
  row.names(mldhigh$dataset) <- 1:(ninstances*2)

  # Decouple majority and minority labels
  mldhigh$dataset[1:ninstances, minIndexes] <- 0
  mldhigh$dataset[(ninstances+1):(ninstances*2), majIndexes] <- 0

  mldbase + mldhigh # Join the instances without changes with the filtered ones
}

#' @title Generates a label concurrence report
#' @description This function produces a label concurrence report, providing the average SCUMBLE, SCUMBLE by label, a list with the minority labels
#' most affected by this problem indicating which majority labels they appear with, and a concurrence plot. The report output is written in the
#' standard output by default, but it could be redirected to a PDF file.
#' @param mld \code{mldr} object to analyze
#' @param pdfOutput Boolean value indicating if the output has to be sent to a PDF file. Defaults to true, so the output is shown in the console.
#' @param file If the \code{pdfOutput} parameter is \code{true} the output will be written in the file specified by this parameter. The default file name is \code{Rocurrence.pdf}.
#' @return None
#' @seealso \code{\link{remedial}}, \code{\link{labelInteractions}}
#' @examples
#'
#' library(mldr)
#'\dontrun{
#' concurrenceReport(birds)
#'}
#' @export
concurrenceReport <- function(mld, pdfOutput = FALSE, file = "Rconcurrence.pdf") {
  textFile <- tempfile()

  if(pdfOutput) sink(textFile)

  # Display dataset scumble
  cat("Dataset ", mld$name, ": Mean SCUMBLE ", mld$measures$scumble, " with CV ", mld$measures$scumble.cv, "\n\n",
      sep = "")

  # Display scumble per label
  cat("SCUMBLE mean values by label:\n")
  values <- mld$labels[order(mld$labels$SCUMBLE, decreasing = TRUE),]
  cat(t(cbind("# ", rownames(values), ": ", signif(values$SCUMBLE, digits = 4), "\n")), sep = "")
  cat("\n")

  # Show interesting label interactions
  lblint <- labelInteractions(mld)
  printInteractions(mld, intList = lblint)

  if(pdfOutput) {
    sink()
    pdf(file)
    par(mar = rep(0, 4))
    plot.new()

    ytop <- par()$usr[4]
    cheight <- par()$cxy[2] / 2
    lines <- readLines(textFile)
    file.remove(textFile)

    text(0, y = seq(from = ytop - cheight, by = -cheight, length.out = length(lines)), lines, cex = .5, adj = 0)
    plot(mld, type = "LC", labelIndices = as.numeric(c(lblint$indexes, unique(unlist(lapply(lblint$interactions, names))))))

    dev.off()
    cat('Concurrence report saved into ', file)
  } else
    plot(mld, type = "LC", labelIndices = as.numeric(c(lblint$indexes, unique(unlist(lapply(lblint$interactions, names))))))
}

#' @title Provides data about interactions between labels
#' @description This function facilitates a list with the minority labels most affected by the problem of concurrence with majority labels,
#' provinding the indexes of the majority labels interacting with each minority and also the number of instances in which they appear together.
#' @param mld \code{mldr} object to analyze
#' @param labelProportion A value in the (0,1] range establishing the proportion of minority labels to be included as result. By default at least
#' 3 or 10\% of minority labels are included, or all of them if there are fewer than 3.
#' @return A list with two slots, \code{indexes} and \code{interactions}. The former contains the indexes of the minority labels, sorted from
#' higher to lower SCUMBLE metric. The latter will provide an element for each of the previous labels, communicating the indexes of the majority
#' labels it interacts with and the number of samples in which they appear together.
#' @seealso \code{\link{remedial}}, \code{\link{concurrenceReport}}
#' @examples
#'
#' library(mldr)
#' labelInteractions(birds)
#' @export
labelInteractions <- function(mld, labelProportion) {
  # Extract minority labels
  minority <- mld$labels[mld$labels$IRLbl > mld$measures$meanIR,]
  majority <- mld$labels[mld$labels$IRLbl <= mld$measures$meanIR,]

  if (missing(labelProportion))
    labelProportion <- if (0.1 * nrow(minority) > 3)
        0.1
      else if (nrow(minority) > 3)
        3 / nrow(minority)
      else
        1

  maxScumble <- minority[order(minority$SCUMBLE, decreasing = TRUE)[1:ceiling(nrow(minority) * labelProportion)],]

  labeldata <- mld$dataset[,mld$labels$index]

  # For each selected label, extracts the rest of labels that interact with it
  interactionsByInstance <- lapply(maxScumble$index, function(i) {
    curInstances <- mld$dataset[mld$dataset[i] == 1,]
    indexes <- mld$labels$index[1] + unlist(apply(curInstances[mld$labels$index] == 1, 1, which)) - 1
    table(indexes[indexes %in% majority$index])
  })

  names(interactionsByInstance) <- rownames(maxScumble)

  list(indexes = maxScumble$index, interactions = interactionsByInstance)
}

printInteractions <- function(mld, labelProportion = 0.1, intList, interactionLimit = 10) {
  if (missing(intList))
    intList <- labelInteractions(mld, labelProportion)

  for (i in 1:length(intList$indexes)) {
    cat("Minority label ", nameOfLabel(mld, intList$indexes[i]), " (", intList$indexes[i], ", SCUMBLE ", mld$labels[mld$labels$index == intList$indexes[i], ]$SCUMBLE, ") interacts with:\n", sep="")

    intTable <- sort(intList$interactions[[i]], decreasing = TRUE)

    if (interactionLimit > 0)
      intTable <- intTable[1:min(interactionLimit, length(intTable))]

    for (l in 1:length(intTable)) {
      cat("# ", nameOfLabel(mld, names(intTable)[l]),
          " (", names(intTable)[l], ", SCUMBLE ", mld$labels[mld$labels$index == names(intTable)[l], ]$SCUMBLE, "): ",
          intTable[l], " interaction", ifelse(intTable[l] == 1, "", "s"), "\n", sep="")
    }

    cat("\n")
  }
}

nameOfLabel <- function(mld, index) {
  rownames(mld$labels[mld$labels$index == index,])
}
