#' hsmm_pred class
#'
#' A single prediction of \code{signalHsmm}.
#'
#' @details Always a named list of five elements
#' \enumerate{
#' \item \code{sp_probability} is a probability of signal peptide presence.
#' \item \code{sp_start} is a start of potential signal peptide (naively 1 aminoacid).
#' \item \code{sp_end} is a position of last amino acid of signal peptide.
#' \item \code{struc} is numeric vector representing predicted structure of input 
#' protein.
#' \item \code{prot} is character vector containing input sequence of amino acids.
#' \item \code{str_approx} has value bigger than 0 if the predicted signal peptide 
#' structure was approximated (usually in case of sequences that have no signal peptides).
#' }
#' @seealso \code{\link{summary.hsmm_pred}} \code{\link{plot.hsmm_pred}}
#' @name hsmm_pred
NULL

#' Plot single signalHsmm prediction
#'
#' Plots objects of class \code{\link{hsmm_pred}}.
#'
#' @param x object of class \code{\link{hsmm_pred}}.
#' @param add_legend \code{logical}, if \code{TRUE}, legend is added to the plot.
#' @param only_sure \code{logical}, if \code{FALSE} does not draw signal peptide structure
#' when probability is smaller than 0.5.
#' @param ... ignored.
#' @return Nothing.
#' @export
#' @keywords hplot methods

plot.hsmm_pred <- function(x, add_legend = TRUE, only_sure = TRUE, ...) {
  plot(c(1, 50), c(1, 5), type="n", axes=F, ylab = "", xlab = "Amino acid index",
       main = x[["name"]])
  axis(1, 1L:50, labels = FALSE)
  axis(1, 1L:25*2 - 1, labels = 1L:25*2 - 1)
  
  #get borders of regions, add 0.5 to have countinous regions for purpose of easy plotting
  cstruc <- cumsum(rle(x[["struc"]])[["lengths"]])
  cstruc05 <- c(1, cstruc + 0.5)
  cstruc <- c(0, cstruc)
  sig_colours <- c("#73DC4AFF", "#30E9E9FF", "#FCC753FF", "red")
  
  #old boring black text
  #text(1L:50, 1, x[["prot"]])
  # 
  if(only_sure & x[["sp_probability"]] < 0.5) {
    text(1L:50, 1, 
         x[["prot"]][1L:50], 
         col = sig_colours[4])
    lines(c(cstruc05[1], cstruc05[5]), c(1.5, 1.5), 
          col = sig_colours[4], lwd = 8)
    
  } else {
    for(i in 1L:4) {
      #print amino acids in color!
      text((cstruc[i] + 1):cstruc[i + 1], 1, 
           x[["prot"]][(cstruc[i] + 1):cstruc[i + 1]], 
           col = sig_colours[i])
      lines(c(cstruc05[i], cstruc05[i + 1]), c(1.5, 1.5) + ifelse(i == 4, 0, 1), 
            col = sig_colours[i], lwd = 8)
      lines(c(cstruc05[4], cstruc05[4]), c(1.5, 2.5), lty = "dashed", lwd = 2)
    }
  }
  
  
  if (add_legend)
    legend("topright", 
           col = rev(c(sig_colours, "black", "white", "white")),
           lwd = rev(c(5, 5, 5, 5, 2, 1, 1)), 
           lty = rev(c(rep("solid", 4), "dashed", "blank", "blank")),
           legend = rev(c("n-region", 
                          "h-region", 
                          "c-region", 
                          "mature protein", 
                          "cleavage site", 
                          paste0("Signal peptide probability: ", 
                                 formatC(x[["sp_probability"]], digits = 2, format = "f")),
                          ifelse(x[["str_approx"]] > 0, 
                                 "Signal peptide structure interpolated",
                                 " "))), 
           bty = "n")
}

#' Summarize single signalHsmm prediction
#'
#' Summarizes objects of class \code{\link{hsmm_pred}}.
#'
#' @param object of class \code{\link{hsmm_pred}}.
#' @param only_sure \code{logical}, if \code{FALSE} does not draw signal peptide structure
#' when probability is smaller than 0.5.
#' @param ... ignored
#' @return Nothing.
#' @export
#' @keywords print methods

summary.hsmm_pred <- function(object, only_sure = TRUE, ...) {
  struc <- rle(object[["struc"]])[["lengths"]]
  cstruc <- cumsum(struc)
  cat(object[["name"]],
      paste0("Probability of signal peptide presence: ", 
             format(object[["sp_probability"]], digits = 4)),
      paste0("Signal peptide", ifelse(object[["sp_probability"]] < 0.5, " not ", " "), 
             "detected."),
      sep = "\n"
  )
  
  if(only_sure & object[["sp_probability"]] > 0.5) {
    cat("\nThe prediction of the regional structure is an experimental feature, use at your own risk.",
        paste0("Start of signal peptide: ", object[["sp_start"]]),
        paste0("End of signal peptide: ", object[["sp_end"]]),
        paste0("n-region (length ", struc[1], "):"),
        paste0(c("         ", object[["prot"]][1L:cstruc[1]]), collapse = ""),
        paste0("\nh-region (length ", struc[2], "):"),
        paste0(c("         ", object[["prot"]][(cstruc[1] + 1):cstruc[2]]), collapse = ""),
        paste0("\nc-region (length ", struc[3], "):"),
        paste0(c("         ", object[["prot"]][(cstruc[2] + 1):cstruc[3]], "\n"), collapse = ""),
        sep = "\n"
    )
    if(object[["str_approx"]] > 0)
      cat("Signal peptide structure interpolated.\n")
  }
}