new.page <- function(workflow.name, i.case = 1, model = "roc", case.main = NULL,
    default.layout = TRUE){
  if(is.null(case.main)){
    case.main <- get.case.main(i.case, model)
  }
  if(default.layout){
    nf <- layout(matrix(c(1, 1, 2, 3, 4, 5, 6, 7),
                        nrow = 4, ncol = 2, byrow = TRUE),
                 c(1, 1), c(2, 8, 8, 8), respect = FALSE)
  }

  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, paste(workflow.name, ", ", case.main, ", NPS", sep = ""))
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))

  invisible()
} # End of new.page().
