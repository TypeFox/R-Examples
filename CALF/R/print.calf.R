#'@method print calf
#'@export
print.calf <- function(x, ...){
  if (x$randomize == TRUE) cat("Randomized Output:", "\n", "\n")
  if (!is.null(x$proportion)) cat("Proportion of Data:", x$proportion, "\n", "\n")
  print.data.frame(x$selection, row.names = FALSE)
  cat("\nAUC:", x$auc)
}

#'@method print calf_randomize
#'@export
print.calf_randomize <- function(x, ...){
  if (x$times == 1) {
    cat("Randomized Output Across", x$times, "Replication:", "\n", "\n")
    print.data.frame(x$selection, row.names = FALSE)
    cat("\nAUC:", x$auc)
  } else {
    cat("Randomized Output Across", x$times, "Replications:", "\n", "\n")
    print.data.frame(x$multiple, row.names = FALSE)
    cat("\n", "\n")
    print.data.frame(x$auc, row.names = F)
  }
}

#'@method print calf_subset
#'@export
print.calf_subset <- function(x, ...){
  if (x$times == 1) {
    cat("Proportion =", x$proportion, "Output Across", x$times, "Replication:", "\n", "\n")
    print.data.frame(x$selection, row.names = FALSE)
    cat("\nAUC:", x$auc)
  } else {
    cat("Proportion =", x$proportion, "Output Across", x$times, "Replications:", "\n", "\n")
    print.data.frame(x$multiple, row.names = FALSE)
    cat("\n", "\n")
    print.data.frame(x$auc, row.names = F)
  }
}

