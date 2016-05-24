auto_hinges <- function(w.length,
                        step.limit = 0.5) {
  # this uses average step.size but is "cheaper" to compute
  ((w.length[length(w.length)] - w.length[1]) / length(w.length)) > (step.limit * 0.75)
  # this was earlier used and searches for the largest step
#  stepsize(w.length)[2] > step.limit
}
