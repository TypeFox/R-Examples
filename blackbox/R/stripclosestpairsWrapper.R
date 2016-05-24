stripclosestpairsWrapper <- function(table, finallength, scales, fixedNbr=NA) {
  if (is.na(fixedNbr)) {
    stop.redef("(!) from stripclosestpairsWrapper : at least one fixed point required")
  }
  if (nrow(table)> blackbox.getOption("kriglength")) {
    greedyMAXMINwithFixed(table, finallength, scales, fixedNbr) ## more economical...
  } else {
    stripclosestpairs(table, finallength, scales, fixedNbr) ## takes memory... opaque code using somme banwidth...
  }
} ##end stripclosestpairsWrapper
