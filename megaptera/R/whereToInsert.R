whereToInsert <- function(phy, tax, tip){
  
  tax <- rbind(tax[which(tax$spec == tip), ],
               tax[tax$spec %in% phy$tip.label, ])
  tax <- tax[ , tax[1, ] != "-"]
  id <- apply(tax, 2, function(x) x[1] %in% x[-1])
  id <- max(which(id))
  rank <- names(tax)[id]
  cat(rank, levels(tax[1, id])[tax[1, id]])
  id <- which(tax[, id] == tax[1, id])
  clade <- tax$spec[id[-1]]
  an <- noi(phy, clade)
  cat(" (node ", an, ") ", sep = "")
  
  ## check monophyly
  oops <- setdiff(descendants(phy, an, labels = TRUE), clade)
  if ( length(oops) > 0 & length(clade) > 1 ){
    if ( length(clade) / length(oops) > .5 ){
      cat("- WARNING: ", rank, "is paraphyletic:", paste(head(oops, 3), collapse = ", "),
          "[", length(oops), "]")
    } else {
      cat("- WARNING: ", rank, "seems polyphyletic:", paste(head(oops, 3), collapse = ", "),
          "[", length(oops), "]")
    }
  } else {
    cat("- ", rank, "is monophyletic")
  }
  an
}


