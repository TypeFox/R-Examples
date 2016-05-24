# make a pdf or eps based on a switch

pdf_eps <- function(filename, width, height, type = "pdf", ...) {
  if(type == "pdf")
    pdf(paste0(filename, ".pdf"), width = width, height = height, ...)
  if(type == "eps") {
    postscript(paste0(filename, ".eps"), horizontal = FALSE, onefile = FALSE, paper =
      "special", width = width, height = height, ...)
  }
}
