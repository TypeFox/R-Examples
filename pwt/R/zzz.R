.onAttach <- function(lib, pkg) {
  packageStartupMessage(paste(" Package", sQuote("pwt"), 
    "provides Penn World Table versions 5.6, 6.x, 7.x.\n",
    "For more recent versions see package", sQuote("pwt8"), "(or subsequent packages)."))
}
