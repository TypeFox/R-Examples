.onAttach <- function(...) {

  meta <- packageDescription("cna")
  authors <- sprintf("%s.", meta$Author)
  year <- sub("-.*", "", meta$Date)
  title <- meta$Title
  note <- sprintf("R Package Version %s.", meta$Version)

  msg <- paste(authors, sprintf("%s.", year), sprintf("%s.", title),
               note, "URL: http://cran.r-project.org/package=cna.")

  msg <- paste(strwrap(msg, indent = 2, exdent = 2), collapse = "\n")

  packageStartupMessage("\nPlease cite the cna package as:\n\n", msg,
   "\n\nA BibTeX entry is provided by:\n  citation(\"cna\")\n")
}
