defSconcordance <- function() cat('
\\newcommand{\\Sconcordance}[1]{%
  \\ifx\\pdfoutput\\undefined%
  \\csname newcount\\endcsname\\pdfoutput\\fi%
  \\ifcase\\pdfoutput\\special{#1}%
  \\else%
   \\begingroup%
     \\pdfcompresslevel=0%
     \\immediate\\pdfobj stream{#1}%
     \\pdfcatalog{/SweaveConcordance \\the\\pdflastobj\\space 0 R}%
   \\endgroup%
  \\fi}
')

useknitr <- function(writeMacro) {
    if (!requireNamespace("knitr"))
        stop("This function requires the knitr package.")
    knitr::opts_knit$set(concordance = TRUE)
    infile <- knitr::current_input()
    if (missing(writeMacro))
        writeMacro <- any(grepl(knitr::knit_patterns$get("header.begin"),
        		    readLines(infile, 100)))   
    if (writeMacro) 
    	defSconcordance()
    cat("\\input{", tools::file_path_sans_ext(infile), "-concordance}",
        sep = "")
}
