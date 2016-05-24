## create errata.html file
## cd DIR; Rscript runErrata.R
library(markdown)
markdownToHTML("errata.md", "errata.html")
