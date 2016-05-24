## ----opts, echo=FALSE----------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 7,
  fig.height = 5
)

## ----setup, echo=TRUE, eval=TRUE-----------------------------------------
library(rodham)

# get list of emails
data("emails")

# equivalent to:
em <- search_emails()

identical(emails, em)

## ----edges, echo=TRUE, eval=TRUE-----------------------------------------
edges <- edges_emails(emails)
knitr::kable(head(edges))

## ----simple network, echo=TRUE, eval=TRUE--------------------------------
g <- igraph::graph.data.frame(edges)
# plot network
plot(g, layout = igraph::layout.fruchterman.reingold(g),
     vertex.label.color = hsv(h = 0, s = 0, v = 0, alpha = 0.0), 
     vertex.size = log1p(igraph::degree(g)) * 2, edge.arrow.size = 0.1, 
     edge.arrow.width = 0.1, edge.width = log1p(igraph::E(g)$freq)/4,
     vertex.frame.color="#FFFFFF")

## ----get xpdf, echo=TRUE, eval=FALSE-------------------------------------
#  xpdf <- get_xpdf(dest = "C:/") # get extractor
#  # or if you downloaded manually point to pdftotext
#  xpdf <- "your/path/xpdfbin-win-3.04/bin64/pdftotext.exe"

## ----get emails, echo=TRUE, eval=FALSE-----------------------------------
#  dir <- "./rodham" # this is where the emails will be saved
#  dir.create(dir) # directory must exist
#  emails_bengh <- get_emails(release = "Benghazi", save.dir = dir, extractor = xpdf)

## ----read emails, echo=TRUE, eval=FALSE----------------------------------
#  files <- list.files(emails_bengh) # list extracted files
#  content <- lapply(1:length(files), function(x){
#    readLines(paste0(emails_bengh, "/", files[[x]])) # read files
#  })

