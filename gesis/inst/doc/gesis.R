## ----setup, eval=FALSE---------------------------------------------------
#  if(!dir.exists("downloads")) dir.create("downloads")
#  gesis_remDr <- setup_gesis(download_dir = "downloads")

## ---- eval=FALSE---------------------------------------------------------
#  login_gesis(gesis_remDr, user = "myusername", pass = "mypassword")

## ---- eval=FALSE---------------------------------------------------------
#  download_dataset(gesis_remDr, doi = 5928, filetype = "dta", purpose = 1)

## ---- eval=FALSE---------------------------------------------------------
#  dir("downloads")
#  gesis_remDr$Close()
#  gesis_remDr$closeServer()

## ---- eval=FALSE---------------------------------------------------------
#  browse_codebook(doi = 5928)

## ------------------------------------------------------------------------
library(xml2)

# Browsing the gesis website, we find the url for the main page for these studies
url <- "https://dbk.gesis.org/dbksearch/GDesc2.asp?no=0074&ll=10&db=d&notabs=1"

page <- read_html(url)
doi_links <- xml_find_all(page, "//a[contains(text(), 'ZA')]")
doi <- substr(xml_text(doi_links), 3, 7)
str(doi)

## ---- eval = FALSE-------------------------------------------------------
#  # Setup preliminaries
#  if(!dir.exists("downloads")) dir.create("downloads")
#  gesis_remDr <- setup_gesis(download_dir = "downloads")
#  
#  # Log in
#  login_gesis(gesis_remDr, user = "myusername", pass = "mypassword")
#  
#  # Loop over DOIs to download
#  lapply(doi, download_dataset, remDr = gesis_remDr)

