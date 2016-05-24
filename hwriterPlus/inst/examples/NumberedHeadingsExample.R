### Example file for including CSS in the html file
### and for automatically numbering headings
###
### DJS, 15/10/2012

date()
options(width=80)

### Set some state variables
opSys <- Sys.info()["sysname"]
if (opSys == "Windows"){
  linux <- FALSE
} else {
  linux <- TRUE
}

require(hwriterPlus)

### Prepare report file
today <- format(Sys.Date(), "%d%B%Y")
reportName <- "NumberedHeadingsExample.html"
reportName

### Include CSS file
cssText <- paste(readLines("numberedheadings.css"), collapse = "\n")


### Open file for writing
pg <- newPage(reportName, css = cssText)

### Start writing
hwrite("<span class = 'title'> Example File for Numbered Headings</span>",
       pg, center = TRUE, br = TRUE)
today <- format(Sys.Date(), "%B %d, %Y")
today
hwrite(paste("<span class = 'subtitle'> Dr David J Scott <br>",
             today, "</span>"),
       pg, center = TRUE, br = TRUE)


hwrite("Level 1 Heading",
       pg, heading = 1, br = TRUE)

hwrite("Level 2 Heading", pg, heading = 2)

hwrite("Another Level 2 Heading", pg, heading = 2)

hwrite("Level 3 Heading", pg, heading = 3)

hwrite("Another Level 1 Heading",
       pg, heading = 1, br = TRUE)


hwrite("Some ordinary text.", pg, br = TRUE)

hwrite("", pg, br = TRUE)
hwrite("", pg, br = TRUE)
hwrite("Source File: NumberedHeadingsExample.R", pg,
       style = 'font-size:8pt', br = TRUE)

closePage(pg)
directory <- getwd()
reportName <- paste("file://", directory,"/", reportName, sep = "")
reportName

browseURL(reportName)

q(save = "no")
