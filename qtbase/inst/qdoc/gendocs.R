## Compile Qt documentation that is particular to the R API.

## This script should be executed within this directory.
## Make sure you have QtHelp support in qtbase.

library(qtbase)

source("typeMap.R")

## If you have ever run 'assistant', this should exist
dataDir <- Qt$QDesktopServices$storageLocation(Qt$QDesktopServices$DataLocation)
version <- Qt$QGlobalSpace$qVersion()
qhcFile <- file.path(dataDir, "Trolltech/Assistant",
                     paste("qthelpcollection_", version, ".qhc", sep = ""))

## Steps:
## - Extract files

engine <- Qt$QHelpEngine(qhcFile)
ns <- grep("^com.trolltech.qt.", engine$registeredDocumentations(),
               value=TRUE)
files <- engine$files(ns, character(), "html")

## delete the irrelevant class documentation
irrelevant <-
  tolower(unique(gsub("<.*", "", grep("^Q", names(typeMapBase), value=TRUE))))
irrelevantFiles <- paste(irrelevant, "html", sep = ".")

names(files) <- basename(unlist(lapply(files, qinvoke, "toString")))
files <- files[setdiff(names(files), irrelevantFiles)]

## this takes a while
rawHtml <- lapply(files, engine$fileData)
html <- lapply(rawHtml, rawToChar)

## - Transform them

## map types
for (type in names(typeMap))
  html <- gsub(type, typeMap[type], html)

## pointers, references => values
html <- gsub(" \\*+", "", html)
html <- gsub("&amp;", "", html, fixed=TRUE)

## delete const (well, maybe leave it in there.. it's informative)
##html <- gsub("const ", "", html, fixed=TRUE)

## false -> FALSE, true -> TRUE
html <- gsub("false", "FALSE", html)
html <- gsub("true", "TRUE", html)

## drop property accessor information
html <- gsub("<p><b>Access functions:</b></p>.*?</table>", "", html)

## drop links to irrelevant classes

linkexprs <- paste('<a href="', irrelevant, '.html">(.*?)</a>', sep = "")
for (linkexpr in linkexprs)
  html <- gsub(linkexpr, "\\1", html)

## - Write it back to disk

mapply(writeLines, html, names(files))

## read it back in
## engine <- Qt$QHelpEngineCore("qtbase.qhc")
## files <- engine$files("org.r-project.qtbase.1", character(), "html")
## fileData <- rawToChar(engine$fileData(files[[1]]))

## Parse the qt.qhp file and patch it. This file needs to be copied
## from the Qt source tree, after making it. It is not present in the
## installed version (only the qhc/qch files are necessary). Also
## necessary are the 'images' directory and the 'classic.css' file.

parse <- xmlInternalTreeParse("qt.qhp")
root <- xmlRoot(parse)
project <- xmlChildren(root)

## Fix a bug in XML package
setMethod("xmlAttrs<-", "XMLNode",
          function(node, append = TRUE,
                   suppressNamespaceWarning =
                   getOption("suppressXMLNamespaceWarning", FALSE), value)
          {
            node <- addAttributes(node, .attrs = value,
                                  suppressNamespaceWarning =
                                  suppressNamespaceWarning, append = append)
            node
          })

xmlValue(project$namespace) <- "org.r-project.qtbase.1"
xmlValue(project$virtualFolder) <- "qtbase"
xmlAttrs(project$customFilter)["name"] <- "qtbase 1.x"
xmlValue(xmlChildren(project$customFilter)$filterAttribute) <- "qtbase"
xmlValue(xmlChildren(project$filterSection)$filterAttribute) <- "qtbase"
xmlAttrs(xmlChildren(xmlChildren(project$filterSection)$toc)[[1]])["title"] <-
  "Qt (qtbase) Reference"

xmlChildren(root) <- project
saveXML(root, "qtbase.qhp")

## - Compile help files
system("qcollectiongenerator qtbase.qhcp -o qtbase.qhc")

## - Delete the generated HTML
unlink(outfiles)

####### NOT WORKING YET STUFF #######

## Attempt to get information from QtHelp, instead of using qt.qhp template

## No efficient way to obtain this. Could get the model index, obtain
## every keyword, and look up the links. Then use the XML package to
## insert the keywords.

index <- engine$indexModel()
index$createIndex("qt")
keywords <- index$stringList()

## this takes an ungodly amount of time
links <- unlist(lapply(keywords, function(keyword) {
  links <- sapply(engine$linksForIdentifier(keyword), qinvoke, "toString")
}))

anchor <- grepl("#", links, fixed=TRUE)
names(links)[anchor] <- paste(names(links)[anchor],
                              sub(".*#", "", links[anchor]), sep = ": ")
links <- links[unique(names(links))]
names(links) <- sub("^Qt 4.6: ", "", names(links))

attrs <- c(sub("\\..*", "", basename(links)),
           XML:::insertEntities(names(links)), basename(links))
names(attrs) <- rep(c("id", "name", "ref"), each = length(links))
attrs <- split(attrs, seq(length(links)))
keywordNodes <- mapply(xmlNode, "keyword", attrs = attrs)

root <- xmlRoot(xmlTreeParse("qtbase.qhp"))
xmlChildren(xmlChildren(root)$filterSection)$keywords <-
  xmlNode("keywords", .children = keywordNodes)

saveXML(root, "qtbase.qhp")

## Similar thing for TOC? Not as useful as index.
##contents <- engine$contentModel()

