# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-11-15 06:59:39 gjw>
#
# 080921 TEXT MINING DATA
#
# Copyright (c) 2009 Togaware Pty Ltd
#
# This file is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.
#
########################################################################
#
# First some notes:
#
#

## > show(corpus)
## A text document collection with 5 text documents

## > summary(corpus)
## A text document collection with 5 text documents

## > inspect(corpus[1])

## tdm <-  TermDocMatrix(corpus)
## findFreqTerms(tdm, 5, Inf)
## findAssocs(tdm, "ads", 0.97)

## ##
## ## Add in the target
## ##

## target <- c(1, 0, 0, 1, 0)
## crs$dataset <- as.data.frame(cbind(tdm@.Data, target))
## set.seed(123)
## crs$sample <- sample(nrow(crs$dataset), 4)

## ##
## ## Ignore 1 (15th), 61 (_is_), 238 (30%) or get error, probably
## ## because of their names.
## ##

## crs$rpart <- rpart(target ~ .,
##                    data=crs$dataset[crs$sample,c(2:60,62:237,239:285)],
##                    method="class")

## crs$rf <- randomForest(as.factor(target) ~ .,
##                        data=crs$dataset[crs$sample,c(2:60,62:237,239:285)],
##                        importance=TRUE, na.action=na.omit)


## crs$glm <- glm(target ~ .,
##                data=crs$dataset[crs$sample,c(2:60,62:237,239:285)],
##                family=binomial(logit))

## ##
## ## The others dont yet work:
## ##


## crs$ada <- ada(target ~ ., data=crs$dataset[crs$sample,c(2:60,62:237,239:285)])

## crs$ksvm <- ksvm(as.factor(target) ~ .,
##                  data=crs$dataset[crs$sample,c(2:60,62:237,239:285)],
##                  prob.model=TRUE)

executeDataCorpus <- function()
{
  # 080921 Load all documents in the specified corpus as a document
  # corpus except target.csv, if there is one. Load .target.csv if
  # there is one as the target for each document in the corpus. The
  # .target.csv file must have two columns, comma separated. The first
  # row should name the columns, but we don't actually use the column
  # names here. The first column is the document id and must be the
  # filename without its extension. The second column is the
  # classification, for example 0 or 1. I use the name ".target.csv"
  # so that the corpus loader will ignore it as a hidden file.

  # 130310 For now, each time we Execute, reload the dataset. Effect
  # this with the following:

  crs$dataset <- NULL
  theWidget("select_treeview")$getModel()$clear()
  
  # Obtain interface information.

  location <- theWidget("data_corpus_location_filechooserbutton")$getFilename()
  strip <- theWidget("data_corpus_strip_checkbutton")$getActive()
  lcase <- theWidget("data_corpus_lowercase_checkbutton")$getActive()
  stopw <- theWidget("data_corpus_stopwords_checkbutton")$getActive()
  stemw <- theWidget("data_corpus_stem_checkbutton")$getActive()

  # Start the log for this task.
  
  startLog("LOAD A CORPUS")

  # Ensure the package is available.

  lib.cmd <- "library(tm, quietly=TRUE)"
  if (! packageIsAvailable("tm", "text mining")) return(FALSE)
  appendLog("Use the tm package to support text mining.", lib.cmd)
  eval(parse(text=lib.cmd))

  # This seems to be avaiable somewhere? library(RStem)
  
  # Load the document corpus.

  corpus.cmd <- sprintf('my.corpus <- Corpus(DirSource("%s"))',
                        gsub("\\\\", "/", location))
  appendLog("Load the document corpus.", corpus.cmd)
  setStatusBar(Rtxt("Loading corpus from the documents found in"), location, "...")
  eval(parse(text=corpus.cmd))

  # Process the documents.

  map.cmd <- ""
  
  if (strip)
    map.cmd <- sprintf("%s\nmy.corpus <- tm_map(my.corpus, stripWhitespace)", map.cmd)
  if (lcase) 
    map.cmd <- sprintf("%s\nmy.corpus <- tm_map(my.corpus, content_transformer(tolower))", map.cmd)
  if (stopw) 
    map.cmd <- sprintf(paste("%s\nmy.corpus <- tm_map(my.corpus,",
                             'removeWords, stopwords("english"))'), map.cmd)
  if (stemw)
  {
    lib.cmd <- "library(SnowballC, quietly=TRUE)"
    if (! packageIsAvailable("SnowballC", "word stemming")) return(FALSE)
    appendLog(packageProvides("SnowballC", "stemDocument"), lib.cmd)
    eval(parse(text=lib.cmd))

    map.cmd <- sprintf("%s\nmy.corpus <- tm_map(my.corpus, stemDocument)", map.cmd)
  }
  

  # 111020 For now, always remove punctuation and numbers.
  
  map.cmd <- sprintf("%s\nmy.corpus <- tm_map(my.corpus, removePunctuation)", map.cmd)
  map.cmd <- sprintf("%s\nmy.corpus <- tm_map(my.corpus, removeNumbers)", map.cmd)

  # 111020 TODO Update and include some more information.

##   Dictionary(TermDocumentMatrix(my.corpus))

## tdm <- TermDocumentMatrix(my.corpus, 
##                           control = list(removePunctuation = TRUE, 
##                                          removeNumbers = TRUE, 
##                                          stopwords = TRUE))

## plot(tdm, corThreshold = 0.8, weighting = TRUE, 
##      attrs = list(graph = list(rankdir = "BT"), 
##                   node = list(shape = "circle"))) 
 

## dissimilarity(my.corpus[[1]], my.corpus[[2]], method = "eJaccard") 
## dissimilarity(tdm, method = "cosine")

## rownames(tdm) 
## colnames(tdm) 
## dimnames(tdm) 
## Docs(tdm) 
## nTerms(tdm) 
## Terms(tdm)

## inspect(my.corpus[1:3]) 
## tdm <- TermDocumentMatrix(my.corpus)[1:10, 1:10] 
## inspect(tdm)

## summary(my.corpus)

## findFreqTerms(tdm, 2, 3 )

## removeSparseTerms(tdm,0.4)

## searchFullText(my.corpus[[3]], "accounts")

## termFreq(my.corpus[[1]])


  
  appendLog("Transform the documents.", sub("^\n", "", map.cmd))
  setStatusBar(Rtxt("Transforming the documents"), "...")
  eval(parse(text=map.cmd))

  # Convert into a keyword count dataset.

  ds.cmd <- "crs$dataset <- as.data.frame(t(as.matrix(TermDocumentMatrix(my.corpus))))"
  appendLog("Convert into a dataset.", ds.cmd)
  eval(parse(text=ds.cmd))

  # Add in targets if they exist.

  target.fname <- paste(location, ".target.csv", sep="/")
  if (file.exists(target.fname))
  {
    read.cmd <- sprintf('target <- read.csv("%s", encoding="%s")',
                        target.fname, crv$csv.encoding)
    appendLog("Read in the targets.", read.cmd)
    eval(parse(text=read.cmd))

    if (nrow(crs$dataset) != nrow(target))
    {
      errorDialog(Rtxt("The number of targets is different to the",
                       "number of documents:"),
                  sprintf("%s %s %s.", nrow(target), Rtxt("versus"), nrow(crs$dataset)),
                  Rtxt("You may need to update the file"),
                  target.fname,
                  Rtxt("to match the number of documents in the corpus."))
      return(FALSE)
    }
    
    target.cmd <- "crs$dataset <- cbind(crs$dataset, TARGET=target[[2]])"
    appendLog("Add the targets to the dataset.", target.cmd)
    eval(parse(text=target.cmd))
  }

  # Set the title and dataname correctly.

  crs$dataname <- basename(location)
  setMainTitle(crs$dataname)

  # For now, always succeed.
  
  setStatusBar(Rtxt("Corpus has been loaded from the documents in"),
               location,
               ifelse(file.exists(target.fname),
                      paste(Rtxt("with targets from"), ".target.csv"),
                      ""))

  return(TRUE)
}

  
