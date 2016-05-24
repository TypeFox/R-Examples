
## ------------------------------------------------------------------------
library (package = "TiddlyWikiR")


## ------------------------------------------------------------------------
newEmptyWiki ("myTemplate.html")
dir ()


## ----, include=FALSE-----------------------------------------------------
dir.create ("files")


## ----, results="hide"----------------------------------------------------
png ("files/myplot.png")
plot (1:10)
dev.off ()


## ------------------------------------------------------------------------
myImage <- twImage (imgf = "files/myplot.png", label = "drag the image with your mouse", ref = "LinkToTiddler", width = "10%+")


## ------------------------------------------------------------------------
tagList <- list ('@@my_tag_for_plot@@' = myImage, '@@my_tag_for_number@@' = 1000000)
writeTags (tagList, infile = "myTemplate.html", outfile = "myReport_1.html")


## ------------------------------------------------------------------------
myLink <- twLink ("follow the link to a tiddler", ref = "TiddlerName")


## ------------------------------------------------------------------------
myList <- twList (elements = c("line 1", "line 2", "line 3","line 4"), 
                  level = c(1,1,2,2), 
                  type = c("o", "o", "u", "u")) 


## ------------------------------------------------------------------------
cat (wikify (myList), sep = "\n")


## ------------------------------------------------------------------------
myData <- as.data.frame (matrix (rnorm (12), ncol = 4))
rownames (myData) <- c("one", "two", "three")
myData 


## ------------------------------------------------------------------------
myTable <- twTable (dat = myData, sortable = TRUE)


## ------------------------------------------------------------------------
ref (myTable)[,"V3"] <- c("tiddlerOne", "tiddlerTwo", "http://tiddlywiki.com/")


## ------------------------------------------------------------------------
color (myTable)[,"V2"] <- c("red", "blue", "green")


## ------------------------------------------------------------------------
x <- 1:100
y <- rnorm (100)
my.stats <- summary (glm (y ~ x))


## ------------------------------------------------------------------------
myVector <- c("This may be some ''bold text''",
              "!This may be a heading")


## ------------------------------------------------------------------------
myTiddler <- newTiddler (title = "MyNewTiddler", content = list (myVector, myLink, myList, myTable, my.stats))


## ------------------------------------------------------------------------
writeTiddlers (tid = myTiddler, infile = "myReport_1.html", outfile = "myReport_2.html")


## ------------------------------------------------------------------------
writeTags (list ('@@currently does not exist@@' = "has been created"), file = "myReport_2.html")


## ------------------------------------------------------------------------
sessionInfo ()


