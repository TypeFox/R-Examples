#' @title readExcelSheet
#' @description Function reads data from an Excel file from a specified sheet
#' @author Lech Madeyski
#' @export readExcelSheet
#' @param path Path to an Excel file, e.g. /User/lma/datasets/MyDataSet.xls
#' @param sheet Name of a sheet within an Excel file we want to read
#' @param colNames If TRUE, first row of data will be used as column names.
#' @examples
#' myPath=system.file("extdata", "DataSet.xlsx", package = "reproducer")
#' Madeyski15SQJ.NDC<-readExcelSheet(path=myPath, sheet="Madeyski15SQJ.NDC", colNames=TRUE)
readExcelSheet <- function(path, sheet, colNames){
  dataset = openxlsx::read.xlsx(xlsxFile=path, sheet=sheet, colNames=colNames)
  return(dataset)
}



#' @title densityCurveOnHistogram
#' @description Density curve overlaid on histogram
#' @author Lech Madeyski
#' @export densityCurveOnHistogram
#' @param df Data frame with data to be displayed
#' @param colName Name of the selected column in a given data frame
#' @param limLow the limit on the lower side of the displayed range
#' @param limHigh the limit on the higher side of the displayed range
#' @return A figure being a density curve overlaid on histogram
#' @examples
#' densityCurveOnHistogram(Madeyski15EISEJ.PropProjects, "STUD", 0, 100)
#' densityCurveOnHistogram(data.frame(x<-rnorm(50, mean=50, sd=5)), "x", 0, 100)
densityCurveOnHistogram <- function(df, colName, limLow, limHigh) {
  p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x=colName), environment = environment()) +
    ggplot2::xlab("") +
    ggplot2::ggtitle(colName) +
    ggplot2::geom_histogram(ggplot2::aes_string(y="..density.."), fill="cornsilk", colour="grey60") +
    ggplot2::geom_density(fill = 'green', alpha = 0.4) +
    ggplot2::xlim(limLow, limHigh)
  return(p1)
}

#' @title boxplotHV
#' @description Box plot
#' @author Lech Madeyski
#' @export boxplotHV
#' @param df Data frame with data to be displayed
#' @param colName Name of the selected column in a given data frame
#' @param limLow the limit on the lower side of the displayed range
#' @param limHigh the limit on the higher side of the displayed range
#' @param isHorizontal Boolean value to control whether the box plot should be horizontal or not (i.e., vertical)
#' @return A box plot
#' @examples
#' boxplotHV(Madeyski15EISEJ.PropProjects, "STUD", 0, 100, TRUE)
#' boxplotHV(Madeyski15EISEJ.PropProjects, "STUD", 0, 100, FALSE)
#' boxplotHV(Madeyski15SQJ.NDC, "simple", 0, 100, FALSE)
#' boxplotHV(Madeyski15SQJ.NDC, "simple", 0, 100, TRUE)
boxplotHV <- function(df, colName, limLow, limHigh, isHorizontal) {
  p2 <- ggplot2::ggplot(df, ggplot2::aes_string(x=1, y = colName)) +
    ggplot2::ylab("") +
    ggplot2::ggtitle(colName) +
    ggplot2::geom_boxplot(fill = 'orange') + ggplot2::theme_bw() + ggplot2::ylim(limLow, limHigh)  +
    #ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, size=4, fill="white") +
    ggplot2::stat_summary(fun.y="mean", geom="point", shape=23, fill="white") +
    ggplot2::scale_x_continuous(breaks=NULL)

  if(isHorizontal) {
    p2 <- p2 + ggplot2::coord_flip()
  }
  return(p2)
}

#' @title boxplotAndDensityCurveOnHistogram
#' @description Boxplot and density curve overlaid on histogram
#' @author Lech Madeyski
#' @export boxplotAndDensityCurveOnHistogram
#' @param df Data frame with data to be displayed
#' @param colName Name of the selected column in a given data frame
#' @param limLow the limit on the lower side of the displayed range
#' @param limHigh the limit on the higher side of the displayed range
#' @return A figure being a density curve overlaid on histogram
#' @examples
#' library(ggplot2)
#' library(grid)
#' library(gridExtra)
#' boxplotAndDensityCurveOnHistogram(Madeyski15EISEJ.PropProjects, "STUD", 0, 100)
#' boxplotAndDensityCurveOnHistogram(Madeyski15SQJ.NDC, "simple", 0, 100)
boxplotAndDensityCurveOnHistogram <- function(df, colName, limLow, limHigh) {

  if (!"package:ggplot2" %in% search()) {
    suppressPackageStartupMessages(attachNamespace("ggplot2"))
    on.exit(detach("package:ggplot2"))
  }

  p1 <- densityCurveOnHistogram(df, colName, limLow, limHigh)
  p2 <- boxplotHV(df, colName, limLow, limHigh, TRUE)

  #arrange the plots together, with appropriate height and width for each row and column
  p1 <- p1 + ggplot2::labs(
    axis.title.x = ggplot2::element_blank(),
    text = ggplot2::element_text(),
    x = ggplot2::element_blank())

  p1 <- p1 + ggplot2::theme(axis.text.y = ggplot2::element_text(size=12), axis.title=ggplot2::element_text(size=12), axis.text=ggplot2::element_text(size=12)) + ggplot2::scale_y_continuous(labels = fmt())

  p2 <- p2 + ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::theme(axis.title=ggplot2::element_text(size=12), axis.text=ggplot2::element_text(size=12))
  p2 <- p2 + ggplot2::theme(title=ggplot2::element_blank()) + ggplot2::ylim(0,100) + ggplot2::scale_x_continuous(limits=c(0.5,1.5), breaks=c(1), labels=c("Box plot"))

  gridExtra::grid.arrange(p1, p2, nrow=2, heights=c(4, 1))
}


#' @title printXTable
#' @description print data table using xtable R package
#' @author Lech Madeyski
#' @export printXTable
#' @param data Data structure including columns to be printed.
#' @param selectedColumns Columns selected to be printed.
#' @param tableType Type of table to produce. Possible values are "latex" or "html". Default value is "latex".
#' @param alignCells Defines how to align data cells.
#' @param digits Defines the number of decimal points in each column.
#' @param caption Caption of the table.
#' @param label Label of the table.
#' @param fontSize Size of the font used to produce a table.
#' @param captionPlacement The caption will be have placed at the bottom of the table if captionPlacement is "bottom" and at the top of the table if it equals "top". Default value is "bottom".
#' @param alignHeader Defines how to align column headers of a table.
#' @return A table generated on the fly on a basis of passed data (data, selectedColumns etc.).
#' @examples
#' d <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
#' reproducer::printXTable(d, '"Study"', "latex", '"cc"', 0, "C", "L", "tiny", "top", "l")
printXTable <- function(data, selectedColumns, tableType="latex", alignCells, digits, caption, label, fontSize, captionPlacement="bottom", alignHeader)
{
  df <- as.data.frame(unclass(data), optional=TRUE)
  sourcedata.xtable <- xtable::xtable(df[eval(parse(text = selectedColumns))], digits = digits, caption=caption, label=label)
  if (exists("alignHeader")){
    names(sourcedata.xtable)=alignHeader
  }
  xtable::align(sourcedata.xtable) <- eval(parse(text = alignCells))
  print(sourcedata.xtable, booktabs = TRUE, NA.string="", include.rownames=FALSE, size=fontSize, caption.placement=captionPlacement, type = tableType, sanitize.text.function = function(x){x})
} #sanitize.text.function = identity



#' @title cloudOfWords
#' @description cloud of words
#' @author Lech Madeyski
#' @export cloudOfWords
#' @param textFile A text file used to produce a word cloud
#' @return A figure being a word cloud
#' @examples
#' myPath=system.file("NAMESPACE", package = "reproducer")
#' cloudOfWords(textFile=myPath)
cloudOfWords <- function(textFile) {
  # Text mining, see http://davetang.org/muse/2013/04/06/using-the-r_twitter-package/

  # Read the text file
  myText <- readLines(textFile) #myText <- Corpus(DirSource("./tmp"))

#  if(require("tm") && require("wordcloud")){
    # Load the data as a corpus
    docs = tm::VCorpus(tm::VectorSource(myText))

    # Clean up
    docs <- tm::tm_map(docs,
                       tm::content_transformer(function(x) iconv(x, to='UTF-8', sub='byte')),
                       mc.cores=1)

      # Text transformation
    toSpace <- tm::content_transformer(function (x , pattern ) gsub(pattern, " ", x))
    docs <- tm::tm_map(docs, toSpace, "/", mc.cores=1)
    docs <- tm::tm_map(docs, toSpace, "@", mc.cores=1)
    docs <- tm::tm_map(docs, toSpace, "\\|", mc.cores=1)
    docs <- tm::tm_map(docs, toSpace, "<", mc.cores=1)
    docs <- tm::tm_map(docs, toSpace, ">", mc.cores=1)
    # Cleaning the text
    # Convert the text to lower case
    docs <- tm::tm_map(docs, tm::content_transformer(tolower), mc.cores=1)

    # Remove numbers
    docs <- tm::tm_map(docs, tm::removeNumbers, mc.cores=1)

    # Remove english common stopwords
    #docs <- tm::tm_map(docs, removeWords, stopwords("english"), lazy=TRUE, mc.cores=1)

    # Remove your own stop word
    # specify your stopwords as a character vector
    #docs <- tm::tm_map(docs, removeWords, c("blabla1", "blabla2"), lazy=TRUE, mc.cores=1)

    # Remove punctuations
    docs <- tm::tm_map(docs, tm::removePunctuation, mc.cores=1)

    # Eliminate extra white spaces
    docs <- tm::tm_map(docs, tm::stripWhitespace, mc.cores=1)

    # Text stemming
    # docs <- tm_map(docs, stemDocument, mc.cores=1)


    # Build a term-document matrix
    dtm <- tm::TermDocumentMatrix(docs)
    m <- as.matrix(dtm)
    v <- sort(rowSums(m),decreasing=TRUE)
    d <- data.frame(word = names(v),freq=v)
    # head(d, 10)


    # Generate the Word cloud
    set.seed(1234)
    wordcloud::wordcloud(words = d$word, freq = d$freq, min.freq = 1,
              max.words=200, random.order=FALSE, rot.per=0.35,
              colors=RColorBrewer::brewer.pal(8, "Dark2"))

  #}
}


#' @title fmt
#' @description Formatting function to set decimal precision in labels
#' @author Lech Madeyski
#' @export fmt
fmt <- function(){
  function(x) format(x,nsmall = 3,scientific = FALSE)
}
