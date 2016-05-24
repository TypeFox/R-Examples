# http://nlp.stanford.edu/software/stanford-corenlp-full-2015-01-29.zip

volatiles = new.env(parent=emptyenv())

#' Initialize the CoreNLP java object
#'
#' This must be run prior to calling any other CoreNLP
#' functions. It may be called multiple times in order
#' to specify a different parameter set, but note that
#' if you use a different configuration during the same
#' R session it must have a unique name.
#'
#' @importFrom           rJava .jinit .jaddClassPath .jcall .jnew
#' @param libLoc         a string giving the location of the CoreNLP java
#'                       files. This should point to a directory which
#'                       contains, for
#'                       example the file "stanford-corenlp-*.jar", where "*" is the
#'                       version number. If missing, the function will try to find the
#'                       library in the environment variable CORENLP_HOME, and otherwise
#'                       will fail.
#' @param parameterFile  the path to a parameter file. See the CoreNLP documentation for
#'                       an extensive list of options. If missing, the package will simply
#'                       specify a list of standard annotators and otherwise only use default
#'                       values.
#' @param mem            a string giving the amount of memory to be assigned to the rJava
#'                       engine. For example, "6g" assigned 6 gigabytes of memory. At least
#'                       2 gigabytes are recommended at a minimum for running the CoreNLP
#'                       package. On a 32bit machine, where this is not possible, setting
#'                       "1800m" may also work. This option will only have an effect the first
#'                       time \code{initCoreNLP} is called, and also will not have an effect if
#'                       the java engine is already started by a seperate process.
#' @param annotators     optional character string. When parameterFile is missing, this
#'                       is taken to be the list of annotators that you want to run. Either
#'                       a length one character vector with annotator names seperated by
#'                       commas, or a character vector with one annotator per element.
#'@examples
#'\dontrun{
#'initCoreNLP()
#'sIn <- "Mother died today. Or, maybe, yesterday; I can't be sure."
#'annoObj <- annotateString(sIn)
#'}
#' @export
initCoreNLP = function(libLoc, parameterFile, mem="4g", annotators) {
  # Find location of the CoreNLP Libraries
  if (missing(libLoc)) {
    libLoc = paste0(system.file("extdata",package="coreNLP"),
                    "/stanford-corenlp-full-2015-04-20")
    if (!file.exists(libLoc))
      stop("Please run downloadCoreNLP() in order to install required jar files.")
  }
  if (!file.exists(libLoc) || !file.info(libLoc)$isdir)
    stop("libLoc does not point to an existing directory path")
  path = Sys.glob(paste0(libLoc,"/*.jar"))

  # Start java engine, if not already done, and add to classpath
  options(java.parameters = paste0("-Xmx", mem))
  rJava::.jinit()
  rJava::.jaddClassPath(path)

  # Determine if the corenlp files have been loaded correctly
  len = length(grep("stanford-corenlp-", basename(rJava::.jclassPath())))
  if (len == 0L)
    stop("The coreNLP jar files are were not found in libLoc.")
  if (len < 4L)
    warning("The set of coreNLP jar files may be incomplete. Proceed with caution")

  # Read parameter file and add to classpath
  if (missing(parameterFile)) {
    path = paste0(system.file("extdata",package="coreNLP"),"/config.properties")
  } else path = Sys.glob(paste0(parameterFile[[1]]))
  rJava::.jaddClassPath(dirname(path))

  if (!is.null(volatiles$cNLP))
    rJava::.jcall(volatiles$cNLP, "V", "clearAnnotatorPool")

  if (!missing(annotators) & missing(parameterFile)) {
    annotators = paste(annotators,collapse=",")
    prop = rJava::.jnew("java.util.Properties")
    rJava::.jcall(prop, "Ljava/lang/Object;", "setProperty", "annotators", annotators)
    volatiles$cNLP = rJava::.jnew("edu.stanford.nlp.pipeline.StanfordCoreNLP", prop)
  } else {
    volatiles$cNLP = rJava::.jnew("edu.stanford.nlp.pipeline.StanfordCoreNLP", basename(path))
  }

  volatiles$xmlOut = rJava::.jnew("edu.stanford.nlp.pipeline.XMLOutputter")
}

#' Annotate a string of text
#'
#' Runs the CoreNLP annotators over a given string of text. The details
#' for which annotators to run and how to run them are specified in the
#' properties file loaded in via the \code{initCoreNLP} function (which
#' must be run prior to any annotation).
#'
#' @importFrom        rJava .jcall
#' @param text        a vector of strings for which an annotation is desired.
#'                    Will be collapsed to length 1 using new line characters
#'                    prior to the annotation.
#' @param format      the desired output format. Option \code{obj}, the default,
#'                    returns an R object of class \code{annotation} and will
#'                    likely be the desired choice for most users. The \code{xml}
#'                    and \code{text} exist primarily for subsequently saving
#'                    to disk.
#' @param outputFile  character string indicating where to put the output. If
#'                    set to NA, the output will be returned by the function.
#' @param includeXSL  boolean. Whether the xml style sheet should be included
#'                    in the output. Only used if format is \code{xml} and
#'                    outputFile is not \code{NA}.
#'@examples
#'\dontrun{
#'initCoreNLP()
#'sIn <- "Mother died today. Or, maybe, yesterday; I can't be sure."
#'annoObj <- annotateString(sIn)
#'}
#' @export
annotateString = function(text, format=c("obj", "xml", "text"), outputFile=NA,
                          includeXSL=FALSE) {
  if (is.null(volatiles$cNLP))
    stop("Must initilize with 'initCoreNLP'!")

  format = match.arg(format)
  if (length(text) > 1L) text = paste(text,collapse="\n")
  if (!is.na(outputFile) & format == "obj") {
    warning("Cannot output 'obj' as text file, saving xml instead.")
    format = "xml"
  }

  anno = rJava::.jcall(volatiles$cNLP, "Ledu/stanford/nlp/pipeline/Annotation;", "process", text)

  # Parse the output
  if (format == "text") {
    rJava::.jcall(volatiles$cNLP, "V", "prettyPrint", anno, .jnew("java.io.PrintWriter", tf <- tempfile()))
    on.exit(file.remove(tf))
    out = readLines(tf)
  } else {
    d = rJava::.jcall(volatiles$xmlOut, "Lnu/xom/Document;", "annotationToDoc", anno, volatiles$cNLP)
    xml = rJava::.jcall(d, "Ljava/lang/String;", "toXML")
    if (format == "xml") out = xml else out = parseAnnoXML(xml)
  }

  # Return the output if asked for; otherwise just write to disk
  if (is.na(outputFile)) return(out)
  outputCon = file(outputFile, "w")
  on.exit(close(outputCon),add=TRUE)
  writeLines(out, outputCon)
  if (includeXSL & format == "xml")
    file.copy(from=paste0(system.file("extdata",package="coreNLP"), "/CoreNLP-to-HTML.xsl"),
      to=dirname(outputFile))

  invisible(outputFile)
}

#' Annotate a text file
#'
#' Runs the CoreNLP annotators for the text contained in a given file.
#' The details for which annotators to run and how to run them are
#' specified in the properties file loaded in via the \code{initCoreNLP}
#' function (which must be run prior to any annotation).
#'
#' @importFrom        rJava .jcall
#' @param file        a string giving the location of the file to be loaded.
#' @param format      the desired output format. Option \code{obj}, the default,
#'                    returns an R object of class \code{annotation} and will
#'                    likely be the desired choice for users loading the output
#'                    into R. The \code{xml} and \code{text} exist primarily for
#'                    saving the files on the disk.
#' @param outputFile  character string indicating where to put the output. If
#'                    set to NA, the output will be returned by the function.
#' @param includeXSL  boolean. Whether the xml style sheet should be included
#'                    in the output. Only used if format is \code{xml} and
#'                    outputFile is not \code{NA}.
#' @export
annotateFile = function(file, format=c("obj", "xml", "text"), outputFile=NA,
                            includeXSL=FALSE) {

  if (is.null(volatiles$cNLP))
    stop("Must initilize with 'initCoreNLP'!")

  # Processing inputs
  format = match.arg(format)
  if (!is.na(outputFile) & format == "obj") {
    warning("Cannot output 'obj' as text file, saving xml instead.")
    format = "xml"
  }
  if (is.character(file)) {
    file = file(file, "r")
    on.exit(close(file))
  }

  # Read input into R and annotate
  #   TODO: Use native java interface instead of reading into R.
  text = readLines(file)
  if (length(text) > 1L) text = paste(text,collapse="\n")
  anno = rJava::.jcall(volatiles$cNLP, "Ledu/stanford/nlp/pipeline/Annotation;",
                "process", text)

  # Parse the output
  if (format == "text") {
    rJava::.jcall(volatiles$cNLP, "V", "prettyPrint", anno, .jnew("java.io.PrintWriter", tf <- tempfile()))
    on.exit(file.remove(tf))
    out = readLines(tf)
  } else {
    d = rJava::.jcall(volatiles$xmlOut, "Lnu/xom/Document;", "annotationToDoc", anno, volatiles$cNLP)
    xml = rJava::.jcall(d, "Ljava/lang/String;", "toXML")
    if (format == "xml") out = xml else out = parseAnnoXML(xml)
  }

  # Return the output if asked for asked for; otherwise just write to disk
  if (is.na(outputFile)) return(out)
  outputCon = file(outputFile, "w")
  on.exit(close(outputCon),add=TRUE)
  writeLines(out, outputCon)
  if (includeXSL && format == "xml")
    file.copy(from=paste0(system.file("extdata",package="coreNLP"), "CoreNLP-to-HTML.xsl"),
      to=dirname(outputFile))
}

#' Print a summary of an annotation object
#'
#' @method print annotation
#' @param x    an annotation object
#' @param ...  other arguments. Currently unused.
#'
#'@examples
#'print(annoEtranger)
#'
#' @export
print.annotation = function(x, ...) {
  cat("\nA CoreNLP Annotation:\n")
  cat("  num. sentences:", x$token[nrow(x$token),1], "\n")
  cat("  num. tokens:", nrow(x[[1]]), "\n")
  cat("\n")
  invisible(x)
}

#' Get tokens as data frame
#'
#' Returns a data frame of the tokens from an annotation
#' object.
#'
#' @param annotation    an annotation object
#'
#'@examples
#'getToken(annoEtranger)
#'
#' @export
getToken = function(annotation) {
  annotation$token
}

#' Get parse tree as character vector
#'
#' Returns a character vector of the parse trees.
#' Mostly use for visualization; the output of
#' \code{\link{getToken}} will generally be more
#' conveniant for manipulating in R.
#'
#' @param annotation    an annotation object
#'
#'@examples
#'getParse(annoEtranger)
#'
#' @export
getParse = function(annotation) {
  annotation$parse
}

#' Get Dependencies
#'
#' Returns a data frame of the coreferences of an annotation
#'
#' @param annotation    an annotation object
#' @param type          the class of coreference desired
#'
#'@examples
#'getDependency(annoEtranger)
#'getDependency(annoHp)
#'
#' @export
getDependency = function(annotation, type=c("CCprocessed","basic","collapsed")) {
  type = match.arg(type)
  dep <- if (type == "basic") annotation$basicDep
  else if (type == "collapsed") annotation$collapsedDep
  else if (type == "CCprocessed") annotation$collapsedProcDep

  dep$govIndex = match(paste0(dep$sentence,"-",dep$governorIdx),
                     paste0(annotation$token$sentence,"-",annotation$token$id))
  dep$depIndex = match(paste0(dep$sentence,"-",dep$dependentIdx),
                     paste0(annotation$token$sentence,"-",annotation$token$id))
  dep
}

#' Get Sentiment scores
#'
#' Returns a data frame of the sentiment scores from an annotation
#'
#' @param annotation    an annotation object
#'
#'@examples
#'getSentiment(annoEtranger)
#'getSentiment(annoHp)
#'
#' @export
getSentiment = function(annotation) {
  annotation$sentiment
}

#' Get Coreference
#'
#' Returns a dataframe containing all coreferences detected in the text.
#'
#' @param annotation    an annotation object
#'
#'@examples
#'getCoreference(annoHp)
#'
#' @export
getCoreference = function(annotation) {
  coref = annotation$coref[,-grep("text",names(annotation$coref))]
  coref$startIndex = match(paste0(coref$sentence,"-",coref$start),
                     paste0(annotation$token$sentence,"-",annotation$token$id))
  coref$endIndex = match(paste0(coref$sentence,"-",as.numeric(coref$end)-1),
                     paste0(annotation$token$sentence,"-",annotation$token$id))
  coref
}

#' Load CoreNLP XML file
#'
#' Loads a properly formated XML file output by the CoreNLP
#' library into an \code{annotation} object in R.
#'
#' @param file      connection or character string giving the file name to load
#' @param encoding  encoding to be assumed for input strings.  It is used to mark
#'                  character strings as known to be in Latin-1 or UTF-8: it is
#'                  not used to re-encode the input. Passed to \code{readLines}.
#' @export
loadXMLAnnotation = function(file, encoding="unknown") {
  if (is.character(file)) {
    file = file(file, "r")
    on.exit(close(file))
  }
  x = readLines(file)
  tryCatch(output <- parseAnnoXML(x),
           error = function(e) stop("Not a valid CoreNLP XML file."))
  output
}

#' Convert Penn TreeBank POS to Universal Tagset
#'
#' Maps a character string of English Penn TreeBank
#' part of speech tags into the universal tagset
#' codes. This provides a reduced set of tags (12), and
#' a better cross-linguist model of speech.
#'
#' @param pennPOS   a character vector of penn tags to match
#'
#'@examples
#'tok <- getToken(annoEtranger)
#'cbind(tok$POS,universalTagset(tok$POS))
#'
#' @export
universalTagset = function(pennPOS) {
  mtab = structure(c("!", "#", "$", "''", "(", ")", ",", "-LRB-",
      "-RRB-",  ".", ":", "?", "CC", "CD", "CD|RB", "DT", "EX", "FW", "IN",
      "IN|RP",  "JJ", "JJR", "JJRJR", "JJS", "JJ|RB", "JJ|VBG", "LS", "MD", "NN",
      "NNP", "NNPS", "NNS", "NN|NNS", "NN|SYM", "NN|VBG", "NP", "PDT",  "POS", "PRP",
      "PRP$", "PRP|VBP", "PRT", "RB", "RBR", "RBS", "RB|RP",  "RB|VBG", "RN", "RP", "SYM",
      "TO", "UH", "VB", "VBD", "VBD|VBN",  "VBG", "VBG|NN", "VBN", "VBP", "VBP|TO",
      "VBZ", "VP", "WDT",  "WH", "WP", "WP$", "WRB", "``", ".", ".", ".", ".", ".",
      ".",  ".", ".", ".", ".", ".", ".", "CONJ", "NUM", "X", "DET", "DET",  "X",
      "ADP", "ADP", "ADJ", "ADJ", "ADJ", "ADJ", "ADJ", "ADJ",  "X", "VERB", "NOUN",
      "NOUN", "NOUN", "NOUN", "NOUN", "NOUN",  "NOUN", "NOUN", "DET", "PRT", "PRON",
      "PRON", "PRON", "PRT",  "ADV", "ADV", "ADV", "ADV", "ADV", "X", "PRT", "X",
      "PRT", "X",  "VERB", "VERB", "VERB", "VERB", "VERB", "VERB", "VERB", "VERB",
      "VERB", "VERB", "DET", "X", "PRON", "PRON", "ADV", "."), .Dim = c(68L,  2L))

  index = match(pennPOS, mtab[,1])
  output = rep("X", length(pennPOS))
  output[!is.na(index)] = mtab[index[!is.na(index)],2]
  output
}

#' Parse annotation xml
#'
#' Returns an annotation object from a character vector containing
#' the xml. Not exported; use \code{loadXMLAnnotation} instead.
#'
#' @importFrom   XML xmlRoot xmlParse xmlToDataFrame xmlChildren xmlAttrs
#' @param xml    character vector containing the xml file from an annotation
parseAnnoXML = function(xml) {
  xml = XML::xmlRoot(XML::xmlParse(xml))[[1]]
  coref = xml[[2]]
  xml = xml[[1]]
  sentences = XML::xmlChildren(xml)

  out = list(token=NULL,parse=NULL,basicDep=NULL,collapsedDep=NULL,
              collapsedProcDep=NULL, coref=NULL)

  if (length(sentences)==0L) {
    class(out) = "annotation"
    return(out)
  }

  tokenNames = c("sentence", "id", "word", "lemma", "CharacterOffsetBegin",
                 "CharacterOffsetEnd", "POS", "NER", "Speaker")
  depNames = c("sentence", "governor", "dependent", "type", "governorIdx",
               "dependentIdx")
  corefNames = c("corefId", "sentence", "start", "end", "head", "text")
  sentNames = c("id", "sentimentValue", "sentiment")

  for (i in 1:length(sentences)) {
    sent = sentences[[i]]
    df = data.frame(sentence=i, id=sapply(XML::xmlChildren(sent[[1L]]),function(v) XML::xmlAttrs(v)[1]),
                    XML::xmlToDataFrame(sent[[1]], stringsAsFactors=FALSE), stringsAsFactors=FALSE)

    index = match(tokenNames, names(df))
    if (length(index) != ncol(df)) df = df[,index[!is.na(index)]]
    if (any(is.na(index))) df = fillDF(df, tokenNames)

    out$token = rbind(out$token, df)

    if (!is.null(sent[[2L]]))
      out$parse = c(out$parse, XML::xmlValue(sent[[2]]))

    if (!is.null(sent[[3L]]) && length(XML::xmlToDataFrame(sent[[3L]]))) {
      df = data.frame(sentence=i, XML::xmlToDataFrame(sent[[3L]],stringsAsFactors=FALSE),
                type=sapply(XML::xmlChildren(sent[[3L]]),function(v) XML::xmlAttrs(v)[[1]]),
                governorIdx=sapply(XML::xmlChildren(sent[[3L]]), function(v) XML::xmlAttrs(v[[1]])[1]),
                dependentIdx=sapply(XML::xmlChildren(sent[[3L]]), function(v) XML::xmlAttrs(v[[2]])[1]),
                stringsAsFactors=FALSE)

      index = match(depNames, names(df))
      if (length(index) != ncol(df)) df = df[,index[!is.na(index)]]
      if (any(is.na(index))) df = fillDF(df, depNames)

      out$basicDep = rbind(out$basicDep, df)
    }

    if (!is.null(sent[[4L]]) && length(XML::xmlToDataFrame(sent[[4L]]))) {
      df = data.frame(sentence=i, XML::xmlToDataFrame(sent[[4L]],stringsAsFactors=FALSE),
                type=sapply(XML::xmlChildren(sent[[4L]]),function(v) XML::xmlAttrs(v)[[1]]),
                governorIdx=sapply(XML::xmlChildren(sent[[4L]]), function(v) XML::xmlAttrs(v[[1]])[1]),
                dependentIdx=sapply(XML::xmlChildren(sent[[4L]]), function(v) XML::xmlAttrs(v[[2]])[1]))

      index = match(depNames, names(df))
      if (length(index) != ncol(df)) df = df[,index[!is.na(index)]]
      if (any(is.na(index))) df = fillDF(df, depNames)

      out$collapsedDep = rbind(out$collapsedDep, df)
    }

    if (!is.null(sent[[5L]]) && length(XML::xmlToDataFrame(sent[[5L]]))) {
      df = data.frame(sentence=i, XML::xmlToDataFrame(sent[[5L]],stringsAsFactors=FALSE),
                type=sapply(XML::xmlChildren(sent[[5L]]),function(v) XML::xmlAttrs(v)[[1]]),
                governorIdx=sapply(XML::xmlChildren(sent[[5L]]), function(v) XML::xmlAttrs(v[[1]])[1]),
                dependentIdx=sapply(XML::xmlChildren(sent[[5L]]), function(v) XML::xmlAttrs(v[[2]])[1]))

      index = match(depNames, names(df))
      if (length(index) != ncol(df)) df = df[,index[!is.na(index)]]
      if (any(is.na(index))) df = fillDF(df, depNames)

      out$collapsedProcDep = rbind(out$collapsedProcDep, df)
    }

    sm = XML::xmlAttrs(sent)
    df = data.frame(matrix(sm,nrow=1),stringsAsFactors=FALSE)
    names(df) = names(sm)

    index = match(sentNames, names(df))
    if (length(index) != ncol(df)) df = df[,index[!is.na(index)],drop=FALSE]
    if (any(is.na(index))) df = fillDF(df, sentNames)

    if (nrow(df)) out$sentiment = rbind(out$sentiment, df)
  }

  if (!is.null(out$token)) {
    if (!is.na(index <- match("word", names(out$token))[1]))
      names(out$token)[index] = "token"
    if (sum(!is.na((index <- match(c("CharacterOffsetBegin","CharacterOffsetEnd"),names(out$token))))))
      out$token[,index] = apply(out$token[,index,drop=FALSE],2,as.integer)
  }

  if (!is.null(out$basicDep)) {
    if (sum(is.na((index <- match(c("governorIdx","dependentIdx"),names(out$basicDep))))))
      out$basicDep[,index] = apply(out$basicDep[,index,drop=FALSE],2,as.integer)
  }

  if (!is.null(out$collapsedDep)) {
    if (sum(!is.na((index <- match(c("governorIdx","dependentIdx"),names(out$collapsedDep))))))
      out$collapsedDep[,index] = apply(out$collapsedDep[,index,drop=FALSE],2,as.integer)
  }

  if (!is.null(out$collapsedProcDep)) {
    if (sum(!is.na((index <- match(c("governorIdx","dependentIdx"),names(out$collapsedProcDep))))))
      out$collapsedProcDep[,index] = apply(out$collapsedProcDep[,index,drop=FALSE],2,as.integer)
  }

  if (!is.null(coref)) {
    coref = XML::xmlChildren(coref)
    for (corefId in 1:length(coref)) {
      df = data.frame(corefId=corefId, XML::xmlToDataFrame(coref[[corefId]], stringsAsFactors=FALSE),
                      stringsAsFactors=FALSE)

      index = match(corefNames, names(df))
      if (length(index) != ncol(df)) df = df[,index[!is.na(index)]]
      if (any(is.na(index))) df = fillDF(df, corefNames)

      out$coref = rbind(out$coref, df)
    }
  }

  if (!is.null(out$sentiment)) {
    if (sum(!is.na(index <- match(c("id","sentimentValue"),names(out$sentiment)))))
      out$sentiment[,index[!is.na(index)]] = apply(out$sentiment[,index[!is.na(index)],drop=FALSE],2,as.integer)
  }

  class(out) = "annotation"
  out
}

fillDF = function(df, nameVec) {
  if (nrow(df) == 0L) return(df)
  index = match(nameVec, names(df))
  df[,nameVec[is.na(index)]] = NA
  return(df)
}
