#' # pyOptions
require(testthat)
require(PythonInR)
invisible(capture.output(pyConnect()))

## PySource
tmpfile <- tempfile(fileext=".R")
writeLines(c("x <- 3", "print(x)", "BEGIN.Python()", 
             "x=3**3", "print(u'Hello R!\\n')", 
             "END.Python"), tmpfile)
expect_that(pySource(tmpfile), prints_text("Hello R!"))

readLines(tmpfile)

x <- c("x=3**3", "print(u'Hello R!\\nHello R!')")
x <- paste(x, collapse="\n")
cat(x)
pyExec(x)
             
pyExec("\nx=3**3\nprint(u'Hello R!\\nHello R!\\n')\n")
pyExecp("3")


    code <- readLines(tmpfile)
    temp <- tempfile()
    code <- paste(code, collapse = "\n")
    m <- unlist(regmatches(code, gregexpr("BEGIN\\.Python\\(\\)(.*?)END\\.Python", 
        code)))
    repl <- gsub("END.Python", "", gsub("BEGIN.Python()", "", 
        m, fixed = T), fixed = T)
    repl <- sprintf('pyExec("""%s""")', repl)
    eval(parse(text=repl))
    for (i in 1:length(m)) {
        code <- sub(m[i], repl[i], code, fixed = TRUE)
    }
    cat(code)
    writeLines(code, temp)
    source(temp, verbose=TRUE, max.deparse.length=500)
    str(source)
    repl

    getOption("encoding")
    cat(readLines(temp), sep="\n")
