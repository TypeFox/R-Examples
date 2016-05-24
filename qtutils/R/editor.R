
##' An editor widget.
##'
##' This function creates an editor widget, optionally initializing it
##' with the contents of a file.
##' 
##' @title An editor widget.
##' @param file Character string giving the name of the file to be
##' displayed.
##' @param readonly Logical flag indicating whether the editor should
##' be in read-only initially.
##' @param richtext Logical flag indicating whether the editor should
##' be cabaple of displaying rich-text (HTML).
##' @param ... Further arguments, passed on to
##' \code{\link{RCodeEditor}} if relevant.
##' @param rsyntax Logical flag indicating whether R syntax
##' highlighting should be applied to the contents of the editor. 
##' 
##' @return If \code{richtext=TRUE}, a QTextEdit instance.  Otherwise
##' a \code{\link{RCodeEditor}} instance.
##' 
##' @author Deepayan Sarkar
qeditor <- function(file = NULL,
                    readonly = FALSE,
                    richtext = FALSE,
                    ...,
                    rsyntax = tail(strsplit(basename(file), ".", fixed = TRUE)[[1]], 1) %in% c("R", "r", "S", "r"))
{
    if (richtext) 
    {
        edit <- Qt$QTextEdit()
        edit$setAcceptRichText(TRUE)
    }
    else
    {
        edit <- RCodeEditor(...)
        ## edit$setAcceptRichText(FALSE)
        ## ## FIXME: edit$setFont(qfont(family = "monospace"))
        ## edit$setFontFamily("monospace")
	## edit$setLineWrapMode(Qt$QTextEdit$NoWrap)
    }
    if (rsyntax)
        .Call(qt_qsetRSyntaxHighlighter, edit)
    
### NOTE NOTE NOTE: Cannot use QFile and QIODevice (Error: Cannot
### handle Moc type 'qint64')
    ## if (!is.null(file))
    ## {
    ##     qfile <- Qt$QFile(file)
    ##     status <-
    ##         if (readonly) qfile$open(Qt$QIODevice$ReadOnly)
    ##         else qfile$open(Qt$QIODevice$ReadWrite)
    ##     if (!status) return(NULL)
    ## }
    ## if (!is.null(file)) 
    ## {
    ##     stream <- Qt$QTextStream()
    ##     stream$setDevice(qfile)
    ##     ## stream <- Qt$QTextStream(qfile)
    ##     txt <- stream$readAll()
    ##     if (!is.null(txt)) edit$setText(txt)
    ##     if (readonly) edit$setReadOnly(TRUE)
    ##     qfile$close()
    ## }

### Instead, handle purely in R
    if (!is.null(file) && nzchar(file)) 
    {
        edit$setPlainText(paste(readLines(file, warn = FALSE),
                                collapse = "\n"))
    }
    if (readonly) edit$setReadOnly(TRUE)
    cursor <- edit$textCursor()
    cursor$setPosition(0)
    edit$setTextCursor(cursor)
    edit$ensureCursorVisible()

    edit$resize(600, 400)
    edit
}




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' 
##' @title A pager replacement
##' @param file The file whose content is to be displayed.
##' @param header Ignored.
##' @param title Character string to be used as the widget title if applicable.
##' @param delete.file Logical flag indicating whether the file should
##' be deleted after being shown. 
##' @param parent An optional QTabWidget instance.  If supplied, the
##' editor widget is added as a child.
##' 
##' @return A QTextEdit or RCodeEditor instance (see \code{\link{qeditor}}).
##' 
##' @author Deepayan Sarkar
qpager <- function(file,
                   header = "",
                   title = "R Information",
                   delete.file = FALSE,
                   parent = NULL)
{
    ans <- qeditor(file = file, readonly = TRUE)
    if (delete.file)
    {
        if (getOption("verbose"))
            warning(sprintf("Deleting file: %s", file))
        unlink(file)
    }
    if (!is.null(parent) && is(parent, "QTabWidget"))
    {
        parent$addTab(ans, title)
        invisible(ans)
    }
    else 
    {
        ans$setWindowTitle(title)
        ans$resize(600, 400)
        ans
    }
}

