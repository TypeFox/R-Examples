
gedit <- function(..., handler = NULL, action = NULL) {
    wid <- gWidgets::gedit(...)
    if (!is.null(handler))
        geditAddGoodHandlers(wid, handler = handler,
                             action = action)
    wid
}

geditAddGoodHandlers <- function(wid, handler, ...) {
    ## gedit event handler only triggered when Enter pressed:
    addHandlerChanged(wid, handler = handler, ...)

    ## need to also detect changes (keystrokes)
    ## and update when lose focus

    ## local state variable (exists in function environment)
    iNeedUpdating <- FALSE
    ## keystroke events trigger the flag
    setNeedUpdate <- function(h, ...)
        iNeedUpdating <<- TRUE
    addHandlerKeystroke(wid, handler = setNeedUpdate)
    ## when the widget loses focus, do the update
    doUpdateIfNeeded <- function(h, ...) {
        if (iNeedUpdating)
            handler(h, ...)
        iNeedUpdating <<- FALSE
    }
    addHandlerBlur(wid, handler = doUpdateIfNeeded, ...)
}


shrinkrange <- function(r, f = 0.1)
{
  stopifnot(length(r) == 2)
  orig.d <- diff(r) / (1 + 2*f)
  orig.r <- r - c(-f, f) * orig.d
  orig.r
}

is.somesortoftime <- function(x) {
  inherits(x, "Date") ||
  inherits(x, "POSIXt") ||
  inherits(x, "yearmon") ||
  inherits(x, "yearqtr")
}

gmessage.error <- function(message, title="Error", icon="error", ...)
    gmessage(message, title=title, icon=icon, ...)

Filters <- matrix(c(
                    "R or S files (*.R,*.q,*.ssc,*.S)", "*.R;*.q;*.ssc;*.S",
                    "Postscript files (*.ps)",          "*.ps",
                    "Encapsulated Postscript (*.eps)",  "*.eps",
                    "PDF files (*.pdf)",                "*.pdf",
                    "Png files (*.png)",                "*.png",
                    "Jpeg files (*.jpeg,*.jpg)",        "*.jpeg;*.jpg",
                    "Text files (*.txt)",               "*.txt",
                    "R images (*.RData,*.rda)",         "*.RData;*.rda",
                    "Zip files (*.zip)",                "*.zip",
                    "SVG files (*.svg)",                "*.svg",
                    "Windows Metafiles (*.wmf,*.emf)",  "*.wmf;*.emf",
                    "xfig files (*.fig)",               "*.fig",
                    "All files (*.*)",                  "*.*"), ncol=2, byrow=T,
                  dimnames=list(c('R','ps','eps','pdf','png','jpeg','txt',
                  'RData','zip','svg','wmf','fig','All'),NULL))

get.extension <- function(path)
{
    ## Extract and return the extension part of a filename

    parts <- strsplit(path, "\\.")[[1]]
    if (length(parts) > 1)
        last <- parts[length(parts)]
    else
        last <- ""
    last
}

## RGtk2 helpers

guiTextView <-
    function(text, title = "Text View",
             wrap.mode=c("none", "char", "word", "word_char"),
             size=c(640, 400))
{
    wrap.mode <- match.arg(wrap.mode)
    win <- gtkWindow(show = FALSE)
    win["title"] <- title
    win$setDefaultSize(size[1], size[2])
    editTV <- gtkTextView()
    setTextviewMonospace(editTV)
    editTV$setWrapMode(GtkWrapMode[wrap.mode])
    setTextview(editTV, text)
    scroller <- gtkScrolledWindow()
    scroller$add(editTV)
    scroller$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"])
    win$add(scroller)
    win$show()
}

        ## possible with gWidgets, but way too slow.
        #txtBox <- gtext(callTxt, font.attr=c(family="monospace"), wrap=FALSE, width=600)
        #gbasicdialog(title="Edit plot call", widget=txtBox,
        #             action=environment(), handler=function(h, ...)
        #             assign("newTxt", svalue(h[[1]]), env=h$action)
        #             )


guiTextInput <-
    function(text="",
             title="Text Input",
             prompt="",
             oneLiner=FALSE,
             accepts.tab=TRUE,
             wrap.mode=c("none", "char", "word", "word_char"),
             size=c(640, 320),
             width.chars=-1,
             focus.on.ok=!oneLiner)
{
    wrap.mode <- match.arg(wrap.mode)
    ## construct dialog
    editBox <- gtkDialog(title=title, NULL, NULL,
                         "OK", GtkResponseType["ok"], "Cancel", GtkResponseType["cancel"],
                         show = FALSE)
    editBox$setDefaultResponse(GtkResponseType["ok"])
    if (nchar(prompt) > 0) {
        editBox$vbox$packStart(gtkLabel(prompt), expand=FALSE, pad=2)
    }
    if (oneLiner) {
        editEntry <- gtkEntry()
        editEntry["activates-default"] <- TRUE
        editEntry["text"] <- text
        editEntry["width-chars"] <- width.chars
        editBox$vbox$packStart(editEntry, pad=10)
    } else {
        editBox$setDefaultSize(size[1], size[2])
        editTV <- gtkTextView()
        setTextviewMonospace(editTV)
        editTV$setWrapMode(GtkWrapMode[wrap.mode])
        editTV$setAcceptsTab(accepts.tab)
        setTextview(editTV, text)
        scroller <- gtkScrolledWindow()
        scroller$add(editTV)
        scroller$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"])
        editBox$vbox$packStart(scroller)
    }
    ## put focus on the OK button
    if (focus.on.ok) editBox$actionArea$getChildren()[[2]]$grabFocus()
    result <- editBox$run() ## make it modal
    newTxt <- if (oneLiner) editEntry["text"] else getTextviewText(editTV)
    editBox$destroy()
    if (result != GtkResponseType["ok"]) return(invisible(NULL))
    newTxt
}

setTextview <- function(tv, ..., sep="")
{
    msg <- paste(sep=sep, ...)
    if (length(msg) == 0) msg <-""
    tv$getBuffer()$setText(msg)
    invisible(NULL)
}

getTextviewText <- function(tv)
{
  ## Extract text content of specified textview
  log.buf <- tv$getBuffer()
  start <- log.buf$getStartIter()$iter
  end <- log.buf$getEndIter()$iter
  return(log.buf$getText(start, end))
}

setTextviewMonospace <- function(tv)
{
    tv$modifyFont(pangoFontDescriptionFromString("monospace 10"))
    invisible(NULL)
}

pangoEscape <- function(x)
{
    x <- gsub('%', '%%', x)
    x <- gsub('&', '&amp;', x)
    x <- gsub('<', '&lt;', x)
    x <- gsub('>', '&gt;', x)
                                        #x <- gsub('&&', '&amp;&amp;', x)
                                        #x <- gsub('& ', '&amp; ', x)
                                        #x <- gsub('<<', '&lt;&lt;', x)
                                        #x <- gsub('<-', '&lt;-', x)
                                        #x <- gsub('< ', '&lt; ', x)
    x
}
