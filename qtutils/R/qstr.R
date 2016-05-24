
qstr <- function(x, ...)
{
    UseMethod("qstr")
}


qstr.default <- function(x, ...)
{
    if (is.list(x)) return(qstr.list(x, ...))
    if (is(x, "standardGeneric")) return(qstr(S4methodsList(x), ...))
    temp <- tempfile()
    ## ostr <- capture.output(str(x, ...))
    ostr <- capture.output(print(x, ...))
    cat(paste(ostr, collapse = "\n"), file = temp)
    ## FIXME: no need to write out file
    w <- qeditor(temp, readonly = TRUE, richtext = FALSE,
                 rsyntax = FALSE, pointsize = 10)
    unlink(temp)
    w
}

## Qt objects should not be shown, especially QWidgets.  Even others
## may occasionally have problems when treated as environments, but
## not sure if it would be more effective to (1) write methods or (2)
## put in some try() statements in the environment method

qstr.QWidget <- function(x, ...) qstr.default(class(x))

## "qstr.try-error" <- ??? default OK?

qstr.data.frame <- function(x, ...)
{
    qdataview(x)
}

qstr.table <- function(x, ...)
{
    if (length(dim(x)) == 2) qdataview(x)
    else qstr.default(x, ...)
}

qstr.matrix <- function(x, ...)
{
    qdataview(x)
}

S4methodsList <- function(x)
{
    generic <- x@generic
    allDefs <- showMethods(generic, includeDefs = TRUE, printTo = FALSE)
    wz <- which(!nzchar(allDefs))[c(FALSE, TRUE)] ## blank lines come in sets of 2
    wz <- head(wz, -1) ## last one not needed
    starts <- 1L + c(1L, wz)
    ends <- c(wz, length(allDefs))
    stopifnot(length(starts) == length(ends))
    ans <- vector(mode = "list", length = length(starts))
    names(ans) <- sprintf("%s (%s)", generic, allDefs[starts])
    for (i in seq_along(ans))
    {
        ans[[i]] <- eval(parse(text = allDefs[(starts[i]+1):ends[i]]))
    }
    ans
}


qstr.function <- function(x, ...)
{
    temp <- tempfile()
    ostr <- capture.output(print(x, ...))
    cat(paste(ostr, collapse = "\n"), file = temp)
    wfun <- qeditor(temp, readonly = TRUE, richtext = FALSE,
                    rsyntax = TRUE, pointsize = 10)
    unlink(temp)
    ## qsetStyleSheet("font-family : monospace", widget = wfun)
    wfun
}

qstr.listOrEnv <- function(x, ...)
{
    isList <- is.list(x)
    objs <- 
        if (isList)
        {
            if (is.null(names(x)))
                names(x) <- as.character(seq_along(x))
            names(x)
        }
        else
            ls(x, all.names = TRUE)
    container <- Qt$QSplitter(1L) ## qsplitter(horizontal = TRUE)
    container$opaqueResize <- FALSE

    wlist <- Qt$QListWidget()
    wlist$addItems(objs)
    for (i in seq_along(objs))
    {
        if (!isList && bindingIsActive(objs[i], x))
        {
            obj.class <- obj.mode <- "Active binding"
        }
        else
        {
            ## obj.class <- class(x[[ objs[i] ]]) ## FIXME: use is()? S3 works?
            obj.class <- try(setdiff(is(x[[ objs[i] ]]), "oldClass"), silent = TRUE)
            obj.mode <- try(mode(x[[ objs[i] ]]), silent = TRUE)
        }
        if (is(obj.class, "try-error"))
            wlist$item(i-1L)$setToolTip(sprintf("<html>%s<br><strong>Error on evaluation: </strong>%s</html>",
                                                objs[i],
                                                as.character(obj.class)))
        else
            wlist$item(i-1L)$setToolTip(sprintf("<html>%s<br><strong>Class: </strong>%s<br><strong>Mode: </strong>%s</html>",
                                                objs[i],
                                                paste(obj.class, collapse = ","),
                                                obj.mode))
    }

    ## If wlist represents an environment, add a context-menu action
    ## to open an evaluation environment.
    if (!isList)
    {
        wlist$setContextMenuPolicy(Qt$Qt$ActionsContextMenu)
        replAct <- Qt$QAction(text = "Start Evaluation Interface",
                              parent = wlist)
        replAct$setShortcutContext(Qt$Qt$WidgetShortcut)
        replAct$setShortcut(Qt$QKeySequence("Ctrl+Return"))
        replHandler <- function(checked) {
            qreplu(x, ...)$show()
        }
        qconnect(replAct, signal = "triggered", handler = replHandler)
        wlist$addAction(replAct)
    }

    preview.container <- Qt$QWidget()
    preview.layout <- Qt$QGridLayout()
    preview.layout$setContentsMargins(0L, 0L, 0L, 0L)
    preview.layout$setSpacing(0L)
    preview.container$setLayout(preview.layout)

    container$addWidget(wlist)
    container$addWidget(preview.container)
    wlist$setSizePolicy(Qt$QSizePolicy$Preferred,
                        Qt$QSizePolicy$Expanding)
    ##qsetExpanding(wlist, horizontal = FALSE)
    preview.container$setSizePolicy(Qt$QSizePolicy$Preferred,
                                    Qt$QSizePolicy$Expanding)
    ## qsetExpanding(preview.container, horizontal = TRUE)
    container$setStretchFactor(0L, 0L)
    container$setStretchFactor(1L, 1L)
    ## qsetStretchFactor(container, 0L, 0L)
    ## qsetStretchFactor(container, 1L, 10L)

    sub.env <- new.env(parent = emptyenv())
    sub.env$preview <- NULL
    sub.env$objects <- objs
    sub.env$wlist <- wlist
    sub.env$preview.layout <- preview.layout
    sub.env$preview.container <- preview.container

    user.data <- list(x = x, sub.env = sub.env)
    handleSelection <- function(item) #user.data
    {
        i <- 1L + user.data$sub.env$wlist$currentRow
        obj <- user.data$sub.env$objects[i]
        new.preview <- qstr(try(user.data$x[[obj]], silent = TRUE))
        if (!is.null(user.data$sub.env$preview))
            user.data$sub.env$preview$close()
        user.data$sub.env$preview.layout$addWidget(new.preview, 0, 0)
        user.data$sub.env$preview <- new.preview
    }

##     qconnect(wlist,
##              user.data = user.data,
##              handler = handleSelection,
##              ## which = "cellClicked_int_int")
##              which = "itemClicked_qlistwidgetitem")
    ## attr(container, "activation.handler") <- 
    qconnect(wlist,
             signal = "itemActivated",
             handler = handleSelection)    
    ## user.data = user.data,
    container$resize(600, 400)
    container
}

qstr.list <- qstr.listOrEnv
qstr.environment <- qstr.listOrEnv




