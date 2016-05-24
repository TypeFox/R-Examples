
## Convenient interfaces to various built-in modal dialogs such as
## file chooser, color picker, integer/text/etc. input.


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Qt-based file chooser
##' @param caption 
##' @param dir 
##' @param parent 
##' @author Deepayan Sarkar
qfile.choose <-
    function(caption = "", dir = "", filter = "",
             allow.new = FALSE, parent = NULL)
{
    FUN <-
        if (allow.new) Qt$QFileDialog$getSaveFileName
        else Qt$QFileDialog$getOpenFileName
    ans <- FUN(parent, as.character(caption),
               path.expand(dir), as.character(filter))
    ans
}

qdir.choose <- function(caption = "", dir = "", parent = NULL)
{
    FUN <- Qt$QFileDialog$getExistingDirectory
    ans <- FUN(parent,
               as.character(caption),
               as.character(dir))
    ans
}

qgetColor <- function(parent = NULL, title = "", alpha = 1)
{
    FUN <- Qt$QColorDialog$getColor
    ans <- FUN(Qt$QColor(255, 255, 255),
               parent,
               as.character(title))
    if (ans$isValid())
        rgb(ans$red(), ans$green(), ans$blue(),
            alpha * 255, maxColorValue = 255)
    else NA_character_
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' 
##' @title qgetDouble
##' 
##' @param title window title
##' @param label Character string giving a label for the selection. 
##' @param value Initial value.
##' @param minValue Minimum value allowed.
##' @param maxValue Maximum value allowed.
##' @param decimals The maximum number of decimals allowed. 
##' @param parent 
##' @return NULL
##' @author Deepayan Sarkar
qgetDouble <-
    function(title = "Choose value",
             label = "Enter a numeric value",
             value = 0.0,
             minValue = -.Machine$double.xmax,
             maxValue = .Machine$double.xmax,
             decimals = 3L, parent = NULL)
{
    FUN <- Qt$QInputDialog$getDouble
    ans <- FUN(parent,
               as.character(title),
               as.character(label),
               as.double(value),
               as.double(minValue),
               as.double(maxValue),
               as.integer(decimals))
    (ans)
}

qgetInteger <-
    function(title = "Choose value", label = "Enter an integer value", value = 0L,
             minValue = -.Machine$integer.max,
             maxValue = .Machine$integer.max,
             step = 1L, parent = NULL)
{
    FUN <- Qt$QInputDialog$getInt
    ans <- FUN(parent,
               as.character(title),
               as.character(label),
               as.integer(value),
               as.integer(minValue),
               as.integer(maxValue),
               as.integer(step))
    ans
}

qgetText <-
    function(title = "Choose value", label = "Enter a text string",
             text = "",
             echomode = c("normal", "noecho", "password", "echoonedit"),
             parent = NULL)
{
    echomode <- match.arg(echomode)
    emode <- c("normal" = 0L, "noecho" = 1L, "password" = 2L, "echoonedit" = 3L)[echomode]
    FUN <- Qt$QInputDialog$getText
    ans <- FUN(parent,
               as.character(title),
               as.character(label),
               emode,
               as.character(text))
    ans
}

## getItem ( QWidget * parent, const QString & title, const QString &
##          label, const QStringList & items, int current = 0, bool
##          editable = true, bool * ok = 0, Qt::WindowFlags flags = 0 )

qgetItem <-
    function(title = "Choose value", label = "Choose a list item",
             items = "",
             current = 1L, editable = TRUE,
             parent = NULL)
{
    FUN <- Qt$QInputDialog$getItem
    ans <- FUN(parent,
               as.character(title),
               as.character(label),
               items,
               as.integer(current - 1L),
               as.logical(editable))
    ans
}

