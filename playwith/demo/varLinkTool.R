
library(playwith)

## A tool to create a linked scatterplot of covariates.

varLink_handler <- function(widget, playState) {
    dat <- getDataArg(playState)
    if (is.null(dat)) {
        gmessage("Can not extract data object from the call.",
                 icon = "error", parent = playState$win)
        return()
    }
    NONE <- "(none)"
    nm <- c(NONE, names(dat))
    var1W <- gdroplist(nm, editable = TRUE)
    var2W <- gdroplist(nm, editable = TRUE)
    lay <- glayout()
    lay[1,1:3] <-
        glabel("<b>Create a linked plot with formula:</b>",
               markup = TRUE)
    lay[2,1] <- "[Response expression]"
    lay[2,2] <- " ~ "
    lay[2,3] <- "[Explanatory expression]"
    lay[3,1] <- var1W
    lay[3,2] <- " ~ "
    lay[3,3] <- var2W
    visible(lay) <- TRUE
    gbasicdialog(title = "Link covariates",
                 widget = lay,
                 handler = function(h, ...) {
                     v1txt <- svalue(var1W)
                     v2txt <- svalue(var2W)
                     dispose(h$obj)
                     if (v1txt %in% c("", NONE))
                         v1txt <- "NULL"
                     if (v2txt %in% c("", NONE))
                         v2txt <- "NULL"
                     var1 <- parse(text = v1txt)[[1]]
                     var2 <- parse(text = v2txt)[[1]]
                     form <- call("~", var1, var2)
                     newCall <- call("xyplot", form)
                     if (is.null(var1) && is.null(var2))
                         return()
                     if (is.null(var1) || is.null(var2)) {
                         if (is.null(var1))
                             form <- call("~", var2)
                         else form <- call("~", var1)
                         newCall <- call("qqmath", form)
                     }
                     newCall$data <- quote(dat)
                     playwith(plot.call = newCall,
                              new = TRUE,
                              link.to = playState,
                              click.mode = "Brush")
                 })
}

varLinkTool <- list("VarLink", "gtk-dnd-multiple", "Link vars",
                   callback = varLink_handler)

playwith(xyplot(temperature ~ wind, data = environmental),
         tools = list(varLinkTool), click.mode = "Brush")
