
###################################################
### code chunk number 76: Controls.Rnw:772-775
###################################################
avail_DFs <- function() {
  c("", ".GlobalEnv", ProgGUIinR:::avail_dfs(.GlobalEnv))
}


###################################################
### code chunk number 77: Controls.Rnw:777-781
###################################################
##' get numeric variables
##'
##' @param where ".GlobalEnv" or name of data frame in global workspace
##' @return vector of numeric variable names


###################################################
### code chunk number 78: getNumeric
###################################################
get_numeric <- function(where) {
  val <- get(where, envir = .GlobalEnv)
  ProgGUIinR:::find_vars(val, is.numeric)
}


###################################################
### code chunk number 79: Controls.Rnw:791-817
###################################################
window <- gwindow("Find the mean", visible = FALSE)
group <- ggroup(cont = window, horizontal = FALSE)
group1 <- ggroup(cont = group)
glabel("Select data frame:", cont = group1)
df_combo_box <- gcombobox(avail_DFs(), cont = group1)
##
frame <- gframe("Arguments:", cont = group, horizontal=FALSE)
enabled(frame) <- FALSE
lyt <- glayout(cont = frame, expand = TRUE)
widget_list <- list() 
##
lyt[1,1] <- "x"
lyt[1,2] <- (widget_list$x <- gcombobox("           ",
                                        cont = lyt))
##
lyt[2,1] <- "trim"
lyt[2,2] <- 
  (widget_list$trim <- gslider(from = 0, to = 0.5, by = 0.01,
                               cont = lyt))
##
lyt[3,1] <- "na.rm"
lyt[3,2] <- 
  (widget_list$na.rm <- gcheckbox("", checked = TRUE, 
                                  cont = lyt))
group2 <- ggroup(cont = group)
compute_button <- gbutton("compute", cont = group2)


###################################################
### code chunk number 80: Controls.Rnw:829-837
###################################################
addHandlerChanged(df_combo_box, handler = function(h,...) {
  val <- svalue(h$obj)
  enabled(frame) <- val !=""
  enabled(compute_button) <- val != ""
  if(val != "") 
    widget_list$x[] <- get_numeric(val)
  svalue(widget_list$x, index = TRUE) <- 0
})


###################################################
### code chunk number 81: computeHandler
###################################################
addHandlerChanged(compute_button, handler = function(h,...) {
  out <- lapply(widget_list, svalue)
  out$x <- get(out$x, get(svalue(df_combo_box),
                          envir = .GlobalEnv))
  print(do.call(mean.default, out))
})


###################################################
### code chunk number 82: visible
###################################################
visible(window) <- TRUE
