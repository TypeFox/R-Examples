##' The Starting of Missing Data GUI.
##'
##' This function starts an open-files GUI, allowing 1) selecting one
##' or more data files; 2)opening the main missing-data GUI for one
##' data file. The missing data GUI consists of two tabs. In the
##' summary tab, there are a list of all variables, a list of
##' variables having missing values to color by, two radios for
##' imputation methods and graph types respectively, a checkbox group
##' for the conditional variables, four buttons and a graphics
##' device. In the help tab, the layout is the same as the summary
##' tab.  But when the users move their mouse on those widgets, or
##' click any of those radios or buttons, the functions of all widgets
##' will be described at the place of the graphics device. The
##' attributes of the variables can be changed. If the user double
##' clicks on any variables in the top left table of missing-data GUI,
##' an attribute window will pop up. Then the name could be edited,
##' and the class could be changed to one of the four classes:
##' integer, numeric, factor, and character. When a numeric variable
##' is changed to a categorical variable, the condtions in the bottom
##' left checkbox group will be updated. If the list of the color by
##' variables is very long, the selector allows text entry to find the
##' variable when this widget is active.
##'
##' If more than one files are listed in the window but no file is
##' focused when clicking the "Watch Missing Values", then the first
##' file is selected for the main missing-data GUI. If more than one
##' files are focused, then the first file of the focused files is
##' selected for the main GUI.
##' @param data A data frame which is shown in the main missing-data
##' GUI. If it is null, then the open-files GUI opens.
##' @param width the width of window. Default to be 1000, and the
##' minimal is 800.
##' @param height the height of window. Default to be 750, and the
##' minimal is 600.
##' @return NULL
##' @author Xiaoyue Cheng <\email{xycheng@@unomaha.edu}>
##' @export
##' @examples
##' if (interactive()) {
##' MissingDataGUI()
##'
##' data(tao)
##' MissingDataGUI(tao)
##'
##' data(brfss)
##' MissingDataGUI(brfss)
##' }
##'
MissingDataGUI = function(data=NULL, width=1000, height=750) {
    if (is.null(data)) {
        combo0 = gwindow("Open A File...", visible = TRUE)
        group = ggroup(horizontal = FALSE, container = combo0)
        f.list = matrix(nrow = 0, ncol = 1, dimnames = list(NULL, "File"))
        gt = gtable(f.list, multiple = T, container = group, expand = T)
		gb1 = gbutton("Open", container = group, handler = function(h,...) gt[,] = union(gt[,],na.omit(gfile(multiple=TRUE))))
        gb2 = gbutton("Watch Missing Values", container = group,handler = function(h, ...) WatchMissingValues(h, data=NULL, gt=gt))
    } else {
        if (is.data.frame(data)) {
            WatchMissingValues(data=data, size.width=width, size.height=height)
        } else {
            gmessage("Please use a data frame.")
            warning("The input needs to be a data frame.")
        }
    }
}
