
#' Select item(s) from a list.
#'
#' Provides the same functionality as choose.files from utils package,
#' but relies on Java and rJava package and therefore is system independent provided Java 1.5 and higher is installed.
#' This brings up a modal dialog box with a (scrollable) list of items, which can be selected by the mouse. 
#' If multiple is true, further items can be selected or deselected by holding the control key down whilst selecting,
#' and shift-clicking can be used to select ranges.
#' Normal termination is via the 'OK' button or by hitting Enter or double-clicking an item. 
#' Selection can be aborted via the 'Cancel' button or pressing Escape.
#' If no graphical widget is available it displays a text list from which the user can choose by number(s). 
#' The multiple = FALSE case uses menu. Preselection is only supported for multiple = TRUE, where it is indicated by a '+' preceding the item.
#' It is an error to use select.list in a non-interactive session.
#'
#' @note jselect.list() is called internally by rchoose.list() if it's appropriate for a given platform/graphics combination.
#' Calling jselect.list() directly forces the package to use Java based dialog regardless of system capabilities and therefore may fail.
#' Use the direct call to jselect.list() only if it seems beneficial to bypass the rchoose.list() decision logic.
#'
#' @name jselect.list
#' @param choices A character vector of items.
#' @param preselect A character vector, or NULL. If non-null and if the string(s) appear in the list, the item(s) are selected initially.
#' @param multiple Logical: can more than one item be selected?
#' @param title Optional character string for window title, or NULL for no title.
#' @param modal Indicates how the modality of the dialog is implemented.
#' return A character vector of selected items. If multiple is false and no item was selected (or Cancel was used), '' is returned.
#' If multiple is true and no item was selected (or Cancel was used) then a character vector of length 0 is returned. 
#' @seealso {\code{\link{rselect.list}}, \code{\link{select.list}}}
#' @examples
#' \dontrun{
#' jselect.list(c("Peter", "Alex", "Roger", "Leah"),title="Select", multiple=TRUE);
#' }
#' @export jselect.list
#' @author  Alex Lisovich, Roger Day

jselect.list<-function(choices, preselect = NULL, multiple = FALSE, title = NULL,modal=canUseJavaModal()){

	jchoices<-.jarray(choices);

	selector<-new(J("rjavautils.rJavaListChooser"),NULL,modal);
	selector$setTitle(title);
	selector$setMultipleSelection(multiple);
	selector$setSelectionValues(jchoices);
	if(!is.null(preselect)){
		indices<-as.integer(which(choices%in% preselect)-1);
		selector$setInitialSelection(.jarray(indices));
	}

	selection<-selector$showDialog();
	if(!modal){
		while (selector$isRunning()){};
		selector$dispose();
		selection<-selector$getSelection();
	}

	return(as.vector(selection));
}

#' Select item(s) from a list.
#'
#' Provides the same functionality as select.list from the utils package, but is intended to broaden
#' the range of systems where selection can be made through the graphical dialog provided Java 1.5 or higher is installed.
#' This brings up a modal dialog box with a (scrollable) list of items, 
#' which can be selected by the mouse. If multiple is true, further items can be selected or deselected by holding
#' the control key down whilst selecting, and shift-clicking can be used to select ranges.
#' Normal termination is via the 'OK' button or by hitting Enter or double-clicking an item. 
#' Selection can be aborted via the 'Cancel' button or pressing Escape.
#' If no graphical widget is available it displays a text list from which the user can choose by number(s). 
#' The multiple = FALSE case uses menu. Preselection is only supported for multiple = TRUE, where it is indicated by a '+' preceding the item.
#' It is an error to use rselect.list in a non-interactive session.
#'
#' @name rselect.list
#' @param choices A character vector of items.
#' @param preselect A character vector, or NULL. If non-null and if the string(s) appear in the list, the item(s) are selected initially.
#' @param multiple Logical: can more than one item be selected?
#' @param title Optional character string for window title, or NULL for no title.
#' @param graphics logical indicating if a graphical widget should be used.
#' return A character vector of selected items. If multiple is false and no item was selected (or Cancel was used), '' is returned.
#' If multiple is true and no item was selected (or Cancel was used) then a character vector of length 0 is returned.
#' @seealso {\code{\link{jselect.list}}, \code{\link{select.list}}}
#' @examples
#' \dontrun{
#' rselect.list(c("Peter", "Alex", "Roger", "Leah"),title="Select", multiple=TRUE);
#' }
#' @export rselect.list
#' @author  Alex Lisovich, Roger Day

rselect.list<-function(choices, preselect = NULL, multiple = FALSE, title = NULL,graphics = getOption("menu.graphics")) {

	if(canUseJava())
		return(jselect.list(choices,preselect,multiple,title));

	if(.Platform$OS.type=="windows" || .Platform$GUI=="AQUA")
		return(select.list(choices,preselect,multiple,title));
		
	if(canUseTclTk())
		return(tcltk::tk_select.list(choices,preselect,multiple,title));
 
	return(select.list(choices,preselect,multiple,title));
}



