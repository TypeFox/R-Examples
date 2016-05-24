# The goal of RGtkDataFrame is to provide a (flat) GtkTreeModel implementation
# that responds to the familiar subset and replacement operators for data frames
# with minimal performance lost relative to a native data frame.

# This allows the R programmer to interact with the GUI without sacrificing
# familiarity and speed.

"[<-.RGtkDataFrame" <-
function(x, i, j, value)
{
	frame <- as.data.frame(x)
	
	old_nrow <- nrow(frame)
	old_ncol <- ncol(frame)
	
	if (missing(i))
		i <- 1:old_nrow
	if (missing(j))
		j <- 1:old_ncol
	
	frame[i, j] <- value
	
	if (is.character(i))
		i <- match(i, rownames(frame))
  else if (is.logical(i))
    i <- which(i)
	if (is.character(j))
		j <- match(j, colnames(frame))
  else if (is.logical(j))
    j <- which(j)
	
	changed <- integer(0)
	if (length(unique(j)) > ncol(frame) - old_ncol)
		changed <- i # existing columns changed, all specified rows "changed"
	else if (nrow(frame) > old_nrow) # otherwise, just add new rows
		changed <- ((old_nrow+1):nrow(frame))
	
	resort <- x$getSortColumnId() %in% j
	
	.RGtkCall("R_r_gtk_data_frame_set", x, frame, as.list(as.integer(changed-1)), resort)
	
	x
}

"[.RGtkDataFrame" <- function(x, i, j, drop = T)
{
	frame <- as.data.frame(x)
	if (!missing(i) && length(i) > 0 && inherits(i[[1]], "GtkTreePath"))
		i <- .RGtkCall("R_gtk_tree_paths_to_indices", i)+1
	frame[i, j, drop=drop]
}

rGtkDataFrame <- rGtkDataFrameNew <- function(frame = data.frame())
{
  sort_closure <- function(frame, col, order) {
    new_order <- order(frame[,col+1],decreasing=order)
    list(frame[new_order,drop=F],
         as.integer((1:length(new_order))[new_order]-1))
  }
  
  w <- .RGtkCall("R_r_gtk_data_frame_new", as.data.frame(frame), sort_closure)
  w
}

as.data.frame.RGtkDataFrame <- function(x, ...) .RGtkCall("R_r_gtk_data_frame_get", x)

rGtkDataFrameAppendRows <- function(x, ...) {
	frame <- as.data.frame(x)
  new_frame <- rbind(frame, ...)
	new_ind <- (nrow(frame)+1):nrow(new_frame)
	if (nrow(new_frame) > nrow(frame))
		.RGtkCall("R_r_gtk_data_frame_set", x, new_frame, as.list(as.integer(new_ind-1)), T)
	x
}
rGtkDataFrameAppendColumns <- function(x, ...) {
	frame <- as.data.frame(x)
	new_frame <- cbind(frame, ...)
	if (ncol(new_frame) > ncol(frame))
		.RGtkCall("R_r_gtk_data_frame_set", x, new_frame, NULL, F)
	x
}

 rGtkDataFrameSetFrame <- function(x, frame = data.frame()) {
  rows <- list()
  if (nrow(frame) > 0)
    rows <- as.list(as.integer(1:nrow(frame)-1))
  .RGtkCall("R_r_gtk_data_frame_set", x, frame, rows, T)
}

dim.RGtkDataFrame <- function(x, ...) {
  dim(as.data.frame(x))
}
dimnames.RGtkDataFrame <- function(x, ...) {
	dimnames(as.data.frame(x))
}
"dimnames<-.RGtkDataFrame" <- function(x, value) {
  frame <- as.data.frame(x)
  dimnames(frame) <- value
  .RGtkCall("R_r_gtk_data_frame_set", x, frame, NULL, F)
  x
}
