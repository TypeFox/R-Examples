f.aggregate <- function(data){
## AGGREGATES A DATA FRAME (USING SOMETHING SIMILAR TO aggregate.data.frame, 
## ONLY BETTER...).
## REMOVES IDENTICAL ROW COPIES AND ADDS A UNIQUE tag PLUS THE freq.
## IN ADDITION, AN ATTRIBUTE IS ADDED, GIVING THE ORIGINAL LINE NUMBERS
## CORRESPONDING TO EACH LINE IN THE AGGREGATED FILE
## THE RETURNED DATA FRAME IS IN THE SAME ORDER AS THE ORIGINAL, ONLY THE 
## NON-UNIQUE ROWS HAVE BEEN REMOVED (FIRST OF EACH ROW TYPE IS KEPT)
##
	.nlines <- dim(data)[1]
	# CREATE A UNIQUE TAG:
	.tag <- f.create.tag(data) 
	# KEEP TRACK OF OLD LINE NUMBERS. .tag CAN LATER ON INDEX .lines IF NECESS.:
	.lines <- tapply(1:.nlines, .tag, function(x)x)
	.freq <- f.groupsum(X = rep(1, .nlines), INDICES = .tag)
	.unique <- !duplicated(.tag)
	.tag.unique <- .tag[.unique]
	.data.agg <- dframe(data, freq = .freq)
	.data.agg <- .data.agg[.unique, , drop = F]
	.lines <- .lines[.tag.unique] # RE-ORDERS .lines TO SAME SEQUENCE AS LINES IN DATA. NOTE THAT NAMES OF .lines ALSO CORRESPONDS TO THE TAG VALUES
#
# ADD ORIGINAL LINE NUMBERS TO OUTPUT:

	attr(.data.agg, "orig.lines") <- .lines
	
	return(.data.agg)
}
