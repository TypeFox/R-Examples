# make a tableplot for collinearity diagnostics


tableplot.colldiag <- function(values,
		cell.specs, 
		prop.col = c("white", "pink", "red"),      # colors for variance proportions
		cond.col = c("#A8F48D", "#DDAB3E", "red"), # colors for condition indices
		cond.max = 100,                            # scale.max for condition indices
		prop.breaks = c(0, 20, 50, 100),           # breaks for variance proportions
		cond.breaks = c(0, 5, 10, 1000),           # breaks for condition indices
		show.rows = nvar:1,                        # order and rows to display
		...
)
{
# TODO: this is cheating! Should we just require a colldiag object?
  if(inherits(values,"lm")) {
  	if(!require("perturb")) stop("requires perturb package for lm objects")
  	values <- colldiag(values, add.intercept=FALSE, center=TRUE)
  } else {
  	stopifnot(inherits(values,"colldiag"))
  }

	collin <- round(100 * values$pi)        # variance proportions
	condind <- round(values$condindx,2)     # condition indices
	vars <- colnames(values$pi)             # variable names
	nvar <- ncol(values$pi)                 # number of variables

  if(missing(cell.specs)) {
    cell.specs <- make.specs(
    	shape=c(rep(0,3), rep(2,3)),           # squares and circles
      cell.fill=c(prop.col, rep("white",3)),
    	back.fill=c(rep("white",3), cond.col),
    	scale.max=rep(c(100,cond.max), each=3),
    	label=1, label.size=1, 
    	ref.lines=rep(c(TRUE,FALSE), each=3) 
    	)
  	}

#browser()

  r.label = paste("#", show.rows, sep="")
  table <-cbind(condind, collin)[show.rows,]
  rownames(table) <- r.label
  vars <- c("CondIndex", vars)
  colnames(table) <- vars

  type.mat <- matrix(cut(collin, breaks=prop.breaks-0.1, labels=FALSE), dim(collin))  
  type.cond <- 3+cut(condind, breaks=cond.breaks-0.1, labels=FALSE)
  types <- cbind(type.cond, type.mat)[show.rows,]

#  h.parts <- as.vector(xtabs(~ cut(condind, breaks=cond.breaks-0.1, labels=FALSE)))
  h.parts <- rev(as.vector(xtabs(~ type.cond)))
  tableplot.default(table, types,
  	cell.specs=cell.specs, 
  	v.parts=c(1,nvar), 
  	h.parts=h.parts,
  	...
  	)
}
