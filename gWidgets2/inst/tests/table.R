w <- gwindow("gtable", visible=TRUE)
g <- ggroup(cont=w, horizontal=FALSE)

x <- data.frame("key"=state.name[1:10], value=state.x77[1:10,'Population'], stringsAsFactors=FALSE)
y <- x
y$icons <- rep("ok", length=10)
y$tooltips <- toupper(state.name[1:10])
##

tbl1 <- gtable(x, cont=g)
tbl2 <- gtable(y, icon.col=3, tooltip.col=4, cont=g)



## test
## svalue
expect_equal(svalue(tbl1), x[[1]][integer(0)])

## svalue<-
svalue(tbl1, index=TRUE) <- 1
expect_equal(svalue(tbl1, index=TRUE), 1)
expect_equal(svalue(tbl1, drop=TRUE), state.name[1])

## replace
tbl3 <- gtable(x[, 1:2], cont=g)
tbl3[] <- x
expect_equal(length(names(tbl3)), length(names(x)))
