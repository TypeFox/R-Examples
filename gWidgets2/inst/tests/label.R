w <- gwindow("test")
g <- ggroup(cont=w, horizontal=FALSE)

l1 <- glabel("one", cont=g)
l2 <- glabel("<b>bold markup</b>", markup=TRUE, cont=g)
l3 <- glabel("bold", cont=g)
l4 <- glabel("red", cont=g)
l5 <- glabel("LARGE", cont=g)

font(l3) <- list(weight="bold")
font(l4) <- list(color="red")
font(l5) <- list(scale="large")

## tests
## svalue
expect_equal(svalue(l1), "one")

## svalue<-
svalue(l1) <- "two"
expect_equal(svalue(l1), "two")
