w <- gwindow("slider spinbutton")
g <- ggroup(cont=w, horizontal=FALSE)
g$set_borderwidth(11)

sl1 <- gslider(from=0, to=100, by=1, cont=g)
sl2 <- gslider(from=0, to=100, by=1, value=50, cont=g, handler=function(h,...) print("call handler"))


## test

## svalue
expect_equal(svalue(sl1), 0)
expect_equal(svalue(sl2), 50)

## svalue<-
svalue(sl2) <- 5
expect_equal(svalue(sl2), 5)


## handler
expect_output(sl2$invoke_change_handler(), "call handler")


## spinbutton

sp1 <- gspinbutton(from=0, to=100, by=1, cont=g)
sp2 <- gspinbutton(from=0, to=100, by=1, value=50, cont=g, handler=function(h,...) print("call handler"))


## test

## svalue
expect_equal(svalue(sp1), 0)
expect_equal(svalue(sp2), 50)

## svalue<-
svalue(sp2) <- 5
expect_equal(svalue(sp2), 5)


## handler
expect_output(sp2$invoke_change_handler(), "call handler")
