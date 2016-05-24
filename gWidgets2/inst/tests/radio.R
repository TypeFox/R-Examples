w <- gwindow("radio buttons")
g <- ggroup(cont=w, horizontal=FALSE)

## vanilla
x <- state.name[1:4]
rb1 <- gradio(x, cont=g)

## with handler
rb2 <- gradio(x, cont=g, handler=function(h,...) print(svalue(rb2)))

## add handler
rb3 <- gradio(x, cont=g)
addHandlerChanged(rb3, handler=function(h,...) print("add handler"))

## tests

## svalue
expect_equal(svalue(rb1), x[1])
expect_equal(svalue(rb1, index=TRUE), 1)

## svalue <-
svalue(rb1) <- x[2]
expect_equal(svalue(rb1), x[2])
## by index
svalue(rb1, index=TRUE) <- 3
expect_equal(svalue(rb1), x[3])

##expect_warning(svalue(rb1, index=TRUE) <- length(rb1) + 1)

## [
expect_equal(rb1[], x)
expect_equal(rb1[2], x[2])

## [<-
rb1[] <- x[1:3]

## length
expect_equal(length(rb1), 3)

## invoke handler
expect_output(rb3$invoke_change_handler(), "add handler")

## numeric
rb2 <- gradio(1:4, cont=g)
svalue(rb2) <- 2
expect_equal(svalue(rb2), "2")

## coerce
rb2$coerce_with <- as.numeric
expect_equal(svalue(rb2), 2)
