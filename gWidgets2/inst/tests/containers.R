w <- gwindow("containers")
g <- ggroup(cont=w, horizontal=FALSE)

##########
## notebook
nb <- gnotebook(cont=g)
glabel("page one", label="page one", cont=nb)
glabel("page two", label="page two", cont=nb)

## svalue
expect_equal(svalue(nb), 2)
## svalue<-
svalue(nb) <- 1
expect_equal(svalue(nb), 1)

## length
expect_equal(length(nb), 2)

## names
expect_equal(names(nb), c("page one", "page two"))

names(nb) <- c("one", "and two")
expect_equal(names(nb), c("one", "and two"))

## dispose
dispose(nb)
expect_equal(length(nb), 1)

##########
## stack widget
cards <- gstackwidget(cont=g)
glabel("card 1", cont=cards)
glabel("card 2", cont=cards)


## svalue
expect_equal(svalue(cards), 2)
## svalue<-
svalue(cards) <- 1
expect_equal(svalue(cards), 1)

## length
expect_equal(length(cards), 2)

## dispose
dispose(cards)
expect_equal(length(cards), 1)

#########
## frame
f <- gframe("label", cont=g)
glabel("in frame", cont=f)

## names
expect_equal(names(f), "label")



#########
## gexpandgroup
expgp <- gexpandgroup("expandgroup", cont=g)
glabel("in expandgroup", cont=expgp)

expect_equal(names(expgp), "expandgroup")
visible(expgp) <- FALSE
expect_equal(visible(expgp), FALSE)

#########
## gpanedgroup
pg <- gpanedgroup(cont=g)
glabel("left", cont=pg)
glabel("right", cont=pg)

#########
## glayout
lyt <- glayout(cont=g)
lyt[1,1] <- "upper left"
lyt[2,2] <- gbutton("lower right", cont=lyt)

