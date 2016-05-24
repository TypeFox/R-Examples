w <- gwindow("notebook ")


nb <- gnotebook(cont=w)

## add pages
mapply(gbutton, state.name[1:4], cont=list(nb), label=paste("tab", 1:4))

## switchpages
svalue(nb) <- 2
expect_equal(svalue(nb), 2)

## names
expect_equal(names(nb), paste("tab", 1:4))

## names<_
names(nb) <- toupper(state.name[1:4])
expect_equal(names(nb), toupper(state.name[1:4]))

## length
expect_equal(length(nb), 4)

## dispose
svalue(nb) <- 4
dispose(nb)
expect_equal(length(nb), 3)

## hcange handler
ctr <- 1
addHandlerChanged(nb, handler=function(h,...) ctr <<- 2)
svalue(nb) <- 1 ## invoke handler
expect_equal(ctr, 2)


