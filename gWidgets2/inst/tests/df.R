w <- gwindow("gdf", visible=FALSE)
g <- ggroup(cont=w, horizontal=FALSE)

x <- mtcars

##
DF <- gdf(x, cont=g)
visible(w) <- TRUE


## test
## [
expect_equal(DF[1,1], x[1,1])
#expect_equal(DF[1,], x[1,])
expect_equal(DF[,1], x[,1])
expect_equal(DF[,], x[,])

## [<-
DF[1,1] <- 22
expect_equal(DF[1,1], 22)

## length
expect_equal(as.vector(length(DF)), length(x))
##dim
expect_equal(as.vector(dim(DF)), dim(x))
## names
expect_equal(names(DF), names(x))



DF[,1] <- 2*mtcars[,1]
expect_equal(DF[,1], 2*mtcars[,1])

# replace
DF[] <- head(mtcars)
expect_equal(as.vector(dim(DF)), dim(head(mtcars)))
