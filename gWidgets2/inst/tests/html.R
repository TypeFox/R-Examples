w <- gwindow("ghtml")
g <- gnotebook(cont=w, horizontal=FALSE)

h1 <- ghtml("no markup", cont=g, fill="x", label="no markup")
h2 <- ghtml("<i>emphasize</i>", cont=g, label="markup")
## from url
f <- "http://www.yahoo.com"
h3 <- ghtml(f, cont=g, width=300, label="url")

## test
expect_equal(svalue(h1), "no markup")


svalue(h1) <- "new text"
expect_equal(svalue(h1), "new text")

svalue(h3) <- (x <- "http://www.google.com")
expect_equal(svalue(h3), x)
               

