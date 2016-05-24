w <- gwindow("statusbar")
sb <- gstatusbar("statusbar", cont=w)

## test
## svsalu
expect_equal(svalue(sb), "statusbar")

## svalue<--
svalue(sb) <- "new"
expect_equal(svalue(sb), "new")
