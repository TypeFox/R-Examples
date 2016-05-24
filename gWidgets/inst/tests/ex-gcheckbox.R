w <- gwindow("checkboxes", visible=FALSE)
g <- ggroup(cont = w, horizontal = FALSE)

## checkbox
cb <- gcheckbox("label", cont = g)

# svalue
print(svalue(cb))

# svalue<-
svalue(cb) <- FALSE


# [ (names)
print(cb[])

# [ (names<-)
cb[1] <- "new label"

# handler
addHandlerChanged(cb, function(h,...) print("clicked"))

## radio

r <- gradio(letters[1:3], cont = g, horizontal=TRUE)

# svalue
print(svalue(r))
print(svalue(r, index=TRUE))


#svalue<-
svalue(r) <- "b"
svalue(r,index = TRUE) <- 3

# [ -- names
print(r[])

# [<-
r[1] <- "A"
r[] <- c("A","B","C")

# handler
addHandlerChanged(r, handler = function(h,...) print("clicked radio"))

## gcheckboxgroup

cbg <- gcheckboxgroup(letters[1:3], cont = g)

## empty?
print(svalue(cbg))

## svalue<-
svalue(cbg) <- c(TRUE,TRUE, FALSE)

## svalue
print(svalue(cbg))
print(svalue(cbg, index=TRUE))

## names
print(cbg[])

## [<-
cbg[1] <- "A"
cbg[] <- c("A","B","C")

# handler
addHandlerChanged(cbg, handler = function(h,...) print(svalue(h$obj)))



visible(w) <- TRUE
