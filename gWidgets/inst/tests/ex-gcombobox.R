w <- gwindow("combobox example", visible=FALSE)
g <- ggroup(cont=w, horizontal = FALSE)

m = data.frame(labels = letters[1:3],
  icons = c("quit","open","file"),
  tips = paste("the letter", letters[1:3])
  )

## vector argument
cb1 = gcombobox(m[,1, drop=TRUE], selected=1, cont =g)

# svalue
print(svalue(cb1))
print(svalue(cb1, index=TRUE))

#savlue<-
svalue(cb1) <- "b"
svalue(cb1, index=TRUE) <- 3

# [
print(cb1[])

# [<-
cb1[] <- toupper(letters[1:3])

## handler
addHandlerChanged(cb1, handler = function(h,...) print(svalue(h$obj)))

##### Width
cb1.5 <- gcombobox(m, width=100, cont = g)


##### editable
cb2 <- gcombobox(m[,1,drop=TRUE], editable = TRUE, cont = g)

## svalue<-
svalue(cb2) <- "editable"

## svalue
print(svalue(cb2))

## icons? toolkit specific, but should handle this gracefully
cb3 <- gcombobox(m, cont = g)


visible(w) <- TRUE
