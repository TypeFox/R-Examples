w <- gwindow("test gcheckboxgroup", visible=FALSE)
items <- letters[1:4]
cbg <- gcheckboxgroup(items, checked= items=="a", cont=w, horizontal=TRUE)

## svalue
print(svalue(cbg))
svalue(cbg) <- c(T,F,T,T)
print(svalue(cbg))


## [
print(cbg[])
cbg[] <- letters[12:15]                 # same length

## handler
addHandlerChanged(cbg, handler=function(h,...) {
  print(svalue(h$obj))
})


visible(w) <- TRUE
