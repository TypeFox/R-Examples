
w <- gwindow("gtext test")
g <- gvbox(cont=w)
lorem <- "
Lorem ipsum dolor sit amet, consectetur adipisicing elit,
sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi
ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit
in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur
sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit
anim id est laborum.
"
txt1 <- gtext(lorem, cont=g)


## replace text
svalue(txt1) <- "new"
expect_equal(svalue(txt1), "new")

## add text
insert(txt1, "more", do.newline=FALSE)
expect_equal(svalue(txt1), "newmore")

insert(txt1, "more", do.newline=TRUE)
expect_equal(svalue(txt1), "newmoremore\n")

## eispose text
dispose(txt1)
expect_equal(svalue(txt1), "")

## fonts
insert(txt1, "blue text", font.attr=list(weight="bold",  color="blue"))

## font for a butter
txt2 <- gtext("monospace font for entire buffer", cont=g,
              font.attr=list(family="monospace"))


## handlers
addHandlerSelectionChanged(txt2, function(h,...) {
  print("You have selected:")
  print(svalue(h$obj, drop=TRUE))
})

addHandlerBlur(txt2, function(h,...) print("focus out"))
addHandlerFocus(txt2, function(h,...) print("focus in"))
addHandlerChanged(txt2, function(h,...) print(h$key))
