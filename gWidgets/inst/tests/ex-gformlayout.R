w <- gwindow("gformlayout test", visible=FALSE)

## layout a collection of widgets to generate a t.test
tTest <- list(type = "ggroup",
              horizontal = FALSE,
              children = list(
                list(type="fieldset",
                     columns = 2,
                     label = "Variable(s)",
                     label.pos = "top",
                     label.font = c(weight="bold"),
                     children = list(
                       list(name = "x",
                            label = "x",
                            type = "gedit",
                            text = ""),
                       list(name = "y",
                            label = "y",
                            type = "gedit",
                            text = "",
                            depends.on = "x",
                            depends.FUN = function(value) nchar(as.character(value)) > 0,
                            depends.signal = "addHandlerBlur"
                            )
                       )
                     ),
                list(type = "fieldset",
                     label = "Hypotheses",
                     columns = 2, 
                     children = list(
                       list(name = "mu",
                            type = "gedit",                            
                            label = "Ho: mu=",
                            text = "0",
                            coerce.with = as.numeric),
                       list(name = "alternative",
                            type="gcombobox",
                            label = "HA: ",
                            items = c("two.sided","less","greater")
                            )
                       )
                     ),
                list(type = "fieldset",
                     label = "two sample test",
                     columns = 2,
                     depends.on = "y",
                     depends.FUN = function(value) nchar(as.character(value)) > 0,
                     depends.signal = "addHandlerBlur",                     
                     children = list(
                       list(name = "paired",
                            label = "paired samples",
                            type = "gcombobox",
                            items = c(FALSE, TRUE)
                            ),
                       list(name = "var.equal",
                            label = "assume equal var",
                            type = "gcombobox",
                            items = c(FALSE, TRUE)
                            )
                       )
                     ),
                list(type = "fieldset",
                     columns = 1,
                     children = list(
                       list(name = "conf.level",
                            label = "confidence level",
                            type = "gedit",
                            text = "0.95",
                            coerce.with = as.numeric)
                       )
                     )
                )
              )


g <- ggroup(horizontal = FALSE, cont = w)
fl <- gformlayout(tTest, cont = g, expand=TRUE)
bg <- ggroup(cont = g)
addSpring(bg)
b <- gbutton("run t.test", cont = bg)
addHandlerChanged(b, function(h,...) {
  out <- svalue(fl)
  out$x <- svalue(out$x) # turn text string into numbers via get()
  if(out$y == "") {
    out$y <- out$paired <- NULL 
  } else {
   out$y <- svalue(out$y)
  }
  print(do.call("t.test",out))
})

visible(w) <- TRUE
