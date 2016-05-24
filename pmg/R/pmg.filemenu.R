
exit.pmg = function(w) {
  ## use pmg.window global by name
  pmg.window$Hide()
  pmg.window$Destroy()
  assignInNamespace("pmg.window",NULL, "pmg")
}


## genericWidget lists
read.table.list = list(
  title = "read.table()",
  help = "read.table",
  action = list(
    beginning = "read.table(",
    ending = ")"
    ),
  type = "text",                        #either text or graphic
  variableType = "fileurl",
  assignto = TRUE,
  arguments = list(
    arguments = list(
      header=FALSE.list,
      sep = list(
        type = "gedit",
        text = "\"\""
        ),
#      quote = list(
#        type = "gedit",
#        text = "\"'"
#        ),
      dec = list(
        type = "gedit",
        text = "\".\""
        ),
      skip = list(
        type = "gedit",
        text = 0
        ),
      check.names = TRUE.list,
      comment.char = list(
        type = "gedit",
        text = "\"#\""
        )
      )
    )
  )


## fwf
read.fwf.list = list(
  title = "read.fwf()",
  help = "read.fwf",
  action = list(
    beginning = "read.fwf(",
    ending = ")"
    ),
  type = "text",                        #either text or graphic
  variableType = "fileurl",
  assignto = TRUE,
  arguments = list(
    arguments = list(
      widths = list(
        type="gedit",
        text = ""
        ),
      header=FALSE.list,
      sep = list(
        type = "gedit",
        text = "\"\""
        ),
      as.is = FALSE.list,
      skip = list(
        type = "gedit",
        text = 0
        )
      )
    )
  )


## csv
read.csv.list = list(
  title = "read.csv()",
  help = "read.csv",
  action = list(
    beginning = "read.csv(",
    ending = ")"
    ),
  type = "text",                        #either text or graphic
  variableType = "fileurl",
  assignto = TRUE,
  arguments = list(
    arguments = list(
      header=FALSE.list,
      sep = list(
        type = "gedit",
        text = "\"\""
        ),
#      quote = list(
#        type = "gedit",
#        text = "\"'"
#        ),
      dec = list(
        type = "gedit",
        text = "\".\""
        )
      )
    )
  )

