## stuff to make generic widgets


## some predefined types for GenericWidget
default.color.list = list(
  type="gdroplist",
  items=c("","\"black\"","\"blue\"","\"red\"","\"yellow\"","\"brown\"","\"green\"","\"pink\"","NULL"),
  editable=TRUE
  )
lty.list = list(
  type = "gdroplist",
  items = c("\"solid\"","\"dashed\"","\"dotted\"","\"dotdash\"","\"longdash\"","\"twodash\"","\"blank\""),
  editable=TRUE
  )
pch.list = list(
  type = "gspinbutton",
  from=0,
  to=26,
  by=1,
  value=1
  )
EMPTY.list = list(
  type = "gedit",
  text = ""
  )
BLANK.list = list(                      # for putting in a space
  type = "glabel",
  text = ""
  )
NULL.list = list(
  type = "gedit",
  text = "NULL"
  )
FALSE.list = list(
  type = "gradio",
  items = c("TRUE","FALSE"),
  index = FALSE,
  selected = 2
  )
TRUE.list = list(
  type = "gradio",
  index = FALSE,
  items = c("TRUE","FALSE")
  )
emptyTRUE.list = list(
  type = "gdroplist",
  items = c("","TRUE","FALSE")
  )
alternative.list = list(
  type="gdroplist",
  items=c("\"two.sided\"","\"less\"","\"greater\"")
  )
conf.level.list = list(
  type = "gdroplist",
  items = c(0.95, 0.99, 0.90, 0.80),
  editable = TRUE
  )
labels.list = list(
    main = EMPTY.list,
    sub = BLANK.list,
    xlab = EMPTY.list,
    ylab = EMPTY.list
  )


## trick to assign in global env. Not idea!!!
assign_global <- function(nm, obj) {
    envir <- eval(parse(text=".GlobalEnv"))
    assign(nm, obj, envir=envir)
}
