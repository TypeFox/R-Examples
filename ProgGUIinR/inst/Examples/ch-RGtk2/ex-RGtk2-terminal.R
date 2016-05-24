### R code from vignette source 'ex-RGtk2-terminal.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-terminal.Rnw:13-15
###################################################
## make a *basic* terminal in RGtk2
library(RGtk2)


###################################################
### code chunk number 2: TextViewWidget
###################################################
view <- gtkTextView()
buffer <- view$getBuffer()
font <- pangoFontDescriptionFromString("Monospace")
view$modifyFont(font)                     # widget wide


###################################################
### code chunk number 3: ex-RGtk2-terminal.Rnw:30-36
###################################################
buffer$createTag(tag.name = "cmdInput")
buffer$createTag(tag.name = "cmdOutput", 
                 weight = PangoWeight["bold"])
buffer$createTag(tag.name = "cmdError", 
       weight = PangoStyle["italic"], foreground = "red")
buffer$createTag(tag.name = "uneditable", editable = FALSE)


###################################################
### code chunk number 4: ex-RGtk2-terminal.Rnw:41-46
###################################################
start_cmd <- buffer$createMark("start_cmd", 
                              buffer$getStartIter()$iter, 
                              left.gravity = TRUE)
bufferEnd <- buffer$createMark("bufferEnd", 
                               buffer$getEndIter()$iter)


###################################################
### code chunk number 5: ex-RGtk2-terminal.Rnw:53-67
###################################################
add_prompt <- function(obj, prompt = c("prompt", "continue"),
                      set_mark = TRUE) 
{
  prompt <- match.arg(prompt)
  prompt <- getOption(prompt)
  
  end_iter <- obj$getEndIter()
  obj$insert(end_iter$iter, prompt)
  if(set_mark)
    obj$moveMarkByName("start_cmd", end_iter$iter)
  obj$applyTagByName("uneditable", obj$getStartIter()$iter, 
                     end_iter$iter)
}
add_prompt(buffer) ## place an initial prompt


###################################################
### code chunk number 6: add_ouput
###################################################
add_ouput <- function(obj, output, tag_name = "cmdOutput") {
  end_iter <- obj$getEndIter()
  if(length(output) > 0)  
    sapply(output, function(i)  {
      obj$insertWithTagsByName(end_iter$iter, i, tag_name)
      obj$insert(end_iter$iter, "\n", len=-1)
    })
}


###################################################
### code chunk number 7: ex-RGtk2-terminal.Rnw:90-98
###################################################
find_cmd <- function(obj) {
  end_iter <- obj$getEndIter()
  start_iter <- obj$getIterAtMark(start_cmd)
  cmd <- obj$getText(start_iter$iter, end_iter$iter, TRUE)
  regex <- paste("\n[", getOption("continue"), "] ", sep = "")
  cmd <- unlist(strsplit(cmd, regex))
  cmd
}


###################################################
### code chunk number 8: evalCmd
###################################################
require(evaluate)
eval_cmd <- function(view, cmd) {
  buffer <- view$getBuffer()
  out <- try(evaluate:::evaluate(cmd, .GlobalEnv), 
             silent = TRUE)

  if(inherits(out, "try-error")) {
    ## parse error
    add_ouput(buffer, out, "cmdError")
  } else if(inherits(out[[2]], "error")) {
    if(grepl("end", out[[2]])) {        # a hack here
      add_prompt(buffer, "continue", set_mark = FALSE)
      return()
    } else {
      add_ouput(buffer, out[[2]]$message, "cmdError")
    }
  } else {
    add_ouput(buffer, out[[2]], "cmdOutput")
  }
  add_prompt(buffer, "prompt", set_mark = TRUE)
}


###################################################
### code chunk number 9: connectBinding
###################################################
gSignalConnect(view, "key-release-event", 
               f=function(view, event) {
                 buffer <- view$getBuffer()
                 keyval <- event$getKeyval()
                 if(keyval == GDK_Return) {
                   cmd <- find_cmd(buffer)
                   if(length(cmd) && nchar(cmd) > 0)
                     eval_cmd(view, cmd)
                 }
               })


###################################################
### code chunk number 10: ex-RGtk2-terminal.Rnw:154-160
###################################################
scroll_viewport <- function(view, ...) {
  view$scrollToMark(bufferEnd, within.margin = 0)
  return(FALSE)
}
gSignalConnect(buffer, "changed", scroll_viewport, data=view, 
               after = TRUE, user.data.first = TRUE)


###################################################
### code chunk number 11: makeGUI
###################################################
## scroll window
sw <- gtkScrolledWindow()
sw$setPolicy("automatic", "automatic")
sw$add(view)

## top-level window
w <- gtkWindow(show=FALSE)
w$setTitle("A terminal")
w$add(sw)
w$setSizeRequest(400,200)
w$showAll()


###################################################
### code chunk number 12: ex-RGtk2-terminal.Rnw:179-251
###################################################
## History features
## This is not illustrated in text, but is added here to illustrate how this might be implemented
## The major issue with this example is we can't trap the return or arrow keys before they move 
## the cursor so any thing ends up looking jerky

## store the stack and a pointer to the current command with the text buffer
buffer$setData("history", list())
buffer$setData("ptr", 0)


## replace cmd with that in str.
replace_cmd <- function(obj, str) {
  end_iter <- obj$getEndIter()
  start_iter <- obj$getIterAtMark(start_cmd)
  obj$delete(start_iter$iter, end_iter$iter)
  end_iter <- obj$getEndIter()
  obj$insertWithTagsByName(end_iter$iter, str[1], "cmdInput")
  if(length(str) > 1) {
    for(i in str[-1]) {
      obj$insert(end_iter$iter, "\n")
      obj$insertWithTagsByName(end_iter$iter, getOption("continue"), "cmdInput")
      obj$insertWithTagsByName(end_iter$iter, i, "cmdInput")
    }
  }
  moveViewport(obj)
}

## This adds the command to the history stack and moves the pointer.
add_history <- function(obj, cmd) {
  history <- obj$GetData("history"); ptr <- obj$GetData("ptr")
  history <- c(history, cmd)
  ptr <- length(history)
  obj$SetData("ptr", ptr)
  obj$SetData("history", history)
}

## these next two functions scroll through the history
scroll_history_up <- function(obj) {
  ## move through history
  ptr <- obj$GetData("ptr") - 1
  if(ptr > 0)
    replace_cmd(obj, obj$GetData("history")[[ptr]])
  obj$SetData("ptr", max(ptr,0))
  obj$PlaceCursor(obj$GetEndIter()$iter)
}

scroll_history_down <- function(obj) {
  ## move through history
  ptr <- obj$GetData("ptr") + 1
  history <- obj$GetData("history")
  if(ptr <= length(history)) 
    replace_cmd(obj, history[[ptr]])
  obj$SetData("ptr", min(ptr,length(history)))
  obj$PlaceCursor(obj$GetEndIter()$iter)
}

## History bindings
## this uses Control-p and Control-n to move
ID <- gSignalConnect(view, "key-release-event", f=function(w, e, data) {
  if(e$GetState() != GdkModifierType['control-mask'])
    return(TRUE)

  obj <- w$GetBuffer()
  keyval <- e$GetKeyval()

  if(keyval == GDK_p) {
    scroll_history_up(obj)
  } else if(keyval == GDK_n) {
    scroll_history_down(obj)
  }
  return(TRUE)
})


