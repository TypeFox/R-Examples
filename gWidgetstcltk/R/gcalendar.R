## add calendar widget: shoule I have gcalendar, gcalendarbrowser?
## no handler function, can add to entry object with addhandler

setClass("gCalendartcltk",
         representation = representation("gComponenttcltk",
            format="character"),
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )


setMethod(".gcalendar",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   text="",
                   format="%Y-%m-%d",
                   handler = NULL, action=NULL,
                   container=NULL,...) {

            force(toolkit)

            theArgs <- list(...)

            if(format != "%Y-%m-%d") {
              message(gettext("The format argument is not employed. Pass in coercion function through the coerce.with argument if year-month-day is not desired."), "\n")
              format <- "%Y-%m-%d"
            }

            ## No initial date
            ##            if(text == "" && format != "")
            ##              text <- format(Sys.Date(), format)


            g <- ggroup(container=container, ...)
            e <- gedit(text, container=g, width=11, expand=TRUE)
            b <- gbutton("date", container=g)


            obj <- new("gCalendartcltk",
                       block= getBlock(g),
                       widget = getWidget(e),
                       format=format,
                       toolkit=toolkit, ID=getNewID(), e = new.env())
            
            
            addHandlerClicked(b, action=obj, handler=function(h,...) {
              text <- svalue(obj)
              year <- as.numeric(format(as.Date(text, tag(obj, "format")), format="%Y"))
              month <- as.numeric(format(as.Date(text, tag(obj, "format")), format="%m"))
              makeCalendar(obj, year, month)
            })

            if(!is.null(theArgs$coerce.with))
              coerce.with <- theArgs$coerce.with
            else
              coerce.with <- function(x, ...) {
                as.Date(x, format=format)
              }
            theArgs <- list(...)
            tag(obj, "..entry") <- e
            tag(obj,"format") <- format
            tag(obj,"coerce.with") <- coerce.with

            return(obj)          # drop down to tcltk widget
          })

setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gCalendartcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            svalue(tag(obj, "..entry"))
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gCalendartcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   widget <-tag(obj, "..entry")
                   svalue(widget) <- value
                 })

## helper
makeCalendar <- function(widget, year, month) {

  toplevel <- tktoplevel()
  f <- ttkframe(toplevel, padding=c(3,3,12,12))
  tkpack(f, expand=TRUE, fill="both", side="top")
  cframe <- ttkframe(f)
  calframe <- ttkframe(f)
  tkpack(cframe, fill="x", side="top")
  tkpack(calframe, expand=TRUE, anchor="n")


  year <- year; month <- month          # function local

  
  ##' from chron with slight change to arguments
  day.of.week <- function (year, month, day) {
    ix <- year + trunc((month - 14)/12)
    jx <- (trunc((13 * (month + 10 - (month + 10)%/%13 * 12) - 
                  1)/5) + day + 77 + (5 * (ix - (ix%/%100) * 100))%/%4 + 
           ix%/%400 - (ix%/%100) * 2)
    jx%%7
  }
  
  
  ## is this a valid date
  validDate <- function(year, month, day) 
    !is.na(as.Date(sprintf("%s-%s-%s", year, month, day), format="%Y-%m-%d"))
  
  ## how many days in a month
  days.in.month <- function(year, month) {
    for(i in c(31, 30, 29, 28)) {
      if(validDate(year, month, i))
        return(i)
    }
  }
  ## 0-based week of month
  week.of.month <- function(year, month, day) {
    first.day <- day.of.week(year, month, 1)
    (first.day + day - 1) %/% 7
  }
  
  makeMonth <- function(w, year, month) {
    ## add headers
    days <- c("S","M","T","W","Th","F","S")
    lapply(1:7, function(i) {
      l <- ttklabel(w, text=days[i])           # color
      tkgrid(l, row=0, column=i-1, sticky="")
    })
    ## add days
    lapply(1:days.in.month(year, month),  function(day) {
      l <- ttklabel(w, text=day)

      ## bind to each day
      ## might be more efficient to bind to toplevel and intercept
      tkbind(l, "<Button-1>", function(W) {
        day <- tclvalue(tkcget(W,"-text"))        
        svalue(widget) <- sprintf("%s-%s-%s", year, month, day)
        tkdestroy(toplevel)
      })

      
      tkgrid(l, row=1 + week.of.month(year, month, day),
             column=day.of.week(year, month, day),
             sticky="e")
    })
  }

  ## controls
  prevb <- ttklabel(cframe, text="<")
  nextb <- ttklabel(cframe, text=">")
  curmo <- ttklabel(cframe)
  
  tkpack(prevb, side="left", anchor="w")
  tkpack(curmo, side="left", anchor="center", expand=TRUE)
  tkpack(nextb, side="left", anchor="e")

  
  setMonth <- function() {
    tkpack("forget", calframe)
    calframe <<- ttkframe(f); tkpack(calframe)
    makeMonth(calframe, year, month)
    tkconfigure(curmo, text=sprintf("%s %s", month.abb[month], year))
  }
  setMonth()                              # initial calendar
  
  tkbind(prevb, "<Button-1>", function() {
    if(month > 1) {
      month <<- month - 1
    } else {
      month <<- 12; year <<- year - 1
    }
    setMonth()
  })
  
  tkbind(nextb, "<Button-1>", function() {
    if(month < 12) {
      month <<- month + 1
  } else {
    month <<- 1; year <<- year + 1
  }
    setMonth()
  })
  
}
