### R code from vignette source 'ex-tcltk-grid-layout.Rnw'

###################################################
### code chunk number 1: ex-tcltk-grid-layout.Rnw:7-34
###################################################
library(tcltk)
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
  !is.na(as.Date(sprintf("%s-%s-%s", year, month, day), format = "%Y-%m-%d"))

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


###################################################
### code chunk number 2: ex-tcltk-grid-layout.Rnw:47-64
###################################################
make_month <- function(parent, year, month) {
  ## add headers
  days <- c("S","M","T","W","Th","F","S")
  sapply(1:7, function(i) {
    label <- ttklabel(parent, text = days[i])       
    tkgrid(label, row = 0, column = i-1, sticky = "")
  })
  ## add days
  sapply(seq_len(ProgGUIinR:::days.in.month(year, month)), 
         function(day) {
           label <- ttklabel(parent, text = day)
           row <- ProgGUIinR:::week.of.month(year, month, day)
           col <- ProgGUIinR:::day.of.week(year, month, day)
           tkgrid(label, row = 1 + row, column = col,
                  sticky = "e")
         })
}


###################################################
### code chunk number 3: ex-tcltk-grid-layout.Rnw:69-70
###################################################
year <- 2000; month <- 1


###################################################
### code chunk number 4: ex-tcltk-grid-layout.Rnw:75-82
###################################################
window <- tktoplevel()
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both", side = "top")
c_frame <- ttkframe(frame)
cal_frame <- ttkframe(frame)
tkpack(c_frame, fill = "x", side = "top")
tkpack(cal_frame, expand = TRUE, anchor = "n")


###################################################
### code chunk number 5: ex-tcltk-grid-layout.Rnw:88-96
###################################################
previous_button <- ttklabel(c_frame, text = "<")
next_button <- ttklabel(c_frame, text = ">")
current_month <- ttklabel(c_frame)
#
tkpack(previous_button, side = "left", anchor = "w")
tkpack(current_month, side = "left", anchor = "center", 
       expand = TRUE)
tkpack(next_button, side = "left", anchor = "e")


###################################################
### code chunk number 6: stackedWidget
###################################################
set_month <- function() {
  tkpack("forget", cal_frame)
  cal_frame <<- ttkframe(frame)
  make_month(cal_frame, year, month)
  tkconfigure(current_month,              # month label
              text = sprintf("%s %s", month.abb[month], year))
  tkpack(cal_frame)
}
set_month()                              # initial calendar


###################################################
### code chunk number 7: connectSignal
###################################################
tkbind(previous_button, "<Button-1>", function() {
  if(month > 1) {
    month <<- month - 1
  } else {
    month <<- 12; year <<- year - 1
  }
  set_month()
})


###################################################
### code chunk number 8: ex-tcltk-grid-layout.Rnw:132-140
###################################################
tkbind(next_button, "<Button-1>", function() {
 if(month < 12) {
    month <<- month + 1
  } else {
    month <<- 1; year <<- year + 1
  }
  set_month()
})


###################################################
### code chunk number 9: ex-tcltk-grid-layout.Rnw:146-151
###################################################
tkbind("TLabel", "<Button-1>", function(W) {
  day <- as.numeric(tkcget(W, "-text"))
  if(!is.na(day))
    print(sprintf("You selected: %s/%s/%s", month, day, year))
})


