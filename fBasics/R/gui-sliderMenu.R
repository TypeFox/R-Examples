
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


###############################################################################
# FUNCTION:                 SLIDER MENU:
#  .sliderMenu               Opens a teching demo slider menu
#  .tdSliderMenu             Opens a teching demo slider and button menu
###############################################################################


.slider.env = new.env()


# ------------------------------------------------------------------------------

.sliderMenu <- function(refresh.code, names, minima, maxima,
                        resolutions, starts,
                        title = "Slider", no = 0, set.no.value = 0)
{
  # A function implemented by Diethelm Wuertz
  
  # Description:
  #   Starts a slider menu
  
  # Source:
  #   Built on code written by Peter Wolf
  
  # FUNCTION:
  
  # Requirement:
  if (!require(tcltk, quietly = TRUE))
    stop("\n -- Package tcltk not available -- \n\n")
  
  # Environment:
  if (!exists(".slider.env")) {
    .slider.env <<- new.env()
  }
  if (no != 0) {
    options(show.error.messages = FALSE)
    ans <- as.numeric(tcltk::tclvalue(get(paste("slider", no, sep = ""),
                                          envir = .slider.env)))
    options(show.error.messages = TRUE)
    return(ans)
  }
  if (set.no.value[1] != 0) {
    try(eval(parse(text =
                     paste("tclvalue(slider", set.no.value[1],
                           ")<-", set.no.value[2], sep = "")),
             envir = .slider.env),
        silent = TRUE)
    return(set.no.value[2])
  }
  
  # Toplevel:
  nt = tcltk::tktoplevel()
  tcltk::tkwm.title(nt, title)
  
  
  # Slider:
  for (i in seq(names)) {
    eval(parse(text = paste("assign(\"slider", i, "\",
            tclVar(starts[i]), envir = .slider.env)", sep = "")))
    tcltk::tkpack(fr<-tcltk::tkframe(nt), anchor = "sw")
    lab = tcltk::tklabel(fr, text = names[i], anchor = "sw")
    sc = tcltk::tkscale(fr, command = refresh.code, from = minima[i],
                        to = maxima[i], showvalue = TRUE, resolution =
                          resolutions[i], orient = "horiz")
    assign("sc", sc, envir = .slider.env)
    tcltk::tkgrid(sc, lab)
    eval(parse(text = paste("tkconfigure(sc, variable = slider", i, ")",
                            sep = "")), envir = .slider.env)
  }
  tcltk::tkpack(fr<-tcltk::tkframe(nt), anchor = "sw")
  
  # Quit:
  quitButton = tcltk::tkbutton(fr, text = "   Quit   ",
                               command = function() {
                                 tcltk::tkdestroy(nt)
                               } )
  
  # Reset:
  resetButton = tcltk::tkbutton(fr, text = "   Start | Reset   ",
                                command = function() {
                                  for (i in seq(starts))
                                    eval(parse(text =
                                                 paste("tclvalue(slider", i, ")<-", starts[i], sep="")),
                                         envir = .slider.env)
                                  refresh.code()
                                }  )
  
  # Compose:
  tcltk::tkgrid(resetButton, quitButton, sticky = "sew")
}


# ------------------------------------------------------------------------------


.tdSliderMenu <-
  function(sl.functions, names, minima, maxima, resolutions, starts,
           but.functions, but.names, no, set.no.value, obj.name, obj.value,
           reset.function, title)
  {
    # A function implemented by Diethelm Wuertz
    
    # Description
    #   Opens a teching demo slider menu
    
    # Notes:
    #   Build on ideas and code from:
    #   R Package: TeachingDemos
    #   Title: Demonstrations for teaching and learning
    #   Version: 1.5
    #   Author: Greg Snow
    #   Description: This package is a set of demonstration functions
    #       that can be used in a classroom to demonstrate statistical
    #       concepts, or on your own to better understand the concepts
    #       or the programming.
    #   Maintainer: Greg Snow <greg.snow@intermountainmail.org>
    #   License: Artistic
    
    # FUNCTION:
    
    # Requirement:
    if (!require(tcltk, quietly = TRUE))
      stop("\n -- Package tcltk not available -- \n\n")
    
    # Setup:
    if(!missing(no)) {
      return(as.numeric(tcltk::tclvalue(get(paste(".tdSlider", no, sep=""),
                                            envir = .slider.env))))
    }
    if(!missing(set.no.value)){
      try(eval(parse(text=paste("tclvalue(.tdSlider", set.no.value[1],")<-",
                                set.no.value[2], sep = "")), envir = .slider.env))
      return(set.no.value[2])
    }
    if(!exists(".slider.env")) {
      .slider.env <<- new.env()
    }
    if(!missing(obj.name)){
      if(!missing(obj.value)) {
        assign(obj.name, obj.value, envir = .slider.env)
      } else {
        obj.value <- get(obj.name, envir = .slider.env)
      }
      return(obj.value)
    }
    if(missing(title)) {
      title = "Control Widget"
    }
    
    # GUI Settings:
    nt <- tcltk::tktoplevel()
    tcltk::tkwm.title(nt, title)
    tcltk::tkwm.geometry(nt, "+0+0")
    
    # Buttons:
    tcltk::tkpack(
      f.but <- tcltk::tkframe(nt), fill = "x")
    
    # Quit Button:
    quitCMD = function() {
      tcltk::tkdestroy(nt)
    }
    tcltk::tkpack(
      tcltk::tkbutton(f.but, text = "Quit", command = quitCMD, anchor = "sw"),
      side = "right",
      fill = "y")
    
    # Reset Button:
    if(missing(reset.function)) {
      reset.function <- function(...) print("relax")
    }
    if(!is.function(reset.function)) {
      reset.function<-eval(parse(text =
                                   paste("function(...){",reset.function,"}")))
    }
    resetCMD = function()
    {
      for(i in seq(names))
        eval(parse(text = paste("tclvalue(.tdSlider",i,")<-",
                                starts[i], sep = "")),
             envir = .slider.env)
      reset.function()
    }
    tcltk::tkpack(
      tcltk::tkbutton(f.but, text = "Reset", command = resetCMD, anchor = "sw"),
      side = "right",
      fill = "y")
    if (missing(but.names)) {
      but.names <- NULL
    }
    for (i in seq(but.names)) {
      but.fun <-
        if (length(but.functions) > 1)
          but.functions[[i]]
      else
        but.functions
      if (!is.function(but.fun)) {
        but.fun <-
          eval(parse(text = paste("function(...){", but.fun, "}")))
      }
      tcltk::tkpack(
        tcltk::tkbutton(f.but, text = but.names[i], command = but.fun,
                        anchor = "nw"),
        # side = "right",
        fill = "x"
      )
    }
    
    # Sliders:
    if(missing(names)) {
      names <- NULL
    }
    if(missing(sl.functions)) {
      sl.functions <- function(...){}
    }
    for(i in seq(names)){
      eval(parse(text = paste("assign('.tdSlider",i,"',
            tclVar(starts[i]), env = .slider.env)", sep = "")))
      tcltk::tkpack(fr <- tcltk::tkframe(nt))
      lab <- tcltk::tklabel(fr,
                            text = names[i],
                            anchor = "sw",
                            width = "35")
      sc <- tcltk::tkscale(fr,
                           from = minima[i],
                           to = maxima[i],
                           showvalue = TRUE,
                           resolution = resolutions[i],
                           orient = "horiz")
      tcltk::tkpack(lab,
                    sc,
                    anchor = "sw",
                    side = "right");
      assign("sc", sc, envir = .slider.env)
      
      eval(parse(text=paste("tkconfigure(sc,variable=.tdSlider",i,")",
                            sep="")), envir = .slider.env)
      sl.fun <-
        if(length(sl.functions)>1)
          sl.functions[[i]]
      else
        sl.functions
      if(!is.function(sl.fun))
        sl.fun<-eval(parse(text=paste("function(...){", sl.fun,"}")))
      
      tcltk::tkconfigure(sc, command = sl.fun)
    }
    assign("slider.values.old", starts, envir = .slider.env)
    
    # Return Value:
    invisible(nt)
  }


################################################################################

