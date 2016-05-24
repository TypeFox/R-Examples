library(gWidgetsRGtk2)

# Example A: a radio button and dialog
CA <- gwindow(title = "Happy or not", visible = TRUE, guiToolkit = "RGtk2")
CB <- ggroup(container = CA, horizontal = FALSE)
wg.title <- glabel(text = "Tell me how you feel ...", container = CB)
wg.button <- gradio(items = c("Haaappy", "Saaad", "Fiiine"), container=CB,
  handler = function(h, ...) {
    gmessage(message = paste("So you are", svalue(h$obj), "!"), 
      icon = "info")  # h$obj refers to wg.button
  }
)    

# Example B: personal income tax
guiTax <- function() { 
  rMain <- gwindow(title = "A taxation example", visible = FALSE,
    width = 700, height = 200, guiToolkit = "RGtk2")
  r0 <- ggroup(container = rMain, horizontal = FALSE)
  wg.tit <- glabel(text = "Calculating Net Monthly Income", container = r0)
  font(wg.tit) <- list(size = 15, weight = "bold", color = "gray30") 
  
  r1 <- ggroup(container = r0, horizontal = TRUE)
  rL <- gframe(text = "Inputs", container = r1, horizontal = FALSE)
  addSpace(obj = r1, value = 20)
  rR <- gframe(text = "Outputs", container = r1)
  rLa <- glayout(container = rL); rLb <- ggroup(container = rL)
  rRa <- glayout(container = rR)
  
  rLa[1, 1] <- glabel(text = "Gross annual income ($)", container = rLa)
  rLa[1, 2] <- ia <- gedit(initial.msg = "a number", container = rLa)
  rLa[2, 1] <- glabel(text = "Number of dependents", container = rLa)
  rLa[2, 2] <- ib <- gspinbutton(from = 0, to = 15, by = 1, container=rLa)
  rRa[1, 1] <- glabel(text = "Monthly taxes ($)", container = rRa)
  rRa[1, 2] <- oa <- gedit(container = rRa)
  rRa[2, 1] <- glabel(text = "Monthly net income ($)", container = rRa)
  rRa[2, 2] <- ob <- gedit(container = rRa)
  
  addSpring(rLb)
  wd.submit <- gbutton(text = "Submit", container = rLb)
  addHandlerClicked(obj = wd.submit, handler = function(h, ...) {
      inc <- as.numeric(svalue(ia))  # input ia from line 30
      nd <- as.numeric(svalue(ib))   # input ib from line 32
      agi <- inc - 5000 * nd         # transformation
      ag <- ifelse(agi < 0, yes = 0, no = agi)      
      tax <- ifelse(test = ag <= 10000, yes = ag * 0.1, no = 
        ifelse(test = ag <= 50000, yes = 1000 + (ag - 10000) * 0.2,
          no = 9000 + (ag - 50000) * 0.3)) 
      m.tax <- tax / 12
      m.inc <- (inc - tax) / 12  
      svalue(oa) <- sprintf(fmt = "%.2f", m.tax)  # output oa to line 34
      svalue(ob) <- sprintf(fmt = "%.2f", m.inc)  # output ob to line 36
    }
  )
  visible(rMain) <- TRUE
  invisible(NULL) 
}
guiTax()