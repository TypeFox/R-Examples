guiApt <- function(pos = 1)
{  
  # -----------------------------------------------------------------------
  # A. Position and size  of containers
  if(!any(pos == 1:4)) stop("'pos' needs to 1, 2, 3, or 4.\n")     
  width <- ifelse(test = pos == 1 || pos == 3, yes = 650, no = 780)

  # -----------------------------------------------------------------------
  # B. Containers: 7 notebook tabs
  rMain <- gwindow(title = "Asymmetric Price Transmission (APT)", 
    visible = FALSE, guiToolkit = "RGtk2", width = width, height = 750)
  rNB <- gnotebook(container = rMain, tab.pos = pos)
  tn <- c("Home", "1. Data", "2. ciTarFit", "3. ciTarThd", 
          "4. ciTarLag", "5. ecmSymFit", "6. ecmAsyFit")
  r0 <- ggroup(container = rNB, label = tn[1], use.scrollwindow = TRUE,
               horizontal = FALSE)
  r1 <- ggroup(container = rNB, label = tn[2], use.scrollwindow = TRUE,
               horizontal = FALSE)
  r2 <- ggroup(container = rNB, label = tn[3], use.scrollwindow = TRUE,
               horizontal = FALSE)
  r3 <- ggroup(container = rNB, label = tn[4], use.scrollwindow = TRUE,
               horizontal = FALSE)
  r4 <- ggroup(container = rNB, label = tn[5], use.scrollwindow = TRUE,
               horizontal = FALSE)
  r5 <- ggroup(container = rNB, label = tn[6], use.scrollwindow = TRUE,
               horizontal = FALSE)
  r6 <- ggroup(container = rNB, label = tn[7], use.scrollwindow = TRUE,
               horizontal = FALSE)
  
  # -----------------------------------------------------------------------
  # C. Menu and toolbar
  hC <- function(h, ...) {dispose(rMain)}
  h0 <- function(h, ...) {svalue(rNB) <- match(x=tn[1], table=names(rNB))}
  h1 <- function(h, ...) {svalue(rNB) <- match(x=tn[2], table=names(rNB))}
  h2 <- function(h, ...) {svalue(rNB) <- match(x=tn[3], table=names(rNB))}
  h3 <- function(h, ...) {svalue(rNB) <- match(x=tn[4], table=names(rNB))}
  h4 <- function(h, ...) {svalue(rNB) <- match(x=tn[5], table=names(rNB))}
  h5 <- function(h, ...) {svalue(rNB) <- match(x=tn[6], table=names(rNB))}
  h6 <- function(h, ...) {svalue(rNB) <- match(x=tn[7], table=names(rNB))}
  h7 <- function(h, ...) {
    gmessage(title = "About this GUI", message = paste(
      "This GUI is developed for demonstration purpose. The apt ",
      "package is chosen for showing how to convert functions ",
      "into an intuitive interface for users. This software ",
      "program can be distributed freely for noncommercial use.\n", 
      "Copyright 2014 by C. Sun", sep = "\n"),     
      icon = "info", parent = rMain)
  }
  
  aSepH <- gseparator(horizontal = TRUE)
  aSepV <- gseparator(horizontal = FALSE)
  aClos <- gaction(handler = hC, label = "Close", icon = "close", 
                   tooltip = "Shut down the whole GUI")
  aHome <- gaction(handler = h0, label = "Home", icon = "home", 
                   tooltip = "Go to tab 'Home'")
  aDat  <- gaction(handler = h1, label = "1 Data property")
  aFit  <- gaction(handler = h2, label = "2 Model Fit", icon = "apply", 
                   tooltip = "Go to tab '2. ciTarFit'")
  aThr  <- gaction(handler = h3, label = "3 Threshold choice")
  aLag  <- gaction(handler = h4, label = "4 Lag choice")
  aSym  <- gaction(handler = h5, label = "5 Symmetric")
  aAsy  <- gaction(handler = h6, label = "6 Asymmetric", icon = "convert", 
                   tooltip = "Go to tab '6. ecmAsymFit'")
  aHel  <- gaction(handler = h0, label = "Welcome page")
  aApt  <- gaction(handler = h7, label = "About")
  
  listMenu <- list(
    File = list(a = aHome, b = aClos), Data = list(a = aDat), 
    Threshold = list(a = aFit, sep = aSepH, b = aThr, c = aLag),
    ECM = list(a = aSym, b = aAsy), 
    Help = list(a = aHel, b = aApt))        
  listTool <- list(a = aHome, sep = aSepV, 
    c = aFit, sep = aSepV, d = aAsy,  sep = aSepV, b = aClos)
  
  gmenu(container = rMain, menulist = listMenu)
  gtoolbar(container = rMain, toolbarlist = listTool, style = "icons")
  
  # -----------------------------------------------------------------------
  # D. Controls for each notebook tab
  
  # Tab 0 = Home ----------------------------------------------------------
  addSpace(obj = r0, value = 20, horizontal = FALSE)
  t1 <- glabel(text = "Welcome to GUI for the apt Package", container = r0)
  font(t1) <- list(weight = "bold", size = 17, color = "purple")
  
  addSpace(obj = r0, value = 15, horizontal = FALSE)
  glabel(container = r0, text = paste(
    "This graphical user interface is designed to help users understand", 
    "the main functionalities of the apt package. The tasks are organized",
    "into six main tabs. All outputs are delivered to the R console.\n",
    "To get started, see help('daVich') and help('guiApt').",
    "For example, run the following codes in R to generate data.",
    "Then use 'prVi' and 'prCh' as inputs for the fuction calls.\n", 
    sep = "\n")
  )
  sample.code <- glabel(container = r0, 
    text = paste(
      "data(daVich)                                ",
      "prVi <- daVich[, 1]", 
      "prCh <- daVich[, 2]", sep = "\n")
  ) 
  font(sample.code) <- list(family = "monospace", color = "blue")
   
  addSpring(obj = r0, horizontal = FALSE)
  glabel(container = r0, text = "Enjoy your play!")
  addSpace(obj = r0, value = 10, horizontal = FALSE)
  glabel(container = r0, text = "Copyright 2014 by C. Sun")
  addSpace(obj = r0, value = 10, horizontal = FALSE)
  
  # Tab 1 = Data  ---------------------------------------------------------
  # Specify data names
  addSpace(obj = r1, value = 10, horizontal = FALSE)
  text11 <- glabel("Supplying two single time series data", container = r1)
  r1a <- gframe(text = "", container = r1)
  glabel(container = r1a, text = "y = ")
  inY <- gedit(container = r1a, text = "prVi")
  glabel(container = r1a, text = "x = ")
  inX <- gedit(container = r1a, initial.msg = "a time series like prCh",
    text = "prCh")   
   
  # Check data consistency
  addSpace(obj = r1, value = 20, horizontal = FALSE)
  text12 <- glabel("Checking data properties", container = r1)
  r1b <- gframe(text = "", container = r1)
  
  gbutton(text = "Check data consistency", container = r1b, handler = 
    function(h, ...) { 
      en <- FALSE; cc <- "red"
      if (!all(svalue(inY) %in% ls(envir=.GlobalEnv), 
               svalue(inX) %in% ls(envir=.GlobalEnv))) {
        mes <- "Error!\nData supplied do not exist."
      } else {   
        y <- get(svalue(inY)); x <- get(svalue(inX))
        if (!inherits(y, what = "ts") || !inherits(x, what = "ts")) {
          mes <- "Error!\nBoth data should be time series."
        } else {
          if (!identical(tsp(y), tsp(x))) {
            mes <- "Error!\nThe two time series have different properties."
          } else {
            if (inherits(y, what = "mts") | inherits(x, what = "mts")) {
              mes <- "Error!\nData cannot be multiple time series."
            } else {        
              mes <- "Congratulation!\nData are good for proposed work."
              en <- TRUE; cc <- "blue"
            }
          }
        }
      }    
      svalue(wg.res) <- mes
      font(wg.res) <- list(color = cc)
      enabled(wg.sum) <- enabled(wg.p1) <- enabled(wg.p2) <- en    
    }
  )  
  rfr <- gframe("Check result", container = r1b)
  wg.res <- gtext(container = rfr, width = 400, height = 70)
   
  # Summary and two plots
  addSpace(obj = r1, value = 20, horizontal = FALSE)
  text13 <- glabel(container =r1, 
    text = "summary() | plot.ts() | ts.plot() = Understanding data")
  r1c <- gframe(text = "", container = r1)  
  wg.sum <- gbutton(text = "Data summary", container = r1c, handler = 
    function(h, ...) {  
      yx <- cbind(get(svalue(inY)), get(svalue(inX)))
      colnames(yx) <- c(svalue(inY), svalue(inX))
      print(summary(yx))
    }
  )
  enabled(wg.sum) <- FALSE
  
  addSpace(obj = r1c, value = 20, horizontal = TRUE)
  wg.p1 <- gbutton(text = "Two separate plots", container = r1c, handler = 
    function(h, ...) {
      yx <- cbind(get(svalue(inY)), get(svalue(inX)))
      colnames(yx) <- c(svalue(inY), svalue(inX))
      dev.new(); plot.ts(yx, main = "Individual data plots")
    }
  )
  enabled(wg.p1) <- FALSE
  
  addSpace(obj = r1c, value = 20, horizontal = TRUE)
  wg.p2 <- gbutton(text = "One integrated plot", container = r1c, handler = 
    function(h, ...) {
      yx <- cbind(get(svalue(inY)), get(svalue(inX)))
      colnames(yx) <- c(svalue(inY), svalue(inX))
      dev.new()
      ts.plot(yx, lty = c(1, 2), main = "One integrated data plot")
      legend("topleft", lty = c(1, 2), legend = colnames(yx)) 
    } 
  )
  enabled(wg.p2) <- FALSE   
  font(text11) <- font(text12) <-  font(text13) <- 
    list(weight = "bold", size = 10, color = "navy")  
   
  # Tab 2 = ciTarFit  -----------------------------------------------------
  arg.wg21 <- svalue(ggenericwidget(lst = "ciTarFit"))             
  arg.wg22 <- svalue(ggenericwidget(lst = "summary.ciTarFit")) 
  arg.wg21$arguments$model$items <- c("'tar'", "'mtar'")
  arg.wg22$arguments$'...' <- NULL    
  text21 <- glabel(container = r2, 
    text = "ciTarFit() = Fitting a threshold model")
  wg21 <- ggenericwidget(lst = arg.wg21, container = r2)
  text22 <- glabel(container = r2,
    text = "summary.ciTarFit() = Summarizing the fitted model")
  wg22 <- ggenericwidget(lst = arg.wg22, container = r2)
  font(text21) <- font(text22) <- 
    list(weight = "bold", size = 10, color = "blue") 
    
  # Tab 3 = ciTarThd  -----------------------------------------------------
  arg.wg31 <- svalue(ggenericwidget(lst = "ciTarThd"))
  arg.wg32 <- svalue(ggenericwidget(lst = "plot.ciTarThd"))             
  arg.wg31$arguments$model$items <- c("'tar'", "'mtar'")
  arg.wg32$arguments$'...' <- NULL
  arg.wg32$assignto <- FALSE   
  text31 <- glabel(container = r3, 
    text = "ciTarThd() = Selecting optimal threshold value")
  ggenericwidget(lst = arg.wg31, container = r3)
  text32 <- glabel(container = r3, 
    text = "plot.ciTarThd() = Plotting threshold value changes") 
  ggenericwidget(lst = arg.wg32, container = r3) 
  font(text31) <- font(text32) <- 
    list(weight = "bold", size = 10, color = "blue")  
    
  # Tab 4 = ciTarLag  -----------------------------------------------------
  arg.wg41 <- svalue(ggenericwidget(lst = "ciTarLag"))
  arg.wg42 <- svalue(ggenericwidget(lst = "plot.ciTarLag"))             
  arg.wg41$arguments$model$items <- c("'tar'", "'mtar'")
  arg.wg42$arguments$'...' <- NULL
  arg.wg42$assignto <- FALSE
  text41 <- glabel(container = r4,
    text = "ciTarLag() = Selecting optimal lag value")
  ggenericwidget(lst = arg.wg41, container = r4)
  text42 <- glabel(container = r4, 
    text = "plot.ciTarLag() = Plotting lag value changes") 
  ggenericwidget(lst = arg.wg42, container = r4) 
  font(text41) <- font(text42) <- 
    list(weight = "bold", size = 10, color = "blue") 
  
  # Tab 5 = ecmSymFit -----------------------------------------------------
  arg.wg53 <- svalue(ggenericwidget(lst = "print.ecm"))
  arg.wg54 <- svalue(ggenericwidget(lst = "summary.ecm"))
  arg.wg53$arguments$'...' <- arg.wg54$arguments$'...' <- NULL
  text51 <- glabel(container = r5, 
    text = "ecmSymFit() = Fitting symmetric ECM")
  ggenericwidget(lst = "ecmSymFit", container = r5)
  text52 <- glabel(container = r5,
    text = "ecmDiag() = Diagnostic statistics for the fitted model") 
  ggenericwidget(lst = "ecmDiag", container = r5)
  text53 <- glabel(container = r5, 
    text = "print.ecm() = Showing the fitted model") 
  ggenericwidget(lst = arg.wg53, container = r5)
  text54 <- glabel(container = r5,
    text = "summary.ecm() = Summarizing the fitted model") 
  ggenericwidget(lst = arg.wg54, container = r5)
  font(text51) <- font(text52) <- font(text53) <- font(text54) <- 
    list(weight = "bold", size = 10, color = "maroon") 
  
  # Tab 6 = ecmAsyFit -----------------------------------------------------
  text61 <- glabel(container = r6, 
    text = "ecmAsyFit() = Fitting asymmetric ECM")
  ggenericwidget(lst = "ecmAsyFit", container = r6) 
  text62 <- glabel(container = r6, 
    text = "ecmAsyTest() = Hypothesis tests on asymmetric ECM")
  ggenericwidget(lst = "ecmAsyTest", container = r6) 
  text63 <- glabel(container = r6, 
    text = "ecmDiag() = Diagnostic statistics for the fitted model")
  ggenericwidget(lst = "ecmDiag", container = r6)
  text64 <- glabel(container = r6, 
    text = "print.ecm() = Showing the fitted model")
  ggenericwidget(lst = arg.wg53, container = r6)
  text65 <- glabel(container = r6,
    text = "summary.ecm() = Summarizing the fitted model")
  ggenericwidget(lst = arg.wg54, container = r6)
  font(text61) <- font(text62) <- font(text63) <- font(text64) <- 
    font(text65) <-   list(weight = "bold", size = 10, color = "maroon") 
  
  svalue(rNB) <- 1  
  visible(rMain) <- TRUE
  invisible(NULL)
}