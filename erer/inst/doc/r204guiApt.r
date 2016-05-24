guiApt <- function(pos = 1)
{
  # -----------------------------------------------------------------------
  # A. Position and size of containers
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
      "... more texts ...", sep = "\n"), icon = "info", parent = rMain)
  }

  aSepH <- gseparator(horizontal = TRUE)
  aSepV <- gseparator(horizontal = FALSE)
  aClos <- gaction(handler = hC, label = "Close", icon = "close",
                  tooltip = "Shut down the whole GUI")
  aHome <- gaction(handler = h0, label = "Home", icon = "home",
                   tooltip = "Go to tab 'Home'")
  aDat <- gaction(handler = h1, label = "1 Data property")
  aFit <- gaction(handler = h2, label = "2 Model Fit", icon = "apply",
                  tooltip = "Go to tab '2. ciTarFit'")
  aThr <- gaction(handler = h3, label = "3 Threshold choice")
  aLag <- gaction(handler = h4, label = "4 Lag choice")
  aSym <- gaction(handler = h5, label = "5 Symmetric")
  aAsy <- gaction(handler = h6, label = "6 Asymmetric", icon = "convert",
                  tooltip = "Go to tab '6. ecmAsymFit'")
  aHel <- gaction(handler = h0, label = "Welcome page")
  aApt <- gaction(handler = h7, label = "About")

  listMenu <- list(
    File = list(a = aHome, b = aClos), Data = list(a = aDat),
    Threshold = list(a = aFit, sep = aSepH, b = aThr, c = aLag),
    ECM = list(a = aSym, b = aAsy),
    Help = list(a = aHel, b = aApt))
  listTool <- list(a = aHome, sep = aSepV,
    c = aFit, sep = aSepV, d = aAsy, sep = aSepV, b = aClos)

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
    "... more texts ...", sep = "\n")
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

  # Tab 1 = Data ---------------------------------------------------------
  # Specify data names
  addSpace(obj = r1, value = 10, horizontal = FALSE)
  text11 <- glabel("Supplying two single time series data", container = r1)
  r1a <- gframe(text = "", container = r1)
  glabel(container = r1a, text = "y = ")
  inY <- gedit(container = r1a, text = "prVi")
  glabel(container = r1a, text = "x = ")
  inX <- gedit(container = r1a, initial.msg = "a time series like prCh",
    text = "prCh")

  # Omitted: Check data consistency; Summary and two plots
  # ... (more codes here) ...

  # Tab 2 = ciTarFit -----------------------------------------------------
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

  # Omitted tabs: 3 = ciTarThd; 4 = ciTarLag; 5 = ecmSymFit; 6 = ecmAsyFit
  # ... (more codes here) ...

  svalue(rNB) <- 1
  visible(rMain) <- TRUE
  invisible(NULL)
}

library(apt); library(RGtk2); library(gWidgetsRGtk2)
data(daVich); prVi <- daVich[, 1] ; prCh <- daVich[, 2]
guiApt()  # tabs at the bottom
guiApt(2) # tabs at the left