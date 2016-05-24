guiCor <- function()
{
  # 1. Overall layout -----------------------------------------------------
  rMain <- gwindow(title = "Correlation Between Two Random Variables",
    width = 660, height = 565, visible = TRUE, guiToolkit = "RGtk2")
  rNB <- gnotebook(container = rMain)
  rSta <- ggroup(container = rNB, horizontal = TRUE, label = "Static")
  rDyn <- ggroup(container = rNB, horizontal = FALSE, label = "Dynamic")

  # 2. first tab ----------------------------------------------------------  
  rL <- ggroup(container = rSta)
  rR <- ggroup(container = rSta, horizontal = FALSE)  
  kk <- c("Number of Points", "Color of Points",
          "Size of Points (1 to 10)", "Parameter (-1 to 1)")
  r1 <- gframe(text = kk[1], container = rR, expand = TRUE)
  r2 <- gframe(text = kk[2], container = rR, expand = TRUE)
  r3 <- gframe(text = kk[3], container = rR, expand = TRUE)
  r4 <- gframe(text = kk[4], container = rR, expand = TRUE)
  
  gsta <- ggraphics(container = rL)
  myplotSta <- function(h, ...) {
    visible(gsta) <- TRUE  # select and set the current graphics device
    par(mar = c(2, 2, 1.5, 0.2), mgp = c(0, 0.8, 0), las = 1, ps = 9)
    m <- as.numeric(gsub(pattern = ",", replacement = "",
      x = svalue(obj.ptNum)))
    sam <- rCopula(n = m, copula = normalCopula(param =svalue(obj.corPar)))
    plot(x = sam, col = svalue(obj.ptCol), xlim = c(0, 1), ylim = c(0, 1),
      type = 'p', pch = '.', cex = as.numeric(svalue(obj.ptSize)),
      xlab = "", ylab = "", main = "Static graph")
  }

  obj.ptNum <- gradio(items = c("100", "1,000", "5,000", "10,000"),
    selected = 2, horizontal = FALSE, container = r1, handler = myplotSta)
  obj.ptCol <- gcombobox(items = c("red", "green", "purple", "black",
    "yellow", "blue", "orange"), selected = 6, container = r2,
    handler = myplotSta)
  obj.ptSize <- gspinbutton(from = 1, to = 10, by = 1, value = 3,
    container = r3, handler = myplotSta)
  obj.corPar <- gslider(from = -1, to = 1, by = 0.002, value = 0.9,
      container = r4, expand = TRUE, handler = myplotSta)
  font(obj.corPar) <- list(weight = "heavy", color = "red")

  # 1.1 initial appearance of the static plot
  addHandlerExpose(obj = obj.corPar , handler = function(h, ...) {
    visible(gsta) <- TRUE
    par(mar = c(2, 2, 1.5, 0.2), mgp = c(0, 0.8, 0), las = 1, ps = 9)
    sam <- rCopula(n = 1000, copula = normalCopula(param = 0.9))
    plot(x = sam, col = "blue", xlim = c(0, 1), ylim = c(0, 1),
      type = 'p', pch = '.', cex = 3, xlab = "", ylab = "",
      main = "Static graph")
    }
  )

  # 3. second tab  --------------------------------------------------------
  rT <- ggroup(container = rDyn)
  rT1 <- gframe(text = "Number of Plots", container = rT, expand = TRUE)
  rB <- ggroup(container = rDyn)  
  gdyn <- ggraphics(container = rB)
  
  myplotDyn <- function(h, ...) {
    visible(gdyn) <- TRUE
    par(mar = c(0, 0, 0, 0))
    num.plot <- as.numeric(svalue(obj.plotNum))
    set.rho <- seq(from = -1, to = 1, length.out = num.plot)
    set.col <- rainbow(num.plot)
    for (i in 1:num.plot) {
      sam <- rCopula(n = 1000, copula = normalCopula(param = set.rho[i]))
      plot(x = sam, col = set.col[i], xlim = c(0, 1), ylim = c(0, 1),
        type = 'p', pch = '.', cex = 4, axes = FALSE, ann = FALSE,
        main = "Dynamic graph")
    }
    loc <- sample(x = seq(from = 0.2, to = 0.8, by = 0.05), size = 1)
    text(x = loc, y = 1 - loc, labels = "Game over!", col="gray", cex=0.7)
  }                                              
  obj.plotNum <- gradio(items = c("10", "50", "300", "1000"),
    selected = 2, horizontal = TRUE, container = rT1, handler = myplotDyn)

  svalue(rNB) <- 1  # The initial appearance of the GUI is tab 1.
  invisible(NULL)
}

library(copula); library(RGtk2); library(cairoDevice)
library(gWidgetsRGtk2)
guiCor()