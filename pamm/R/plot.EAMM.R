`plot.EAMM` <-
  function (x, graphtype = "both", vi, vs, fun3D="wireframe",...) 
  {
    if (!inherits(x, "EAMM")) 
      stop("use only with \"EAMM\" objects")
    if (graphtype == "VI") {
      par(mfrow = c(2, 2))
      plot(x$VI[x$VS == vs], x$int.pval[x$VS == vs], type = "l", 
           xlab = "Var.Intercept", ylab = "Intercept P-value", 
           main = "Intercept P-value", ylim = c(0, 1))
      abline(h = 0.05)
      lines(x$VI[x$VS == vs], x$CIup.ipv[x$VS == vs], lty = 2)
      lines(x$VI[x$VS == vs], x$CIlow.ipv[x$VS == vs], lty = 2)
      plot(x$VI[x$VS == vs], x$sl.pval[x$VS == vs], type = "l", 
           xlab = "Var.Intercept", ylab = "Slope P-value", main = "Slope P-value", 
           ylim = c(0, 1))
      abline(h = 0.05)
      lines(x$VI[x$VS == vs], x$CIup.slpv[x$VS == vs], lty = 2)
      lines(x$VI[x$VS == vs], x$CIlow.slpv[x$VS == vs], lty = 2)
      plot(x$VI[x$VS == vs], x$int.power[x$VS == vs], type = "l", 
           ylim = c(0, 1), xlab = "Var.Intercept", ylab = "Intercept Power", 
           main = "Intercept Power Calculations")
      lines(x$VI[x$VS == vs], x$CIup.ipo[x$VS == vs], lty = 2)
      lines(x$VI[x$VS == vs], x$CIlow.ipo[x$VS == vs], lty = 2)
      plot(x$VI[x$VS == vs], x$sl.power[x$VS == vs], type = "l", 
           ylim = c(0, 1), xlab = "Var.Intercept", ylab = "Slope Power", 
           main = "Slope Power Calculations")
      lines(x$VI[x$VS == vs], x$CIup.slpo[x$VS == vs], lty = 2)
      lines(x$VI[x$VS == vs], x$CIlow.slpo[x$VS == vs], lty = 2)
    }
    if (graphtype == "VS") {
      par(mfrow = c(2, 2))
      plot(x$VS[x$VI == vi], x$int.pval[x$VI == vi], type = "l", 
           ylim = c(0, 1), xlab = "Var.Slope", ylab = "P-value", 
           main = "Intercept P-value")
      abline(h = 0.05)
      lines(x$VS[x$VI == vi], x$CIup.ipv[x$VI == vi], lty = 2)
      lines(x$VS[x$VI == vi], x$CIlow.ipv[x$VI == vi], lty = 2)
      plot(x$VS[x$VI == vi], x$sl.pval[x$VI == vi], type = "l", 
           ylim = c(0, 1), xlab = "Var.Slope", ylab = "P-value", 
           main = "Slope P-value")
      abline(h = 0.05)
      lines(x$VS[x$VI == vi], x$CIup.slpv[x$VI == vi], lty = 2)
      lines(x$VS[x$VI == vi], x$CIlow.slpv[x$VI == vi], lty = 2)
      plot(x$VS[x$VI == vi], x$int.power[x$VI == vi], type = "l", 
           ylim = c(0, 1), xlab = "Var.Slope", ylab = "Power", 
           main = "Intercept Power Calculations")
      lines(x$nb.ID[x$VI == vi], x$CIup.ipo[x$VI == vi], 
            lty = 2)
      lines(x$nb.ID[x$VI == vi], x$CIlow.ipo[x$VI == vi], 
            lty = 2)
      plot(x$VS[x$VI == vi], x$sl.power[x$VI == vi], type = "l", 
           ylim = c(0, 1), xlab = "Var.Slope", ylab = "Power", 
           main = "Slope Power Calculations")
      lines(x$VS[x$VI == vi], x$CIup.slpo[x$VI == vi], lty = 2)
      lines(x$VS[x$VI == vi], x$CIlow.slpo[x$VI == vi], lty = 2)
    }
    if (graphtype == "both") {
      if (fun3D == "wireframe"){
        par.set <- list(axis.line = list(col = "transparent"), clip = list(panel = "off"))
        p1 <- wireframe(int.pval ~ VI + VS, x, colorkey = FALSE, 
                        drape = TRUE, scales = list(arrows = FALSE, distance = c(2,2,2)), xlab = "Var.Intercept", 
                        ylab = "Var.Slope", main = "Intercept P-value",
                        zlab = list ("P-value", rot =90), 
                        screen = list(z = -50, x = -70, y = 0), par.settings = par.set)
        p2 <- wireframe(int.power ~ VI + VS, x, colorkey = FALSE, 
                        drape = TRUE, scales = list(arrows = FALSE, distance = c(2,2,2)), xlab = "Var.Intercept", 
                        ylab = "Var.Slope", main = "Intercept Power", 
                        zlab = list ("Power", rot =90), 
                        screen = list(z = -50, x = -70, y = 0), par.settings = par.set)
        p3 <- wireframe(sl.pval ~ VI + VS, x, colorkey = FALSE, 
                        drape = TRUE, scales = list(arrows = FALSE, distance = c(2,2,2)), xlab = "Var.Intercept", 
                        ylab = "Var.Slope", main = "Slope P-value", 
                        zlab = list ("P-value", rot =90), 
                        screen = list(z = -50, x = -70, y = 0), par.settings = par.set)
        p4 <- wireframe(sl.power ~ VI + VS, x, colorkey = FALSE, 
                        drape = TRUE, scales = list(arrows = FALSE, distance = c(2,2,2)), xlab = "Var.Intercept", 
                        ylab = "Var.Slope", main = "Slope Power", 
                        zlab = list ("power", rot =90), 
                        screen = list(z = -50, x = -70, y = 0), par.settings = par.set)
        print(p1, split=c(1,1,2,2),more =TRUE)
        print(p2, split=c(1,2,2,2),more =TRUE)
        print(p3, split=c(2,1,2,2),more =TRUE)
        print(p4, split=c(2,2,2,2))
      }
      if (fun3D == "open3d") {
        if (!requireNamespace("rgl", quietly = TRUE)) {
          warning("rgl package needed for this function to work. Please install it. ",
                  call. = FALSE)
        }
        else{
          rgl::open3d()
          rgl::bg3d("white")
          rgl::material3d(col = "white")
          rgl::persp3d(unique(x$VI), unique(x$VS), matrix(x$int.pval, 
                                                     nrow = length(unique(x$nb.ID)), ncol = length(unique(x$nb.repl))), 
                  col = rainbow(10), box = FALSE, zlim = c(0, 1), xlab = "VI", 
                  ylab = "VS", zlab = "int.p-value")
          rgl::open3d()
          rgl::bg3d("white")
          rgl::material3d(col = "black")
          rgl::persp3d(unique(x$VI), unique(x$VS), matrix(x$int.power, 
                                                     nrow = length(unique(x$nb.ID)), ncol = length(unique(x$nb.repl))), 
                  col = rainbow(10), box = FALSE, zlim = c(0, 1), xlab = "VI", 
                  ylab = "VS", zlab = "int.power")
          rgl::open3d()
          rgl::bg3d("white")
          rgl::material3d(col = "black")
          rgl::persp3d(unique(x$VI), unique(x$VS), matrix(x$sl.pval, 
                                                     nrow = length(unique(x$nb.ID)), ncol = length(unique(x$nb.repl))), 
                  col = rainbow(10), box = FALSE, zlim = c(0, 1), xlab = "VI", 
                  ylab = "VS", zlab = "slope.p-value")
          rgl::open3d()
          rgl::bg3d("white")
          rgl::material3d(col = "black")
          rgl::persp3d(unique(x$VI), unique(x$VS), matrix(x$sl.power, 
                                                     nrow = length(unique(x$nb.ID)), ncol = length(unique(x$nb.repl))), 
                  col = rainbow(10), box = FALSE, zlim = c(0, 1), xlab = "VI", 
                  ylab = "VS", zlab = "slope.power")
        }
      }
    }
  }
