##
## Factor analysis with lasso: GTK object
##
##  file name: plot.fanc.R
##  file contents:
##  name;YAMAMOTO, Michio
##  date created:2012.02.02.
##  date modified:2013.12.29.
##  comments: deleted the display of  "X", "F" (031012)
##            created the plot function that can be used (031012)
##           minus loading is red(031012)
##           gamma bar is added, plot for odd number of factors is modified(031912)
##           added the heatmap(100613)
##           resolved the scroll bar in the heatmap(111113)
##           added output button, and several indices for heatmap(131229)
##           modified the figures(131229)
##

##-----define the callback function-----##
##callback function for expose-event signal
##function for depicting figures


##-----variables for info.fanc-----##

##*Rad.Ellipse: Major axis of the ellipse of factor 
##  -> Minor axis is calculated by major axis
##*Len.Rec: Length of rectangular of variables
##*Window.Height: window height
##  -> window width is calculated by the number of variables and the number of factors
##*num.lambda: initial number of Lambda(Ordinal number of Lambda)
##*L: 3-dimensional array of Lambda
##  -> the dimension of Lambda corresponds to (variable, factor, lambda).
##*lambdas, *gammas: vectors of lambda, gamma.
##*lambda.current, *gamma.current: current values of lambda, gamma
##*N.var, *N.fac, *N.lambda, *N.gamma: the number of variables, the number of factors, and the number of uning parameters
##  -> calculated from L

##type in c("heatmap", "path")


##plot.fanc
##modified the cancidates of lambda as gamma varies
plot.fanc <- function (x, Window.Height=500, type=NULL, df.method="reparametrization", ...){

  ##path and heatmap is chose according to the number of variables
  if (identical(type, "path") || (is.null(type) && dim(x$loadings[[1]][[1]])[1] < 50)) {

    if(dim(x$loadings[[1]][[1]])[1] > 100) stop("The number of variables must be less than or equal to 100 to plot the solution path.")
    if(Window.Height<250 || Window.Height>2000) stop("'Window.Height' must be in [250,2000].")
    if(nchar(system.file(package="RGtk2")) == 0) stop("Package 'RGtk2' is required to plot the solution path.")
    requireNamespace("RGtk2", quietly=TRUE)
    ##if(sum(ls()=="info.fanc")!=0) stop('Object "info.fanc" must be deleted')

    ##extract the factor pattern
    L <- x$loadings
    lambdas <- x$rho
    gammas <- x$gamma
    if(df.method=="reparametrization"){
        GFIs <- x$GFI
        AGFIs <- x$AGFI
        AICs <- x$AIC
        BICs <- x$BIC
        CAICs <- x$CAIC
    }
    if(df.method=="active"){
        GFIs <- x$GFI
        AGFIs <- x$AGFI_dfnonzero
        AICs <- x$AIC_dfnonzero
        BICs <- x$BIC_dfnonzero
        CAICs <- x$CAIC_dfnonzero
    }
    info.fanc <- list("Rad.Ellipse"=50, "Len.Rec"=50, "Window.Height"=Window.Height,
                      "N.var"=NULL, "N.fac"=NULL, "N.lambda"=NULL,
                      "L"=NULL, "lambdas"=NULL, "num.lambda"=1, "num.gamma"=1,"num.GFI"=1)
                                        #if(x$factors==1) L <- array(L,dim=c(nrow(x$rho),1,ncol(x$rho)))
                                        #assign("info,fanc", info.fanc, envir=infofanc)
    info.fanc$lambda.current <- lambdas[1,1]
    info.fanc$gamma.current <- gammas[1]
    info.fanc$GFI.current <- GFIs[1]
    info.fanc$AGFI.current <- AGFIs[1]
    info.fanc$AIC.current <- AICs[1]
    info.fanc$BIC.current <- BICs[1]
    info.fanc$CAIC.current <- CAICs[1]
    info.fanc$L <- L
    info.fanc$lambdas <- lambdas
    info.fanc$gammas <- gammas
    info.fanc$num.lambda <- 1
    info.fanc$num.gamma <- 1
    info.fanc$GFIs <- GFIs
    info.fanc$AGFIs <- AGFIs
    info.fanc$AICs <- AICs
    info.fanc$BICs <- BICs
    info.fanc$CAICs <- CAICs
    info.fanc$N.var <- dim(info.fanc$L[[1]][[1]])[1]
    info.fanc$N.fac <- dim(info.fanc$L[[1]][[1]])[2]
    info.fanc$N.lambda <- length(info.fanc$L[[1]])
    info.fanc$N.gamma <- length(info.fanc$L)
    info.fanc$Window.Width <- max(((info.fanc$N.var+1) * info.fanc$Len.Rec + (info.fanc$N.var+2) * info.fanc$Len.Rec / 2), ((info.fanc$N.fac) * 1.5 * info.fanc$Rad.Ellipse),650)

    cbExposeCanvas <- function (gtk.widget, data)
      {
        N.var <- info.fanc$N.var
        N.fac <- info.fanc$N.fac
        Window.Width <- info.fanc$Window.Width
        Window.Height <- info.fanc$Window.Height
        Rad.Ellipse <- info.fanc$Rad.Ellipse
        Len.Rec <- info.fanc$Len.Rec

        drawable <- gtk.widget$window
        HN.var <- N.var / 2
        HN.fac <- N.fac / 2

        cr <- RGtk2::gdkCairoCreate (drawable)
        cr.t <- RGtk2::gdkCairoCreate (drawable)
        Mat <- RGtk2::cairoMatrixInit (0, 0, 0, 0, 0, 0)$matrix

        ##-------------------
        ##   depict the Ellipse of factor
        ##-------------------
        Rem.N.fac <- N.fac %% 2
        RGtk2::cairoTranslate (cr, Window.Width / 2, Rad.Ellipse + 20) ##movement of origin
        RGtk2::cairoTranslate (cr.t, Window.Width / 2, Rad.Ellipse + 20) ##movement of origin

        RGtk2::cairoSetLineWidth (cr, 2.5)

        RGtk2::cairoScale (cr, 1.0, 0.5) ##scale of ellipse

        if (Rem.N.fac == 0) { ##for odd number of factors
          for (i in 1:HN.fac) {
            ii <- i - 1

            ##depict the ellipse
            RGtk2::cairoArc (cr, -Rad.Ellipse * 5 / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 0,
                      Rad.Ellipse, 0.0, 2.0 * pi)
            RGtk2::cairoStroke (cr)
            RGtk2::cairoArc (cr,  Rad.Ellipse * 5 / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 0,
                      Rad.Ellipse, 0.0, 2.0 * pi)
            RGtk2::cairoStroke (cr)

            RGtk2::cairoSetFontSize (cr.t, 20)

            ##depict the label
            text <- sprintf ("f%d", HN.fac - ii)
            RGtk2::cairoMoveTo (cr.t, -Rad.Ellipse * 5 / 4 - Rad.Ellipse / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 5)
            RGtk2::cairoShowText (cr.t, text)

            text <- sprintf ("f%d", HN.fac + 1 + ii)
            RGtk2::cairoMoveTo (cr.t,  Rad.Ellipse * 5 / 4 - Rad.Ellipse / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 5)
            RGtk2::cairoShowText (cr.t, text)
          }
        }
        else if (Rem.N.fac != 0) { ##for even number of factors
          RGtk2::cairoArc (cr, 0.0, 0.0, Rad.Ellipse, 0.0, 2.0 * pi)
          RGtk2::cairoStroke (cr)
          RGtk2::cairoSetFontSize (cr.t, 20)
          text <- sprintf ("f%d", floor(HN.fac) + 1)
          RGtk2::cairoMoveTo (cr.t, 0 - Rad.Ellipse / 4, 5)
          RGtk2::cairoShowText (cr.t, text)

          if (floor(HN.fac) != 0) {
            for (i in 1:floor(HN.fac)) {
              ii <- i - 1

            ##depict the ellipse
              RGtk2::cairoArc (cr, -Rad.Ellipse * 10 / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 0,
                        Rad.Ellipse, 0.0, 2.0 * pi)
              RGtk2::cairoStroke (cr)
              RGtk2::cairoArc (cr,  Rad.Ellipse * 10 / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 0,
                        Rad.Ellipse, 0.0, 2.0 * pi)
              RGtk2::cairoStroke (cr)

              RGtk2::cairoSetFontSize (cr.t, 20)

            ##depict the label
              text <- sprintf ("f%d", floor(HN.fac) - ii)
              RGtk2::cairoMoveTo (cr.t, -Rad.Ellipse * 10 / 4 - Rad.Ellipse / 4 - ii * 2 * Rad.Ellipse * 5 / 4, 5)
              RGtk2::cairoShowText (cr.t, text)

              text <- sprintf ("f%d", floor(HN.fac) + ii + 2)
              RGtk2::cairoMoveTo (cr.t,  Rad.Ellipse * 10 / 4 - Rad.Ellipse / 4 + ii * 2 * Rad.Ellipse * 5 / 4, 5)
              RGtk2::cairoShowText (cr.t, text)
            }
          }
        } ##finish the depict of ellipse

        ##change back to the coordinates
        RGtk2::cairoGetMatrix (cr, Mat)
        RGtk2::cairoMatrixInvert (Mat)
        RGtk2::cairoTransform (cr, Mat)
        RGtk2::cairoGetMatrix (cr.t, Mat)
        RGtk2::cairoMatrixInvert (Mat)
        RGtk2::cairoTransform (cr.t, Mat)


        ##-------------------
        ##  depict rectangular of variables 
        ##-------------------
        ##Change the pattern when the number of variables is even (odd)
        Rem.N.var <- N.var %% 2
        RGtk2::cairoTranslate (  cr, Window.Width / 2, Window.Height - Len.Rec * 5 / 2)
        RGtk2::cairoTranslate (cr.t, Window.Width / 2, Window.Height - Len.Rec * 5 / 2)

        RGtk2::cairoSetLineWidth (cr, 2.5)

        if (Rem.N.var == 0) {##for even number of variables
          for (i in 1:HN.var) {
            ii <- i - 1
            RGtk2::cairoRectangle (cr, -Len.Rec * 5 / 4 - ii * Len.Rec * 3 / 2, 0, Len.Rec, Len.Rec)
            RGtk2::cairoStroke (cr)
            RGtk2::cairoRectangle (cr, Len.Rec / 4 + ii * Len.Rec * 3 / 2, 0, Len.Rec, Len.Rec)
            RGtk2::cairoStroke (cr)

            RGtk2::cairoSetFontSize (cr.t, 20)

            text <- sprintf ("x%d", floor(HN.var) - ii)
            RGtk2::cairoMoveTo (cr.t, -Len.Rec * 5 / 4 - ii * Len.Rec * 3 / 2 + Len.Rec / 4, Len.Rec * 2 / 3)
            RGtk2::cairoShowText (cr.t, text)

            text <- sprintf ("x%d", floor(HN.var) + 1 + ii)
            RGtk2::cairoMoveTo (cr.t, Len.Rec / 4 + ii * Len.Rec * 3 / 2 + Len.Rec / 4, Len.Rec * 2 / 3)
            RGtk2::cairoShowText (cr.t, text)
          }
        }
        else if (Rem.N.var != 0) { ##for odd number of variables
          RGtk2::cairoRectangle (cr, -Len.Rec / 2, 0, Len.Rec, Len.Rec)
          RGtk2::cairoStroke (cr)
          RGtk2::cairoSetFontSize (cr.t, 20)
          text <- sprintf("x%d", floor(HN.var) + 1)
          RGtk2::cairoMoveTo (cr.t, 0 - Len.Rec / 2 + Len.Rec / 3, Len.Rec * 2 / 3)
          RGtk2::cairoShowText (cr.t, text)

          for (i in 1:HN.var) {
            ii <- i - 1
            RGtk2::cairoRectangle (cr, -Len.Rec * 2 - ii * Len.Rec * 3 / 2, 0, Len.Rec, Len.Rec)
            RGtk2::cairoStroke (cr)
            RGtk2::cairoRectangle (cr, Len.Rec      + ii * Len.Rec * 3 / 2, 0, Len.Rec, Len.Rec)
            RGtk2::cairoStroke (cr)

            RGtk2::cairoSetFontSize (cr.t, 20)

            text <- sprintf ("x%d", floor(HN.var) - ii)
            RGtk2::cairoMoveTo (cr.t, -Len.Rec * 2 - ii * Len.Rec * 3 / 2 + Len.Rec / 3 - 0.5, Len.Rec * 2 / 3)
            RGtk2::cairoShowText (cr.t, text)

            text <- sprintf ("x%d", floor(HN.var) + 2 + ii)
            RGtk2::cairoMoveTo (cr.t, Len.Rec      + ii * Len.Rec * 3 / 2 + Len.Rec / 3 - 0.5, Len.Rec * 2 / 3)
            RGtk2::cairoShowText (cr.t, text)
          }
        }

        ##change back to the coordinates
        RGtk2::cairoGetMatrix (cr, Mat)
        RGtk2::cairoMatrixInvert (Mat)
        RGtk2::cairoTransform (cr, Mat)
        RGtk2::cairoGetMatrix (cr.t, Mat)
        RGtk2::cairoMatrixInvert (Mat)
        RGtk2::cairoTransform (cr.t, Mat)



        ##---------------------------
        ##  depict the line between variables and factors
        ##---------------------------
        if (Rem.N.fac == 0) { ##for even number of factors
          cr.Fac.X <- -Rad.Ellipse * 5 / 4 - (HN.fac - 1) * 2 * Rad.Ellipse * 5 / 4 + Window.Width / 2
          cr.Fac.Y <- Rad.Ellipse + 20 + Rad.Ellipse / 2
        }
        else if (Rem.N.fac != 0) { ##for odd number of factors
          cr.Fac.X <- -Rad.Ellipse * 10 / 4 - (floor(HN.fac) - 1) * 2 * Rad.Ellipse * 5 / 4 + Window.Width / 2
          cr.Fac.Y <- Rad.Ellipse + 20 + Rad.Ellipse / 2
        }

        if (Rem.N.var == 0) { ##for even number of variables
          cr.Var.X <- -Len.Rec * 5 / 4 - (HN.var - 1) * Len.Rec * 3 / 2 + Window.Width / 2 + Len.Rec / 2
          cr.Var.Y <- Window.Height - Len.Rec * 5 / 2
        }
        else if (Rem.N.var != 0) { ##for odd number of variables
          cr.Var.X <- -Len.Rec * 2 - (floor(HN.var) - 1) * Len.Rec * 3 / 2 + Window.Width / 2 + Len.Rec / 2
          cr.Var.Y <- Window.Height - Len.Rec * 5 / 2
        }

        ##depict all of the paths individually
        RGtk2::cairoSetLineJoin (cr, 'CAIRO_LINE_JOIN_BEVEL')
        ##  RGtk2::cairoSetLineJoin (cr, 'CAIRO_LINE_JOIN_MITER')
        for (i in 1:N.var) {
          ii <- i - 1
          for (j in 1:N.fac) {
            jj <- j - 1
            ## RGtk2::cairoSetLineWidth (cr, abs(info.fanc$L[i, j, info.fanc$num.lambda, info.fanc$num.gamma] * 10))
            RGtk2::cairoSetLineWidth (cr, abs(info.fanc$L[[info.fanc$num.gamma]][[info.fanc$num.lambda]][i, j]  * 10))

            ## if (info.fanc$L[i, j, info.fanc$num.lambda, info.fanc$num.gamma] < 0) {
            if (info.fanc$L[[info.fanc$num.gamma]][[info.fanc$num.lambda]][i, j] < 0) {
              RGtk2::cairoSetSourceRgb (cr, 255.0, 0.0, 0.0)
            }
            else {
              RGtk2::cairoSetSourceRgb (cr, 0.0, 0.0, 0.0)
            }
            RGtk2::cairoMoveTo (cr, cr.Fac.X + jj * 2 * Rad.Ellipse * 5 / 4, cr.Fac.Y)
            RGtk2::cairoLineTo (cr, cr.Var.X + ii * Len.Rec * 3 / 2, cr.Var.Y)
            RGtk2::cairoStroke (cr)
          }
        }
      }



    ##callback function for expose-event signal of GtkWidget label
    ##-----for 'rho'-----
    cbExposeLabelLambda <- function (gtk.widget, data)
      {
        drawable <- gtk.widget$window

        cr <- RGtk2::gdkCairoCreate (drawable)

        text <- sprintf ("rho : %f", info.fanc$lambda.current)
        RGtk2::cairoSetFontSize (cr, 15)
        RGtk2::cairoMoveTo (cr, 0, 30)
        RGtk2::cairoShowText (cr, text)
      }

    ##-----for 'GFI'-----
    cbExposeLabelGFI <- function (gtk.widget, data)
      {
        drawable <- gtk.widget$window

        cr <- RGtk2::gdkCairoCreate (drawable)

        text <- paste("GFI: ", signif(info.fanc$GFI.current, digits=6), sep="")
        RGtk2::cairoSetFontSize (cr, 15)
        RGtk2::cairoMoveTo (cr, 0, 30)
        RGtk2::cairoShowText (cr, text)
      }

    ##-----for 'AGFI'-----
    cbExposeLabelAGFI <- function (gtk.widget, data)
      {
        drawable <- gtk.widget$window

        cr <- RGtk2::gdkCairoCreate (drawable)

        text <- paste("AGFI: ", signif(info.fanc$AGFI.current, digits=6), sep="")
        RGtk2::cairoSetFontSize (cr, 15)
        RGtk2::cairoMoveTo (cr, 0, 30)
        RGtk2::cairoShowText (cr, text)
      }

    ##-----for 'AIC'-----
    cbExposeLabelAIC <- function (gtk.widget, data)
      {
        drawable <- gtk.widget$window

        cr <- RGtk2::gdkCairoCreate (drawable)

        text <- paste("AIC: ", signif(info.fanc$AIC.current, digits=6), sep="")
        RGtk2::cairoSetFontSize (cr, 15)
        RGtk2::cairoMoveTo (cr, 0, 30)
        RGtk2::cairoShowText (cr, text)
      }

    ##-----for 'BIC'-----
    cbExposeLabelBIC <- function (gtk.widget, data)
      {
        drawable <- gtk.widget$window

        cr <- RGtk2::gdkCairoCreate (drawable)

        text <- paste("BIC: ", signif(info.fanc$BIC.current, digits=6), sep="")
        RGtk2::cairoSetFontSize (cr, 15)
        RGtk2::cairoMoveTo (cr, 0, 30)
        RGtk2::cairoShowText (cr, text)
      }

    ##-----for 'CAIC'-----
    cbExposeLabelCAIC <- function (gtk.widget, data)
      {
        drawable <- gtk.widget$window

        cr <- RGtk2::gdkCairoCreate (drawable)

        text <- paste("CAIC: ", signif(info.fanc$CAIC.current, digits=6), sep="")
        RGtk2::cairoSetFontSize (cr, 15)
        RGtk2::cairoMoveTo (cr, 0, 30)
        RGtk2::cairoShowText (cr, text)
      }


    ##-----for 'gamma'-----
    cbExposeLabelGamma <- function (gtk.widget, data)
      {
        drawable <- gtk.widget$window

        cr <- RGtk2::gdkCairoCreate (drawable)

        text <- sprintf ("gam : %f", info.fanc$gamma.current)
        RGtk2::cairoSetFontSize (cr, 15)
        RGtk2::cairoMoveTo (cr, 0, 30)
        RGtk2::cairoShowText (cr, text)
      }


    ##callback function for value-changed signal of GtkScale widget ("canvas")
    ##-----for 'rho'-----
    cbValueChangedLambda <- function (gtk.scale, data)
      {
        if (RGtk2::gtkRangeGetValue (gtk.scale) == 0) {
          info.fanc$num.lambda <<- 1
          ##    assign("info.fanc$num.lambda" ,1, pos=globalenv())
        }
        else {
          info.fanc$num.lambda <<- ceiling(RGtk2::gtkRangeGetValue (gtk.scale) * info.fanc$N.lambda) ##modify error
        }
        RGtk2::gtkWidgetQueueDraw (data) ##make expose-event signal for canvas
      }

    ##-----for 'gamma'-----
    cbValueChangedGamma <- function (gtk.scale, data)
      {
        if (RGtk2::gtkRangeGetValue (gtk.scale) == 0) {
          info.fanc$num.gamma <<- 1
        }
        else {
          info.fanc$num.gamma <<- ceiling(RGtk2::gtkRangeGetValue (gtk.scale) * info.fanc$N.gamma) ##modify error
        }

        info.fanc$lambda.current <<- info.fanc$lambdas[info.fanc$num.lambda, info.fanc$num.gamma]

        RGtk2::gtkWidgetQueueDraw (data) ##make expose-event signal for canvas
      }



    ##callback function for value-changed signal of GtkScale widget ("label")
    ##-----for 'rho'-----
    cbValueChangedLabelLambda <- function (gtk.scale, data)
      {
                                        #if (RGtk2::gtkRangeGetValue (gtk.scale) == 0) {
                                        #  info.fanc$lambda.current <<- info.fanc$lambdas[1, info.fanc$num.gamma]
                                        #}
                                        #else {
                                        #  info.fanc$lambda.current <<- info.fanc$lambdas[ceiling(RGtk2::gtkRangeGetValue (gtk.scale) * info.fanc$N.lambda), info.fanc$num.gamma] ##modify error
                                        #}
        info.fanc$lambda.current <<- info.fanc$lambdas[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$GFI.current <<- info.fanc$GFIs[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$AGFI.current <<- info.fanc$AGFIs[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$AIC.current <<- info.fanc$AICs[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$BIC.current <<- info.fanc$BICs[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$CAIC.current <<- info.fanc$CAICs[info.fanc$num.lambda, info.fanc$num.gamma]

        RGtk2::gtkWidgetQueueDraw (data[[1]]) ##make an expose-event signal for label.lambda
        RGtk2::gtkWidgetQueueDraw (data[[2]]) ##make an expose-event signal for label.GFI
        RGtk2::gtkWidgetQueueDraw (data[[3]]) ##make an expose-event signal for label.AGFI
        RGtk2::gtkWidgetQueueDraw (data[[4]]) ##make an expose-event signal for label.AIC
        RGtk2::gtkWidgetQueueDraw (data[[5]]) ##make an expose-event signal for label.BIC
        RGtk2::gtkWidgetQueueDraw (data[[6]]) ##make an expose-event signal for label.CAIC
      }

    ##-----for 'gamma'-----
    cbValueChangedLabelGamma <- function (gtk.scale, data)
      {
        if (RGtk2::gtkRangeGetValue (gtk.scale) == 0) {
          info.fanc$gamma.current <<- info.fanc$gammas[1]
        }
        else {
          info.fanc$gamma.current <<- info.fanc$gammas[ceiling(RGtk2::gtkRangeGetValue (gtk.scale) * info.fanc$N.gamma)] ##modify error
        }

        info.fanc$lambda.current <<- info.fanc$lambdas[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$GFI.current <<- info.fanc$GFIs[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$AGFI.current <<- info.fanc$AGFIs[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$AIC.current <<- info.fanc$AICs[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$BIC.current <<- info.fanc$BICs[info.fanc$num.lambda, info.fanc$num.gamma]
        info.fanc$CAIC.current <<- info.fanc$CAICs[info.fanc$num.lambda, info.fanc$num.gamma]

        RGtk2::gtkWidgetQueueDraw (data[[1]]) ##make an expose-event signal for label.gamma
        RGtk2::gtkWidgetQueueDraw (data[[2]]) ##make an expose-event signal forlabel.lambda
        RGtk2::gtkWidgetQueueDraw (data[[3]]) ##make an expose-event signal forlabel.GFI
        RGtk2::gtkWidgetQueueDraw (data[[4]]) ##make an expose-event signal forlabel.AGFI
        RGtk2::gtkWidgetQueueDraw (data[[5]]) ##make an expose-event signal forlabel.AIC
        RGtk2::gtkWidgetQueueDraw (data[[6]]) ##make an expose-event signal forlabel.BIC
        RGtk2::gtkWidgetQueueDraw (data[[7]]) ##make an expose-event signal forlabel.CAIC
      }


    ## callback function for clicked signal of GtkScale widget ("button.loadings")
	cbButtonClicked <- function (gtk.widget, data)
      {
        print(info.fanc$L[[info.fanc$num.gamma]][[info.fanc$num.lambda]])
      }


    ##-----main window to depict-----##
    MakeInterface <- function (gchar.title)
      {
        ##"info" includes 'L' and 'lambdas'.
        ##'L' is a p*m*r array of pattern matrices
        ##'lambdas' is a r*1 vector of tuning parameter

        lambdas <- info.fanc$lambdas
        gammas <- info.fanc$gammas

        ##the range of lambda, gamma are restricted in (0, 1)
        ##for 'rho'
        N.lambda <- length(lambdas[,1])
        Min.lambda <- 0
        Max.lambda <- 1
        Step.lambda <- 1 / (N.lambda - 1)

        ##for 'gamma'
        N.gamma <- length(gammas)
        Min.gamma <- 0
        Max.gamma <- 1
        Step.gamma <- 1 / (N.gamma - 1)

        window <- RGtk2::gtkWindowNew (show=FALSE)
        RGtk2::gtkWindowSetTitle (window, gchar.title)
        ##  RGtk2::gtkWidgetSetSizeRequest (window, info.fanc$Window.Width, info.fanc$Window.Height)
        RGtk2::gtkWidgetSetSizeRequest (window, info.fanc$Window.Width, -1)
        RGtk2::gtkContainerSetBorderWidth (window, 5)
        ##  RGtk2::gSignalConnect (window, "destroy", gtkMainQuit, NULL) ##maybe not needed


        ##-------------------
        ##    depict the path diagram
        ##-------------------
        vbox <- RGtk2::gtkVBoxNew (FALSE, 3)
        RGtk2::gtkContainerAdd (window, vbox)
        RGtk2::gtkContainerSetBorderWidth (vbox, 5)

        alignment <- RGtk2::gtkAlignmentNew (0.5, 0.5, 0, 0)
        RGtk2::gtkBoxPackStart (vbox, alignment, TRUE, TRUE, 0)

        canvas <- RGtk2::gtkDrawingAreaNew ()
        RGtk2::gtkContainerAdd (alignment, canvas)
        RGtk2::gtkWidgetSetSizeRequest (canvas, info.fanc$Window.Width, info.fanc$Window.Height)
        RGtk2::gSignalConnect (canvas, "expose-event", cbExposeCanvas, NULL)


        ##-----------------------
        ##  depict the goodness of fit indices
        ##-----------------------
        ##-----for goodness of fit indecies-----
        hbox.GFI <- RGtk2::gtkHBoxNew (FALSE, 5)
        RGtk2::gtkBoxPackStart (vbox, hbox.GFI, FALSE, FALSE, 0)

        ##depict the value of GFI
        alignment.GFI <- RGtk2::gtkAlignmentNew (0.5, 0.5, 0, 0)
        RGtk2::gtkBoxPackStart (hbox.GFI, alignment.GFI, TRUE, FALSE, 0)

        label.GFI <- RGtk2::gtkDrawingAreaNew ()
        RGtk2::gtkContainerAdd (alignment.GFI, label.GFI)
        RGtk2::gtkWidgetSetSizeRequest (label.GFI, 130, 50)
        RGtk2::gSignalConnect (label.GFI, "expose-event", cbExposeLabelGFI, NULL)

        ## RGtk2::gtkBoxPackStart (hbox.GFI, scale.lambda, TRUE, TRUE, 0)


        ##depict the value of AGFI
        ##  vbox.AGFI <- RGtk2::gtkVBoxNew (FALSE, 5)
        ##  RGtk2::gtkBoxPackStart (vbox.AGFI, hbox.GFI, FALSE, FALSE, 0)
        alignment.AGFI <- RGtk2::gtkAlignmentNew (0.5, 0.5, 0, 0)
        RGtk2::gtkBoxPackStart (hbox.GFI, alignment.AGFI, TRUE, FALSE, 0)

        label.AGFI <- RGtk2::gtkDrawingAreaNew ()
        RGtk2::gtkContainerAdd (alignment.AGFI, label.AGFI)
        RGtk2::gtkWidgetSetSizeRequest (label.AGFI, 130, 50)
        RGtk2::gSignalConnect (label.AGFI, "expose-event", cbExposeLabelAGFI, NULL)

        ##depict the value of AIC
        ##  vbox.AIC <- RGtk2::gtkVBoxNew (FALSE, 5)
        ##  RGtk2::gtkBoxPackStart (vbox.AIC, hbox.GFI, FALSE, FALSE, 0)
        alignment.AIC <- RGtk2::gtkAlignmentNew (0.5, 0.5, 0, 0)
        RGtk2::gtkBoxPackStart (hbox.GFI, alignment.AIC, TRUE, FALSE, 0)

        label.AIC <- RGtk2::gtkDrawingAreaNew ()
        RGtk2::gtkContainerAdd (alignment.AIC, label.AIC)
        RGtk2::gtkWidgetSetSizeRequest (label.AIC, 130, 50)
        RGtk2::gSignalConnect (label.AIC, "expose-event", cbExposeLabelAIC, NULL)


        ##depict the value of BIC
        ##  vbox.BIC <- RGtk2::gtkVBoxNew (FALSE, 5)
        ##  RGtk2::gtkBoxPackStart (vbox.BIC, hbox.GFI, FALSE, FALSE, 0)
        alignment.BIC <- RGtk2::gtkAlignmentNew (0.5, 0.5, 0, 0)
        RGtk2::gtkBoxPackStart (hbox.GFI, alignment.BIC, TRUE, FALSE, 0)

        label.BIC <- RGtk2::gtkDrawingAreaNew ()
        RGtk2::gtkContainerAdd (alignment.BIC, label.BIC)
        RGtk2::gtkWidgetSetSizeRequest (label.BIC, 130, 50)
        RGtk2::gSignalConnect (label.BIC, "expose-event", cbExposeLabelBIC, NULL)


        ##depict the value of CAIC
        ##  vbox.CAIC <- RGtk2::gtkVBoxNew (FALSE, 5)
        ##  RGtk2::gtkBoxPackStart (vbox.CAIC, hbox.GFI, FALSE, FALSE, 0)
        alignment.CAIC <- RGtk2::gtkAlignmentNew (0.5, 0.5, 0, 0)
        RGtk2::gtkBoxPackStart (hbox.GFI, alignment.CAIC, TRUE, FALSE, 0)

        label.CAIC <- RGtk2::gtkDrawingAreaNew ()
        RGtk2::gtkContainerAdd (alignment.CAIC, label.CAIC)
        RGtk2::gtkWidgetSetSizeRequest (label.CAIC, 130, 50)
        RGtk2::gSignalConnect (label.CAIC, "expose-event", cbExposeLabelCAIC, NULL)


        ##-----------------------
        ##   make the scale bar
        ##-----------------------
        ##-----for 'rho'-----
        hbox.lambda <- RGtk2::gtkHBoxNew (FALSE, 5)
        RGtk2::gtkBoxPackStart (vbox, hbox.lambda, FALSE, FALSE, 0)

        scale.lambda <- RGtk2::gtkHScaleNewWithRange (Min.lambda, Max.lambda, Step.lambda)
        RGtk2::gtkScaleSetDigits (scale.lambda, 2)
        RGtk2::gtkScaleSetDrawValue (scale.lambda, FALSE) ##The presence or absence of display of value under the bar

        ##depict the value of rho
        alignment.lambda <- RGtk2::gtkAlignmentNew (0.5, 0.5, 0, 0)
        RGtk2::gtkBoxPackStart (hbox.lambda, alignment.lambda, FALSE, FALSE, 0)

        label.lambda <- RGtk2::gtkDrawingAreaNew ()
        RGtk2::gtkContainerAdd (alignment.lambda, label.lambda)
        RGtk2::gtkWidgetSetSizeRequest (label.lambda, 120, 50)
                                        #RGtk2::gtkWidgetSetSizeRequest (label.lambda, 100, 50)
        RGtk2::gSignalConnect (label.lambda, "expose-event", cbExposeLabelLambda, NULL)

        RGtk2::gtkBoxPackStart (hbox.lambda, scale.lambda, TRUE, TRUE, 0)



        ##-----for 'gamma'-----
        hbox.gamma <- RGtk2::gtkHBoxNew (FALSE, 5)
        RGtk2::gtkBoxPackStart (vbox, hbox.gamma, FALSE, FALSE, 0)

        scale.gamma <- RGtk2::gtkHScaleNewWithRange (Min.gamma, Max.gamma, Step.gamma)
        RGtk2::gtkScaleSetDigits (scale.gamma, 2)
        RGtk2::gtkScaleSetDrawValue (scale.gamma, FALSE) ##The presence or absence of display of value under the bar

        ##depict the value of gamma
        alignment.gamma <- RGtk2::gtkAlignmentNew (0.5, 0.5, 0, 0)
        RGtk2::gtkBoxPackStart (hbox.gamma, alignment.gamma, FALSE, FALSE, 0)

        label.gamma <- RGtk2::gtkDrawingAreaNew ()
        RGtk2::gtkContainerAdd (alignment.gamma, label.gamma)
        RGtk2::gtkWidgetSetSizeRequest (label.gamma, 120, 50)
                                        #RGtk2::gtkWidgetSetSizeRequest (label.gamma, 100, 50)
        RGtk2::gSignalConnect (label.gamma, "expose-event", cbExposeLabelGamma, NULL)

        RGtk2::gtkBoxPackStart (hbox.gamma, scale.gamma, TRUE, TRUE, 0)


        ##-----------------------
        ##    display the output button
        ##-----------------------
        hbox.button <- RGtk2::gtkHBoxNew (FALSE, 5)
        RGtk2::gtkBoxPackStart (vbox, hbox.button, FALSE, FALSE, 0)
        button.loadings <- RGtk2::gtkButtonNewWithLabel ("Output loadings")
        alignment.button <- RGtk2::gtkAlignmentNew (0.5, 0.5, 0, 0)
        RGtk2::gtkBoxPackEnd (hbox.button, alignment.button, FALSE, FALSE, 0)
        RGtk2::gtkContainerAdd (alignment.button, button.loadings)
        RGtk2::gtkWidgetSetSizeRequest (button.loadings, 120, 50)
        RGtk2::gSignalConnect (button.loadings, "clicked", cbButtonClicked, NULL)


        ##---------------------------------
        ##   make value-changed signal
        ##--------------------------------
        RGtk2::gSignalConnect (scale.lambda, "value_changed", cbValueChangedLambda, canvas)
        RGtk2::gSignalConnect (scale.lambda, "value_changed", cbValueChangedLabelLambda, list(label.lambda, label.GFI, label.AGFI, label.AIC, label.BIC, label.CAIC))
        RGtk2::gSignalConnect (scale.gamma, "value_changed", cbValueChangedGamma, canvas)
        RGtk2::gSignalConnect (scale.gamma, "value_changed", cbValueChangedLabelGamma, list(label.gamma, label.lambda, label.GFI, label.AGFI, label.AIC, label.BIC, label.CAIC))


        ##---------------------------------
        ##        display Window
        ##--------------------------------
        RGtk2::gtkWidgetShowAll(window)
      }

    MakeInterface ("Factor analysis with MC+")


  } else if (identical(type, "heatmap") || (is.null(type) && dim(x$x)[2] >= 50)) {

    if(nchar(system.file(package="matlab")) == 0){
      answer <- readline("The package 'matlab' is required to plot the solution path. \nDo you want to install 'matlab' now?  (y/n)")
      if(answer=="y"){
        install.packages("matlab")
        if (nchar(system.file(package="matlab")) == 0) stop('The package "matlab" was not able to be installed')
      } else {
        stop("The plot was terminated.")
      }
    }

    if(nchar(system.file(package="tcltk")) == 0){
      answer <- readline("The package 'tcltk' is required to plot the solution path. \nDo you want to install 'tcltk' now?  (y/n)")
      if(answer=="y"){
          install.packages("tcltk")
          if (nchar(system.file(package="tcltk")) == 0) stop('The package "tcltk" was not able to be installed')
      }else{
          stop("The plot was terminated.")
      }
    }

    requireNamespace("matlab", quietly=TRUE)
    requireNamespace("tcltk", quietly=TRUE)

    ##For image function
    n.col <- 256
    col.red <- rgb(red=1, green = (0:n.col)/n.col, blue = (0:n.col)/n.col, names = paste("red", 0:n.col, sep = "."))
    col.black <- rgb(red=(0:n.col)/n.col, green = (0:n.col)/n.col, blue = (0:n.col)/n.col, names = paste("black", 0:n.col, sep = "."))
    col.all <- c(col.red,rev(col.black))
    max.col <- 0
    for(i in 1:length(x$loadings)){
      for(j in 1:length(x$loadings[[1]])){
        max.col <- max(max.col,max(abs(x$loadings[[i]][[j]])))
      }
    }


    N.rho <- nrow(fit$rho); N.gamma <- ncol(fit$rho); sl.names <- c("rho", "gamma"); sl.mins <- c(1, 1); sl.maxs <- c(N.rho, N.gamma); sl.deltas <- c(1, 1); sl.defaults <- c(1, 1); fit <- fit
    slider.env <- new.env()
    slider.fanc <- function (sl.names, sl.mins, sl.maxs, sl.deltas,
                             sl.defaults, but.functions, but.names, no, set.no.value,
                             obj.name, obj.value, reset.function, title, fit) {
      N.rho <- sl.maxs[1]
      N.gamma <- sl.maxs[2]
      if (!missing(no))
        return(as.numeric(tcltk::tclvalue(get(paste("slider", no, sep = ""), envir = slider.env))))
      if (!missing(set.no.value)) {
        try(eval(parse(text = paste("tcltk::tclvalue(slider", set.no.value[1], ")<-", set.no.value[2], sep = "")), envir = slider.env))
        return(set.no.value[2])
      }
      ## if (!exists("slider.env"))
      ##   slider.env <<- new.env()
      if (!missing(obj.name)) {
        if (!missing(obj.value)) {
          assign(obj.name, obj.value, envir = slider.env)
        } else {
          obj.value <- get(obj.name, envir = slider.env)
        }
        return(obj.value)
      }
      if (missing(title))
        title <- "slider control widget"

      nt <- tcltk::tktoplevel()
      tcltk::tkwm.title(nt, title)
      tcltk::tkwm.geometry(nt, "+0+0")
      if (missing(sl.names))
        sl.names <- NULL

      ##output rho
      i <- 1
      eval(parse(text = paste("assign('slider", i, "',tcltk::tclVar(sl.defaults[i]),envir=slider.env)", sep = "")))
      tcltk::tkpack(fr.rho <- tcltk::tkframe(nt))
      sc1 <- tcltk::tkscale(fr.rho, from = sl.mins[i], to = sl.maxs[i], showvalue = F, resolution = sl.deltas[i], orient = "horiz", length = 500)
      tcltk::tkpack(fr.rho, sc1)
      assign("sc1", sc1, envir = slider.env)
      eval(parse(text = paste("tcltk::tkconfigure(sc1,variable=slider", i, ")", sep = "")), envir = slider.env)
      label.temp <- as.character(fit$rho[N.rho, N.gamma])
      my.label1 <- tcltk::tclVar(paste("rho: ", label.temp, sep=""))
      lab.rho <- tcltk::tklabel(fr.rho, textvariable=my.label1, width = "25")
      tcltk::tkpack(lab.rho, sc1, side = "bottom", anchor="w")

      ##insert space under the scale bar of rho 
      tcltk::tkpack(fr.space1 <- tcltk::tkframe(nt))
      lab.space1 <- tcltk::tklabel(fr.space1, textvariable=tcltk::tclVar(""))
      tcltk::tkpack(fr.space1, lab.space1)

      ##output gamma
      i <- 2
      eval(parse(text = paste("assign('slider", i, "',tcltk::tclVar(sl.defaults[i]),envir=slider.env)", sep = "")))
      tcltk::tkpack(fr.gamma <- tcltk::tkframe(nt))
      sc2 <- tcltk::tkscale(fr.gamma, from = sl.mins[i], to = sl.maxs[i], showvalue = F, resolution = sl.deltas[i], orient = "horiz", length = 500)
      tcltk::tkpack(fr.gamma, sc2)
      assign("sc2", sc2, envir = slider.env)
      eval(parse(text = paste("tcltk::tkconfigure(sc2,variable=slider", i, ")", sep = "")), envir = slider.env)
      label.temp <- as.character(fit$gamma[N.gamma])
      my.label2 <- tcltk::tclVar(paste("gamma: ", label.temp, sep=""))
      lab.gamma <- tcltk::tklabel(fr.gamma, textvariable=my.label2, width = "25")
      tcltk::tkpack(lab.gamma, sc2, side = "bottom", anchor='w')

      ##insert space under the scale bar of gamma
      tcltk::tkpack(fr.space2 <- tcltk::tkframe(nt))
      lab.space2 <- tcltk::tklabel(fr.space2, textvariable=tcltk::tclVar(""))
      tcltk::tkpack(fr.space2, lab.space2)


      ##-----------------------
      ##    output goodness of fit indices
      ##-----------------------
      tcltk::tkpack(fr.criteria <- tcltk::tkframe(nt))
      ##GFI
      label.temp <- as.character(signif(fit$GFI[N.rho, N.gamma], digits=6))
      label.GFI <- tcltk::tclVar(paste("GFI: ", label.temp, sep=""))
      lab.GFI <- tcltk::tklabel(fr.criteria, textvariable=label.GFI, width = "14")
      tcltk::tkpack(fr.criteria, lab.GFI, side="left", fill="none", anchor='w')
      ##AGFI
      label.temp <- as.character(signif(fit$AGFI[N.rho, N.gamma], digits=6))
      label.AGFI <- tcltk::tclVar(paste("AGFI: ", label.temp, sep=""))
      lab.AGFI <- tcltk::tklabel(fr.criteria, textvariable=label.AGFI, width = "14")
      tcltk::tkpack(lab.GFI, lab.AGFI, side = "left", anchor='w')
      ##AIC
      label.temp <- as.character(signif(fit$AIC[N.rho, N.gamma], digits=6))
      label.AIC <- tcltk::tclVar(paste("AIC: ", label.temp, sep=""))
      lab.AIC <- tcltk::tklabel(fr.criteria, textvariable=label.AIC, width = "14")
      tcltk::tkpack(lab.AGFI, lab.AIC, side="left", fill="none", anchor='w')
      ##BIC
      label.temp <- as.character(signif(fit$BIC[N.rho, N.gamma], digits=6))
      label.BIC <- tcltk::tclVar(paste("BIC: ", label.temp, sep=""))
      lab.BIC <- tcltk::tklabel(fr.criteria, textvariable=label.BIC, width = "14")
      tcltk::tkpack(lab.AIC, lab.BIC, side="left", anchor='w')
      ##CAIC
      label.temp <- as.character(signif(fit$CAIC[N.rho, N.gamma], digits=6))
      label.CAIC <- tcltk::tclVar(paste("CAIC: ", label.temp, sep=""))
      lab.CAIC <- tcltk::tklabel(fr.criteria, textvariable=label.CAIC, width = "14")
      tcltk::tkpack(lab.BIC, lab.CAIC, side="left", anchor='w')
      ##EBIC
      label.temp <- as.character(signif(fit$EBIC[N.rho, N.gamma], digits=6))
      label.EBIC <- tcltk::tclVar(paste("EBIC: ", label.temp, sep=""))
      lab.EBIC <- tcltk::tklabel(fr.criteria, textvariable=label.EBIC, width = "14")
      tcltk::tkpack(lab.CAIC, lab.EBIC, side="left", anchor='w')


      ##-----------------------
      ##    output factor loadings
      ##-----------------------
      ##insert space under the goodness of fit indices
      tcltk::tkpack(fr.space3 <- tcltk::tkframe(nt))
      lab.space3 <- tcltk::tklabel(fr.space3, textvariable=tcltk::tclVar(""))
      tcltk::tkpack(fr.space3, lab.space3)

      fun.output.loadings <- function (...) {
        rho <-   as.numeric(tcltk::tclvalue(get(paste("slider", 1, sep = ""), envir = slider.env)))
        gamma <- as.numeric(tcltk::tclvalue(get(paste("slider", 2, sep = ""), envir = slider.env)))
        print(fit$loadings[[gamma]][[rho]])
      }
      tcltk::tkpack(fr.button <- tcltk::tkframe(nt))
      button.loadings <- tcltk::tkbutton(fr.button, text = "Output loadings", command=fun.output.loadings)
      tcltk::tkpack(fr.criteria, fr.button, side="top")
      tcltk::tkpack(fr.button, button.loadings, side="top", anchor='e')

      my.fun <- function (...) {
        rho <-   as.numeric(tcltk::tclvalue(get(paste("slider", 1, sep = ""), envir = slider.env)))
        gamma <- as.numeric(tcltk::tclvalue(get(paste("slider", 2, sep = ""), envir = slider.env)))
        label.temp <- as.character(fit$rho[rho, gamma])
        tcltk::tclvalue(my.label1) <- paste("rho: ", label.temp, sep="")
        label.temp <- as.character(fit$gamma[gamma])
        tcltk::tclvalue(my.label2) <- paste("gamma: ", label.temp, sep="")

        ##GFI
        label.temp <- as.character(signif(fit$GFI[rho, gamma], digits=6))
        tcltk::tclvalue(label.GFI) <- paste("GFI: ", label.temp, sep="")
        ##AGFI
        label.temp <- as.character(signif(fit$AGFI[rho, gamma], digits=6))
        tcltk::tclvalue(label.AGFI) <- paste("AGFI: ", label.temp, sep="")
        ##AIC
        label.temp <- as.character(signif(fit$AIC[rho, gamma], digits=6))
        tcltk::tclvalue(label.AIC) <- paste("AIC: ", label.temp, sep="")
        ##BIC
        label.temp <- as.character(signif(fit$BIC[rho, gamma], digits=6))
        tcltk::tclvalue(label.BIC) <- paste("BIC: ", label.temp, sep="")
        ##CAIC
        label.temp <- as.character(signif(fit$CAIC[rho, gamma], digits=6))
        tcltk::tclvalue(label.CAIC) <- paste("CAIC: ", label.temp, sep="")
        ##EBIC
        label.temp <- as.character(signif(fit$EBIC[rho, gamma], digits=6))
        tcltk::tclvalue(label.EBIC) <- paste("EBIC: ", label.temp, sep="")

        loadings <- fit$loadings[[gamma]][[rho]]
        loadings <- as.matrix(loadings)
        loadings  <- t(loadings)
        loadings <- matlab::fliplr(loadings)
        image(loadings, col=col.all, zlim=c(-max.col, max.col))
      }

      tcltk::tkconfigure(sc1, command=my.fun)
      tcltk::tkconfigure(sc2, command=my.fun)

      assign("slider.values.old", sl.defaults, envir = slider.env)

      invisible(nt)
    }

    heatmap.fanc <- function (fit) {
      N.rho <- nrow(fit$rho)
      N.gamma <- ncol(fit$rho)
      slider.fanc(sl.names=c("rho", "gamma"), sl.mins=c(1, 1), sl.maxs=c(N.rho, N.gamma), sl.deltas=c(1, 1), sl.defaults=c(1, 1), fit=fit)
    }

    heatmap.fanc(x)

  } else {
    stop("Only 'path' and 'heatmap' are available for 'type'")
  }
}
