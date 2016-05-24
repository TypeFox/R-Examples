#' Plot function for an RLum.Data.Spectrum S4 class object
#'
#' The function provides a standardised plot output for spectrum data of an
#' RLum.Data.Spectrum S4 class object
#'
#' \bold{Matrix structure} \cr (cf. \code{\linkS4class{RLum.Data.Spectrum}})
#'
#' \itemize{ \item \code{rows} (x-values): wavelengths/channels (xlim, xlab)
#' \item \code{columns} (y-values): time/temperature (ylim, ylab) \item
#' \code{cells} (z-values): count values (zlim, zlab) }
#'
#' \emph{Note: This nomenclature is valid for all plot types of this
#' function!}\cr
#'
#' \bold{Nomenclature for value limiting}
#'
#' \code{xlim}: Limits values along the wavelength axis\cr \code{ylim}: Limits
#' values along the time/temperature axis\cr \code{zlim}: Limits values along
#' the count value axis\cr
#'
#' \bold{Energy axis re-calculation}
#'
#' If the argument \code{xaxis.energy = TRUE} is chosen, instead intensity vs.
#' wavelength the spectrum is plotted as intensiyt vs. energy. Therefore the
#' entire spectrum is re-recaluated (e.g., Appendix 4 in Blasse and Grabmeier,
#' 1994):
#'
#' The intensity of the spectrum (z-values) is re-calcualted using the
#' following equation:
#'
#' \deqn{\phi_{E} = \phi_{\lambda} * \lambda^2 / (hc)}
#'
#' with \eqn{\phi_{E}} the intensity per interval of energy \eqn{E} (eV),
#' \eqn{\phi_{\lambda}} the intensity per interval of wavelength \eqn{\lambda}
#' (nm) and \eqn{h} (eV/s) the Planck constant and \eqn{c} (m/s) the velocity
#' of light.
#'
#' For transforming the wavelength axis (x-values) the equation
#'
#' \deqn{E = hc/\lambda}
#'
#' is used. For further details please see the cited the literature.\cr
#'
#' \bold{Details on the plot functions}
#'
#' Spectrum is visualised as 3D or 2D plot. Both plot types are based on
#' internal R plot functions. \cr
#'
#' \bold{\code{plot.type = "persp"}}
#'
#' Arguments that will be passed to \code{\link{persp}}: \itemize{ \item
#' \code{shade}: default is \code{0.4} \item \code{phi}: default is \code{15}
#' \item \code{theta}: default is \code{-30} \item \code{expand}: default is
#' \code{1} \item \code{ticktype}: default is \code{detailed}, \code{r}: default is \code{10}}
#'
#' \emph{Note: Further parameters can be adjusted via \code{par}. For example
#' to set the background transparent and reduce the thickness of the lines use:
#' \code{par(bg = NA, lwd = 0.7)} previous the function call.}
#'
#' \bold{\code{plot.type = "single"}}\cr
#'
#' Per frame a single curve is returned. Frames are time or temperature
#' steps.\cr
#'
#' \bold{\code{plot.type = "multiple.lines"}}\cr
#'
#' All frames plotted in one frame.\cr
#'
#' \bold{\code{plot.type = "transect"}}\cr
#'
#' Depending on the selected wavelength/channel range a transect over the
#' time/temperature (y-axis) will be plotted along the wavelength/channels
#' (x-axis). If the range contains more than one channel, values (z-values) are
#' summed up. To select a transect use the \code{xlim} argument, e.g.
#' \code{xlim = c(300,310)} plot along the summed up count values of channel
#' 300 to 310.\cr
#'
#' \bold{Further arguments that will be passed (depending on the plot type)}
#'
#' \code{xlab}, \code{ylab}, \code{zlab}, \code{xlim}, \code{ylim},
#' \code{zlim}, \code{main}, \code{mtext}, \code{pch}, \code{type}, \code{col},
#' \code{border}, \code{box} \code{lwd}, \code{bty} \cr
#'
#' @param object \code{\linkS4class{RLum.Data.Spectrum}} (\bold{required}): S4
#' object of class \code{RLum.Data.Spectrum}
#' @param par.local \code{\link{logical}} (with default): use local graphical
#' parameters for plotting, e.g. the plot is shown in one column and one row.
#' If \code{par.local = FALSE} global parameters are inherited.
#' @param plot.type \code{\link{character}} (with default): plot type, for
#' 3D-plot use \code{persp}, or \code{persp3d}, for a 2D-plot \code{contour},
#' \code{single} or \code{multiple.lines} (along the time or temperature axis)
#' or \code{transect} (along the wavelength axis) \cr
#'
#' Note: The use of \code{\link[rgl]{persp3d}} will produce a dynamic 3D surface plot on
#' the screen.
#'
#' @param optical.wavelength.colours \code{\link{logical}} (with default): use
#' optical wavelength colour palette. Note: For this, the spectrum range is
#' limited: \code{c(350,750)}. Own colours can be set with the argument
#' \code{col}.
#'
#' @param bg.channels \code{\link{vector}} (optional): defines channel for
#' background subtraction If a vector is provided the mean of the channels is
#' used for subtraction. Note: Background subtraction is applied prior to
#' channel binning
#'
#' @param bin.rows \code{\link{integer}} (with defaul): allow summing-up
#' wavelength channels (horizontal binning), e.g. \code{bin.rows = 2} two
#' channels are summed up
#'
#' @param bin.cols \code{\link{integer}} (with default): allow summing-up
#' channel counts (vertical binning) for plotting, e.g. \code{bin.cols = 2} two
#' channels are summed up
#'
#' @param rug \code{\link{logical}} (with default): enables or disables colour
#' rug. Currently only implemented for plot type \code{multiple.lines} and
#' \code{single}
#'
#' @param xaxis.energy \code{\link{logical}} (with default): enables or
#' disables energy instead of wavelength axis. Note: This option means not only
#' simnply redrawing the axis, insteadly the spectrum in terms of intensity is
#' recalculated, s. details.
#'
#' @param legend.text \code{\link{character}} (with default): possiblity to
#' provide own legend text. This argument is only considered for plot types
#' providing a legend, e.g. \code{plot.type="transect"}
#'
#' @param \dots further arguments and graphical parameters that will be passed
#' to the \code{plot} function.
#'
#' @return Returns a plot.
#'
#' @note Not all additional arguments (\code{...}) will be passed similarly!
#'
#' @section Function version: 0.4.2
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso \code{\linkS4class{RLum.Data.Spectrum}}, \code{\link{plot}},
#' \code{\link{plot_RLum}}, \code{\link{persp}}, \code{\link[rgl]{persp3d}},
#' \code{\link{contour}}
#'
#' @references Blasse, G., Grabmaier, B.C., 1994. Luminescent Materials.
#' Springer.
#'
#' @keywords aplot
#'
#' @examples
#'
#'
#' ##load example data
#' data(ExampleData.XSYG, envir = environment())
#'
#' ##(1)plot simple spectrum (2D) - contour
#' plot_RLum.Data.Spectrum(TL.Spectrum,
#'                         plot.type="contour",
#'                         xlim = c(310,750),
#'                         ylim = c(0,300),
#'                         bin.rows=10,
#'                         bin.cols = 1)
#'
#' ##(2) plot spectrum (3D)
#' plot_RLum.Data.Spectrum(TL.Spectrum,
#'                         plot.type="persp",
#'                         xlim = c(310,750),
#'                         ylim = c(0,100),
#'                         bin.rows=10,
#'                         bin.cols = 1)
#'
#' ##(3) plot multiple lines (2D) - multiple.lines (with ylim)
#' plot_RLum.Data.Spectrum(TL.Spectrum,
#'                         plot.type="multiple.lines",
#'                         xlim = c(310,750),
#'                         ylim = c(0,100),
#'                         bin.rows=10,
#'                         bin.cols = 1)
#'
#' \dontrun{
#'  ##(4) plot real 3d spectrum using rgl
#'  plot_RLum.Data.Spectrum(TL.Spectrum, plot.type="persp3d",
#'  xlim = c(310,750), ylim = c(0,300), bin.rows=10,
#'  bin.cols = 1)
#' }
#'
#' @export
plot_RLum.Data.Spectrum <- function(
  object,
  par.local = TRUE,
  plot.type = "contour",
  optical.wavelength.colours = TRUE,
  bg.channels,
  bin.rows = 1,
  bin.cols = 1,
  rug = TRUE,
  xaxis.energy = FALSE,
  legend.text,
  ...
){


  # Integrity check -----------------------------------------------------------

  ##check if object is of class RLum.Data.Spectrum
  if(class(object) != "RLum.Data.Spectrum"){

    stop("[plot_RLum.Data.Spectrum()] Input object is not of type RLum.Data.Spectrum")

  }

  ##XSYG
  ##check for curveDescripter
  if("curveDescripter" %in% names(object@info) == TRUE){

    temp.lab <- strsplit(object@info$curveDescripter, split = ";")[[1]]
    xlab <- if(xaxis.energy == FALSE){
      temp.lab[2]}else{"Energy [eV]"}
    ylab <- temp.lab[1]
    zlab <- temp.lab[3]

  }else{

    xlab <- if(xaxis.energy == FALSE){
      "Row values [a.u.]"}else{"Energy [eV]"}
    ylab <- "Column values [a.u.]"
    zlab <- "Cell values [a.u.]"

  }

  # Do energy axis conversion -------------------------------------------------------------------
  if (xaxis.energy) {
    temp.object.data <- sapply(1:ncol(object@data), function(x) {
      object@data[,x] * x ^ 2 / (4.13566733e-015 * 299792458e+09)
    })

    ##preserve column and rownames
    colnames(temp.object.data) <- colnames(object@data)
    rownames(temp.object.data) <-
      4.13566733e-015 * 299792458e+09 / as.numeric(rownames(object@data))

    ##write back to original data
    object@data <-
      temp.object.data[order(as.numeric(rownames(temp.object.data))),]

  }



  ##deal with addition arguments
  extraArgs <- list(...)

  main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
  {"RLum.Data.Spectrum"}

  zlab <- if("zlab" %in% names(extraArgs)) {extraArgs$zlab} else
  {ifelse(plot.type == "multiple.lines", ylab, zlab)}

  xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
  {xlab}

  ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
  {ifelse(plot.type == "single" | plot.type == "multiple.lines",
          "Luminescence [cts/channel]", ylab)}

  xlim <- if("xlim" %in% names(extraArgs)) {extraArgs$xlim} else
  {c(min(as.numeric(rownames(object@data))),
     max(as.numeric(rownames(object@data))))}

  ylim <- if("ylim" %in% names(extraArgs)) {extraArgs$ylim} else
  {c(min(as.numeric(colnames(object@data))),
     max(as.numeric(colnames(object@data))))}

  #for zlim see below

  mtext <- if("mtext" %in% names(extraArgs)) {extraArgs$mtext} else
  {""}

  cex <- if("cex" %in% names(extraArgs)) {extraArgs$cex} else
  {1}

  phi <- if("phi" %in% names(extraArgs)) {extraArgs$phi} else
  {15}

  theta <- if("theta" %in% names(extraArgs)) {extraArgs$theta} else
  {-30}

  r <- if("r" %in% names(extraArgs)) {extraArgs$r} else
  {10}

  shade <- if("shade" %in% names(extraArgs)) {extraArgs$shade} else
  {0.4}

  expand <- if("expand" %in% names(extraArgs)) {extraArgs$expand} else
  {1}

  border <- if("border" %in% names(extraArgs)) {extraArgs$border} else
  {NULL}

  box <- if("box" %in% names(extraArgs)) {extraArgs$box} else
  {TRUE}

  ticktype <- if("ticktype" %in% names(extraArgs)) {extraArgs$ticktype} else
  {"detailed"}

  log<- if("log" %in% names(extraArgs)) {extraArgs$log} else
  {""}

  type<- if("type" %in% names(extraArgs)) {extraArgs$type} else
  {"l"}

  pch<- if("pch" %in% names(extraArgs)) {extraArgs$pch} else
  {1}

  lwd<- if("lwd" %in% names(extraArgs)) {extraArgs$lwd} else
  {1}

  bty <- if("bty" %in% names(extraArgs)) {extraArgs$bty} else
  {NULL}

  sub<- if("sub" %in% names(extraArgs)) {extraArgs$sub} else
  {""}


  # prepare values for plot ---------------------------------------------------
  temp.xyz <- get_RLum(object)

  ##check for the case of a single column matrix
  if(ncol(temp.xyz)>1){

    ##reduce for xlim
    temp.xyz <- temp.xyz[as.numeric(rownames(temp.xyz)) >= xlim[1] &
                           as.numeric(rownames(temp.xyz)) <= xlim[2],]

    ##reduce for ylim
    temp.xyz <- temp.xyz[, as.numeric(colnames(temp.xyz)) >= ylim[1] &
                           as.numeric(colnames(temp.xyz)) <= ylim[2]]

  }

  ## wavelength
  x <- as.numeric(rownames(temp.xyz))

  ## time/temp
  y <- as.numeric(colnames(temp.xyz))


  # Background subtraction ---------------------------------------------------

  if(missing(bg.channels) == FALSE){

    if(length(bg.channels) > 1){

      temp.bg.signal <- rowMeans(temp.xyz[,bg.channels])
      temp.xyz <- temp.xyz[,1:ncol(temp.xyz)] - temp.bg.signal

    }else{

      temp.xyz <- temp.xyz[,1:ncol(temp.xyz)] - temp.xyz[,bg.channels]
      temp.xyz <- ifelse(temp.xyz < 0, mean(temp.xyz[,bg.channels]), temp.xyz)

    }

    ##set values < 0 to 0
    temp.xyz <- ifelse(temp.xyz < 0, mean(temp.xyz[,bg.channels[1]]), temp.xyz)

  }


  # Channel binning ---------------------------------------------------------

  if(missing(bin.rows) == FALSE){

    ##calculate n.rows
    n.rows <- nrow(temp.xyz)

    ##modulo operation for the number of groups
    bin.group.rest <- n.rows%%bin.rows

    ##define groups for binning
    bin.group <- rep(1:(n.rows/bin.rows), 1, each = bin.rows)

    ##add last group
    bin.group <- c(bin.group, rep(n.rows/bin.rows + 1, 1, each = bin.group.rest))

    ##sum up rows
    temp.xyz <- rowsum(temp.xyz, bin.group)

    ##correct labeling
    x <- x[seq(1, n.rows, bin.rows)]

    ## to avoid odd plots remove last group if bin.rows is not a multiple
    ## of the row number
    if(bin.group.rest != 0){

      temp.xyz <- temp.xyz[-nrow(temp.xyz),]
      x <- x[-length(x)]

      warning("Last wavelength channel has been removed due to binning.")

    }


    rm(bin.group.rest)

  }


  if(missing(bin.cols) == FALSE){

    ##calculate n.cols
    n.cols <- ncol(temp.xyz)

    ##check for validity
    if(bin.cols > n.cols){

      bin.cols <- n.cols

      warning("bin.cols > the number of columns. Value reduced to number of cols.")

    }

    ##modulo operation for the number of groups
    bin.group.rest <- n.cols%%bin.cols

    ##define groups for binning
    bin.group <- rep(1:(n.cols/bin.cols), 1, each = bin.cols)

    ##add last group
    bin.group <- c(bin.group, rep(n.cols/bin.cols + 1, 1, each = bin.group.rest))

    ##sum up cols
    temp.xyz <- rowsum(t(temp.xyz), bin.group)
    temp.xyz <- t(temp.xyz)

    ##correct labeling
    y <- y[seq(1, n.cols, bin.cols)]

    ## to avoid odd plots remove last group if bin.cols is not a multiple
    ## of the col number
    if(bin.group.rest != 0){

      temp.xyz <- temp.xyz[,-ncol(temp.xyz)]
      y <- y[-length(y)]

      warning("Last count channel has been removed due to column binning.")

    }

  }

  ##check for zlim
  zlim <- if("zlim" %in% names(extraArgs)) {extraArgs$zlim} else
  {range(temp.xyz)}


  # set color values --------------------------------------------------------

  if("col" %in% names(extraArgs) == FALSE | plot.type == "single" | plot.type == "multiple.lines"){

    if(optical.wavelength.colours == TRUE | rug == TRUE){

      ##make different colour palette for energy valuesw
      if (xaxis.energy) {
        col.violet <- c(2.76, ifelse(max(xlim) <= 4.13, max(xlim), 4.13))
        col.blue <- c(2.52, 2.76)
        col.green <- c(2.18, 2.52)
        col.yellow <- c(2.10, 2.18)
        col.orange <- c(2.00, 2.10)
        col.red <- c(1.57, 2.00)
        col.infrared <-
          c(1.55, ifelse(min(xlim) >= 1.55, min(xlim), 1.57))


        #set colour palette
        col <- unlist(sapply(1:length(x), function(i){

          if(x[i] >= col.violet[1] & x[i] < col.violet[2]){"#EE82EE"}
          else if(x[i] >= col.blue[1] & x[i] < col.blue[2]){"#0000FF"}
          else if(x[i] >= col.green[1] & x[i] < col.green[2]){"#00FF00"}
          else if(x[i] >= col.yellow[1] & x[i] < col.yellow[2]){"#FFFF00"}
          else if(x[i] >= col.orange[1] & x[i] < col.orange[2]){"#FFA500"}
          else if(x[i] >= col.red[1] & x[i] < col.red[2]){"#FF0000"}
          else if(x[i] <= col.infrared[2]){"#BEBEBE"}

        }))


      }else{
        col.violet <- c(ifelse(min(xlim) <= 300, min(xlim), 300),450)
        col.blue <- c(450,495)
        col.green <- c(495,570)
        col.yellow <- c(570,590)
        col.orange <- c(590,620)
        col.red <- c(620,790)
        col.infrared <-
          c(790, ifelse(max(xlim) >= 800, max(xlim), 800))


        #set colour palette
        col <- unlist(sapply(1:length(x), function(i){

          if(x[i] >= col.violet[1] & x[i] < col.violet[2]){"#EE82EE"}
          else if(x[i] >= col.blue[1] & x[i] < col.blue[2]){"#0000FF"}
          else if(x[i] >= col.green[1] & x[i] < col.green[2]){"#00FF00"}
          else if(x[i] >= col.yellow[1] & x[i] < col.yellow[2]){"#FFFF00"}
          else if(x[i] >= col.orange[1] & x[i] < col.orange[2]){"#FFA500"}
          else if(x[i] >= col.red[1] & x[i] < col.red[2]){"#FF0000"}
          else if(x[i] >= col.infrared[1]){"#BEBEBE"}

        }))


      }



      ##find unique colours
      col.unique <- unique(col)

      ##if only one colour value, then skip gradient calculation as it causes
      ## an error

      if(length(col.unique) > 1){

        ##set colour function for replacement
        colfunc <- colorRampPalette(col.unique)

        ##get index for colour values to be cut from the current palette
        col.unique.index <-
          sapply(1:length(col.unique), function(i) {
            max(which(col == col.unique[i]))

          })


        ##remove last index (no colour gradient needed), for energy axis use the first value
        col.unique.index <- col.unique.index[-length(col.unique.index)]


        ##set borders for colour gradient recalculation
        col.unique.index.min <- col.unique.index - (50)/bin.rows
        col.unique.index.max <- col.unique.index + (50)/bin.rows

        ##set negative values to the lowest index
        col.unique.index.min[col.unique.index.min<=0] <- 1



        ##build up new index sequence (might be better)
        col.gradient.index <- as.vector(unlist((
          sapply(1:length(col.unique.index.min), function(j){

            seq(col.unique.index.min[j],col.unique.index.max[j], by = 1)

          }))))


        ##generate colour ramp and replace values
        col.new <- colfunc(length(col.gradient.index))
        col[col.gradient.index] <- col.new

        ##correct for overcharged colour values (causes zebra colour pattern)
        if (diff(c(length(col), nrow(temp.xyz))) < 0) {
          col <- col[1:c(length(col) - diff(c(length(col), nrow(temp.xyz))))]

        }else if(diff(c(length(col), nrow(temp.xyz))) > 0){
          col <- col[1:c(length(col) + diff(c(length(col), nrow(temp.xyz))))]


        }


      }


    }else{

      col <- "black"

    }

  }else{

    col <- extraArgs$col

  }


  # Do log scaling if needed -------------------------------------------------

  ##x
  if(grepl("x", log)==TRUE){x <- log10(x)}

  ##y
  if(grepl("y", log)==TRUE){y <- log10(y)}

  ##z
  if(grepl("z", log)==TRUE){temp.xyz <- log10(temp.xyz)}


  # PLOT --------------------------------------------------------------------

  ##par setting for possible combination with plot method for RLum.Analysis objects
  if(par.local == TRUE){par(mfrow=c(1,1), cex = cex)}

  ##rest plot type for 1 column matrix
  if(ncol(temp.xyz) == 1){
    plot.type = "single"
    warning("[plot_RLum.Data.Spectrum()] Single column matrix: plot.type has been automatically reset to 'single'")
  }

  if(plot.type == "persp3d" && ncol(temp.xyz) > 1){

    ## ==========================================================================#
    ##perspective plot 3D screen (package rgl)
    ## ==========================================================================#

    ##check whether rgl is available
    ##code snippet taken from
    ##http://r-pkgs.had.co.nz/description.html
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("[plot_RLum.Data.Spectrum()] Package 'rgl' needed for this plot type. Please install it.",
           call. = FALSE)
    }

    rgl::persp3d(x, y, temp.xyz,
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            zlim = zlim,
            col = col,
            main = main)

  }else if(plot.type == "persp" && ncol(temp.xyz) > 1){
    ## ==========================================================================#
    ##perspective plot
    ## ==========================================================================#

    persp(x, y, temp.xyz,
          shade = shade,
          phi = phi,
          theta = theta,
          xlab = xlab,
          ylab = ylab,
          zlab = zlab,
          zlim = zlim,
          scale = TRUE,
          col = col[1:(length(col)-1)], ##needed due to recycling of the colours
          main = main,
          expand = expand,
          border = border,
          box = box,
          r = r,
          ticktype = ticktype)


    ##plot additional mtext
    mtext(mtext, side = 3, cex = cex*0.8)


  }else if(plot.type == "contour" && ncol(temp.xyz) > 1) {
    ## ==========================================================================#
    ##contour plot
    ## ==========================================================================#
    contour(x,y,temp.xyz,
            xlab = xlab,
            ylab = ylab,
            main = main,
            col = "black"
    )

    ##plot additional mtext
    mtext(mtext, side = 3, cex = cex*0.8)


  } else if(plot.type == "single") {
    ## ==========================================================================#
    ## single plot
    ## ==========================================================================#

    col.rug <- col

    col<- if("col" %in% names(extraArgs)) {extraArgs$col} else
    {"black"}



    for(i in 1:length(y)){

      if("zlim" %in% names(extraArgs) == FALSE){zlim <- range(temp.xyz[,i])}

      plot(x, temp.xyz[,i],
           xlab = xlab,
           ylab = ylab,
           main = main,
           xlim = xlim,
           ylim = zlim,
           col = col,
           sub = paste(
             "(frame ",i, " | ",
             ifelse(i==1,
                    paste("0.0 :", round(y[i], digits = 1)),
                    paste(round(y[i-1], digits = 1),":",
                          round(y[i], digits =1))),")",
             sep = ""),
           type = type,
           pch = pch)

      if(rug == TRUE){
        ##rug als continous polygons
        for(i in 1:length(x)){
          polygon(x = c(x[i],x[i+1],x[i+1],x[i]),
                  y = c(min(zlim),min(zlim), par("usr")[3], par("usr")[3]),
                  border = col.rug[i], col = col.rug[i])
        }
      }

    }

    ##plot additional mtext
    mtext(mtext, side = 3, cex = cex*0.8)


  }else if(plot.type == "multiple.lines" && ncol(temp.xyz) > 1) {
    ## ========================================================================#
    ## multiple.lines plot
    ## ========================================================================#

    col.rug <- col

    col<- if("col" %in% names(extraArgs)) {extraArgs$col} else
    {"black"}

    ##change graphic settings
    par.default <- par()[c("mfrow", "mar", "xpd")]
    par(mfrow = c(1,1), mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)

    ##grep zlim
    if("zlim" %in% names(extraArgs) == FALSE){zlim <- range(temp.xyz)}

    ##open plot area
    plot(NA, NA,
         xlab = xlab,
         ylab = ylab,
         main = main,
         xlim = xlim,
         ylim = zlim,
         sub = sub,
         bty = bty)

    if(rug == TRUE){
      ##rug als continous polygons
      for(i in 1:length(x)){
        polygon(x = c(x[i],x[i+1],x[i+1],x[i]),
                y = c(min(zlim),min(zlim), par("usr")[3], par("usr")[3]),
                border = col.rug[i], col = col.rug[i])
      }
    }

    ##add lines
    for(i in 1:length(y)){

      lines(x,
            temp.xyz[,i],
            lty = i,
            lwd = lwd,
            type = type,
            col = col)
    }

    ##for missing values - legend.text
    if(missing(legend.text)){

      legend.text <- as.character(paste(round(y,digits=1), zlab))

    }

    ##legend
    legend(x = par()$usr[2],
           y = par()$usr[4],

           legend = legend.text,

           lwd= lwd,
           lty = 1:length(y),
           bty = "n",
           cex = 0.6*cex)

    ##plot additional mtext
    mtext(mtext, side = 3, cex = cex*0.8)

    ##reset graphic settings
    par(par.default)
    rm(par.default)

  }else if(plot.type == "transect" && ncol(temp.xyz) > 1) {
    ## ========================================================================#
    ## transect plot
    ## ========================================================================#

    ##sum up rows (column sum)
    temp.xyz <- colSums(temp.xyz)

    ##consider differences within the arguments
    #check for zlim
    zlim <- if("zlim" %in% names(extraArgs)) {extraArgs$zlim} else
    {c(0,max(temp.xyz))}

    #check for zlim
    zlab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
    {paste("Counts [1/summed channels]")}

    plot(y, temp.xyz,
         xlab = ylab,
         ylab = zlab,
         main = main,
         xlim = ylim,
         ylim = zlim,
         col = col,
         sub = paste("(channel range: ", min(xlim), " : ", max(xlim), ")", sep=""),
         type = type,
         pch = pch)

    ##plot additional mtext
    mtext(mtext, side = 3, cex = cex*0.8)


  }else{

    stop("[plot_RLum.Data.Spectrum()] Unknown plot type.")

  }

}
