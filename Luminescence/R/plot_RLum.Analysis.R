#' Plot function for an RLum.Analysis S4 class object
#'
#' The function provides a standardised plot output for curve data of an
#' RLum.Analysis S4 class object
#'
#' The function produces a multiple plot output. A file output is recommended
#' (e.g., \code{\link{pdf}}).
#'
#' \bold{curve.transformation}\cr
#'
#' This argument allows transforming continuous wave (CW) curves to pseudo
#' (linear) modulated curves. For the transformation, the functions of the
#' package are used. Currently, it is not possible to pass further arguments to
#' the transformation functions. The argument works only for \code{ltype}
#' \code{OSL} and \code{IRSL}.\cr
#'
#' Please note: The curve transformation within this functions works roughly,
#' i.e. every IRSL or OSL curve is transformed, without considerung whether it
#' is measured with the PMT or not! However, for a fast look it might be
#' helpful.\cr
#'
#'
#' @param object \code{\linkS4class{RLum.Analysis}} (\bold{required}): S4
#' object of class \code{RLum.Analysis}
#'
#' @param subset named \code{\link{list}} (optional): subsets elements for plotting. The
#' arguments in the named \code{\link{list}} will be directly passed to the function \code{\link{get_RLum}}
#' (e.g., \code{subset = list(curveType = "measured")})
#'
#' @param nrows \code{\link{integer}} (optional): sets number of rows for
#' plot output, if nothing is set the function tries to find a value.
#'
#' @param ncols \code{\link{integer}} (optional): sets number of columns
#' for plot output, if nothing is set the function tries to find a value.
#'
#' @param abline \code{\link{list}} (optional): allows to set similar ablines
#' in each plot. This option uses the function \code{\link{do.call}}, meaning
#' that every argument in the \code{list} has to be provided as \code{list},
#' e.g. \code{abline = list(list(v = 120), list(v = 350))} produces two
#' vertical ablines: One at 150 and another one at 350. Within the call all
#' arguments supported by \code{\link{abline}} are fully supported,
#'
#' @param combine \code{\link{logical}} (with default): allows to combine all
#' \code{\linkS4class{RLum.Data.Curve}} objects in one single plot.
#'
#' @param curve.transformation \code{\link{character}} (optional): allows
#' transforming CW-OSL and CW-IRSL curves to pseudo-LM curves via
#' transformation functions. Allowed values are: \code{CW2pLM}, \code{CW2pLMi},
#' \code{CW2pHMi} and \code{CW2pPMi}. See details.
#'
#' @param plot.single \code{\link{logical}} (with default): global par settings are
#' considered, normally this should end in one plot per page
#'
#' @param \dots further arguments and graphical parameters will be passed to
#' the \code{plot} function. Supported arguments: \code{main} (can be provided as
#' vector for \code{combine = TRUE}), \code{mtext},
#' \code{log}, \code{lwd}, \code{lty} \code{type}, \code{pch}, \code{col},
#' \code{norm}, \code{ylim}, \code{xlab} ... and for \code{combine = TRUE}
#' also: \code{xlim}, \code{ylab}, \code{sub}, \code{legend.text},
#' \code{legend.pos} (typical plus 'outside'), \code{legend.col}
#'
#' @return Returns multiple plots.
#'
#' @note Not all arguments available for \code{\link{plot}} will be passed!
#' Only plotting of \code{RLum.Data.Curve} and \code{RLum.Data.Spectrum}
#' objects are currently supported.
#'
#' @section Function version: 0.2.9
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso \code{\link{plot}}, \code{\link{plot_RLum}},
#' \code{\link{plot_RLum.Data.Curve}}
#'
#' @references #
#'
#' @keywords aplot
#'
#' @examples
#'
#'
#' ###load data
#' data(ExampleData.BINfileData, envir = environment())
#'
#' ##convert values for position 1
#' temp <- Risoe.BINfileData2RLum.Analysis(CWOSL.SAR.Data, pos=1)
#'
#' ##plot all values
#' plot_RLum.Analysis(temp)
#'
#' ##plot (combine) TL curves in one plot
#' temp.sel <- get_RLum(temp, recordType = "TL", drop = FALSE)
#' plot_RLum.Analysis(temp.sel, combine = TRUE, norm = TRUE, main = "TL combined")
#'
#'
#' @export
plot_RLum.Analysis <- function(
  object,
  subset,
  nrows,
  ncols,
  abline,
  combine = FALSE,
  curve.transformation,
  plot.single = FALSE,
  ...
){

  # Integrity check ----------------------------------------------------------------------------

  ##check if object is of class RLum.Data.Curve
  if(is(object,"RLum.Analysis") == FALSE){

    stop("[plot_RLum.Analysis()] Input object is not of type 'RLum.Analysis'")

  }

  # Make selection if wanted  -------------------------------------------------------------------
  if(!missing(subset)){

    ##check whether the user set the drop option ...
    subset <- subset[!sapply(names(subset), function(x){"drop" %in% x})]
    object <- do.call(get_RLum, c(object = object, subset, drop = FALSE))

  }


  ##deal with addition arguments
  extraArgs <- list(...)

  ##main
  main <- if ("main" %in% names(extraArgs)) {

      ##main - allow to set different mains
      if(length(extraArgs$main) == 1 | length(extraArgs$main) < length(object)){
        rep(x =  extraArgs$main, length(object))

      } else{
        extraArgs$main

      }
    } else{
      NULL
    }

  ##mtext
  mtext <- if("mtext" %in% names(extraArgs)) {extraArgs$mtext} else
  {NULL}

  ##log
  log <- if("log" %in% names(extraArgs)) {extraArgs$log} else
  {""}

  ##lwd
  lwd <- if("lwd" %in% names(extraArgs)) {extraArgs$lwd} else
  {1}

  ##lty
  lty <- if("lty" %in% names(extraArgs)) {extraArgs$lty} else
  {1}

  ##type
  type <- if("type" %in% names(extraArgs)) {extraArgs$type} else
  {"l"}

  ##xlim
  xlim <- if("xlim" %in% names(extraArgs)) {extraArgs$xlim} else
  {NULL}

  ##ylim
  ylim <- if("ylim" %in% names(extraArgs)) {extraArgs$ylim} else
  {NULL}

  ##pch
  pch <- if("pch" %in% names(extraArgs)) {extraArgs$pch} else
  {1}

  ##col
  col <- if("col" %in% names(extraArgs)) {extraArgs$col} else
  {"black"}

  ##norm (for RLum.Data.Curve)
  norm <- if("norm" %in% names(extraArgs)) {extraArgs$norm} else
  {FALSE}

  ##cex
  cex <- if("cex" %in% names(extraArgs)) {extraArgs$cex} else
  {1}

  ##try to find optimal parameters, this is however, a little bit stupid, but
  ##better than without any

  if(combine){
    n.plots <- length(unique(as.character(structure_RLum(object)$recordType)))

  }else{
    n.plots <- length_RLum(object)

  }

  if (missing(ncols) | missing(nrows)) {
    if (missing(ncols) & !missing(nrows)) {
      if (n.plots  == 1) {
        ncols <- 1

      } else{
        ncols <- 2

      }

    }
    else if (!missing(ncols) & missing(nrows)) {
      if (n.plots  == 1) {
        nrows <- 1

      }
      else if (n.plots  > 1 & n.plots <= 4) {
        nrows <- 2

      } else{
        nrows <- 3

      }


    } else{
      if (n.plots  == 1) {
        nrows <- 1
        ncols <- 1

      }
      else if (n.plots  > 1 & n.plots  <= 2) {
        nrows <- 1
        ncols <- 2

      } else if (n.plots  > 2 & n.plots <= 4) {
        nrows <- 2
        ncols <- 2

      }
      else{
        nrows <- 3
        ncols <- 2

      }

    }

  }




  # Plotting ------------------------------------------------------------------

  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ##(1) NORMAL (combine == FALSE)
  if(combine == FALSE || length(object@records) == 1){

    ##show warning message
    if(combine == TRUE & length(object@records) == 1){

      warning("Nothing to combine, object contains a single curve.")

    }

    ##grep RLum.Data.Curve or RLum.Data.Spectrum objects
    temp <- lapply(1:length(object@records), function(x){

      if(is(object@records[[x]], "RLum.Data.Curve") == TRUE ||
           is(object@records[[x]], "RLum.Data.Spectrum") == TRUE){

        object@records[[x]]

      }})

    ##calculate number of pages for mtext
    if(length(temp)%%(nrows*ncols)>0){

      n.pages <- round(length(temp)/(nrows*ncols), digits=0)+1

    }else{

      n.pages <- length(temp)/(nrows*ncols)

    }

    ##set par
    par.default <- par("mfrow")
    if(!plot.single){par(mfrow=c(nrows,ncols))}

    ##plot curves
    for(i in 1:length(temp)){

      if(is(temp[[i]], "RLum.Data.Curve") == TRUE){

        ##set curve transformation if wanted
        if((grepl("IRSL", temp[[i]]@recordType) | grepl("OSL", temp[[i]]@recordType)) &
             !missing(curve.transformation)){

          if(curve.transformation=="CW2pLM"){

            temp[[i]] <- CW2pLM(temp[[i]])

          }else if(curve.transformation=="CW2pLMi"){

            temp[[i]] <- CW2pLMi(temp[[i]])

          }else if(curve.transformation=="CW2pHMi"){

            temp[[i]]<- CW2pHMi(temp[[i]])

          }else if(curve.transformation=="CW2pPMi"){

            temp[[i]] <- CW2pPMi(temp[[i]])

          }else{

            warning("Function for 'curve.transformation' is unknown. No transformation is performed.")

          }

        }


        ##check xlim and ylim values and adjust if where necessary
        ##xlim
        if (!is.null(xlim)) {
          xlim.set <- xlim
          if (xlim[1] < min(temp[[i]]@data[,1])) {
            xlim.set[1] <- min(temp[[i]]@data[,1])
          }
          if (xlim[2] > max(temp[[i]]@data[,1])) {
            xlim.set[2] <- max(temp[[i]]@data[,1])
          }

        }else{
          xlim.set <- xlim

        }

        ##ylim
        if (!is.null(ylim)) {
          ylim.set <- ylim
          if (ylim[1] < min(temp[[i]]@data[,2])) {
            ylim.set[1] <- min(temp[[i]]@data[,2])
          }
          if (ylim[2] > max(temp[[i]]@data[,2])) {
            ylim.set[2] <- max(temp[[i]]@data[,2])
          }

        }else{
          ylim.set <- ylim

        }

        plot_RLum.Data.Curve(temp[[i]],
                             col = if(unique(col) != "black"){col} else{
                               if(grepl("IRSL", temp[[i]]@recordType) == TRUE){"red"} else
                                 if(grepl("OSL", temp[[i]]@recordType) == TRUE){"blue"} else
                                 {col}
                             },
                             mtext = paste("#",i,sep=""),
                             par.local = FALSE,
                             main = if(is.null(main)){temp[[i]]@recordType}else{main[i]},
                             log = log,
                             lwd = lwd,
                             type = type,
                             lty = lty,
                             xlim = xlim.set,
                             ylim = ylim.set,
                             pch = pch,
                             norm = norm,
                             cex = cex,
                             ...)

        ##add abline
        if(!missing(abline)){

          for(k in 1:length(abline)){

            do.call("abline", abline[[k]])

          }

        }


      } else if(is(temp[[i]], "RLum.Data.Spectrum") == TRUE) {

        plot_RLum.Data.Spectrum(temp[[i]],

                                mtext = paste("#",i,sep=""),
                                par.local = FALSE,
                                main = if(main==""){temp[[i]]@recordType}else{main[i]})

      }

      if(i%%(nrows*ncols)==0){
        mtext(mtext, outer = TRUE, side=3, line=-2)
      }

    }#end for loop


    ##reset par
    if(!plot.single){par(mfrow = par.default)}

  }else{

    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##(2) NORMAL (combine == TRUE)

    ##(1) check RLum objects in the set
    object.list <- get_RLum(object)


    sapply(1:length(object.list), function(x){

      if(is(object.list[[x]])[1] != "RLum.Data.Curve"){

        stop("[plot_RLum.Analysis()] Using 'combine' is limited to 'RLum.Data.Curve' objects.")

      }

    })


    ##account for different curve types, combine similar
    temp.object.structure  <- structure_RLum(object)
    temp.recordType <- as.character(unique(temp.object.structure$recordType))


    ##change graphic settings
    if(!plot.single){
      par.default <- par()[c("cex", "mfrow")]

      if(!missing(ncols) & !missing(nrows)){
        par(mfrow = c(nrows, ncols))

      }


      ##this 2nd par request is needed as seeting mfrow resets the par settings ... this might
      ##not be wanted
      par(cex = cex)

    }else{
      par.default <- par()[c("cex")]
      par(cex = cex)

    }

    ##(2) PLOT values

    for(k in 1:length(temp.recordType)) {

      ##main 2
      main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
      {paste0(temp.recordType[[k]], " combined")}

      ###get type of curves
      temp.object <-
        get_RLum(object, recordType = temp.recordType[k], drop = FALSE)

      ##get structure
      object.structure  <- structure_RLum(temp.object)


      ##now get the real list object (note the argument recursive = FALSE)
      object.list <-
        get_RLum(object, recordType = temp.recordType[k], recursive = FALSE)

      ##prevent problems for non set argument
      if (missing(curve.transformation)) {
        curve.transformation <- "None"
      }

      ##transform values to data.frame and norm values
      temp.data.list <- lapply(1:length(object.list), function(x) {
        ##set curve transformation if wanted

        if (grepl("IRSL", object.list[[x]]@recordType) |
              grepl("OSL", object.list[[x]]@recordType)) {
          if (curve.transformation == "CW2pLM") {
            object.list[[x]] <- CW2pLM(object.list[[x]])

          }else if (curve.transformation == "CW2pLMi") {
            object.list[[x]] <- CW2pLMi(object.list[[x]])

          }else if (curve.transformation == "CW2pHMi") {
            object.list[[x]] <- CW2pHMi(object.list[[x]])

          }else if (curve.transformation == "CW2pPMi") {
            object.list[[x]] <- CW2pPMi(object.list[[x]])

          }

        }


        temp.data <- as(object.list[[x]], "data.frame")

        ##normalise curves if argument has been set
        if (norm == TRUE) {
          temp.data[,2] <- temp.data[,2] / max(temp.data[,2])

        }

        return(temp.data)

      })

      ##get some extra arguments, here new, as the values have to considered differently
      sub <- if ("sub" %in% names(extraArgs)) {
        extraArgs$sub
      } else
      {
        ""
      }


      ##xlab
      xlab <- if ("xlab" %in% names(extraArgs)) {
        extraArgs$xlab
      } else
      {
        "x"
      }

      ##ylab
      ylab <- if ("ylab" %in% names(extraArgs)) {
        extraArgs$ylab
      } else
      {
        "y"
      }

      ##xlim
      xlim <- if ("xlim" %in% names(extraArgs)) {
        extraArgs$xlim
      } else
      {
        c(min(object.structure$x.min), max(object.structure$x.max))
      }

      ##ylim
      ylim <- if ("ylim" %in% names(extraArgs)) {
        extraArgs$ylim
      } else
      {
        temp.ylim  <- t(sapply(1:length(temp.data.list), function(x) {
          temp.data <- temp.data.list[[x]]
          range(temp.data[temp.data[,1] >= min(xlim) &
                            temp.data[,1] <= max(xlim),2])

        }))

        c(min(temp.ylim), max(temp.ylim))

      }

      ##col (again)
      col <- if ("col" %in% names(extraArgs)) {
        extraArgs$col
      } else
      {
        get("col", pos = .LuminescenceEnv)
      }

      ##if length of provided colours is < the number of objects, just one colour is supported
      if (length(col) < length(object.list)) {
        col <- rep_len(col, length(object.list))

      }

      ##col (again)
      lty <- if ("lty" %in% names(extraArgs)) {
        extraArgs$lty
      } else
      {
        1
      }

      ##if length of provided lty values is < the number of objects, just the first supported
      if (length(lty) < length(object.list)) {
        lty <- rep(lty[1], times = length(object.list))

      }

      ##legend.text
      legend.text <-
        if ("legend.text" %in% names(extraArgs)) {
          extraArgs$legend.text
        } else
        {
          paste("Curve", 1:length(object.list))
        }

      ##legend.col
      legend.col <-
        if ("legend.col" %in% names(extraArgs)) {
          extraArgs$legend.col
        } else
        {
          NULL
        }


      ##legend.pos
      legend.pos <-
        if ("legend.pos" %in% names(extraArgs)) {
          extraArgs$legend.pos
        } else
        {
          "topright"
        }

      if (legend.pos == "outside") {
        par.default.outside <- par()[c("mar", "xpd")]
        par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
      }

      ##main - allow to set different mains
      if(length(main) == 1 | length(main) < length(temp.recordType)){
        main <- rep(x = main, length(temp.recordType))

      }

      ##open plot area
      plot(
        NA,NA,
        xlim = xlim,
        ylim = ylim,
        main = main[k],
        xlab = xlab,
        ylab = ylab,
        log = log,
        sub = sub
      )

      ##loop over all records
      for (i in 1:length(object.list)) {
        lines(temp.data.list[[i]],
              col = col[i],
              lty = lty,
              lwd = lwd)

      }


      ##mtext
      mtext(mtext, side = 3, cex = .8 * cex)

      ##legend
      legend(
        x = ifelse(legend.pos == "outside", par()$usr[2], legend.pos),
        y = ifelse(legend.pos == "outside", par()$usr[4], NULL),
        legend = legend.text,
        lwd = lwd,
        lty = lty,
        col = if(is.null(legend.col)){col[1:length(object.list)]}else{legend.col},
        bty = "n",
        cex = 0.8 * cex
      )


    }

    ##reset graphic settings
    if (exists("par.default.outside")) {
      par(par.default.outside)
      rm(par.default.outside)
    }
    par(par.default)
    rm(par.default)

  }

}
