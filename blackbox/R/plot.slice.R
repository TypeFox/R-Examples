plotSlice <- function(xvar, yvar, fixedlist) { ## fixedlist typically a rosglobal subset, potentially in log scale
  ##type="p" 3D plot; "c" contour plot; "b" a two panel figure with
  ##perspective and contour plots; "I" image plot with legend strip;
  ##"C" image plot with contours overlaid.
  ## NB les plots 'pastel' avec les axes mal legend['e]s sont ceux completement hors du convex hull. Je laisse tel quel
  ## ===================
  ## grids
  fittedNames <- blackbox.getOption("fittedNames")
  FONKgNames <- blackbox.getOption("FONKgNames")
  plotOptions <- blackbox.getOption("plotOptions")
  #x11bool <- (blackbox.getOption("interactiveGraphics") && (blackbox.getOption("fittedparamnbr")<4))
  FixedParam <- names(fixedlist)
  EssaiGridList <- list()
  for (st in FixedParam) EssaiGridList[st] <- as.list(fixedlist[st]) ##, original, potentially logscale, values
  EssaiGridList[xvar] <- list(gridfn(varname=xvar, varnameS=c(xvar, yvar))) ## handles logscale
  EssaiGridList[yvar] <- list(gridfn(varname=yvar, varnameS=c(xvar, yvar)))
  EssaiGridList <- EssaiGridList[fittedNames] ## simply reorder EssaiGridList elements according to fittedNames order
  ## plot main label
  FixedValue <- unlist(fixedlist)
  logs <- islogscale(FixedParam) ## vector of T/F of the length of FixedParam
  FixedValue[logs] <- exp(FixedValue[logs]) ## UNLOGS values
  main <- prettyPlotMain(FixedValue)
  if (inherits(main,"expression")) { ## FR modif 10/2015
    main <- parse(text=paste("\"logL for\"",main,sep="~")) ## 02/2016 sep="~" !!
  } else main <- paste("Likelihood for ", main, sep="")
  if("3Dplots" %innc% plotOptions) { ## %in% -> %innc% 03/2014; and FR->... is this the most appropriate test ?
    ##Color 3-D plot of the likelihood surface
    makeplottypes(xvar, yvar, types=c("p"), #x11bool=x11bool,
                  ## 'main' title will show UNLOG'ed values
                  main=main,
                  grid.list= EssaiGridList)
  }
  ## color 2D contour plot ("b" is both drape and contour plots on the same graphic -- mais pas implem dans makeplottypes)
  makeplottypes(xvar, yvar, types=c("C"), #x11bool=x11bool,
                ## 'main' title will show UNLOG'ed values
                main=main,
                grid.list= EssaiGridList)
  if ("BWPlots" %innc% plotOptions) {
    ## Black & White 2D contour plot
    makeplottypes(xvar, yvar, types="C", # x11bool=x11bool,
                  main=main,
                  grid.list=EssaiGridList,
                  bw=TRUE)## bw is surfaceKrig argument
  }
} ## end def plotSlice
