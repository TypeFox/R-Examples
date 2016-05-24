#' Produce simple gridded maps of paleofire data
#' 
#' Produce gridded map graphics representing spatial variability in charcoal 
#' data from the Global Charcoal Database.
#' 
#' Takes any pfTransform object as input, and allows any set of one or more
#' time bins to be specified for plotting (one plot per bin). Time bins are 
#' specified as for pfCompositeLF (which is called by pfSimpleGrid. The extent, 
#' resolution, and projection of the desired grid are also user-specified. 
#' 
#' Records are first composited, and then aggregated with other sites falling 
#' in the same grid cell according to the specified function 'fun' (defauts to 
#' mean). This is a considerably simpler approach than the distance-based spatial
#' binning used by pfDotMap, although it has its own tradeoffs (e.g. grid cells
#' are unlikely to represent equal area). 
#' 
#' A flexible bootstrapped significance test is implemented. Within each time 
#' bin X grid cell combination, composite z-score values are randomly sampled 
#' (with replacement) from sites within the grid cell. The function is applied 
#' to the sampled values. Quantiles of all bootstrap function evaluations are 
#' computed, and significance is reported if a user-specified test value is 
#' outside of these bootstrap CI. Note that bootstrap CI calculated here reflect 
#' only spatial variability, as no temporal resampling is performed. 
#' 
#' @param TR An object returned by \code{\link{pfTransform}}
#' @param tarAge Numeric, the target ages for prebinning given in years (e.g.
#' tarAge = seq(0, 10000, 20)). If unspecified the sequence is defined as
#' tarAge=seq(from=min age, to=max Age, by=median resolution).
#' @param hw Numeric, the half window width for the locfit procedure (in
#' years).
#' @param binhw Numeric, bin half width for the prebinning procedure (use the
#' same value as tarAge intervals for overlapping bins or tarAge intervals/2
#' for non-overlapping bins).
#' @param fun Function to be used for aggregating across sites.
#' @param n.boot Number of bootstrap replicates to use when creating confidence
#' intervals around each grid-cell value. 
#' @param prob.CI Vector of two quantiles to define the bootstrap CI for 
#' significance testing 
#' @param test.val Test value for bootstrap significance test.
#' @param proj4 proj.4 string representing the desired projection for plotted
#' maps. Default is unprojected. See \url{http://www.spatialreference.org} to
#' look up the string for your favorite projections.
#' @param res,ext Desired grid resolution and extent. If
#' \code{grd.res} is a single number, the grid will be defined with equal
#' x/y resolution; a two-element vector (x,y) can also be supplied for
#' unequal resolution. \code{grd.ext} is specified as a vector, matrix, or Extent
#' object, as for the function raster::extent.
#' @param fig.file.name Character sequence representing the file name for the
#' output figures. Can be preceded by a path as long as all directories in the path
#' exist. The file will be a PDF with one figure per time bin, each on a separate page.
#' @param show.plots Logical indicating whether plots will be printed to the screen.
#' @param title.text Character sequence for labeling figures. Time bin bounds will
#' be added automatically. 
#' @param cols, cuts Vectors of color specifications and values defining the plot 
#' legend. Grid-cell values will be binned by \code{cuts} and assigned the colors in
#' \code{cols}. If either are NULL, the function tries to guess at a good scheme. 
#' \code{cuts} may also be a single value specifying the number of bins.
#' @param cuts  Defines range and resolution of color scale
#' @param zlim Two-element vector representing the bounds of the color scale. Ignored
#' if \code{cuts} is fully specified, but otherwise used in defining the color bins.
#' @param base.map Currently, either \code{'coasts'} or \code{'countries'} to
#' choose which base map (from required library \code{'rworldmap'}) to be
#' plotted as the base map for all plots. Could easily be modified to accept
#' any SpatialPolygons object.
#' @param base.map.col,  Color specifications for plotting
#' the basemap.
#' @param  base.map.lwd Line width specifications for plotting
#' the basemap.
#' @seealso \code{\link[paleofire]{pfGridding}}
#' 
#' @return Plots are produced on the current device and/or in pdf files according to
#' input arguments. In addition, a named list of useful objects is
#' returned:
#' 
#' \item{COMP}{ The binned composite generated for plotting.  } 
#' \item{tarAge}{ The list of target ages used for temporal binning.  }
#' \item{sg.rast}{ A \code{\link[raster]{Raster-class}} object containing the gridded 
#' output data } 
#' \item{sg.plots}{ A list of trellis objects representing the composed plots. 
#' Note that these objects can be edited to some degree with the 
#' \code{\link[lattice]{update.trellis}} function, and plotted or used in layouts as 
#' any other trellis graphics can. }
#'
#' @author R. Kelly
#' @references Power, M., J. Marlon, N. Ortiz, P. Bartlein, S. Harrison, F.
#' Mayle, A. Ballouche, R. Bradshaw, C. Carcaillet, C. Cordova, S. Mooney, P.
#' Moreno, I. Prentice, K. Thonicke, W. Tinner, C. Whitlock, Y. Zhang, Y. Zhao,
#' A. Ali, R. Anderson, R. Beer, H. Behling, C. Briles, K. Brown, A. Brunelle,
#' M. Bush, P. Camill, G. Chu, J. Clark, D. Colombaroli, S. Connor, A. L.
#' Daniau, M. Daniels, J. Dodson, E. Doughty, M. Edwards, W. Finsinger, D.
#' Foster, J. Frechette, M. J. Gaillard, D. Gavin, E. Gobet, S. Haberle, D.
#' Hallett, P. Higuera, G. Hope, S. Horn, J. Inoue, P. Kaltenrieder, L.
#' Kennedy, Z. Kong, C. Larsen, C. Long, J. Lynch, E. Lynch, M. McGlone, S.
#' Meeks, S. Mensing, G. Meyer, T. Minckley, J. Mohr, D. Nelson, J. New, R.
#' Newnham, R. Noti, W. Oswald, J. Pierce, P. Richard, C. Rowe, M. Sanchez
#' Goni, B. Shuman, H. Takahara, J. Toney, C. Turney, D. Urrego-Sanchez, C.
#' Umbanhowar, M. Vandergoes, B. Vanniere, E. Vescovi, M. Walsh, X. Wang, N.
#' Williams, J. Wilmshurst, and J. Zhang. 2008. Changes in fire regimes since
#' the Last Glacial Maximum: an assessment based on a global synthesis and
#' analysis of charcoal data. Climate Dynamics 30:887-907.
#' @examples
#' 
#' \dontrun{
#' ID=pfSiteSel(id_region==c("WNA0"), l12==1 & long<(-130))
#' plot(ID)
#' 
#' ## Transform data
#' res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' ## Plot maps for 1000-yr bins spanning 3-0 kBP
#' # dev.new(width=10,height=10) # A big plot area helps. 
#' gridmap = pfSimpleGrid( TR=res3, tarAge=seq(0,2000,1000), hw=500, ext=c(-170,-80,40,80))
#' summary(gridmap)
#' 
#' # Plot the mean map from the first time bin
#' newmap = update(gridmap$sg.plots[[1]], main="A relabeled map")
#' newmap
#' }
#' 

pfSimpleGrid = function(TR, tarAge, hw, binhw=0.5*mean(diff(tarAge)), fun=mean,
                    n.boot=0, prob.CI=c(0.025,0.975), test.val=0, 
                    proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                    res=5, ext=c(-180,180,-90,90), 
                    fig.file.name=NULL, show.plots=TRUE, title.text="",
                    cols = NULL, cuts = NULL, zlim=NULL,
                    base.map='coasts', base.map.col=grey(0.7), base.map.lwd=0.5
                    
) {
 
# ---------------- TEST BLOCK
# 
# # Easier to test without running the code as a function. Comment everything above here (function definition) and 
# # uncomment this code to run the whole thing as a standard script. The load() line needs to point 
# # to where a pre-made pfTransform object has been saved. You will also need to load the two additional 
# # functions at the end of this file (just run them once at the beginning of each session). 
# # And technically, don't run the '}' that closes the main function definition (though I think if you do it will 
# # run everything and just give a harmless error at the end.
# rm(list=ls())
# 
# load('/Work/Research/GPWG/pfDotMap working/All_GPCD_Transformed_v2.rdata')
# bins              = seq(-500,21500,1000)
# show.plots        = T
# fig.file.name     = '/Work/Research/GPWG/pfDotMap working/SimpleGrid Maps vTest7.pdf' # NULL for no saved plot
# base.map          = 'coasts'
# res           = 5
# ext           = c(-180,180,-90,90)
# proj4         = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
# n.boot            = 100    # too small, but OK for a teest
# zlim              = c(-2.5,2.5)
# fun               = mean
# cols              = c("#0571B0","#92C5DE",grey(0.9),"#F4A582","#CA0020")  # Mean
# cuts = seq(-2.5,2.5,by=1) # Defines range and resolution of color scale
# 
# # zlim              = NULL
# # fun               = "count"
# # cols              = brewer.pal(9, 'YlOrRd')  # others
# # cuts              = NULL
# 
# prob.CI           = c(0.025,0.975)
# test.val          = 0
# 
# base.map.col      = grey(0.7)
# base.map.lwd      = 0.5
# title.text        = ""
# ---------------- END TEST BLOCK
  
  # ----- Constants
  # proj4 string for the GCD data
    proj4.GCD = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "

  # ----- Load base map
  countriesCoarse<-coastsCoarse<-NULL
  rm(countriesCoarse);rm(coastsCoarse)
  
  data(countriesCoarse,envir = environment(), package="rworldmap")  # A dataset in rworldmap used in the plots below
  data(coastsCoarse,envir = environment(), package="rworldmap")     # An alternative base map. Needs one fix:
  ind = which(coastsCoarse@lines[[94]]@Lines[[1]]@coords[,1] > 180)
  coastsCoarse@lines[[94]]@Lines[[1]]@coords[ind,1] = 180
  
  # Select and transform base map based on inputs (only two options currently)
  if(base.map=='countries') {
    base.map = countriesCoarse
  } else {
    base.map = coastsCoarse
  }
  base.map = sp::spTransform(base.map, sp::CRS(proj4))
  
  
  # ----- Create composite
  if(class(TR)=="pfTransform") {
    cat("Creating composite...")
    # Run pfComposite. Not interested in the composite, but this will do the binning for us. 
    # (I've checked "by hand" and it is accurate and efficient.) 
    # No need for bootstraps--we're only going to use the actual composite for now. Will bootstrap by site later. 
    COMP = pfCompositeLF(TR, tarAge=tarAge, hw=hw, binhw=binhw, nboot=1)
    CHAR = t(COMP$BinnedData)
    n.bin = length(tarAge)
    cat("done!\n")
  } else {
    stop("Input PF must be a pfTransform object.")
  }
  
  
  # ----- Get site locations from GCD
  cat("Retrieveing site coordinates from GCD...")
  # Lookup info from the database
    paleofiresites=NULL; rm(paleofiresites)  
    data(paleofiresites,envir = environment())
    site.dat=paleofiresites[paleofiresites$id_site %in% TR$params$ID$id_site, ]

  # Convert data to spatial data into a SpatialPointsDataFrame
    # Add binned + transformed charcoal data to original data frame. 
      dat = data.frame(lon=site.dat$long, lat=site.dat$lat, CHAR)

  # Specify the columns that are coordinates, and then the coordinate ref. system. This converts the data frame to class 'SpatialPointsDataFrame', a class in the 'sp' library. 
    sp::coordinates(dat) = c('lon','lat')
    sp::proj4string(dat) = proj4.GCD
    
  # NA locations
    na.mat = is.na(CHAR)
  cat("done!\n")
  
  
  # ----- Define prediction grid
  grd = raster(extent(ext))     # Extent
    res(grd) = res              # Resolution
    sp::proj4string(grd) = proj4    # Projection


  # ----- Calculate desired stat for each time bin
  sg.rast = rasterize(dat, grd, field=names(dat), fun=fun, na.rm=T)
    # If don't specify "names(dat)" then get an extra layer "ID". 

    # Convert to stack (will be layer if only one time bin)
    sg.rast = stack(sg.rast)

  # ----- Bootstrap CI
  if(n.boot>0) {
    # Bootstrap CI function.
    # Actually should work for various functions and CI probs but this code currently only relies 
    # on the defaults (mean, 95% CI).
    bootTest = function(x, ...) {
      x = na.omit(x)
      nx = length(x)
  
      if(nx<2) {
        return( NA )
      } else {
        xboot = matrix(x[sample(1:nx, n.boot*nx, replace=T)], nrow=n.boot, ncol=nx)
        quants = quantile(apply(xboot, 1, fun, ...), probs=prob.CI)
        return(ifelse( quants[1]>max(test.val) | quants[2]<min(test.val), 1, 0))
#         return(quants)
      }
    }
  
    isSingle = function(x, ...) ifelse(length(na.omit(x))==1,1,NA)
#     boot.lCI = function(x,...) bootTest(x, prob=min(prob.CI))
#     boot.uCI = function(x,...) bootTest(x, prob=max(prob.CI))

    # Even though it doesn't need to be done in a loop, do so in order to print progress
    sg.rast.lCI = sg.rast.uCI = sg.rast.1 = sg.rast.sig = NA*sg.rast
    timer = Sys.time()
    for(i in 1:n.bin) {
      cat(i,"...")
#       sg.rast.lCI[[i]] = rasterize(dat, grd, field=names(dat)[i], fun=boot.lCI)
#       sg.rast.uCI[[i]] = rasterize(dat, grd, field=names(dat)[i], fun=boot.uCI)
#       sg.rast.sig[[i]] = (sg.rast.lCI[[i]]>=max(test.val)) | (sg.rast.uCI[[i]]<=min(test.val))
      sg.rast.sig[[i]] = rasterize(dat, grd, field=names(dat)[i], fun=bootTest)
      sg.rast.1[[i]] = rasterize(dat, grd, field=names(dat)[i], fun=isSingle)
      print(summary(as.factor(values(sg.rast.sig[[i]]))))
    }
    sg.rast.sig[ sg.rast.sig<1 ] = NA
    Sys.time()-timer
  }


  # ----- Handle default colors
  if(is.null(zlim)) zlim = range(values(sg.rast), na.rm=T)
  if(length(cuts)==1) cuts = seq(min(zlim), max(zlim), length=cuts+1)
  
  if(is.null(cols) & is.null(cuts)) {
    cuts = seq(min(zlim), max(zlim), length=10)
    cols = rev(topo.colors(length(cuts)-1))
  } else if(is.null(cols)) {
    cols = rev(topo.colors(length(cuts)-1))
  } else if(is.null(cuts)) {
    cuts = seq(min(zlim), max(zlim), length=length(cols)+1)
  }


  # ----- Make Plots
  cat("\nPlotting ", n.bin, " figures...", sep="")
  
  # Can call spplots on the whole raster brick, which creates a separate panel for each layer (in this case, one per time bin). E.g.

    # sg.plots = spplot(sg.rast, at=cuts, col.regions=cols, scales=list(draw=T),layout=c(1,1))
    # 
    # # Then can update the "proper" way for trellis objects. E.g.:
    #   sg.plots = update(sg.plots, 
    #               strip=strip.custom(factor.levels=
    #                 paste0(title.text, bins[1:n.bin], "-", bins[2:(n.bin+1)], " BP")),
    #               sp.layout=list("sp.lines",base.map,col=base.map.col,lwd=base.map.lwd))
      
  # Pretty snazzy. But a little confusing if like me you're not very familiar with lattice graphics. Also, not the way we set up pfDotMaps. 
  
  # So for now, making a separate trellis plot for each time bin, storing them in a list. 
  sg.plots = list()
  base.map.layout = list("sp.lines",base.map,col=base.map.col,lwd=base.map.lwd)
  pb = txtProgressBar(0,n.bin,style=2) # progress bar
  for(j in 1:n.bin) {  #     j=1
    setTxtProgressBar(pb, j) # update progress bar

    site.dots.layout = list("sp.points", dat[ !na.mat[,j],j], col=1, pch=16, cex=0.3)

    sg.plots[[j]] = sp::spplot(sg.rast[[j]], at=cuts, col.regions=cols, scales=list(draw=T),
      main = as.list(paste0(title.text, tarAge[j]-binhw, "-", tarAge[j]+binhw, " BP")))
    
    # Should be better way, but can't figure out how to add multiple extra layers in separate steps. This works...
    if(n.boot>0) {
      if(all(is.na(values(sg.rast.sig[[j]])))) {
        sg.sig.layout = list("sp.text", "(no sig vals)", loc=c(0,0))
      } else {
        sg.sig.layout = list("sp.polygons", rasterToPolygons(sg.rast.sig[[j]]), first=F)
      }
      
      if(all(is.na(values(sg.rast.1[[j]])))) {
        sg.1.layout = list("sp.text", "(no single-site cells)", loc=c(0,0))
      } else {
        sg.1.layout = list("sp.polygons", rasterToPolygons(sg.rast.1[[j]]), col=2, first=F)
      }
      sg.plots[[j]] = update(sg.plots[[j]], 
        sp.layout=list(base.map.layout,sg.1.layout,sg.sig.layout,site.dots.layout))
    } else {
      sg.plots[[j]] = update(sg.plots[[j]], sp.layout=list(base.map.layout, site.dots.layout))
    }
  }


  # ----- Print plots to screen and/or file
  # Plot to screen
    if(show.plots) print(sg.plots)

  # Save to file
    if(!is.null(fig.file.name)) {
      pdf(fig.file.name, width=12, height=12)
      print(sg.plots)
      dev.off()
    }
  cat("done!\n")


  # ----- Return
    output = list(COMP=COMP, tarAge=tarAge, sg.rast=sg.rast, sg.plots=sg.plots, site.dat=site.dat)
    return(output)
  
cat("\nAll done!\n\n")
} # End main function definition

