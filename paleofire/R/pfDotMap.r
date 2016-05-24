#' Produce maps of paleofire data
#' 
#' Produce map graphics representing spatial variability in charcoal data from
#' the Global Charcoal Database.
#' 
#' Takes any pfTransform object as input, and allows any set of one or more
#' time bins to be specified for plotting (one plot per bin). Time bins are 
#' specified as for pfCompositeLF (which is called by pfSimpleGrid. The extent, 
#' resolution, and projection of the desired grid are also user-specified. 
#' 
#' Results will be plotted on a regular lon/lat grid. To determine which sites
#' contribute to each grid cell value, the code searches within a specified
#' great circle distance (i.e. on the surface of the globe) around each grid
#' cell center. To avoid missing any sites, the distance is set equal to the
#' greatest distance from a grid cell center to its most distant corner, which
#' occurs at the equator where grid cells are largest. This conservative
#' approach will result in many sites falling within multiple grid boxes. At
#' all latitudes, the defined radii will overlap near the edges of the grid
#' boxes. At higher latitudes, the lon/lat grid cells are physically much
#' smaller, so overlap will be considerably greater. There are alternatives,
#' like using a grid that is irregular in terms of lon/laton, or changing the
#' area of grid cells depending on latitude. But all have their tradeoffs, and
#' this one is simple.
#' 
#' Current version produces plots of mean CHAR, number of sites per grid cell,
#' and number of grid cells contributed to by each site (due to overlapping
#' radii described above). The mean plot additionally shows points in two
#' sizes, representing those mean values whose 95"\%" confidence intervals do
#' (small dots) or do not (large dots) contain zero. Finally, a time series is
#' plotted in each figure with the current time bin highlighted.
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
#' @param fig.base.name Character sequence representing the base name for the
#' figures. Can be preceded by a path as long as all directories in the path
#' exist. One figure will be produced for each time bin, with years (and file
#' suffix) appended to the base name automatically. A value of \code{NULL}
#' (default) causes figures to be plotted to the current device in sequence.
#' @param grd.res,grd.ext Desired grid resolution and extent in degrees. If
#' \code{grd.res} is a single number, the grid will be defined with equal
#' lon/lat resolution; a two-element vector (lon,lat) can also be supplied for
#' unequal resolution. \code{grd.ext} is specified as a vector of the form
#' \code{c(min-lon,max-lon,min-lat,max-lat)}.
#' @param grd.lonlat A data frame of coordinates for every grid cell center, to
#' be used in cases where an irregular grid is desired. Columns must be named
#' 'lon' and 'lat'. If specified, grd.res and grd.ext are ignored. Note that
#' this option could have undesirable results for unusual grid definitions. In
#' particular, the maximum radius for including sites in a grid cell is always
#' calculated at the equator. For a regular lon/lat grid, this guarantees all
#' sites will be included in at least one cell, because equatorial cells are
#' largest at the equator. If an irregular grid is specified such that this is
#' not true, the maximum radius calculated could lead to sites excluded from
#' all cells. In this case a warning is printed but the function proceeds
#' anyway.
#' @param base.map Currently, either \code{'coasts'} or \code{'countries'} to
#' choose which base map (from required library \code{'rworldmap'}) to be
#' plotted as the base map for all plots. Could easily be modified to accept
#' any SpatialPolygons object.
#' @param proj4 proj.4 string representing the desired projection for plotted
#' maps. Default is unprojected. See \url{http://www.spatialreference.org} to
#' look up the string for your favorite projections.
#' @param n.boot Number of bootstrap replicates to use when creating confidence
#' intervals around each grid-cell mean. In each time bin X grid cell
#' combination, replicates consist of composite z-score values for that bin,
#' randomly sampled (with replacement) from sites within the grid cell (see
#' 'Details' for precise description of sites included in each cell). I.e., no
#' temporal bootstrapping is done here, so that bootstrap CI reflect only
#' spatial variability.
#' @param cx.minsize,cx.mult Parameters that crudely adjust plotted dot size.
#' cx.minsize defines the minimum cex applied to any point in any map, cx.mult
#' scales all points by an equivalent factor.
#' @return Plots are produced on the current device or in pdf files defined by
#' \code{fig.base.name}. In addition, a named list of useful objects is
#' returned:
#' 
#' \item{COMP}{ The binned composite generated for plotting.  } \item{bins}{
#' The list of bin endpoints.  } \item{sp.grd}{ A
#' \code{\link[sp]{SpatialPointsDataFrame-class}} object containing all the
#' grid-level statistics produced and plotted (mean influx value, bootstrap
#' confidence interval, and number of sites per grid cell).  } \item{sp.sites}{
#' A \code{\link[sp]{SpatialPointsDataFrame-class}} object representing the
#' number of grid cells influenced by each site.  } \item{plots}{ A list with
#' one element for each bin. These elements are themselves named lists of
#' trellis objects representing each of the plots produced ("mean",
#' "sitesPerCell", "cellsPerSite", "timeSeries"). Note that these objects can
#' be edited to some degree with the \code{\link[lattice]{update.trellis}}
#' function, and plotted or used in layouts as any other trellis graphics can.
#' }
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
#' ## Composite charcoal record for North America:
#' ID=pfSiteSel(id_region==c("WNA0"), l12==1 & long<(-130))
#' plot(ID)
#' 
#' ## Transform data
#' res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' ## Plot maps for 1000-yr bins spanning 3-0 kBP
#' # dev.new(width=10,height=10) # A big plot area helps. 
#' dotmap = pfDotMap( TR=res3, tarAge=seq(0,2000,1000), hw=500, grd.ext=c(-170,-80,40,80), 
#'                    cx.minsize=2,cx.mult=3)
#' summary(dotmap)
#' 
#' # Plot the mean map from the first time bin
#' # newmap = update(dotmap$plots[[1]]$mean, main="A relabeled map")
#' # newmap
#' }
#' 
pfDotMap = function(TR, tarAge, hw, binhw=0.5*mean(diff(tarAge)), 
                    fig.base.name=NULL, base.map='coasts',
                    grd.res=5, grd.ext=c(-180,180,-90,90), grd.lonlat=NULL, 
                    proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", n.boot=1000,
                    cx.minsize=0.3, cx.mult=1) {
  
  if (!requireNamespace("rworldmap", quietly = TRUE)) {
    install.packages("rworldmap")
  }

  
  # ---------------- TEST BLOCK
  # Easier to test without running the code as a function. Comment everything above here (function definition) and 
  # uncomment this code to run the whole thing as a standard script. The load() line needs to point 
  # to where a pre-made pfTransform object has been saved. You will also need to load the two additional 
  # functions at the end of this file (just run them once at the beginning of each session). 
  # And technically, don't run the '}' that closes the main function definition (though I think if you do it will 
  # run everything and just give a harmless error at the end.
#   # 
#   rm(list=ls())
#   library(lattice)
#   TR                = readRDS('/Work/Research/GPWG/GCD v3.0 Paper figures/Data/All_GCDv1.1_Transformed_v02.rds')
#   tarAge            = seq(0,2000,1000)
#   hw                = 250
#   binhw             = 500
#   fig.base.name     = '~/Desktop/'
#   base.map          = 'coasts'
#   grd.res           = 5
#   grd.ext           = c(-180,180,-90,90)
#   grd.lonlat        = NULL
#   proj4             = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
#   n.boot            = 10    # too small, but OK for a teest
#   cx.minsize        = 0.3   # minimum dot size
#   cx.mult           = 1     # multiplicative factor for scaling all dots
#   
  # ---------------- END TEST BLOCK
  
  # ----- Load base map
  countriesCoarse<-coastsCoarse<-NULL
  rm(countriesCoarse);rm(coastsCoarse)
  
  data(countriesCoarse,envir = environment(),package="rworldmap")  # A dataset in rworldmap used in the plots below
  data(coastsCoarse,envir = environment(),package="rworldmap")     # An alternative base map. Needs one fix:
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
  
  
  # ----- Get lat/lon from GCD
  cat("Retrieveing site coordinates from GCD...")
  # A little awkward, but seems good to rely on existing paleofire functions for retrieving site coordinates,
  # and as far as I can tell this is the way to do it at present.

  # Lookup info from the database
  paleofiresites=NULL; rm(paleofiresites)  
  data(paleofiresites,envir = environment())
  site.dat=paleofiresites[paleofiresites$id_site %in% TR$params$ID$id_site, ]
  
  
  # Extract coordinates, including in radians for use in the distance function
  sites.lon = site.dat$long
  sites.lonrad = sites.lon*pi/180
  sites.lat = site.dat$lat
  sites.latrad = sites.lat*pi/180
  
  # Define n.site
  n.site = nrow(site.dat)
  
  cat("done!\n")
  
  
  # ----- Define prediction grid
  if(!is.null(grd.lonlat)) {
    # If the grid is already defined via the grd.lonlat input, conservatively set grd.res
    # to the maximum lat/lon gap between adjacent rows/cols. This isn't strictly right
    # and could be way off for unusual grids. However it's fine for something like T31, etc.
    # where the grid is only slightly irregular. 
    grd.res = c(max(diff(sort(grd.lonlat$lon))), max(diff(sort(grd.lonlat$lat))))

    # Also standardize longitude to (-180,180), in case it's not already
    grd.lonlat$lon = ((grd.lonlat$lon+360) %% 360) # First put in (0,360) range
      grd.lonlat$lon[ grd.lonlat$lon>180 ] = grd.lonlat$lon[grd.lonlat$lon>180] - 360
  } else {
    # Otherwise, use grd.res and grd.ext to define a regular grid
    if(length(grd.res)==1)    # Assume equal x/y resolution if single number given
      grd.res     = rep(grd.res,2)
  
    # Find grid cell centers
    grd.lon = seq( grd.ext[1]+grd.res[1]/2, grd.ext[2], grd.res[1])
    grd.lat = seq( grd.ext[3]+grd.res[2]/2, grd.ext[4], grd.res[2])
  
    # Expand to obtain every combination of lon/lat
    grd.lonlat = expand.grid(grd.lon,grd.lat)
    names(grd.lonlat) = c("lon", "lat") # column names for convenience
  }
  
  # Now define derived variables
  n.grd = nrow(grd.lonlat)
  grd.lonlat.rad = grd.lonlat*pi/180   # Will need lat/lon in radians later  
  

  # ----- Figure out radius for including sites in cell-level stats
  #  As discussed at AGU (Marlon, Bartlein, Higuera, Kelly), to avoid missing any sites 
  # we want to search within a radius equal to the greatest distance from a grid cell 
  # center to its most distant corner. 
  #
  # For a regular lat/lon grid, this should occur at the equator, where grid cells are largest. 
  # Note that this conservative approach will result in many sites falling within multiple grid 
  # boxes--even at the equator,  the defined circles will overlap near the edges of the grid boxes. 
  # At higher latitudes, the grid cells are much smaller, so overlap will be considerably greater. 
  # There are alternatives, like using a grid that is irregular in terms of lat/lon, or changing
  # the area of grid cells depending on latitude. But all have their tradeoffs (we thought), 
  # and this one is simple. 
  #
  # For an arbitrary grid (input as grd.lonlat), it is possible that this definition of max.dist will
  # fail us (e.g. if polar grid cells are defined larger than equatorial). For now we will use 
  # it anyway, and just post a warning below if it results in any site being omitted from
  # all grid cells.
  
  # Find the max distance just discussed, i.e. the center-to-corner distance for a cell at the origin:
  max.dist = haverdist(0,0, (grd.res[1]/2)*(pi/180), (grd.res[2]/2)*(pi/180))


  # ------ Calculate distances from sites to grid cell centers
  cat("\nCalculating distances...\n")
  # Output space and progress bar definition
  dists = matrix(NA, nrow=n.grd, ncol=n.site)
  pb = txtProgressBar(0,n.grd,style=2)
  
  # Loop over every grid cell, save distance from the grid cell center to every site
  for(i in 1:n.grd) {   # i=1
    dists[i,] = haverdist(grd.lonlat.rad$lon[i], grd.lonlat.rad$lat[i],sites.lonrad,sites.latrad)
    setTxtProgressBar(pb, i)
  } 
  close(pb) # Close progress bar

  # Warning if somehow we've failed to include sites in at least one grid cell
  ind = which(apply( (dists<=max.dist), 2, sum ) < 1)
  if(length(ind>0)) {
    warning(paste0("Some sites aren't contributing to any grid cells. Check max.dist! (Sites are ", paste(ind, collapse=","), ")."))
  }

  
  # ----- Compute stats for each grid cell
  cat("\nComputing stats for each grid cell X time slice...\n")
  
  # Set aside space for some grid-cell statistics
  grd.n   = matrix(0, nrow=n.grd, ncol=n.bin)
  grd.mean  = matrix(NA, nrow=n.grd, ncol=n.bin)
  grd.lCI = matrix(NA, nrow=n.grd, ncol=n.bin)
  grd.uCI = matrix(NA, nrow=n.grd, ncol=n.bin)
  dat.n.contributions = matrix(0, nrow=n.site, ncol=n.bin) 
  
  # Now loop over each time bin X grid cell and compute desired stats.
  pb = txtProgressBar(0,n.grd*n.bin,style=2) # progress bar
  for(j in 1:n.bin) {  # i=1887; j=1
    for(i in 1:n.grd) {
      ind = which(dists[i,]<max.dist)         # all sites within range
      ind = ind[which(!is.na(CHAR[ind,j]))] # remove NA (sites that don't span bin j)
      if(length(ind)>0) {
        grd.n[i,j]    = length(ind)            # number of sites contributing
        grd.mean[i,j] = mean(CHAR[ind,j])      # mean CHAR
        CI = bootCI(CHAR[ind,j], nboot=n.boot) # CI
        grd.lCI[i,j] = CI[1] # split into lower/upper to keep variables 2-D
        grd.uCI[i,j] = CI[2]
        
        # For each site that contributed to this grid cell, increment dat.n.contributions
        dat.n.contributions[ind,j] = dat.n.contributions[ind,j] + 1
      }
      setTxtProgressBar(pb, (j-1)*n.grd+i) # update progress bar
    }
  }
  close(pb) # close progress bar

  grd.site.ind = list()
  for(i in 1:n.grd) {
      grd.site.ind[[i]] = which(dists[i,]<max.dist)         # all sites within range
  }


# ---------- Plot
  cat("\nPlotting ", n.bin, " figures...\n", sep="")
  pb = txtProgressBar(0,n.bin,style=2) # progress bar
  plotlist = spgrdlist = spsitelist = list()

  
  # ----- Loop over all bins, one plot for each
  for(j in 1:n.bin) {  #     j=1
    # ----- Setup
    #       cat(j,"...",sep="")
    setTxtProgressBar(pb, j) # update progress bar
    
    # Open file connection, if saving the plot
    if(!is.null(fig.base.name))
      pdf(paste(fig.base.name, "_", tarAge[j]-binhw, "-", tarAge[j]+binhw, "bp.pdf", sep=""), 
          width=12, height=12)
    
    # Convert stats to spatial data frames. Might not be necessary but sure makes it easy to use spplot() below.
    sp.grd = cbind(grd.lonlat, grd.n[,j], grd.mean[,j], grd.lCI[,j], grd.uCI[,j])
    names(sp.grd) = c("lon","lat","sitesPerCell","mean.CHAR","CI.lower","CI.upper")
    sp::coordinates(sp.grd) = c('lon','lat')
    sp::proj4string(sp.grd) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
    
    sp.sites = data.frame(sites.lon, sites.lat, dat.n.contributions[,j])
    names(sp.sites) = c("lon","lat","cellsPerSite")
    sp::coordinates(sp.sites) = c('lon','lat')
    sp::proj4string(sp.sites) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
    
    # Now project the data to be plotted
    sp.grd     = sp::spTransform(sp.grd, sp::CRS(proj4))
    sp.sites   = sp::spTransform(sp.sites, sp::CRS(proj4))
    
    # Get x/y lims (could add option to override later). Only purpose currently is that sites vs. grid have 
    # different bounding boxes, which is the default extent for spplot(). I.e. without setting all the same, 
    # the bottom-right plot will have different extent than other two.
    x.lim = sp::bbox(sp.grd)[1,]
    y.lim = sp::bbox(sp.grd)[2,]

    # ----- Create mean plot
      # Define colors and cut locations 
#       cols = c("#CA0020","#F4A582","#F7F7F7","#92C5DE","#0571B0") # from colorbrewer
#       cols = c("#CA0020","#F4A582",grey(0.9),"#92C5DE","#0571B0") # modified from colorbrewer
      cols = c("#0571B0","#92C5DE",grey(0.9),"#F4A582","#CA0020") # modified from colorbrewer
      cuts = seq(-2.5,2.5,by=1) # Defines range and resolution of color scale
#   cols = rev(c("#b2182b","#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac"))
#   zlim = c(-1.75,1.75)
#   cuts = c(-1.75,-1.25,-0.75,-0.25,0.25,0.75,1.25,1.75)

        
      # Determine cuts for sizing point.
        # Specify symbol sizes for the two classes
        cx.sizes = cx.mult*c(0.5,1) 
        
        # Assign symbol size based on whether CI contain 0 
        cx = ifelse(sp.grd$CI.lower>0 | sp.grd$CI.upper<0, max(cx.sizes), min(cx.sizes))
        
        # The previous line will produce NA for cells with n=1 since CI are undefined. Give these "non-significant" symbol size by default.
        cx[which(sp.grd$sitesPerCell==1)] = min(cx.sizes)
          
      # Create plot object (actually plotted later)
      mean.plot = 
        sp::spplot(sp.grd, 'mean.CHAR', xlim=x.lim, ylim=y.lim,
          cuts=cuts, colorkey=T, col.regions=cols, cex=cx, edge.col=grey(0.3), lwd=0.3,
          scales=list(draw=T), sp.layout=list("sp.lines",base.map,col=grey(0.8)),
          main=paste("Charcoal Influx z-Scores: ", tarAge[j]-binhw, "-", tarAge[j]+binhw, " BP", sep="")) 


    # ----- Plot Number of sites per grid cell
      # Generate dot sizes/colors and corresponding key. 
        # Specify scale and legend. Would be good to automate this but it's pretty tricky to produce a good general algorithm. So, for now, hard-coding symbol sizes / labels that work well for global map at 5?? resolution. 
        cuts      = c(0,1,5,10,20,1000) # Where to divide symbol sizes
        cols      = grey(0.2)   # Can be replaced by a vector if different colors are desired
        cx.legend = c("1", "2-5", "6-10", "11-20",">20") # legend text
          n.cx = length(cuts)-1   # number of bins represented
      
        # Define sizes of data points and legend entries by dividing data into bins and scale to range [cx.minsize,1]*cx.mult. This is a pretty good range for symbol 'cex' sizes, although can be modified with the cx.mult and cx.minsize arguments. 'cut(..., labels=F)' returns integer classes from 1:n.cx. These are the same sizes to use in the key. 
        cx.key = ( ((1:n.cx)-1)*(cx.mult-cx.minsize)/(n.cx-1) + cx.minsize )
        cx = cx.key[ cut(sp.grd$sitesPerCell, cuts, labels=F) ]


      # Adjust scale so that the low end corresponds to specified minimum symbol size
      ind.non0 = which(cx>0) # Don't want to change size 0 (== not plotted)
      cx.key = cx.key + cx.minsize - min(cx[ind.non0]) 
      cx[ind.non0] = cx[ind.non0] + cx.minsize - min(cx[ind.non0])


      # Create plot object (actually plotted later)
      sitesPerCell.plot = 
        sp::spplot(sp.grd, 'sitesPerCell', xlim=x.lim, ylim=y.lim, scales=list(draw=F), 
          cex=cx, cex.key=cx.key, legendEntries=cx.legend, cuts=cuts, 
          col.regions=cols, edge.col="white", lwd=0.3,
          sp.layout=list("sp.lines",base.map,col=grey(0.8)), key.space="right",
          main="Number of sites per grid cell")

 
    # ----- Plot Number of grid cells contributed per site
      # Generate dot sizes/colors and corresponding key. 
        # Specify scale and legend. Would be good to automate this but it's pretty tricky to produce a good general algorithm. So, for now, hard-coding symbol sizes / labels that work well for global map at 5?? resolution. 
        cuts   = c(0,1,2,3,4,100) # Where to divide symbol sizes
        cols      = grey(0.2)   # Can be replaced by a vector if different colors are desired
        cx.legend = c("1", "2", "3", "4",">4") # legend text
          n.cx = length(cuts)-1   # number of bins represented
      
        # Define sizes of data points and legend entries. 
        cx.key = ( ((1:n.cx)-1)*(cx.mult-cx.minsize)/(n.cx-1) + cx.minsize )
        cx = cx.key[ cut(sp.sites$cellsPerSite, cuts, labels=F) ]

      # Adjust scale so that the low end corresponds to specified minimum symbol size
      ind.non0 = which(cx>0) # Don't want to change size 0 (== not plotted)
      cx.key = cx.key + cx.minsize - min(cx[ind.non0]) 
      cx[ind.non0] = cx[ind.non0] + cx.minsize - min(cx[ind.non0])
        
      # Create plot object (actually plotted later)
      cellsPerSite.plot = 
        sp::spplot(sp.sites, 'cellsPerSite', xlim=x.lim, ylim=y.lim,
          cex=cx, cex.key=cx.key, legendEntries=cx.legend, cuts=cuts, 
          scales=list(draw=F), col.regions=cols, edge.col="white", lwd=0.3,
          sp.layout=list("sp.lines",base.map,col=grey(0.8)), key.space="right",
          main="Number of grid cells influenced by each site")

  
    # ----- Create time series plot
    if(n.bin>1) {
      timeSeries.dat = data.frame(
                age  = as.numeric(rbind(tarAge-binhw, tarAge+binhw)),
                char = rep(COMP$Result$LocFit, each=2),
                lCI  = rep(COMP$Result[,4], each=2),
                uCI  = rep(COMP$Result[,5], each=2) )


      timeSeries.plot =       
        xyplot( char~age, data=timeSeries.dat, 
          ylim=c(-0.05,0.05)*diff(range(timeSeries.dat[,3:4]))+range(timeSeries.dat[,3:4]),
          panel = function(x,y, ...) {
            ind.j = (2*j-1):(2*j)
            x.j   = x[ind.j]
            y.j   = y[ind.j]
            lCI.j = timeSeries.dat$lCI[ind.j]
            uCI.j = timeSeries.dat$uCI[ind.j]
            
            panel.polygon(c(x,rev(x)), c(timeSeries.dat$lCI, rev(timeSeries.dat$uCI)), 
              col=grey(0.8), border=FALSE)
            panel.polygon(c(x.j,rev(x.j)), c(lCI.j, rev(uCI.j)), 
              col=rgb(1,0,0,0.5), border=FALSE)

            panel.xyplot(x,y, type='l', col=1)
            panel.xyplot(x.j,y.j, type='l', col=2)
          })
      } else {
        timeSeries.plot = NULL
      }

    # ----- Produce plots
      print(mean.plot, position=c(0,0.5,1,1), more=T)
      if(n.bin>1) print(timeSeries.plot, position=c(0,0.3,1,0.5), more=T)
      print(sitesPerCell.plot, position=c(0,0,0.5,0.3), more=T)
      print(cellsPerSite.plot, position=c(0.5,0,1,0.3))


    # Close connection to external figure
    if(!is.null(fig.base.name))  dev.off()

    # ----- Store plot objects
      plotlist[[j]] = list("mean"=mean.plot,"sitesPerCell"=sitesPerCell.plot,"cellsPerSite"=cellsPerSite.plot,"timeSeries"=timeSeries.plot)
      spgrdlist[[j]] = sp.grd
      spsitelist[[j]] = sp.sites
      

  } # End loop over all bins

  # ----- Return
    output = list(COMP=COMP, tarAge=tarAge, sp.grd=spgrdlist, sp.sites=spsitelist, grd.site.ind=grd.site.ind, site.dat=site.dat, plots=plotlist)
    return(output)
  
cat("\nAll done!\n\n")
  

} # End main function definition



# ---------- Auxilliary functions
# Haversine distance function (http://en.wikipedia.org/wiki/Haversine_formula)
# R is the approximate radius of the earth, defaulting to 6371 km
# Note: at least on pair of inputs (i.e., either (lon1,lat1) or (lon2/lat2)) must represent a single location. 
# The other may be vectors representing many locations, but the function does not work to find all pairwise 
# distances between two sets (n>1) of locations. It wouldn't be difficult to code the pairwise version, 
# but I think it wouldn would require a loop so it wouldn't save any run time anyway.

haverdist = function(lon1, lat1, lon2, lat2, R=6371) {
  return( 
    2*R*asin(sqrt( (sin((lat2-lat1)/2))^2 + cos(lat1)*cos(lat2)*((sin((lon2-lon1)/2))^2) ))
  )
}

# Bootstrap CI function.
# Actually should work for various functions and CI probs but this code currently only relies 
# on the defaults (mean, 95% CI).
bootCI = function(x, nboot=1000, fun=mean, ..., probs=c(0.025,0.975)) {
  nx = length(x)
  
  if(nx<2) {
    return( rep(NA, length(probs)) )
  } else {
    xboot = matrix(x[sample(1:nx, nboot*nx, replace=T)], nrow=nboot, ncol=nx)
    return( quantile(apply(xboot, 1, fun, ...), probs=probs, na.rm=T) )
  }
}
