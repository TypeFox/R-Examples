#' Process MODIS cloud mask product files to TIF, and then extract data
#'
#' \tabular{ll}{
#' Package: \tab modiscloud\cr
#' Type: \tab Package\cr
#' Version: \tab 0.14\cr
#' Date: \tab 2013-02-08\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' This package helps the user process downloaded MODIS cloud product HDF files to TIF, and then
#' extract data. Specifically, MOD35_L2 cloud product files, and the associated MOD03 geolocation
#' files (for MODIS-TERRA); and MYD35_L2 cloud product files, and the associated MYD03 geolocation
#' files (for MODIS-AQUA).
#'
#' The package will be most effective if the user installs MRTSwath (MODIS
#' Reprojection Tool for swath products; \url{https://lpdaac.usgs.gov/tools/modis_reprojection_tool_swath}),
#' and adds the directory with the MRTSwath executable to the default R PATH by editing \code{~/.Rprofile}.
#'
#' Each MOD35_L2/MYD35_L2 file requires a corresponding MOD03/MYD03 geolocation
#' file to be successfully processed with the MRTSwath tool.
#'
#' MRTSwath is the MRT (MODIS Reprojection Tool) for the MODIS
#' level 1 and level 2 products (cloud mask is level 2, I think).
#' 
#' A few example MODIS Cloud Product files, and derived TIFs, are found in the data-only package \code{modiscdata}. 
#' These were too big to put in the main package, according to CRAN repository policies (\url{http://cran.r-project.org/web/packages/policies.html}).
#'
#' Note: This code was developed for the following publication. Please cite if used: Goldsmith, Gregory; Matzke, Nicholas J.; 
#' Dawson, Todd (2013). "The incidence and implications of clouds for cloud forest plant water relations."
#' Ecology Letters, 16(3), 207-314. DOI: \url{http://dx.doi.org/10.1111/ele.12039}
#'
#' @name modiscloud-package
#' @aliases modiscloud
#' @docType package
#' @title Process MODIS cloud mask product files to TIF
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @references
#' \url{https://lpdaac.usgs.gov/get_data/}
#' @bibliography /Dropbox/_njm/__packages/modiscloud_setup/modiscloud_refs.bib
#'   @cite NASA2001
#'   @cite Ackerman2010
#'   @cite GoldsmithMatzkeDawson2013
#' @keywords package
#' @seealso \code{\link{check_for_matching_geolocation_files}}
#' @examples
#' # Test function for checking roxygen2, roxygenize package documentation building
#' is.pseudoprime(13, 4)
#'
#' # Some MODIS files are stored in this package's "extdata/" directory 
#' # Here are some example MODIS files in modiscloud/extdata/
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' moddir = system.file("extdata/2002raw/", package="modiscdata")
#'
#' # This directory actually has MYD files (from the MODIS-AQUA platform)
#' # (*will* work with the default files stored in modiscloud/extdata/2002raw/)
#' list.files(path=moddir, pattern="MYD")
#'
#' # Check for matches (for MODIS-AQUA platform)
#' # (*will* work with the default files stored in modiscloud/extdata/2002raw/)
#' fns_df = check_for_matching_geolocation_files(moddir=moddir, modtxt="MYD35_L2", geoloctxt="MYD03", return_geoloc=FALSE, return_product=FALSE)
#'
#' }
#'
#'
#' #######################################################
#' # Run MRTSwath tool "swath2grid"
#' #######################################################
#' 
#' # Source MODIS files (both data and geolocation)
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' # install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' moddir = system.file("extdata/2002raw/", package="modiscdata")
#' 
#' # Get the matching data/geolocation file pairs
#' fns_df = check_for_matching_geolocation_files(moddir, modtxt="MYD35_L2", geoloctxt="MYD03")
#' fns_df
#' 
#' # Resulting TIF files go in this directory
#' tifsdir = getwd()
#' 
#' 
#' # Box to subset
#' ul_lat = 13
#' ul_lon = -87
#' lr_lat = 8
#' lr_lon = -82
#' 
#' for (i in 1:nrow(fns_df))
#' 	{
#' 	
#'	prmfn = write_MRTSwath_param_file(prmfn="tmpMRTparams.prm", tifsdir=tifsdir, modfn=fns_df$mod35_L2_fns[i], geoloc_fn=fns_df$mod03_fns[i], ul_lon=ul_lon, ul_lat=ul_lat, lr_lon=lr_lon, lr_lat=lr_lat)
#'	print(scan(file=prmfn, what="character", sep="\n"))
#' 	
#'	run_swath2grid(mrtpath="swath2grid", prmfn="tmpMRTparams.prm", tifsdir=tifsdir, modfn=fns_df$mod35_L2_fns[i], geoloc_fn=fns_df$mod03_fns[i], ul_lon=ul_lon, ul_lat=ul_lat, lr_lon=lr_lon, lr_lat=lr_lat)
#'	
#' 	}
#'
#' tiffns = list.files(tifsdir, pattern=".tif", full.names=TRUE)
#' tiffns
#'
#'
#' # For some unit testing etc., swath2grid may not be available.  If so, use the default TIFs:
#' if (length(tiffns) == 0)
#' 	{
#'	library(modiscdata)
#'	tifsdir = system.file("extdata/2002tif/", package="modiscdata")
#' 	tiffns = list.files(tifsdir, pattern=".tif", full.names=TRUE)
#'	}
#' 
#' #######################################################
#' # Load a TIF
#' #######################################################
#' library(rgdal)	# for readGDAL
#' 
#' # numpixels in subset
#' xdim = 538
#' ydim = 538
#' 
#' 
#' # Read the grid and the grid metadata		
#' coarsen_amount = 1
#' xdim_new = xdim / floor(coarsen_amount)
#' ydim_new = ydim / floor(coarsen_amount)
#' 		
#' fn = tiffns[1]
#' grd = readGDAL(fn, output.dim=c(ydim_new, xdim_new))
#' 
#' grdproj = CRS(proj4string(grd))
#' grdproj
#' grdbbox = attr(grd, "bbox")
#' grdbbox
#' 
#' 
#' 
#' 
#' 
#' ###########################
#' # Extract values from a particular pixel
#' ###########################
#' # Greg's field site
#' greglat = 10.2971
#' greglon = -84.79282
#' 
#' grdr = raster(grd)
#' 
#' # Input the points x (longitude), then y (latitude)
#' point_to_sample = c(greglon, greglat)
#' xycoords = adf(matrix(data=point_to_sample, nrow=1, ncol=2))
#' names(xycoords) = c("x", "y")
#' 
#' xy = SpatialPoints(coords=xycoords, proj4string=grdproj)
#' #xy = spsample(x=grd, n=10, type="random")
#' pixelval = extract(grdr, xy)
#' 
#' # Have to convert to 8-bit binary string, and reverse to get the count correct
#' # (also reverse the 2-bit strings in the MODIS Cloud Mask table)
#' pixelval = rev(t(digitsBase(pixelval, base= 2, 8)))
#' print(pixelval)
#' 
#' }
#'
#' 
#' #######################################################
#' # Load a TIF
#' #######################################################
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' # install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' tifsdir = system.file("extdata/2002tif/", package="modiscdata")
#' tiffns = list.files(tifsdir, pattern=".tif", full.names=TRUE)
#' tiffns
#' 
#' library(rgdal)	# for readGDAL
#' 
#' # numpixels in subset
#' xdim = 538
#' ydim = 538
#' 
#' 
#' # Read the grid and the grid metadata		
#' coarsen_amount = 1
#' xdim_new = xdim / floor(coarsen_amount)
#' ydim_new = ydim / floor(coarsen_amount)
#' 		
#' fn = tiffns[1]
#' grd = readGDAL(fn, output.dim=c(ydim_new, xdim_new))
#' 
#' grdproj = CRS(proj4string(grd))
#' grdproj
#' grdbbox = attr(grd, "bbox")
#' grdbbox
#' 
#' #######################################################
#' # Extract a particular bit for all the pixels in the grid
#' #######################################################
#' bitnum = 2
#' grdr_vals_bits = get_bitgrid(grd, bitnum)
#' length(grdr_vals_bits)
#' grdr_vals_bits[1:50]
#'
#' #######################################################
#' # Extract a particular pair of bits for all the pixels in the grid
#' #######################################################
#' bitnum = 2
#' grdr_vals_bitstrings = get_bitgrid_2bits(grd, bitnum)
#' length(grdr_vals_bitstrings)
#' grdr_vals_bitstrings[1:50]
#' 
#' }
#'
#' #######################################################
#' # Load some bit TIFs (TIFs with just the cloud indicators extracted)
#' # and add up the number of cloudy days, out of the total
#' # number of observation attempts
#' #######################################################
#' # Code excluded from CRAN check because it depends on modiscdata
#' \dontrun{
#' library(devtools)
#' # The modiscdata (MODIS c=cloud data=data) package is too big for CRAN (60 MB); so it is available on github:
#' # https://github.com/nmatzke/modiscdata
#' # If we can't get install_github() to work, try install_url():
#' # install_github(repo="modiscdata", username="nnmatzke")
#' # install_url(url="https://github.com/nmatzke/modiscdata/archive/master.zip")
#' library(modiscdata)
#' tifsdir = system.file("extdata/2002bit/", package="modiscdata")
#' tiffns = list.files(tifsdir, pattern=".tif", full.names=TRUE)
#' tiffns
#' 
#' library(rgdal)	# for readGDAL
#' 
#' # numpixels in subset
#' xdim = 538
#' ydim = 538
#' 
#' # Read the grid and the grid metadata		
#' coarsen_amount = 1
#' xdim_new = xdim / floor(coarsen_amount)
#' ydim_new = ydim / floor(coarsen_amount)
#' 
#' 
#' sum_nums = NULL
#' for (j in 1:length(tiffns))
#' 	{
#' 	fn = tiffns[j]
#' 	
#' 	grd = readGDAL(fn, output.dim=c(ydim_new, xdim_new))
#' 	
#' 	grdr = raster(grd)
#' 	pointscount_on_SGDF_points = coordinates(grd)
#' 	grdr_vals = extract(grdr, pointscount_on_SGDF_points)
#' 	
#' 	# Convert to 1/0 cloudy/not
#' 	data_grdr = grdr_vals
#' 	data_grdr[grdr_vals > 0] = 1
#' 	
#' 	grdr_cloudy = grdr_vals
#' 	grdr_cloudy[grdr_vals < 4] = 0
#' 	grdr_cloudy[grdr_vals == 4] = 1
#' 	
#' 	# Note: Don't run the double-commented lines unless you want to collapse different bit values.
#' 	# grdr_clear = grdr_vals
#' 	# grdr_clear[grdr_vals == 4] = 0
#' 	# grdr_clear[grdr_vals == 3] = 1
#' 	# grdr_clear[grdr_vals == 2] = 1
#' 	# grdr_clear[grdr_vals == 1] = 1
#' 	# grdr_clear[grdr_vals == 0] = 0
#' 	# 
#' 	
#' 	if (j == 1)
#' 		{
#' 		sum_cloudy = grdr_cloudy
#' 		#sum_not_cloudy = grdr_clear
#' 		sum_data = data_grdr
#' 		} else {
#' 		sum_cloudy = sum_cloudy + grdr_cloudy
#' 		sum_data = sum_data + data_grdr
#' 		}
#' 			
#' 	}
#' 
#'  
#' # Calculate percentage cloudy
#' sum_nums = sum_cloudy / sum_data
#' 
#' grd_final = numslist_to_grd(numslist=sum_nums, grd=grd, ydim_new=ydim_new, xdim_new=xdim_new)
#' 
#' # Display the image (this is just the sum of a few images)
#' image(grd_final)
#' 
#' }
#' 
 
