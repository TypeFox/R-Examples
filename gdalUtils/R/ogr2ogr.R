#' ogr2ogr
#' 
#' R wrapper for ogr2ogr: converts simple features data between file formats
#' 
#' @param src_datasource_name Character. Input vector file.
#' @param dst_datasource_name Character. Output vector file.
#' @param layer Character. Layer to use.
#' @param f Character. output file format name (default is ESRI Shapefile), some possible values are: "ESRI Shapefile", "TIGER", "MapInfo File", "GML", "PostgreSQL"
#' @param append Logical. Append to existing layer instead of creating new
#' @param overwrite Logical. Delete the output layer and recreate it empty.
#' @param update Logical. Open existing output datasource in update mode rather than trying to create a new one
#' @param select Character. Comma-delimited list of fields from input layer to copy to the new layer. A field is skipped if mentioned previously in the list even if the input layer has duplicate field names. (Defaults to all; any field is skipped if a subsequent field with same name is found.) Starting with OGR 2.0, geometry fields can also be specified in the list.
#' @param progress Logical. (starting with GDAL 1.7.0) Display progress on terminal. Only works if input layers have the "fast feature count" capability.
#' @param sql Character. SQL statement to execute. The resulting table/layer will be saved to the output.
#' @param dialect Character. SQL dialect. In some cases can be used to use (unoptimized) OGR SQL instead of the native SQL of an RDBMS by passing OGRSQL. Starting with GDAL 1.10, the "SQLITE" dialect can also be used with any datasource.
#' @param where Character. Attribute query (like SQL WHERE).
#' @param skipfailures Logical. Continue after a failure, skipping the failed feature.
#' @param spat Numeric. c(xmin,ymin,xmax,ymax) spatial query extents. Only features whose geometry intersects the extents will be selected. The geometries will not be clipped unless -clipsrc is specified
#' @param spat_srs Character. srs_def. (OGR >= 2.0) Override spatial filter SRS.
#' @param geomfield Character. (OGR >= 1.11) Name of the geometry field on which the spatial filter operates on.
#' @param dsco Character. Dataset creation option (format specific).
#' @param lco Character. Layer creation option (format specific).
#' @param nln Character. Assign an alternate name to the new layer.
#' @param nlt Character. Define the geometry type for the created layer. One of NONE, GEOMETRY, POINT, LINESTRING, POLYGON, GEOMETRYCOLLECTION, MULTIPOINT, MULTIPOLYGON or MULTILINESTRING. Add "25D" to the name to get 2.5D versions. Starting with GDAL 1.10, PROMOTE_TO_MULTI can be used to automatically promote layers that mix polygon or multipolygons to multipolygons, and layers that mix linestrings or multilinestrings to multilinestrings. Can be usefull when converting shapefiles to PostGIS (and other target drivers) that implements strict checks for geometry type.
#' @param dim Numeric. (starting with GDAL 1.10) Force the coordinate dimension to val (valid values are 2 or 3). This affects both the layer geometry type, and feature geometries. Starting with GDAL 2.0, the value can be set to "layer_dim" to instruct feature geometries to be promoted to the coordinate dimension declared by the layer.
#' @param a_srs Character. Assign an output SRS.
#' @param t_srs Character. Reproject/transform to this SRS on output.
#' @param s_srs Character. Override source SRS.
#' @param preserve_fid Logical. Use the FID of the source features instead of letting the output driver to automatically assign a new one.
#' @param fid Character. If provided, only the feature with this feature id will be reported. Operates exclusive of the spatial or attribute queries. Note: if you want to select several features based on their feature id, you can also use the fact the 'fid' is a special field recognized by OGR SQL. So, '-where "fid in (1,3,5)"' would select features 1, 3 and 5.

#' @param oo Character. "NAME=VALUE". (starting with GDAL 2.0) Input dataset open option (format specific).
#' @param doo Character. "NAME=VALUE". (starting with GDAL 2.0) Destination dataset open option (format specific), only valid in -update mode.
#' @param gt Numeric. group n features per transaction (default 200). Increase the value for better performance when writing into DBMS drivers that have transaction support.
#' @param ds_transaction Logical. (starting with GDAL 2.0) Force the use of a dataset level transaction (for drivers that support such mechanism), especially for drivers such as FileGDB that only support dataset level transaction in emulation mode.
#' @param clipsrc Character.  [xmin ymin xmax ymax]|WKT|datasource|spat_extent: (starting with GDAL 1.7.0) clip geometries to the specified bounding box (expressed in source SRS), WKT geometry (POLYGON or MULTIPOLYGON), from a datasource or to the spatial extent of the -spat option if you use the spat_extent keyword. When specifying a datasource, you will generally want to use it in combination of the -clipsrclayer, -clipsrcwhere or -clipsrcsql options
#' @param clipsrcsql Character. Select desired geometries using an SQL query instead.
#' @param clipsrclayer Character. Select the named layer from the source clip datasource.
#' @param clipsrcwhere Character. Restrict desired geometries based on attribute query.
#' @param clipdst Character. (starting with GDAL 1.7.0) clip geometries after reprojection to the specified bounding box (expressed in dest SRS), WKT geometry (POLYGON or MULTIPOLYGON) or from a datasource. When specifying a datasource, you will generally want to use it in combination of the -clipdstlayer, -clipdstwhere or -clipdstsql options
#' @param clipdstsql Character. Select desired geometries using an SQL query instead.
#' @param clipdstlayer Character. Select the named layer from the destination clip datasource.
#' @param clipdstwhere Character. Restrict desired geometries based on attribute query.
#' @param wrapdateline Logical. (starting with GDAL 1.7.0) split geometries crossing the dateline meridian (long. = +/- 180deg).
#' @param datelineoffset Logical. (starting with GDAL 1.10) offset from dateline in degrees (default long. = +/- 10deg, geometries within 170deg to -170deg will be splited)
#' @param simplify Numeric. (starting with GDAL 1.9.0) distance tolerance for simplification. Note: the algorithm used preserves topology per feature, in particular for polygon geometries, but not for a whole layer.
#' @param segmentize Numeric. (starting with GDAL 1.6.0) maximum distance between 2 nodes. Used to create intermediate points
#' @param fieldTypeToString Character. (starting with GDAL 1.7.0) converts any field of the specified type to a field of type string in the destination layer. Valid types are : Integer, Real, String, Date, Time, DateTime, Binary, IntegerList, RealList, StringList. Special value All can be used to convert all fields to strings. This is an alternate way to using the CAST operator of OGR SQL, that may avoid typing a long SQL query.
#' @param mapFieldType Character. srctype|All=dsttype,... (starting with GDAL 2.0) converts any field of the specified type to another type. Valid types are : Integer, Integer64, Real, String, Date, Time, DateTime, Binary, IntegerList, Integer64List, RealList, StringList. Types can also include subtype between parenthesis, such as Integer(Boolean), Real(Float32), ... Special value All can be used to convert all fields to another type. This is an alternate way to using the CAST operator of OGR SQL, that may avoid typing a long SQL query. This is a generalization of -fieldTypeToString. Note that this does not influence the field types used by the source driver, and is only an afterwards conversion.
#' @param unsetFieldWidth Logical. (starting with GDAL 2.0) set field width and precision to 0.		
#' @param splitlistfields Logical. (starting with GDAL 1.8.0) split fields of type StringList, RealList or IntegerList into as many fields of type String, Real or Integer as necessary.
#' @param maxsubfields Numeric. To be combined with -splitlistfields to limit the number of subfields created for each split field.
#' @param explodecollections Logical. (starting with GDAL 1.8.0) produce one feature for each geometry in any kind of geometry collection in the source file.
#' @param zfield Character. (starting with GDAL 1.8.0) Uses the specified field to fill the Z coordinate of geometries.
#' @param gcp Numeric. c(ungeoref_x,ungeoref_y,georef_x georef_y,elevation) (starting with GDAL 1.10.0) Add the indicated ground control point. This option may be provided multiple times to provide a set of GCPs.
#' @param order Numeric. (starting with GDAL 1.10.0) order of polynomial used for warping (1 to 3). The default is to select a polynomial order based on the number of GCPs.
#' @param tps Logical. (starting with GDAL 1.10.0) Force use of thin plate spline transformer based on available GCPs.
#' @param fieldmap Character. (starting with GDAL 1.10.0) Specifies the list of field indexes to be copied from the source to the destination. The (n)th value specified in the list is the index of the field in the target layer definition in which the n(th) field of the source layer must be copied. Index count starts at zero. There must be exactly as many values in the list as the count of the fields in the source layer. We can use the 'identity' setting to specify that the fields should be transferred by using the same order. This setting should be used along with the -append setting.
#' @param addfields Logical. (starting with GDAL 2.0) This is a specialized version of -append. Contrary to -append, -addfields has the effect of adding, to existing target layers, the new fields found in source layers. This option is usefull when merging files that have non-strictly identical structures. This might not work for output formats that don't support adding fields to existing non-empty layers.
#' @param relaxedFieldNameMatch Logical. (starting with GDAL 1.11) Do field name matching between source and existing target layer in a more relaxed way if the target driver has an implementation for it. [-relaxedFieldNameMatch] [-forceNullable]
#' @param forceNullable Logical. (starting with GDAL 2.0) Do not propagate not-nullable constraints to target layer if they exist in source layer.
#' @param unsetDefault Logical. (starting with GDAL 2.0) Do not propagate default field values to target layer if they exist in source layer.
#' @param unsetFid Logical. (starting with GDAL 2.0) Can be specify to prevent the new default behaviour that consists in, if the output driver has a FID layer creation option and we are not in append mode, to preserve the name of the source FID column and source feature IDs.
#' @param nomd Logical. (starting with GDAL 2.0) To disable copying of metadata from source dataset and layers into target dataset and layers, when supported by output driver.
#' @param mo Character. "META-TAG=VALUE". (starting with GDAL 2.0) Passes a metadata key and value to set on the output dataset, when supported by output driver.
## @param additional_commands Character. Additional commands to pass directly to ogrinfo.
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#'  
#' @return character
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'ogr2ogr' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://gdal.org/ogrinfo.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' PERFORMANCE HINTS
#' 
#' When writing into transactional DBMS (SQLite/PostgreSQL,MySQL, etc...), it might be 
#' beneficial to increase the number of INSERT statements executed between BEGIN TRANSACTION 
#' and COMMIT TRANSACTION statements. This number is specified with the -gt option. For 
#' example, for SQLite, explicitly defining -gt 65536 ensures optimal performance while 
#' populating some table containing many hundredth thousand or million rows. However, note 
#' that if there are failed insertions, the scope of -skipfailures is a whole transaction.
#' 
#' For PostgreSQL, the PG_USE_COPY config option can be set to YES for significantly insertion
#' performance boot. See the PG driver documentation page.
#' 
#' More generally, consult the documentation page of the input and output drivers for performance hints.
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/ogr2ogr.html}
#' 
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(valid_install)
#' {
#' src_datasource_name <- system.file("external/tahoe_highrez_training.shp", package="gdalUtils")
#' dst_datasource_name <- paste(tempfile(),".shp",sep="")
#' ogrinfo(src_datasource_name,"tahoe_highrez_training")
#' # reproject the input to mercator
#' ogr2ogr(src_datasource_name,dst_datasource_name,t_srs="EPSG:3395",verbose=TRUE)
#' ogrinfo(dirname(dst_datasource_name),layer=remove_file_extension(basename(dst_datasource_name)))
#' }
#' @export

ogr2ogr <- function(src_datasource_name,dst_datasource_name,
		layer,
		f,append,overwrite,update,select,progress,sql,dialect,
		where,skipfailures,spat,spat_srs,geomfield,dsco,lco,
		nln,nlt,dim,a_srs,t_srs,s_srs,preserve_fid,fid,
		oo,doo,gt,ds_transaction,
		clipsrc,clipsrcsql,clipsrclayer,
		clipsrcwhere,clipdst,clipdstsql,clipdstlayer,
		clipdstwhere,wrapdateline,datelineoffset,
		simplify,segmentize,
		fieldTypeToString,mapFieldType,unsetFieldWidth,
		splitlistfields,maxsubfields,
		explodecollections,zfield,gcp,order,
		tps,fieldmap,addfields,relaxedFieldNameMatch,
		forceNullable,unsetDefault,unsetFid,
		nomd,mo,
#		additional_commands,
		ignore.full_scan=TRUE,
		verbose=FALSE)
{
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("append","overwrite","update",
							"progress","skipfailures",
							"preserve_fid",
							"ds_transaction",
							"wrapdateline",
							"datelineoffset","unsetFieldWidth",
							"splitlistfields","explodecollections",
							"tps","addfields","relaxedFieldNameMatch",
							"forceNullable","unsetDefault","unsetFid",
							"nomd","mo")),
			vector = list(
					varnames <- c("spat","clipdst")),
			scalar = list(
					varnames <- c("dim","gt","simplify","segmentize",
							"maxsubfields","order")),
			character = list(
					varnames <- c("f","select","sql","dialect",
							"where","spat_srs","geomfield","dsco","lco","nln","nlt",
							"a_srs","t_srs","s_srs",
							"fid","oo","doo",
	#						"clipsrc",
							"clipsrcsql","clipsrclayer",
							"clipsrcwhere","clipdst","clipdstsql",
							"clipdstlayer","clipdstwhere","fieldTypeToString",
							"mapFieldType",
							"zfield","fieldmap",
							"dst_datasource_name","src_datasource_name",
							"layer")),
			repeatable = list(
					varnames <- c("gcp"))
	)
	
	# Fix for clipsrc bug reported by Alex Zvoleff, 5/11/2015
	if(!missing(clipsrc))
	{
		if(is.numeric(clipsrc))
		{
			parameter_variables$vector[[1]] <- c(parameter_variables$vector[[1]],"clipsrc")
		} else
		{
			parameter_variables$character[[1]] <- c(parameter_variables$character[[1]],"clipsrc")
		}
	}
	
#	browser()
	
	parameter_order <- c(
			"append","overwrite","update",
			"progress","skipfailures",
			"preserve_fid",
			"ds_transaction",
			"wrapdateline",
			"datelineoffset","unsetFieldWidth",
			"splitlistfields","explodecollections",
			"tps","addfields","relaxedFieldNameMatch",
			"forceNullable","unsetDefault","unsetFid",
			"nomd","mo",
			"spat","clipdst",
			"dim","gt","simplify","segmentize",
			"maxsubfields","order",
			"f","select","sql","dialect",
			"where","spat_srs","geomfield","dsco","lco","nln","nlt",
			"a_srs","t_srs","s_srs",
			"fid","oo","doo","clipsrc","clipsrcsql","clipsrclayer",
			"clipsrcwhere","clipdst","clipdstsql",
			"clipdstlayer","clipdstwhere","fieldTypeToString",
			"mapFieldType",
			"zfield","fieldmap",
			"gcp",			
			"dst_datasource_name","src_datasource_name",
			"layer"
			)
	
	parameter_noflags <- c("dst_datasource_name","src_datasource_name","layer")
	
	parameter_doubledash <- NULL
	
	parameter_noquotes <- c("spat","clipsrc")
	
	executable <- "ogr2ogr"
	# End ogr2ogr setup
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_doubledash=parameter_doubledash,
			parameter_noquotes=parameter_noquotes)
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 

		return(cmd_output)
}
