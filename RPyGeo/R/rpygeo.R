# 2011-09-07
# - rpygeo.build.env and rpygeo.geoprocessor
#   snapraster now supported
# 2008-08-15
# - rpygeo.geoprocessor
#   - Now using os.chdir to set OS working directory in
#     Python prior to calling the geoprocessor; this seems to
#     help with some(!) relative paths in ArcGIS geoprocessor
#     calls...
#   - Mask did not work because of double double quotes;
#     now fixed
# - New functions rpygeo.ASCIIToRaster.conversion and
#   rpygeo.RasterToASCII.conversion
#   rpygeo.FlowDirection.sa, rpygeo.FlowAccumulation.sa,
#   rpygeo.FlowLength.sa, rpygeo.Sink.sa
# - New function rpygeo.Spline.sa, rpygeo.TopoToRaster.sa

# to do:
# - Adapt rpygeo.Delete.management to be able to delete multiple
#   files in one call (and to use Exist calls to prevent errors?).
# - Better: write a more flexible rpygeo.bundle function
#   that aggregates several calls using the same geoprocessing
#   environment!

library(RSAGA)
# RPyGeo uses default.file.extension, match.arg.ext, get.file.extension
# functions from RSAGA package.

##############################################
# Helper functions
##############################################

rpygeo.extent.to.character = function(x) {
    if (is.null(x)) return(x)
    if (is.na(x)) return(NULL)
    if (is.character(x)) return(x)
    x = paste(x$x[1], x$y[1], x$x[2], x$y[2], sep = " ")
}

rpygeo.build.env = function( 
    modules = "arcgisscripting",
    init = "gp = arcgisscripting.create()",
    workspace = NULL, 
    cellsize = NULL,
    extent = NULL,
    mask = NULL,
    snapraster = NULL,
    overwriteoutput = 0,
    extensions = NULL,
    python.path = "C:\\software\\Python24",
    python.command = "python.exe" )
{
    return( list(
        modules = modules,
        init = init,
        workspace = workspace,
        cellsize = cellsize,
        extent = rpygeo.extent.to.character(extent),
        mask = mask,
        snapraster = snapraster,
        overwriteoutput = overwriteoutput,
        extensions = extensions,
        python.path = python.path,
        python.command = python.command
    ) )
}


write.point.shapefile = function(d, file, x.field = "x", y.field = "y", id.field = NULL)
{
    # Prepare data for shapefile - see example in ?write.shapefile
    # in package shapefiles:
    require(shapefiles)
    if (is.null(id.field)) {
        id.field = "id"
        create.id.field = TRUE
    } else if (all(colnames(d) != id.field))
        create.id.field = TRUE

    if (create.id.field)
        d[[id.field]] = c(1:nrow(d))
    dd = data.frame( id = d[,id.field], 
                        x = d[,x.field],
                        y = d[,y.field] )
    ddShapefile = convert.to.shapefile(dd, d, id.field, 1)

    # Write shapefile:
    write.shapefile(ddShapefile, file, arcgis = TRUE)
}

write.temp.point.shapefile = function(d, pattern = "file",
    tmpdir = tempdir(), ...)
{

    # Write shapefile:
    tmpfile = tempfile(pattern = pattern, tmpdir = tmpdir)
    tmpfile = gsub("\\", "/", tmpfile, fixed = TRUE)
    write.point.shapefile(d, file = tmpfile, ...)

    # Create an expression that will allow the caller to
    # delete the temporary shapefile:    
    exit.expr = bquote( {
        unlink(.(paste(tmpfile,".shp",sep="")));
        unlink(.(paste(tmpfile,".shx",sep="")));
        unlink(.(paste(tmpfile,".dbf",sep=""))) } )
    tmpfile = paste(tmpfile, ".shp", sep="")

    return( list( tempfile = tmpfile, exit.expression = exit.expr ) )
}



##############################################
# Geoprocessing environment
##############################################

#rpygeo.env = rpygeo.build.env()

rpygeo.env = list(
    modules = "arcgisscripting",
    init = "gp = arcgisscripting.create()",
    workspace = NULL, 
    cellsize = NULL,
    extent = NULL,
    mask = NULL,
    overwriteoutput = 0,
    extensions = NULL,
    python.path = "C:\\software\\Python24",
    python.command = "python.exe" )
    

rpygeo.required.extensions = function(expr) {
    # See ArcGIS help on the CheckOutExtension method:
    rpygeo.match.extensions = c("sa","3d","stats","na","di")
    names(rpygeo.match.extensions) = c("Spatial","3d","geostats","network","datainteroperability")
    ext = c()
    for (s in expr) {
        sub.s = strsplit(s,"(",fixed=TRUE)[[1]]
        for (t in sub.s) {
            t = gsub(" ","",t)
            if (substring(t,nchar(t)) == ")")  next
            for (i in 1:length(rpygeo.match.extensions)) { 
                the.match = paste("_",tolower(rpygeo.match.extensions[i]),sep="")
                if ( tolower(substring(t,nchar(t)+1-nchar(the.match))) == the.match )
                    ext = c( ext, names(rpygeo.match.extensions)[i] )
            }
        }
    }
    return(unique(ext))
}



##############################################
# RPyGeo Geoprocessor - the workhorse
##############################################

rpygeo.geoprocessor = function(
        fun, args=NULL,
        py.file="rpygeo.py", msg.file="rpygeo.msg",
        env = rpygeo.env, extensions = NULL, working.directory = getwd(),
        quote.args = TRUE, add.gp = TRUE, wait = TRUE,
        clean.up = wait,
        detect.required.extensions = TRUE )
{
    if (is.logical(clean.up)) {
        if (clean.up) {
            clean.up = c("py","msg")
        } else
            clean.up = c()
    }

    # Convert to character string,
    # and add quotation marks if input was already a string:
    convert = function(x) {
        if (is.numeric(x)) {
            return( as.character(x) )
        } else
            return( paste('"', x, '"', sep="" ) )
    }
    # Consistent indentation is important in Python:
    indent = "    "

    # Expecting a list of arguments, not a vector,
    # because arguments may have different data types:
    if (is.vector(args)) args = as.list(args)
    if ((length(fun) > 1)  & (!is.null(args))) {
        warning("Multiple function calls only allowed if args is NULL. Using only first `fun' element.\n")
        fun = fun[1]
    }
    if (!is.null(args)) if (length(quote.args)==1) quote.args = rep(quote.args,length(args))

    # Create list of required ArcGIS extensions:
    extensions = c( env$extensions, extensions )
    if (detect.required.extensions)
        extensions = c( extensions, rpygeo.required.extensions(fun) )
    extensions = unique( extensions )

    # Have to distinguish between an R version and a Windows version of file names.
    to.windows.filename = function(x) gsub("/","\\",x,fixed=TRUE)
    to.R.filename       = function(x) gsub("\\","/",x,fixed=TRUE)
    py.file = to.windows.filename( paste(working.directory,"/",py.file,sep="") )
    R.py.file = to.R.filename(py.file)
    msg.file = to.windows.filename( paste(working.directory,"/",msg.file,sep="") )
    R.msg.file = to.R.filename(msg.file)

    #---------------------------------------
    # Build the Python geoprocessing script:
    expr = ""
    # Added 2008-08-15
    # This may help with some relative paths in ArcGIS:
    if (!is.null(env$workspace)) {
        expr = paste( expr, "import os\n", sep="" )
        expr = paste( expr, "os.chdir(", convert(to.R.filename(env$workspace)), ")\n", sep="" )
    }
    # End added
    for (mod in env$modules)
        expr = paste( expr, "import ", mod, "\n", sep="" )
    for (init in env$init)
        expr = paste( expr, init, "\n", sep="" )
    if (!is.null(env$workspace))
        expr = paste( expr, "gp.Workspace = ", convert(to.R.filename(env$workspace)), "\n", sep="" )
    if (!is.null(env$cellsize))
        expr = paste( expr, "gp.Cellsize = ", convert(env$cellsize), "\n", sep="" )
    if (!is.null(env$extent))
        expr = paste( expr, "gp.Extent = ", convert(env$extent), "\n", sep="" )
    if (!is.null(env$mask))
        expr = paste( expr, "gp.Mask = ", convert(env$mask), "\n", sep="" )
    if (!is.null(env$snapraster))
        expr = paste( expr, "gp.snapRaster = ", convert(env$snapraster), "\n", sep="" )
    if (!is.null(env$overwriteoutput))
        expr = paste( expr, "gp.Overwriteoutput = ", convert(env$overwriteoutput), "\n", sep="" )
    if (!is.null(extensions))
        for (ext in extensions)
            expr = paste( expr, "gp.CheckOutExtension(", convert(ext), ")\n", sep="" )
    expr = paste( expr, 'rpygeoresult = ""', "\n", sep="" )
    expr = paste( expr, "\n", sep="" )
    
    expr = paste( expr, "try:\n", sep="" )
    
    if (is.null(args)) {
        # Simplest case - each element of `fun' is a complete Python expression:
        for (the.fun in fun)
            expr = paste( expr, indent, ifelse(add.gp,"gp.",""), the.fun, "\n", sep="" )
    } else {
        # More complicated:
        # Only one `fun' call, but many arguments need to be concatenated:
        expr = paste( expr, indent, ifelse(add.gp,"gp.",""), fun, "( ", sep="" )
        for (i.arg in 1:length(args)) {
            if (i.arg > 1)  expr = paste( expr, ", ", sep="" )
            # Character string arguments will usually have to be decorated with quotes,
            # assuming that they are string constants, not variable names or expressions:
            expr = paste( expr, 
                ifelse(quote.args[i.arg], convert(args[[i.arg]]), as.character(args[[i.arg]])), 
                sep="" )
            # to do: use intelligent line breaks in expressions??
        }
        expr = paste( expr, " )\n", sep="" )
    }

    # Catch exceptions: get error messages
    expr = paste( expr, "except:\n", sep="" )
    expr = paste( expr, indent, "rpygeoresult = gp.GetMessages()\n", sep="" )
    
    #expr = paste( expr, "\n", sep="")

    # If an error occurred, write the error message to the `msg.file':
    expr = paste( expr, 'if rpygeoresult != "":\n', sep="")
    expr = paste( expr, indent, 'f = open("', R.msg.file, '", "w")\n', sep="")
    expr = paste( expr, indent, "f.write(rpygeoresult)\n", sep="")
    expr = paste( expr, indent, "f.close()\n", sep="")
    #-----------------------------------------

    
    # Write the Python geoprocessing script to the `py.file':
    py.script = file( R.py.file, open="wt" )
    write(expr, file=py.script)
    close(py.script)
    rm(py.script)

    # Delete message file;
    # otherwise msg file would only be overwritten if an geoprocessing error occurred.
    if (file.exists(R.msg.file))
        unlink(R.msg.file)

    # Set up the system call expression:    
    py.call = ""
    if (!is.null(env$python.path))
        py.call = paste( py.call, env$python.path, "\\", sep="" )
    py.call = paste( py.call, env$python.command, " ", py.file, sep="" )
    py.call = to.windows.filename(py.call)
    
    ######## Run Python:
    system(py.call, invisible = TRUE, minimize = TRUE, wait = wait)

    # Read error messages from the `msg.file', if available:
    res = NULL
    if (file.exists(R.msg.file) & wait) {
        f.msg = file(R.msg.file,"rt")
        res = readLines(f.msg, warn=FALSE)
        close(f.msg)
        if ("msg" %in% clean.up)
            unlink(R.msg.file)
    }

    # Delete the `py.file' script:
    if ("py" %in% clean.up)
        if (file.exists(R.py.file))
            unlink(R.py.file)

    # Return error message or NULL:
    return( res )
}




##############################################
# Spatial Analyst tools
##############################################

# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Hillshade
rpygeo.Hillshade.sa = function( in.raster, out.raster,
    azimuth = 315, altitude = 45, 
    model.shadows = c("NO_SHADOWS","SHADOWS"), z.factor = 1, ...)
{
    model.shadows = match.arg.ext(model.shadows, ignore.case = TRUE)

    rpygeo.geoprocessor( fun="Hillshade_sa",
        args=list(in.raster,out.raster,azimuth,altitude,model.shadows,z.factor),
        quote.args=c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE), ... )
}

rpygeo.Slope.sa = function( in.raster, out.raster,
    unit = c("DEGREE","PERCENT_RISE"), z.factor = 1, ... )
{
    unit = match.arg.ext(unit, ignore.case = TRUE)

    rpygeo.geoprocessor( fun="Slope_sa",
        args=list(in.raster, out.raster, unit, z.factor),
        quote.args=c(TRUE,TRUE,TRUE,FALSE), ... )
}

# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Aspect
rpygeo.Aspect.sa = function( in.raster, out.raster, ... )
{
    rpygeo.geoprocessor( fun = "Aspect_sa",
        args = list(in.raster, out.raster), ... )
}

# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=EucDistance
rpygeo.EucDistance.sa = function( in.data, out.raster,
    maxdist=NULL, cellsize=NULL, out.direction.raster=NULL, 
    env = rpygeo.env, ... )
{
    if (!is.null(maxdist)) if (maxdist==Inf) maxdist = NULL
    args = list(in.data, out.raster)
    args = c(args, maxdist, cellsize, out.direction.raster)
    rpygeo.geoprocessor( fun = "EucDistance_sa",
        args=args, ... )
}


# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Curvature
rpygeo.Curvature.sa = function(in.raster, out.curvature.raster, 
    z.factor = 1, out.profile.curve.raster = NULL,
    out.plan.curve.raster = NULL, ...)
{
    args = list(in.raster, out.curvature.raster, z.factor, 
        out.profile.curve.raster, out.plan.curve.raster)

    rpygeo.geoprocessor(fun = "Curvature_sa", args = args, ...)
}

# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Determining_flow_direction
rpygeo.FlowDirection.sa = function(in.surface.raster,
    out.flow.direction.raster, force.flow = c("NORMAL","FORCE"),
    out.drop.raster = NULL, ...)
{
    force.flow = match.arg.ext(force.flow, ignore.case = TRUE)
    args = list(in.surface.raster, out.flow.direction.raster,
        force.flow, out.drop.raster)
    rpygeo.geoprocessor( fun = "FlowDirection_sa", args = args, ...)
}



# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Calculating_flow_accumulation
rpygeo.FlowAccumulation.sa = function(in.flow.direction.raster,
    out.accumulation.raster, in.weight.raster = NULL,
    data.type = c("FLOAT","INTEGER"), ...)
{
    data.type = match.arg.ext(data.type, ignore.case = TRUE)
    args = list(in.flow.direction.raster, out.accumulation.raster,
        in.weight.raster, data.type)
    rpygeo.geoprocessor( fun = "FlowAccumulation_sa", args = args, ...)
}

rpygeo.FlowLength.sa = function(in.flow.direction.raster,
    out.raster, direction.measurement = c("DOWNSTREAM","UPSTREAM"),
    in.weight.raster = NULL, ...)
{
    direction.measurement = match.arg.ext(direction.measurement, ignore.case = TRUE)
    args = list(in.flow.direction.raster, out.raster,
        direction.measurement, in.weight.raster)
    rpygeo.geoprocessor( fun = "FlowLength_sa", args = args, ...)
}

rpygeo.Sink.sa = function(in.flow.direction.raster, out.raster, ...)
{
    args = list(in.flow.direction.raster, out.raster)
    rpygeo.geoprocessor( fun = "Sink_sa", args = args, ...)
}


# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Viewshed
rpygeo.Viewshed.sa = function(in.raster, in.observer.features,
    out.raster, z.factor = 1, 
    curvature.correction = c("FLAT_EARTH", "CURVED_EARTH"),
    refractivity.coefficient = 0.13,
    x.field = "x", y.field = "y", tmpdir = tempdir(), ...)
{
    if (is.data.frame(in.observer.features)) {
        res = write.temp.point.shapefile(in.observer.features, x.field = x.field, y.field = y.field,
            pattern = "rpygeo", tmpdir = tmpdir)
        in.observer.features = res$tempfile
        on.exit(res$exit.expression)
    }
    curvature.correction = match.arg.ext(curvature.correction, ignore.case = TRUE)
    args = list(in.raster, in.observer.features, out.raster, z.factor, 
            curvature.correction, refractivity.coefficient)
    rpygeo.geoprocessor("Viewshed_sa", args = args, ...)
}

# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Area_Solar_Radiation
# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Solar%20radiation%20analysis%20equations
rpygeo.AreaSolarRadiation.sa = function(in.surface.raster, 
        out.global.radiation.raster, latitude = 45, sky.size = 200, 
        time.configuration, day.interval = 14, hour.interval = 0.5, 
        each.interval = c("NOINTERVAL", "INTERVAL"), 
        z.factor = NULL, slope.aspect.input.type = c("FROM_DEM", "FLAT_SURFACE"), 
        calculation.directions = 32, zenith.divisions = 8, 
        azimuth.divisions = 8, 
        diffuse.model.type = c("UNIFORM_SKY", "STANDARD_OVERCAST_SKY"), 
        diffuse.proportion = 0.3, transmittivity = 0.5, 
        out.direct.radiation.raster = NULL, out.diffuse.radiation.raster = NULL, 
        out.direct.duration.raster = NULL, ...)
{
    if (is.list(time.configuration)) {
        if (time.configuration[[1]] == "WithinDay") {
            stopifnot(time.configuration[[2]] > 0 & time.configuration[[2]] <= 366)
            stopifnot(time.configuration[[3]] >= 0)
            stopifnot(time.configuration[[4]] >= 0)
            stopifnot(time.configuration[[3]] <= 24)
            stopifnot(time.configuration[[4]] <= 24)
            stopifnot(time.configuration[[3]] < time.configuration[[4]])
        } else if (time.configuration[[1]] == "Year") {
            if (length(time.configuration) > 1)
                warning("calculating solar radiation for whole year\nignoring additional parameters in 'time.configuration'")
            time.configuration = "Year"
        }
        time.configuration = paste(time.configuration, collapse = " ")
    }

    args = list(in.surface.raster, out.global.radiation.raster, latitude, 
                sky.size, time.configuration, day.interval, hour.interval, 
                each.interval, z.factor, slope.aspect.input.type, 
                calculation.directions, zenith.divisions, azimuth.divisions, 
                diffuse.model.type, diffuse.proportion, transmittivity, 
                out.direct.radiation.raster, out.diffuse.radiation.raster, 
                out.direct.duration.raster)
    
    rpygeo.geoprocessor( fun = "AreaSolarRadiation_sa", args = args, ...)
}



# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Single_Output_Map_Algebra
rpygeo.SingleOutputMapAlgebra.sa = function(expression.string, out.raster, 
    in.data = NULL, ...)
{
    args = list(expression.string, out.raster)
    if (!is.null(in.data)) args = c(args, in.data)
    rpygeo.geoprocessor(fun = "SingleOutputMapAlgebra_sa",
        args = args, ...)
}


##############################################
# Data Management functions
##############################################

# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Delete_(Data_Management)
rpygeo.Delete.management = function( in.data, data.type = NULL, ... )
{
    rpygeo.geoprocessor( fun = "Delete_management",
        args = c(in.data, data.type), ... )
}





##############################################
# Conversion functions
##############################################


# http://webhelp.esri.com/arcgisdesktop/9.2/body.cfm?tocVisable=1&ID=1309&TopicName=ASCII%20to%20Raster%20(Conversion)
rpygeo.ASCIIToRaster.conversion = function( in.ascii.file, out.raster,
    data.type = c("FLOAT","INTEGER"), ... )
{
    in.ascii.file = default.file.extension(in.ascii.file, ".asc")
    if (!(tolower(get.file.extension(in.ascii.file)) %in% c(".asc",".txt")))
        stop("'in.ascii.file' must have extension '.asc' or '.txt'.\n")
    data.type = match.arg.ext(data.type, ignore.case = TRUE)
    args = list(in.ascii.file, out.raster, data.type)
    rpygeo.geoprocessor( fun = "ASCIIToRaster_conversion",
        args = args, ... )
}

# http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=Raster%20to%20ASCII%20(Conversion)
rpygeo.RasterToASCII.conversion = function( in.raster, out.ascii.file, ... )
{
    out.ascii.file = default.file.extension(out.ascii.file, ".asc")
    if (!(tolower(get.file.extension(out.ascii.file)) %in% c(".asc",".txt")))
        stop("'out.ascii.file' must have extension '.asc' or '.txt'.\n")
    args = list(in.raster, out.ascii.file)
    rpygeo.geoprocessor( fun = "RasterToASCII_conversion",
        args = args, ... )
}
