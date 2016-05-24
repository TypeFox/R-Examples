
# G Main Loop timeouts and idles

gTimeoutAdd <-
function(interval, f, data = NULL)
{
 useData <- !missing(data)

 if(!is.function(f))
   stop("Timeout handlers must be function objects")

 .RGtkCall("R_addGTimeoutHandler", as.integer(interval),  f, data, useData, PACKAGE = "RGtk2")
}


gSourceRemove <-
function(id)
{
  checkPtrType(id, "GTimeoutId")
  .Call("R_removeGSource", as.integer(id), PACKAGE = "RGtk2")
}


gIdleAdd <-
function(f, data = NULL)
{
 useData <- !missing(data)

 if(!is.function(f)) {
  stop("Idle functions must be functions!")
 }

  .RGtkCall("R_addGIdleHandler", f, data, useData, PACKAGE = "RGtk2")
}

# transparent coercion

as.GQuark <-
function(x)
{
	if (is.character(x))
		x <- .Call("R_gQuarkFromString", x, PACKAGE = "RGtk2")
	else {
		x <- as.integer(x)
		class(x) <- "GQuark"
	}
	x
}
as.GList <-
function(x)
{
	x <- as.list(x)
	class(x) <- "GList"
	x
}
as.GSList <-
function(x)
{
	x <- as.list(x)
	class(x) <- "GSList"
	x
}
as.GString <-
function(x)
{
	x <- as.character(x)
	class(x) <- "GString"
	x
}
as.GTimeVal <-
function(x)
{
	x <- as.struct(x, "GTimeVal", c("tv_sec", "tv_usec"))
	x[[1]] <- as.numeric(x[[1]])
	x[[2]] <- as.numeric(x[[2]])
	
	return(x)
}

as.GError <-
function(x)
{
  x <- as.struct(x, "GError", c("domain", "code", "message"))
  x[[1]] <- as.GQuark(x[[1]])
  x[[2]] <- as.integer(x[[2]])
  x[[3]] <- as.character(x[[3]])
  x
}

# Provide the G_FILE_ERROR domain
G_FILE_ERROR <- gFileErrorQuark <-
function()
{
	
	w <- .RGtkCall("S_g_file_error_quark", PACKAGE = "RGtk2")

	return(w)
}

# The GFileError enumeration

GFileError <- c(
"exist" = 0,
"isdir" = 1,
"acces" = 2,
"nametoolong" = 3,
"noent" = 4,
"notdir" = 5,
"nxio" = 6,
"nodev" = 7,
"rofs" = 8,
"txtbsy" = 9,
"fault" = 10,
"loop" = 11,
"nospc" = 12,
"nomem" = 13,
"mfile" = 14,
"nfile" = 15,
"badf" = 16,
"inval" = 17,
"pipe" = 18,
"again" = 19,
"intr" = 20,
"io" = 21,
"perm" = 22,
"nosys" = 23,
"failed" = 24
)
storage.mode(GFileError) <- 'integer'
class(GFileError) <- 'enums'

.GFileError <- c(
"G_FILE_ERROR_EXIST" = 0,
"G_FILE_ERROR_ISDIR" = 1,
"G_FILE_ERROR_ACCES" = 2,
"G_FILE_ERROR_NAMETOOLONG" = 3,
"G_FILE_ERROR_NOENT" = 4,
"G_FILE_ERROR_NOTDIR" = 5,
"G_FILE_ERROR_NXIO" = 6,
"G_FILE_ERROR_NODEV" = 7,
"G_FILE_ERROR_ROFS" = 8,
"G_FILE_ERROR_TXTBSY" = 9,
"G_FILE_ERROR_FAULT" = 10,
"G_FILE_ERROR_LOOP" = 11,
"G_FILE_ERROR_NOSPC" = 12,
"G_FILE_ERROR_NOMEM" = 13,
"G_FILE_ERROR_MFILE" = 14,
"G_FILE_ERROR_NFILE" = 15,
"G_FILE_ERROR_BADF" = 16,
"G_FILE_ERROR_INVAL" = 17,
"G_FILE_ERROR_PIPE" = 18,
"G_FILE_ERROR_AGAIN" = 19,
"G_FILE_ERROR_INTR" = 20,
"G_FILE_ERROR_IO" = 21,
"G_FILE_ERROR_PERM" = 22,
"G_FILE_ERROR_NOSYS" = 23,
"G_FILE_ERROR_FAILED" = 24
)
storage.mode(.GFileError) <- 'integer'
class(.GFileError) <- 'enums'
