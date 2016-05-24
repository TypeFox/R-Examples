#' Construct templatebrain object for an image registration template
#'
#' \code{templatebrain} objects encapsulate key information for the reference
#' brain in an image registration. Usually this will be a standard template
#' brain used for many registrations. \strong{It will normally be much more
#' convenient to use  \code{\link{as.templatebrain}} methods to convert an image
#' file or an im3d object into a \code{templatebrain}}.
#'
#' A variety of methods are available to work on \code{templatebrain} objects. See
#' \code{\link{templatebrain-meths}} for basic methods. The two main functions
#' that are availavle for using template brains are \code{\link{xform_brain}} and
#' \code{\link{mirror_brain}}.
#'
#' \code{templatebrain} objects are only useful for transformation processes
#' when the \code{BoundingBox} is specified to define the physical extent of the
#' volume. We use the definition of the Amira 3D visualisation and analysis
#' software. This corresponds to the \strong{node} centers option in the
#' \href{http://teem.sourceforge.net/nrrd/format.html}{NRRD format}. The
#' bounding box can be obtained from NRRD or AmiraMesh format files. See
#' \code{\link[nat]{boundingbox}} for details.
#'
#' @param name the full name of the template.
#' @param regName the short name. This will be the stem used to prefix
#'   registrations (e.g. JFRC2_someimage.list) for this template brain and
#'   likely also the stem of the template brain image (e.g. JFRC2.nrrd).
#' @param type one of \code{c('single brain', 'average')}, indicating whether
#'   the template brain has been created from just one image, or is the average
#'   of multiple images.
#' @param sex the sex of the template brain. For templates with
#'   \code{type=='average'}, the possibility of \code{sex='intersex'} exists.
#' @param dims dimensions of the image (number of voxels).
#' @param BoundingBox physical dimensions of the image (see
#'   \code{\link[nat]{boundingbox}}).
#' @param voxdims physical spacing between voxels.
#' @param units units of physical measurements (e.g. microns).
#' @param description details of the template.
#' @param doi a DOI for the original template brain image.
#' @rdname templatebrain
#' @param ... additional named arguments that will be added as fields to the
#'   \code{templatebrain} object.
#' @return A list with class \code{templatebrain}.
#' @export
#' @seealso \code{\link{as.templatebrain}}, \code{\link{templatebrain-meths}},
#'   \code{\link{xform_brain}}, \code{\link{mirror_brain}}.
templatebrain<-function(name, regName=name, type=NULL, sex=NULL, dims=NULL,
                        BoundingBox=NULL, voxdims=NULL, units=NULL,
                        description=NULL, doi=NULL, ...) {
  template <- structure(list(name=name, regName=regName, type=type, sex=sex,
                             dims=dims, voxdims=voxdims, origin=origin,
                             BoundingBox=BoundingBox, units=units,
                             description=description, doi=doi),
                        class="templatebrain")
  if(!missing(...)) {
    apl=pairlist(...)
    template[names(apl)]=apl
  }
  template
}

#' Use image file or other object to initialise template brain
#'
#' @param x object used to construct the templatebrain, either a character
#'   vector with the path to a file or an \code{im3d} object.
#' @param ... additional named arguments passed to methods and then on to
#'   \code{templatebrain} that will be added as fields to the
#'   \code{templatebrain} object.
#' @return A list with class \code{templatebrain}.
#' @export
#' @seealso  \code{\link[nat]{im3d}}
#' @rdname as.templatebrain
as.templatebrain <- function(x, ...) UseMethod("as.templatebrain")

#' @rdname as.templatebrain
#' @importFrom nat read.im3d voxdims boundingbox origin
#' @export
#' @examples
#' # Make templatebrain object using image info from the template brain NRRD file
#' nhdr=system.file('images','FCWB.nhdr', package='nat.templatebrains')
#' as.templatebrain(nhdr, name = "FlyCircuit Whole Brain")
as.templatebrain.character <- function(x, ...) {
  if(!file.exists(x)) stop("x does not specify a valid path!")
  im3d <- read.im3d(x, ReadData=FALSE)
  as.templatebrain(im3d, ...)
}

#' @param name,regName name and short name of the template brain. Will use the
#'   filename (minus final extension) by default for both fields.
#' @rdname as.templatebrain
#' @export
as.templatebrain.im3d <- function(x, name=NULL, regName=NULL, ...) {
  # This will be incorrect if the directions are not rectilinear
  file_stem = sub("\\.[^.]+$", "", basename(attr(x, 'file')))
  if(is.null(name)) name = file_stem
  if(is.null(regName)) regName = file_stem
  units <- attr(x, 'header')$'space units'
  templatebrain(name=name, regName=regName, dims=dim(x), voxdims=voxdims(x),
                origin=origin(x), BoundingBox=boundingbox(x),
                units=units, ...)
}

#' Template brain methods
#' @description \code{is.templatebrain} tests if object is of class templatebrain
#' @param x an object (usually a \code{templatebrain}).
#' @return A logical indicating whether or not the object is a \code{templatebrain}.
#' @export
#' @name templatebrain-meths
#' @aliases is.templatebrain
#' @examples
#' data(FCWB.demo)
#' is.templatebrain(FCWB.demo)
#' origin(FCWB.demo)
#' dim(FCWB.demo)
#' voxdims(FCWB.demo)
#' boundingbox(FCWB.demo)
is.templatebrain<-function(x) inherits(x, 'templatebrain')

#' @description \code{as.character.templatebrain} converts template brain to
#'   character vector representation (normally used to extract the short name
#'   i.e. \code{regName}).
#'
#' @param field which field to use (defaults to \code{'regName'}).
#' @param ... additional arguments for methods.
#' @return Character vector.
#' @export
#' @rdname templatebrain-meths
as.character.templatebrain<-function(x, field=c('regName','name'), ...){
  field=match.arg(field)
  x[[field]]
}

#' @description \code{print.templatebrain} prints templatebrain information in
#'   human-readable form
#' @export
#' @rdname templatebrain-meths
print.templatebrain <- function(x, ...) {
  cat("=== Template Brain ===", "\n")
  cat("Name:", x$name, "\n")
  cat("Short Name:", x$regName, "\n")
  cat("Type:", x$type, "\n")
  cat("Sex: ", x$sex, "\n")
  cat(paste0("Voxel size:\n"))
  cat("  x =", paste0(x$voxdims[1], " ", x$units[1], "\n"))
  cat("  y =", paste0(x$voxdims[2], " ", x$units[2], "\n"))
  cat("  z =", paste0(x$voxdims[3], " ", x$units[3], "\n"))
  cat(paste0("Bounding box (", x$units[1], "):\n"))
  cat("  x =", paste0(x$BoundingBox[1, 1], ", y = ", x$BoundingBox[1, 2], ", z = ", x$BoundingBox[1, 3], ",\n"))
  cat("  x =", paste0(x$BoundingBox[2, 1], ", y = ", x$BoundingBox[2, 2], ", z = ", x$BoundingBox[2, 3], ".\n"))
  cat("Description:", x$description, "\n")
  cat("DOI:", x$doi)
  if(exists('...')) {
    cat("\n")
    cat(...)
  }
}

#' @export
#' @description \code{as.im3d} converts a template brain to a \code{nat::im3d}
#'   object; this is probably useful for developers.
#' @method as.im3d templatebrain
#' @importFrom nat as.im3d
#' @rdname templatebrain-meths
#' @seealso \code{\link[nat]{im3d}}
as.im3d.templatebrain <- function(x, ...) {
  newim3d <- nat::im3d(NA, dims=x$dims, voxdims=x$voxdims, origin=x$origin)
  newim3d
}

#' @export
#' @description \code{origin} extracts the space origin of a \code{templatebrain}
#'   object.
#' @method origin templatebrain
#' @importFrom nat origin
#' @rdname templatebrain-meths
#' @seealso \code{\link[nat]{origin}}
origin.templatebrain <- function(x, ...) {
  origin(nat::as.im3d(x))
}

#' @description \code{dim} extracts the dimensions (in number of pixels) of the
#'   image associated with a \code{templatebrain} object.
#' @export
#' @method dim templatebrain
#' @rdname templatebrain-meths
dim.templatebrain <- function(x, ...) {
  dim(nat::as.im3d(x))
}

#' @description \code{voxdims} extracts the dimensions (in number of pixels) of
#'   the image associated with a \code{templatebrain} object.
#' @export
#' @method voxdims templatebrain
#' @importFrom nat voxdims
#' @rdname templatebrain-meths
#' @seealso \code{\link[nat]{voxdims}}
voxdims.templatebrain <- function(x, ...) {
  voxdims(nat::as.im3d(x))
}

#' @description \code{boundingbox} extracts the boundingbox (in calibrated
#'   spatial units, typically microns) of the image associated with a
#'   templatebrain object. See \code{\link[nat]{boundingbox}} for details.
#' @export
#' @method boundingbox templatebrain
#' @rdname templatebrain-meths
#' @importFrom nat boundingbox
#' @seealso \code{\link[nat]{boundingbox}}
boundingbox.templatebrain <- function(x, ...) {
  boundingbox(nat::as.im3d(x))
}
