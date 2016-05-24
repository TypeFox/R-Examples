#' Save ffdf data.frames in a directory
#'
#' \code{save.ffdf} saves all ffdf data.frames in the given \code{dir}. Each column
#' is stored as with filename <ffdfname>$<colname>.ff. All variables given in "..." are stored in ".RData" in the same directory.
#' The data can be reloaded by starting a R session in the directory or by using \code{\link{load.ffdf}}.
#' Note that calling \code{save.ffdf} multiple times for the same directory 
#' will only store the ffdf's that were given in the last call. 
#' 
#' Using \code{save.ffdf} automagically sets the \code{\link{finalizer}}s of the \code{ff}
#' vectors to \code{"close"}. This means that the data will be preserved on disk when the 
#' object is removed or the R sessions is closed. Data can be deleted either using
#' \code{\link{delete}} or by removing the directory where the object were saved 
#' (\code{dir}).
#' @note
#' When saving in the temporary directory pointed at by getOption("fftempdir"), \code{ff} assumes that the
#' resulting files are to be deleted. Be sure to change the finalizers of the 
#' ff vectors when saving in the temporary directory.
#' @example ../examples/save_ffdf.R
#' @param ... \code{ffdf} data.frames, \code{ff} vectors, or other variables to be saved in the directory
#' @param dir path where .RData file will be saved and all columns of supplied \code{ffdf}'s. It will be created if it doesn't exist.
#' @param clone should the ff vectors be \code{\link{clone}}'d, creating a snapshot of the supplied ffdf or ff objects?
#' This should only be necessary if you still need the ff vectors in their current storage location.
#' @param relativepath \code{logical} if \code{TRUE} the stored ff vectors will have relative paths, making moving the data to another storage a simple
#' copy operation.
#' @param overwrite \code{logical} If \code{TRUE} \code{save.ffdf} will overwrite 
#' an previous stored \code{ffdf}, \code{.Rdata} file.
#' @seealso \code{\link{load.ffdf}} 
#' @export
save.ffdf <- function(..., dir="./ffdb", clone=FALSE, relativepath=TRUE, overwrite=FALSE){
   names <- as.character(substitute(list(...)))[-1L]
   dir.create(dir, showWarnings=FALSE, recursive=TRUE)
   
   # TODO store individual object as *.rds files instead of .Rdata: this make it easier to add different
   # ffdf data.frames to an existing directory.
   # However, it makes moving/deleting ffdf's more obscure.
   
   oldwd <- setwd(dir)
   on.exit(setwd(oldwd))
   if (!isTRUE(overwrite) && file.exists(".Rdata")){
     stop("Directory '",dir,"' contains existing '.Rdata' file. 
          To force saving use 'overwrite=TRUE'")
   }
   
   # TODO make this fail safe: when one of the 'n' is non-existing
   existing <- sapply(names, exists, envir=parent.frame())
   names <- names[existing]
   if (any(!existing)){
     warning(names[!existing], " were not saved, because not found")
   }
   for (n in names){
     x = get(n, envir=parent.frame())
     if (is.ffdf(x)) {
       if (isTRUE(clone)){
         x <- ff::clone(x)
       }
       assign(n, move.ffdf(x, dir=".", name=n, relativepath=relativepath))
     }
   }
   
   save(list=names, file=".RData")
   
   rp <- file(".Rprofile", "wt")
   writeLines(".First<-", rp)
   writeLines(deparse(first), rp)
   close(rp)
   
   if (relativepath && !clone){
     for (n in names){
       x = get(n, envir = parent.frame())
       if (is.ffdf(x)){
         for (i in bit::physical(x)){
           filename(i) <- filename(i)
         }
         close(x)
       } else if (is.ff(x)){
          filename(x) <- filename(x)
          close(x)
       }
     }
   }
}

#' Moves all the columns of ffdf data.frames into a directory
#'
#' \code{move.ffdf} saves all columns into the given \code{dir}. Each column
#' is stored as with filename <ffdfname>$<colname>.ff. 
#' If you want to store the data for an other session please use \code{\link{save.ffdf}} or \code{\link{pack.ffdf}}
#' @example ../examples/save_ffdf.R
#' @param x \code{ffdf} data.frame to be moved
#' @param dir path were all of supplied \code{ffdf}'s, will be saved. It will be created if it doesn't exist.
#' @param name name to be used as data.frame name
#' @param relativepath If \code{TRUE} the \code{ffdf} will contain relativepaths. Use with care...
#' @seealso \code{\link{load.ffdf}} \code{\link{save.ffdf}}
#' @export
move.ffdf <- function(x, dir=".", name=as.character(substitute(x)), relativepath=FALSE){  
  gc()
  dir.create(dir, showWarnings=FALSE, recursive=TRUE)
  for (colname in names(x)){
    ffcol <- x[[colname]]
    
    ffcolname <- file.path(dir, paste(name, "$", colname, ".ff", sep=""))
    
    # move file to right directory
    filename(ffcol) <- ffcolname
    
    # set path to relative path, BEWARE if wd is changed this should be reset!
    if (isTRUE(relativepath)){
      bit::physical(ffcol)$filename <- ffcolname
    }
  }
  close(x)
  x
}

#' Loads ffdf data.frames from a directory
#'
#' \code{load.ffdf} loads ffdf data.frames from the given \code{dir}, that were stored using \code{\link{save.ffdf}}. Each column
#' is stored as with filename <ffdfname>$<colname>.ff. All variables are stored in .RData in the same directory.
#' The data can be loaded by starting a R session in the directory or by using \code{\link{load.ffdf}}.
#' @example ../examples/save_ffdf.R
#' @param dir path from where the data should be loaded
#' @param envir environment where the stored variables will be loaded into.
#' @seealso \code{\link{load.ffdf}} 
#' @importFrom tools file_path_as_absolute
#' @export
load.ffdf <- function(dir, envir=parent.frame()){
  if (!isTRUE(file.exists(dir))){
    stop("directory '", dir ,"' does not exist")
  }
  oldwd <- setwd(dir)
  on.exit(setwd(oldwd))
  
  env <- new.env()
  
  load(".RData", envir=env)
  names <- ls(envir=env, all.names=TRUE)
  
  for (n in names){
    x = get(n, envir=env)
    if (is.ffdf(x)){
      for (i in bit::physical(x)){
        bit::physical(i)$filename <- file_path_as_absolute(bit::physical(i)$filename)
      }
      close(x)
    } else if (is.ff(x)){
      bit::physical(x)$filename <- file_path_as_absolute(bit::physical(x)$filename)
      close(x)
    }
    assign(n, x, envir=envir)
  }
  invisible(env)
}

#' Packs ffdf data.frames into a compressed file
#'
#' \code{pack.ffdf} stores ffdf data.frames into the given \code{file} for easy archiving and movement of data.
#' The file can be restored using \code{\link{unpack.ffdf}}. If \code{file} ends with ".zip", the package will be zipped
#' otherwise it will be tar.gz-ed.
#' @example ../examples/save_ffdf.R
#' @param file packaged file, zipped or tar.gz.
#' @param ... ff objects to be packed
#' @seealso \code{\link{save.ffdf}} \code{\link{unpack.ffdf}}
#' @importFrom tools file_ext 
#' @export
pack.ffdf <- function(file, ...){
  td <- tempfile("pack")
  save.ffdf(..., dir=td, clone=TRUE, relativepath=TRUE)
  
  file.create(file)
  file <- file_path_as_absolute(file)
  file.remove(file)
  
  oldwd <- setwd(td)
  on.exit(setwd(oldwd))
  
  d <- c(".Rprofile", ".RData", dir(td))
  
  # if file extension is zip, zip it otherwise tar.gz it
  switch( file_ext(file)
        , zip = zip(zipfile=file, files=d)
        , tar(tarfile=file, ".", compression="gzip")
        )
}

#' Unpacks previously stored ffdf data.frame into a directory
#'
#' \code{unpack.ffdf} restores ffdf data.frames into the given \code{dir}, that were stored using \code{\link{pack.ffdf}}.
#' If \code{dir} is \code{NULL} (the default) the data.frames will restored in a temporary directory.
#' if
#' @example ../examples/save_ffdf.R
#' @param file packaged file, zipped or tar.gz.
#' @param dir path where the data will be saved and all columns of supplied \code{ffdf}'s. It will be created if it doesn't exist.
#' @param envir the environment where the stored variables should be loaded into.
#' @seealso \code{\link{load.ffdf}} \code{\link{pack.ffdf}} 
#' @export
unpack.ffdf <- function(file, dir=NULL, envir=parent.frame()){
  if (is.null(dir)){ 
    dir <- tempfile("unpack")
  }
  
  switch( file_ext(file)
        , zip = unzip(zipfile=file, exdir=dir)
        , untar(tarfile=file, exdir=dir)
  )
  
  env <- load.ffdf(dir, envir=envir)
  invisible(env)
}

first <- function(){
  if (!requireNamespace("ffbase")){
    stop("Please install package ffbase, otherwise the files cannot be loaded.")
  }
  env <- load.ffdf(".", parent.frame())
}

# x <- as.ffdf(iris)
# pack.ffdf("test.zip", x)
