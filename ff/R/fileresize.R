# R file resize
# (c) 2009 Daniel Adler
# Licence: ISC

#! \name{file.resize}
#! \alias{file.resize}
#! \alias{file.move}
#! \title{ Change size of move an existing file }
#! \description{
#!   Change size of an existing file (on some platforms sparse files are
#!   used) or move file to other name and/or location.
#! }
#! \usage{
#!   file.resize(path, size)
#!   file.move(from, to)
#! }
#! \arguments{
#!    \item{path}{ file path (on windows it uses a 'windows' backslash path!)  }
#!    \item{size}{ new filesize in bytes as double }
#!    \item{from}{ old file path  }
#!    \item{to}{ new file path }
#! }
#! \value{
#!   logical scalar repesenting the success of this operation
#! }
#! \details{
#!   \code{file.resize} can enlarge or shrink the file. When enlarged, the file
#!   is filled up with zeros. Some platform implementations feature
#!   sparse files, so that this operation is very fast. We have tested:
#!   \itemize{
#!    \item Ubuntu Linux 8, i386
#!    \item FreeBSD 7, i386
#!    \item Gentoo Linux Virtual-Server, i386
#!    \item Gentoo Linux, x86_64
#!    \item Windows XP
#!   }
#!   The following work but do not support sparse files
#!   \itemize{
#!    \item Mac OS X 10.5, i386
#!    \item Mac OS X 10.4, PPC
#!   }
#!  \code{file.move} tries  to \code{\link{file.rename}}, 
#!  if this fails (e.g. across file systems) the file is copied to the new location and the old file is removed, 
#!  see  \code{\link{file.copy}} and \code{\link{file.remove}}.
#! }
#! \author{ Daniel Adler }
#! \seealso{ \code{\link[base]{file.create}}, \code{\link[base]{file.rename}}, \code{\link[base]{file.info}}, \code{\link{file.copy}}, \code{\link{file.remove}} }
#! \examples{
#!  x <- tempfile()
#!  newsize <- 23       # resize and size to 23 bytes.
#!  file.resize(x, newsize)
#!  file.info(x)$size == newsize
#!  \dontrun{
#!    newsize <- 8*(2^30) # create new file and size to 8 GB.
#!    file.resize(x, newsize)
#!    file.info(x)$size == newsize
#!  }
#!  y <- tempfile()
#!  file.move(x,y)
#!  file.remove(y)
#! }
#! \keyword{ IO }
#! \keyword{ data }

file.resize <- function(path, size)
  .Call("r_file_resize", as.character(path), as.double(size), PACKAGE="ff" )

file.move <- function(from, to){
	oldopt <- options("warn")
	on.exit(options(oldopt))
	options(warn=-1)
	if(file.rename(from, to)){
		TRUE
  }else{
		options(oldopt)
		on.exit()
		file.copy(from, to) && file.remove(from)
	}
}