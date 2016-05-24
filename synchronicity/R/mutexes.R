# MUTEXES

#' @import Rcpp
#' @useDynLib synchronicity

#' @export
setClass('mutex')

#' @export
setGeneric('lock', function(m, ...) standardGeneric('lock'))

#' @export
setGeneric('lock.shared', function(m, ...) standardGeneric('lock.shared'))

#' @export
setGeneric('unlock', function(m, ...) standardGeneric('unlock'))

#' @export
setClass('boost.mutex', contains='mutex', 
  representation(isRead='logical', mutexInfoAddr='externalptr'))

#' @title Is it a read (shared) mutex?
#' 
#' @description Tells the user if a mutex is a read (shared) mutex. If it is 
#' not then it must be a write (exclusive) mutex.
#' @docType methods
#' @rdname read-methods
#' @param m the mutex 
#' @return TRUE if the mutex is read (shared), FALSE otherwise.
#' @export
setGeneric('read', function(m) standardGeneric('read'))

#' @rdname read-methods
#' @aliases read,boost.mutex-method
#' @export
setMethod('read', signature(m='boost.mutex'), function(m) 
  IsRead(m@mutexInfoAddr))

#' @export
setMethod('lock', signature(m='boost.mutex'),
  function(m, ...)
  {
    block = match.call()[['block']]
    if (is.null(block)) block=TRUE
    if (!is.logical(block)) stop('The block argument should be logical')
    block_call = ifelse(block, boost_lock, boost_try_lock)
    block_call(m@mutexInfoAddr)
  })

#' @export
setMethod('lock.shared', signature(m='boost.mutex'),
  function(m, ...)
  {
    block = match.call()[['block']]
    if (is.null(block)) block=TRUE
    if (!is.logical(block)) stop('The block argument should be logical')
    block_call = ifelse(block, boost_lock_shared, boost_try_lock_shared)
    block_call(m@mutexInfoAddr)
  })

#' @export
setMethod('unlock', signature(m='boost.mutex'),
  function(m, ...)
  {
    block_call = ifelse(read(m), boost_unlock_shared, boost_unlock)
    block_call(m@mutexInfoAddr)
  })

#' @export
setGeneric('shared.name', function(m) standardGeneric('shared.name'))

#' @export
setMethod('shared.name', signature(m='boost.mutex'), 
  function(m) GetResourceName(m@mutexInfoAddr) )

#' @export
setGeneric('timeout', function(m) standardGeneric('timeout'))

#' @export
setMethod('timeout', signature(m='boost.mutex'),
  function(m) GetTimeout(m@mutexInfoAddr))

#' @export
setGeneric('is.timed', function(m) standardGeneric('is.timed'))

#' @export
setMethod('is.timed', signature(m='boost.mutex'),
  function(m)
  {
    return(!is.null(timeout(m)))
  })

#' @export
boost.mutex=function(sharedName=NULL, timeout=NULL)
{
  isRead = TRUE
  if (is.null(sharedName)) 
  {
    sharedName = uuid()
  }
  if (!is.null(timeout) && !is.numeric(timeout))
  {
    stop("The timeout parameter must be numeric.")
  }
  if (is.numeric(timeout) && timeout <= 0)
  {
    stop("You must specify a timeout greater than zero.")
  }
  mutexInfoAddr=CreateBoostMutexInfo(sharedName, as.double(timeout))
  return(new('boost.mutex', isRead=isRead, mutexInfoAddr=mutexInfoAddr))
}


#' @title An S4 class holding mutex description information.
#' 
#' @description Objects of class description allow users to ``attach'' 
#' to existing mutexes within or across processes.
#' @slot description the list of description information.
#' @export
setClass('descriptor', representation(description='list'))

#' @title Accessor for descriptor objects
#' 
#' @description Retrieve the list of description information from a 
#' descriptor object.
#' @docType methods
#' @rdname description-methods
#' @param x the descriptor object.
#' @return a list of description information.
#' @export
setGeneric('description', function(x) standardGeneric('description'))

#' @rdname description-methods
#' @aliases description,descriptor-method
#' @export
setMethod('description', signature(x='descriptor'),
  function(x) return(x@description))

#' @title An S4 class holding boost.mutex description information.
#' 
#' @description Objects of class description allow users to ``attach'' 
#' to existing mutexes within or across processes.
#' @slot description the list of description information.
#' @export
setClass('boost.mutex.descriptor', contains='descriptor')

#' @importFrom bigmemory.sri describe
#' @import methods
#' @export
setMethod('describe', signature(x='boost.mutex'),
  function(x)
  {
    return(new('boost.mutex.descriptor',
      description=list(shared.name=shared.name(x),
      timeout=timeout(x))))
  })

#' @title Attach to an existing mutex.
#' 
#' @description Attach to an existing mutex using either a file or 
#' description object
#' @docType methods
#' @rdname attach.mutex-methods
#' @param obj the descriptor object.
#' @param ... other arguments needed by attach.
#' @return A mutex.
#' @export
setGeneric('attach.mutex', function(obj, ...) 
  standardGeneric('attach.mutex'))

#' @rdname attach.mutex-methods
#' @aliases attach.mutex,character-method
#' @export
setMethod('attach.mutex', signature(obj='character'),
  function(obj, ...)
  {
    path = match.call()[['path']]
    if (is.null(path))
    {
      path <- '.'
    }
    path <- path.expand(path)
    if (basename(obj) != obj)
    {

        warning(paste("Two paths were specified in attach.mutex",
          "The one associated with the file will be used.", sep="  "))
      path <- dirname(obj)
      obj <- basename(obj)
    }
    fileWithPath <- file.path(path, obj)
    fi = file.info(fileWithPath)
    print(dir())
    if (is.na(fi$isdir))
      stop( paste("The file", fileWithPath, "could not be found") )
    if (fi$isdir)
      stop( fileWithPath, "is a directory" )
    info <- dget(fileWithPath)
    return(attach.mutex(info, path=path))
  })

#' @rdname attach.mutex-methods
#' @aliases attach.mutex,boost.mutex.descriptor-method
#' @export
setMethod('attach.mutex', signature(obj='boost.mutex.descriptor'),
  function(obj, ...)
  {
    desc = description(obj)
    return(boost.mutex(sharedName = desc$shared.name,
      timeout = desc$timeout))
  })


