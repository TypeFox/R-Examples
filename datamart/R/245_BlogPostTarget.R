#' A Target representing a blog post
#' 
#' This is an internal class representing a blog post. Use
#' MdReport instead.
#'
#' @seealso \code{\link{blogpost}}
#'  
#' @examples
#' getSlots("BlogPostTarget")
#'
#' @name BlogPostTarget-class
#' @rdname BlogPostTarget-class
setClass(
    Class="BlogPostTarget", 
    representation=representation(
        subject="character",
        content="character",
        label="character",
        draft="logical",
        overwrite="logical"
    ),
    contains="Target"
)

#' Constructor for BlogPostTarget objects
#'
#' For internal use only
#'
#' @param name        short (file) name of the blogpost
#' @param subject     title of the blogpost
#' @param content     content of the post
#' @param label       character vector of keywords
#' @param draft       draft or not? default=TRUE
#' @param overwrite   overwrite or not? overwrite=TRUE
#' @param clss        class of the object, default 'BlogPostTarget'
#'
#' @return BlogPostTarget
#' @rdname BlogPostTarget-class
blogpost <- function(name, subject, content, label="", draft=TRUE, overwrite=TRUE, clss="BlogPostTarget") {
  new(clss, name=name, subject=subject, content=content, label=label, draft=draft, overwrite=overwrite)
}


#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,BlogPostTarget,DirectoryLocation-method
setMethod(
    f="put",
    signature=c(target="BlogPostTarget", where="DirectoryLocation"),
    definition=function(target, where, ...) {
        fn <- file.path(where@path, paste(target@name, ".html", sep=""))
        con <- file(fn, encoding="UTF-8")
        writeLines(target@content, con) # FIXME: make sure file is  UTF-8 encoded
        close(con)
        return(fn)
    }
)

#' @rdname put-methods
#' @name put
#' @export
#' @docType methods
#' @aliases put put,BlogPostTarget,SftpLocation-method
setMethod(
    f="put",
    signature=c(target="BlogPostTarget", where="SftpLocation"),
    definition=function(target, where, ...) {
        fn <- put(target, tempdir()); on.exit(unlink(fn))
        ft <- filetarget(name=basename(fn), filename=fn)
        put(ft, where)
    }
)

