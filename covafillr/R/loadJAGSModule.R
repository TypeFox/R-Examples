##' load JAGS module.
##'
##' Calls rjags::load.module with appropriate arguments to load the covafillr module.
##' @return Nothing
##' @seealso \code{\link[rjags]{load.module}}
##' @author Christoffer Moesgaard Albertsen
##' @examples
##' if(require("rjags"))
##'    loadJAGSModule()
##' @export
loadJAGSModule <- function(){
    if(requireNamespace("rjags", quietly = TRUE)){
        if(!.installed_with_jags)
            warning("covafillr was not installed using JAGS. Please re-install covafillr.")

        r_arch <- .Platform$r_arch
        rjags::load.module('covafillr',
                           system.file('libs',r_arch,
                                       package='covafillr'))
        
    }else{
        if(.installed_with_jags){
            stop("rjags must be installed to load the JAGS module.")
        }else{
            stop("rjags must be installed to load the JAGS module. Please re-install covafillr after installing rjags.")
        }
    }
}
