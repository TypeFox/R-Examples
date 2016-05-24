#' @title Calculate marginal correlations with response
#' @description Calculates marginal correlations between a functional covariate and a scalar response.
#' @param object An object of type \code{funeigen} or \code{funreg}.  One or the other of these is 
#' needed in order to provide a smoothed reconstructed curves for the functional covariate
#' for each subject.
#' @param id The vector of subject id's. These tell which responses in \code{response} 
#' correspond to which curves in \code{object}.
#' @param response The vector of responses
#' @param alpha The alpha level for confidence intervals (one minus the two-sided coverage)
#' @return Returns a list with one component for each functional 
#' covariate. Each such component contains the between-subjects correlations
#' between the fitted smoothed latent values of the functional covariate,
#' and the response variable. We call this a marginal correlation because
#' it simply ignores the other functional covariates (rather than trying to
#' adjust or control for them).  Both the functional regression coefficient
#' and the marginal correlation can be useful, although they have different
#' substantive interpretations.
#'@export
marginal.cor <- function(object,          
                         id=NULL,        
                         response=NULL,  
                         alpha=.05       
                          ) {
    if ((class(object)=="funeigen")) {
        return(marginal.cor.funeigen(object=object,
                                     id=id,
                                     response=response,
                                     alpha=alpha));
    }
    if ((class(object)=="funreg")) {
        if (!is.null(id)) {
            stop(paste("A separate id vector should not be specified in",
                       "this context because it is already included in",
                       "the funreg object."));
        }
        if (!is.null(response)) {
            stop(paste("A separate response vector should not be specified in",
                       "this context because it is already included in",
                       "the funreg object."));
        }
        answer <- list();
        for (j in 1:length(object$object.for.x.functions)) {
            answer[[j]] <- marginal.cor.funeigen(object=object$object.for.x.functions[[j]],
                                                 id=object$data$id,
                                                 response=object$data$response,
                                                 alpha=alpha)
        }
        return(answer);
    }
}