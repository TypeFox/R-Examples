`profConfint` <-
function (prof, ...) 
UseMethod("profConfint")
`profConfint.profileModel` <-
function (prof, method = "smooth", endpoint.tolerance = 0.001, 
    max.zoom = 100, n.interpolations = 100, verbose = FALSE, 
    ...) 
{
    if (is.null(prof$quantile)) 
        stop("The profiling object does not have a non-NULL quantile.")
    switch(method, zoom = ci <- profZoom.profileModel(prof = prof, 
        endpoint.tolerance = endpoint.tolerance, max.zoom = max.zoom, 
        verbose = verbose), smooth = ci <- profSmooth.profileModel(prof = prof, 
        n.interpolations = n.interpolations), stop("Invalid method. The supported methods are 'smooth' and 'zoom'"))
    attr(ci, "profileModel object") <- match.call()[["prof"]]
    ci
}
