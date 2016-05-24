ceplot <- 
function (data, model, response = NULL, S = NULL, C = NULL, sigma = NULL, 
    distance = "euclidean", type = "default", cex.axis = NULL, cex.lab = NULL, 
    tck = NULL, view3d = FALSE, Corder = "default", selectortype = "minimal", 
    conf = FALSE, select.colour = "blue", select.cex = 1)
{
    data <- na.omit(data)
    model <- if (!identical(class(model), "list"))
        list(model)
    else model
    model.name <- if (!is.null(names(model)))
        names(model)
    else NULL
    varnamestry <- try(getvarnames(model[[1]]), silent = TRUE)
    response <- if (is.null(response))
        if (class(varnamestry) != "try-error")
           which(colnames(data) == varnamestry$response[1])
        else stop("could not extract response from 'model'.")
    else if (is.character(response))
            which(colnames(data) == response)
        else response
    S <- if(is.null(S)){
         (1:ncol(data))[-response][1L]
        } else if (is.character(S))
            vapply(S, function(x) which(colnames(data) == x), numeric(1))
            else S       
    C <- if (is.null(C))
        arrangeC(data[, -c(response, S)])
    else C
    try(
        if (class(varnamestry) != "try-error"){
            possibleC <- unique(unlist(lapply(lapply(model, getvarnames), `[[`, 
                2)))
            possibleC <- possibleC[possibleC %in% colnames(data)]
            C <- arrangeC(data[, possibleC[!(possibleC %in% colnames(data)[S])], 
                drop = FALSE], method = Corder)
        }     
    , silent = TRUE) 
    C <- if (all(vapply(C, is.numeric, logical(1))))
        as.list(C)
    else if (all(vapply(C, is.character, logical(1))))
            lapply(C, match, table = colnames(data))
        else
            stop("'C' should be a vector or list (containing vectors of length",
                 " 1 or 2) with integer column indices or character variable",
                 " names from 'data'.")
    uniqC <- unique(unlist(C))
    if (any(response %in% uniqC))
        stop("cannot have 'response' variable in 'C'")
    if (any(response %in% S))
        stop("cannot have 'response' variable in 'S'")
    if (!identical(length(intersect(S, uniqC)), 0L))
        stop("cannot have variables common to both 'S' and 'C'")
        
    if (identical(type, "default")){
        ceplot.interactive(data = data, model = model, response = response, 
            S = S, C = C, sigma = sigma, distance = distance, cex.axis = 
            cex.axis, cex.lab = cex.lab, tck = tck, view3d = view3d, Corder = 
            Corder, conf = conf, separate = FALSE, select.colour = 
            select.colour, select.cex = select.cex)
    } else if (identical(type, "separate") & identical(selectortype, 
        "minimal")){
        ceplot.interactive(data = data, model = model, response = response, 
            S = S, C = C, sigma = sigma, distance = distance, cex.axis = 
            cex.axis, cex.lab = cex.lab, tck = tck, view3d = view3d, Corder = 
            Corder, conf = conf, separate = TRUE, select.colour = select.colour, 
            select.cex = select.cex)
    } else if (identical(type, "separate")){
        ceplot.separate(data = data, model = model, response = response, S = S, 
            C = C, sigma = sigma, distance = distance, cex.axis = cex.axis, 
            cex.lab = cex.lab, tck = tck, view3d = view3d, Corder = Corder, 
            selectortype = selectortype, select.colour = select.colour, 
            select.cex = select.cex)
    } else if (identical(type, "shiny")){
        ceplot.shiny(data = data, model = model, response = response, S = S, 
            C = C, cex.axis = cex.axis, cex.lab = cex.lab, tck = tck, 
            Corder = Corder)
    }
}