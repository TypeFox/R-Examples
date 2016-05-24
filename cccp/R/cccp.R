##
## Main function for defining and solving linear and quadratic programs with cone constraints
cccp <- function(P = NULL, q = NULL, A = NULL, b = NULL, cList = list(),
                 x0 = NULL, f0 = NULL, g0 = NULL, h0 = NULL,
                 nlfList = list(), nlgList = list(), nlhList = list(),
                 optctrl = ctrl()){

    if(!is.null(x0)){
        if(!is.null(f0) && !is.null(g0) && !is.null(h0)){
            cpd <- dcp(x0 = x0, f0 = f0, g0 = g0, h0 = h0, cList = cList,
                       nlfList = nlfList, nlgList = nlgList, nlhList = nlhList,
                       A = A, b = b)
        } else if((length(nlfList) > 0) && (length(nlgList) > 0) && (length(nlhList) > 0)){
            cpd <- dnl(q = q, A = A, b = b, cList = cList,
                       x0 = x0, nlfList = nlfList, nlgList = nlgList, nlhList = nlhList)
        } else {
            warning("x0 provided, but missing arguments for either:\nf0, g0, h0 and/or nlfList, nlgList, nlhList.\nDiscarding x0.\n")
            if(is.null(P) && is.null(q)){
                stop("Ill-defined program formulation: At least P or q must be provided.\n")
            }
            if(is.null(P) && !is.null(q)){
                cpd <- dlp(q = q, A = A, b = b, cList = cList)
            } else {
                cpd <- dqp(P = P, q = q, A = A, b = b, cList = cList)
            }
        }
    } else if(is.null(P)){
        cpd <- dlp(q, A, b, cList)
    } else {
        cpd <- dqp(P, q, A, b, cList)
    }

    cps(cpd, optctrl)
}
