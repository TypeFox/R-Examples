# note, in this formulation pos should effectively be a new map with inserted positions, labelled as locx where x is the cM position value
calc.genoprob2 <- function (cross, pos, error.prob = 1e-04, map.function = c("haldane", "kosambi", "c-f", "morgan")) 
{
    ## fixed at 0 instead of an input argument
    off.end <- 0
    if (!any(class(cross) == "cross")) 
        stop("Input should have class \"cross\".")
    map.function <- match.arg(map.function)
    if (map.function == "kosambi") 
        mf <- mf.k
    else if (map.function == "c-f") 
        mf <- mf.cf
    else if (map.function == "morgan") 
        mf <- mf.m
    else mf <- mf.h

    ## removed to allow for more general formulation of positions
#    stepwidth <- match.arg(stepwidth)
    if (error.prob < 1e-50) 
        error.prob <- 1e-50
    if (error.prob > 1) {
        error.prob <- 1 - 1e-50
        warning("error.prob shouldn't be > 1!")
    }
    n.ind <- nind(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)
    type <- class(cross)[1]
    for (i in 1:n.chr) {
        if (n.mar[i] == 1) 
            temp.offend <- max(c(off.end, 5))
        else temp.offend <- off.end
        chrtype <- class(cross$geno[[i]])
        if (chrtype == "X") 
            xchr <- TRUE
        else xchr <- FALSE
        if (type == "riself" || type == "risib" || type == 
            "dh") {
            cfunc <- "calc_genoprob_bc"
            n.gen <- 2
            gen.names <- getgenonames(type, "A", cross.attr = attributes(cross))
        }
        else if (type == "ri8sib" || type == "ri4sib" || type == 
            "ri8self" || type == "ri4self") {
            cfunc <- paste("calc_genoprob_", type, sep = "")
            n.gen <- as.numeric(substr(type, 3, 3))
            gen.names <- LETTERS[1:n.gen]
            if (xchr) 
                warning("calc.genoprob not working properly for the X chromosome for 4- or 8-way RIL.")
        }
        else stop("calc.genoprob not available for cross type ", 
            type, ".")
        gen <- cross$geno[[i]]$data
        gen[is.na(gen)] <- 0
# removed possibility of onemap == FALSE
#            map <- create.map(cross$geno[[i]]$map, step, temp.offend, 
#                stepwidth)
	    map <- pos[[i]]
            rf <- mf(diff(map))
            if (type == "risib" || type == "riself") 
                rf <- qtl:::adjust.rf.ri(rf, substr(type, 3, nchar(type)), 
                  chrtype)
            rf[rf < 1e-14] <- 1e-14
            newgen <- matrix(ncol = length(map), nrow = nrow(gen))
            dimnames(newgen) <- list(NULL, names(map))
            newgen[, colnames(gen)] <- gen
            newgen[is.na(newgen)] <- 0
            n.pos <- ncol(newgen)
            marnames <- names(map)
            z <- .C(cfunc, as.integer(n.ind), as.integer(n.pos), 
                as.integer(newgen), as.double(rf), as.double(error.prob), 
                genoprob = as.double(rep(0, n.gen * n.ind * n.pos)), 
                PACKAGE = "qtl")
        cross$geno[[i]]$prob <- array(z$genoprob, dim = c(n.ind, 
            n.pos, n.gen))
        dimnames(cross$geno[[i]]$prob) <- list(NULL, marnames, 
            gen.names)
        attr(cross$geno[[i]]$prob, "map") <- map
        attr(cross$geno[[i]]$prob, "error.prob") <- error.prob
        attr(cross$geno[[i]]$prob, "step") <- step
        attr(cross$geno[[i]]$prob, "off.end") <- temp.offend
        attr(cross$geno[[i]]$prob, "map.function") <- map.function
    }
    if (type == "ri4self" || type == "ri4sib" || type == "ri8self" || 
        type == "ri8sib") 
        cross <- reorgRIgenoprob(cross)
    cross
}

