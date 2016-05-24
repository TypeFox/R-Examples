# note, in this formulation pos should effectively be a new map with inserted positions, labelled as locx where x is the cM position value
calc.genoprob2 <- 
function(cross, pos, error.prob = 1e-04, map.function = c("haldane", "kosambi", "c-f", "morgan")) 
{
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
        if (type == "f2") {
            if (!xchr) {
                cfunc <- "calc_genoprob_f2"
                n.gen <- 3
                gen.names <- getgenonames("f2", "A", cross.attr = attributes(cross))
            }
            else {
                cfunc <- "calc_genoprob_bc"
                n.gen <- 2
                gen.names <- c("g1", "g2")
            }
        }
        else if (type == "bc") {
            cfunc <- "calc_genoprob_bc"
            n.gen <- 2
            if (!xchr) 
                gen.names <- getgenonames("bc", "A", cross.attr = attributes(cross))
            else gen.names <- c("g1", "g2")
        }
        else if (type == "riself" || type == "risib" || type == 
            "dh") {
            cfunc <- "calc_genoprob_bc"
            n.gen <- 2
            gen.names <- getgenonames(type, "A", cross.attr = attributes(cross))
        }
        else stop("calc.genoprob not available for cross type ", 
            type, ".")
        gen <- cross$geno[[i]]$data
        gen[is.na(gen)] <- 0

        rf <- mf(diff(pos[[i]]))
        if (type == "risib" || type == "riself") 
           rf <- qtl:::adjust.rf.ri(rf, substr(type, 3, nchar(type)), chrtype)
        rf[rf < 1e-14] <- 1e-14
        newgen <- matrix(ncol = length(pos[[i]]), nrow = nrow(gen))
        dimnames(newgen) <- list(NULL, names(pos[[i]]))
        newgen[, colnames(gen)] <- gen
        newgen[is.na(newgen)] <- 0
        n.pos <- ncol(newgen)
        marnames <- names(pos[[i]])

        z <- .C(cfunc, as.integer(n.ind), as.integer(n.pos), 
                as.integer(newgen), as.double(rf), as.double(error.prob), 
                genoprob = as.double(rep(0, n.gen * n.ind * n.pos)), 
                PACKAGE = "qtl")
        cross$geno[[i]]$prob <- array(z$genoprob, dim = c(n.ind, 
            n.pos, n.gen))
        dimnames(cross$geno[[i]]$prob) <- list(NULL, marnames, 
            gen.names)
        attr(cross$geno[[i]]$prob, "map") <- pos[[i]]
        attr(cross$geno[[i]]$prob, "error.prob") <- error.prob
        attr(cross$geno[[i]]$prob, "map.function") <- map.function
    }
    cross
}

