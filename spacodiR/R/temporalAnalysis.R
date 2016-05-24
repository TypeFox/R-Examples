#####
spacodi.all.nodes <- function(sp.plot, phy, sp.parm, return.all = TRUE) {
	sub = subtrees(phy)
	n.nodes <- phy$Nnode
	out <- array(dim = c(length(sub), 4))
	for (i in 1:n.nodes) {
		node = min(sub[[i]]$node)
		time = max(branching.times(sub[[i]]))
		tips = length(sub[[i]]$tip.label)
		ss <- try(spacodi.calc(sp.plot = sp.plot, phy = sub[[i]], prune=TRUE),silent=TRUE)
		ss <- spacodi.calc(sp.plot = sp.plot, phy = sub[[i]], prune=TRUE)

		if(!inherits(ss,"try-error")) {
			if(sp.parm%in%names(ss)) out[i, 1] <- ss[[sp.parm]] else stop("Cannot locate specified parameter")
		} 
		out[i, 2] <- tips
		out[i, 3] <- time
		out[i, 4] <- node
	}
	array.cols <- c(sp.parm, "tips", "node.time", "node.ID")
	dimnames(out) = list(NULL, array.cols)
	if (!return.all) {
		out = out[which(!is.na(out[, 1])), ]
	}
	return(out)
}

#####
spacodi.exp.nodes <- function(sp.plot, phy, n.rep = 10, method = "1s", parm = NULL, dmat = NULL, sp.parm) {
	sub = subtrees(phy)
	out <- array(dim = c(length(sub), 4, n.rep))
	if (method == "1s") {
		for (j in 1:n.rep) {
			out[, , j] <- foo <- spacodi.all.nodes(phy = phy, sp.plot = resamp.1s(sp.plot), sp.parm=sp.parm)
		}
	}
	else if (method == "1a" && !is.null(parm$abund.class.ratio)) {
		for (j in 1:n.rep) {
			out[, , j] <- foo <- spacodi.all.nodes(phy = phy, sp.plot = resamp.1a(sp.plot, parm$abund.class.ratio), sp.parm=sp.parm)
		}
	}
	else if (method == "2s") {
		for (j in 1:n.rep) {
			out[, , j] <- foo <- spacodi.all.nodes(phy = phy, sp.plot = resamp.2s(sp.plot), sp.parm=sp.parm)
		}
	}
	else if (method == "2x" && !is.null(parm$level)) {
		for (j in 1:n.rep) {
			out[, , j] <- foo <- spacodi.all.nodes(phy = phy, sp.plot = resamp.2x(sp.plot, parm$level), sp.parm=sp.parm)
		}
	}
	else if (method == "3i") {
		for (j in 1:n.rep) {
			out[, , j] <- foo <- spacodi.all.nodes(phy = phy, sp.plot = resamp.3i(sp.plot), sp.parm=sp.parm)
		}
	}
	else if (method == "3t") {
		if (!is.null(dmat)) 
		dmat = dmat
		else dmat = NULL
		for (j in 1:n.rep) {
			out[, , j] <- foo <- spacodi.all.nodes(phy = phy, sp.plot = resamp.3t(sp.plot, dmat), sp.parm=sp.parm)
		}
	}
	else if (method == "3x" && !is.null(parm$level)) {
		for (j in 1:n.rep) {
			out[, , j] <- foo <- spacodi.all.nodes(phy = phy, sp.plot = resamp.3x(sp.plot, parm$level), sp.parm=sp.parm)
		}
	}
	else stop(cat("Unrecognized method declaration or insufficient information supplied:\n\tIs 'parm' set to a non-null value?\n\n"))
	dimnames(out) = list(NULL, names(as.data.frame(foo)), paste("iter", seq(1:n.rep), sep = "."))
	return(out)
}

#####
spacodi.by.nodes<-function (sp.plot, phy, sp.parm="Bst", obs.only = FALSE, return.all = TRUE, 
n.rep = 10, method = "1s", parm = NULL, dmat = NULL, rand.test = TRUE, 
r.rep = 10000, ...) 
{
    sp.data = match.spacodi.data(sp.plot = sp.plot, phy = phy, ...)
    sp.plot = sp.data$sp.plot
    phy = sp.data$sp.tree
    if (is.null(sp.plot) && is.null(phy)) {
        stop("Must supply both an sp.plot and tree.")
    }
    if (is.null(r.rep) && rand.test) 
	stop("Must specify 'r.rep' if using the randomization test.")
    if (obs.only == TRUE) {
        rand.test = FALSE
        exp.Bst = FALSE
    } else {
        exp.Bst = TRUE
    }
    
    o.foo = spacodi.all.nodes(sp.plot = sp.plot, phy = phy, sp.parm = sp.parm, return.all = ifelse(obs.only == TRUE, TRUE, FALSE))
    o = o.foo[, sp.parm]
    names(o) = o.foo[, "node.ID"]
    o = o[order(o)]
    n.o = names(o)
    o.orig = as.data.frame(o.foo)
    row.names(o.orig) = o.orig$node.ID
    if (exp.Bst == TRUE) {
        ticker = c(1:35) * ceiling(n.rep/35)
        cat("\nSPACoDi calculations in progress:\n")
        e.out = array(dim = c(length(o), n.rep))
        for (bb in 1:n.rep) {
            b.foo = as.data.frame(spacodi.exp.nodes(sp.plot = sp.plot, phy = phy, n.rep = 1, method = method, parm = parm, sp.parm = sp.parm)[, , 1])
            b.match = match(as.numeric(n.o), b.foo$node.ID)
            e.out[, bb] = as.numeric(b.foo[b.match, sp.parm])
            if (bb %in% ticker) 
			cat(".")
        }
        cat("\n")
        e.out = as.data.frame(e.out)
        row.names(e.out) = as.numeric(n.o)
        names(e.out) = paste("iter", seq(1:n.rep), sep = "")
        orig.match <- match(as.numeric(o.orig$node.ID), as.numeric(row.names(e.out)))
        e.out = e.out[orig.match, ]
        if (rand.test == TRUE) {
            rand.array = array(dim = c(nrow(o.orig), 5))
            for (node in 1:nrow(o.orig)) {
                obs = o.orig[[sp.parm]][node]
                nn = o.orig$node.ID[node]
                exp = e.out[node, which(!is.na(e.out[node, ]))]
                rand.array[node, 1] = ifelse(length(exp) != 0, 
											 resamp.test(obs = obs, exp = exp, 
																   iter = r.rep, two.tailed = TRUE)[[2]], 
											 NA)
                rand.array[node, 2] = nn
                rand.array[node, 3] = obs
                rand.array[node, 4] = ifelse(length(exp) != 0, 
											 mean(unlist(exp)), NA)
                rand.array[node, 5] = length(exp)
            }
            r.test = as.data.frame(rand.array)
            names(r.test) = c("p.value", "node.ID", paste("obs",sp.parm,sep="."), paste("m.exp",sp.parm,sep="."), "valid.comparisons")
            row.names(r.test) = r.test$node.ID
            if (return.all == TRUE) {
                o.orig = as.data.frame(spacodi.all.nodes(sp.plot = sp.plot, phy = phy, return.all = TRUE, sp.parm = sp.parm))
                row.names(o.orig) = o.orig$node.ID
            }
            else {
                o.orig = o.orig
            }
            Bst.out = list(o.orig, e.out, r.test)
            names(Bst.out) = c(paste("observed",sp.parm,sep="."), paste("expected",sp.parm,sep="."), "randomization.test")
            return(Bst.out)
        }
        else {
            if (return.all == TRUE) {
                o.orig = as.data.frame(spacodi.all.nodes(sp.plot = sp.plot, phy = phy, return.all = TRUE, sp.parm = sp.parm))
                row.names(o.orig) = o.orig$node.ID
            }
            else {
                o.orig = o.orig
            }
            Bst.out = list(o.orig, e.out)
            names(Bst.out) = c(paste("observed",sp.parm,sep="."), paste("expected",sp.parm,sep="."))
            return(Bst.out)
        }
    }
    else {
        Bst.out = list(o.orig)
        names(Bst.out) = paste("observed",sp.parm,sep=".")
        return(Bst.out)
    }
}




