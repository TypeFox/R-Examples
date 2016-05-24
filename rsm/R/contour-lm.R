### Reconstructs the data set used in a linear model, and
### returns it as a data.frame

model.data = function (lmobj, lhs = FALSE) {
    form = lmobj$call$formula
    if (is.name(form)) {
        lmobj$call$data = form
        form = formula(lmobj)
    }
    if (lhs) 
        nm = all.vars(form)
    else nm = all.vars(form[[3]])
    if (inherits(lmobj, "rsm") && !is.null(lmobj$data))
        lmobj$data[ ,nm]
    else {
        form = as.formula(paste("~", paste(nm, collapse = "+")))
        envir = attr(lmobj$terms, ".Environment")
        model.frame(form, eval(lmobj$call$data, envir=envir), 
            subset = eval(lmobj$call$subset, envir=envir))
    }
}


### contour plot(s) for a lm
contour.lm = function(x, form, at, bounds, zlim, 
                      xlabs, hook, plot.it=TRUE, atpos = 1, decode=TRUE,
                      image=FALSE, img.col=terrain.colors(50), ...) 
{
    if (!missing(at)) {
        if (is.null(names(at)))
            stop("'at' must be a NAMED list of values used for surface slices")
    }
    if (!missing(bounds)) {
        if (is.null(names(bounds)))
            stop("'bounds' must be a NAMED list of bounds for each variable")
    }
    
    # return decoded values from vector cval for coded variable named cname
    # made to be transparent to uncoded situations
    .decode.value = function(cname, cval) {
        if (decode && !is.null(forms)) {
            if (is.null(form <- forms[[cname]]))
                cval
            else {
                inf = .parse.coding(forms[[cname]])
                inf$const[1] + cval * inf$const[2]
            }
        }
        else
            cval
        
    }
    
    lmobj = x   # generic wants it named 'x', I don't!
    
    data = model.data(lmobj)
    
    # make list of formulas if not already
    if(inherits(form,"formula")) {
        if (length(form)==2) { # rhs only
            vars = all.vars(form[[2]])
            n = length(vars)
            if (n < 2) stop("Need at least two variables")
            form = list()
            elt = 1
            for (i in 1:(n-1))
                for (j in (i+1):n) {
                    form[[elt]] = formula(paste(vars[j],vars[i],sep="~"))
                    elt = elt + 1
                }
        }
        else {
            yvars = all.vars(form[[2]])
            xvars = all.vars(form[[3]])
            form = list()
            elt = 1
            for (i in 1:length(xvars))
                for (j in 1:length(yvars)) {
                    form[[elt]] = formula(paste(yvars[j],xvars[i],sep="~"))
                    elt = elt + 1
                }
        }
    }
    vars = sort(unique(as.character(sapply(form, all.vars))))
    
    dots = list(...)
    if (!is.null(dots$ylab))
        message("'ylab' ignored. Specify axis labels using 'xlabs'")
    if(!missing(xlabs) && length(xlabs) < length(vars))
        stop("'xlabs' does not contain enough labels" )
    
    forms = NULL  # for all non-rsm objects
    if (inherits(lmobj, "rsm")) {
        forms = codings(lmobj)
        if (missing(xlabs) && !is.null(forms)) {
            if(!decode)
                xlabs = sapply(vars, function(v) 
                    paste(as.character(forms[[v]][2:3]), collapse=" = "))
            else {
                xlabs = sapply(vars, function(v) all.vars(forms[[v]][[3]])[1])
            }
        }
        else if (missing(xlabs))
            xlabs = vars
    }   
    else if (missing(xlabs)) 
        xlabs = vars
    
    
    # gather 'at' info
    tmp = lapply(data, function(var) {
        if (is.factor(var)) factor(levels(var)) ####NEW [1], levels=levels(var))
        else mean(var)
    })
    # remember original at list
    orig.atnm = NULL
    if (!missing(at)) {
      orig.atnm = names(at)
      for (nm in orig.atnm) {
        numflag = is.numeric(tmp[[nm]])
        if (numflag) tmp[[nm]] = as.numeric(at[[nm]])
        else         tmp[[nm]] = at[[nm]]
      }
    }
    at = tmp
    
    # gather 'bounds' info -- elts can be vectors of length 2, 3, or n
    tmp = lapply(data, function(x) if (is.numeric(x)) range(x))
    if (!missing(bounds))
      for (nm in names(bounds))
        if (length(bounds[[nm]]) > 1)
          tmp[[nm]] = bounds[[nm]]
    bounds = lapply(tmp, function(x) {
      if (length(x) == 2) seq(x[1], x[2], length=26)
      else if (length(x) == 3) seq(x[1], x[2], length=x[3])
      else x
    })
    
    # get names to use in slice labels
     isnum = sapply(at, is.numeric) # was commented-out
     allnum = names(at)[isnum]      # ditto
#     allfac = names(at)[!isnum]
    
    ### Accumulate the z values
    plot.data = list()
    lbls = rep("", length(form))
    z.rng = NULL
    for (i in 1:length(form)) {
        AT = at
        v = all.vars(form[[i]])
        if (length(unique(v)) == 1) next
        y = AT[[v[1]]] = bounds[[v[1]]]
        x = AT[[v[2]]] = bounds[[v[2]]]
        newdata = do.call(expand.grid, AT)
        ord = order(newdata[[v[1]]], newdata[[v[2]]])
        newdata = newdata[ord, ]
        z = predict (lmobj, newdata = newdata)
# NEW: average over factor levels...
        rep = length(z) / length(x) / length(y)  # copies at each (x,y) pair
        if (rep > 1) 
            z = apply(matrix(z, nrow=rep), 2, mean)
        
        
        z.rng = range (c (z.rng, z))
        if (!missing(zlim)) {
            z[z > zlim[2]] = NA
            z[z < zlim[1]] = NA
        }
        vnames = c(x=v[2], y=v[1])
        labs = c(xlabs[sapply(vnames, charmatch, vars)], vnames, "")
        lbls[i] = paste(labs[3], labs[4], sep=" ~ ")
        
        # figure out slice labels
        if (atpos != 0) {
            #atidx = - sapply(vnames, grep, allnum)
            #atidx = allnum[atidx] 
            # REPLACEMENT:
            atidx = setdiff(allnum, vnames)  ## was vnames in rsm2.0 -- NOT correct!
            if (length(atidx) > 0) { ###### || length(allfac) > 0) {
                atlabs = NULL
                if (length(atidx) > 0) {
                    atvals = round(sapply(atidx, function(v) .decode.value(v, at[[v]])), 2)
                    if (decode && !is.null(forms)) {
                        dclabs = sapply(atidx, function(x) {
                            f = forms[[x]]
                            if (is.null(f)) x
                            else {
                                info = .parse.coding(f)
                                info[[1]][2]
                            }
                        })
                        atlabs = paste(dclabs, atvals, sep=" = ")
                    }
                    else
                        atlabs = paste(atidx, atvals, sep = " = ")
                }
                # added  for factors in 'at'
                facidx = setdiff(orig.atnm, atidx)
                faclabs = paste(facidx, unlist(at[facidx]), sep = " = ")
                atlabs = c(atlabs, faclabs)
#                 faclabs = paste(allfac, sapply(allfac, function(v) at[[v]]), sep=" = ")
#                 atlabs = paste(c(atlabs, faclabs), collapse = ", ")
                atlabs = paste(atlabs, collapse = ", ") ### NEW
                labs[5] = paste("Slice at", atlabs)
                if (atpos < 3)
                    labs[atpos] = paste(labs[atpos], "\n", labs[5], sep = "")
            }
        }
        y = .decode.value(v[1], y) 
        x = .decode.value(v[2], x)
        
        plot.data[[i]] = list(x=x, y=y, 
                              z=matrix(z, nrow=length(x)), labs=labs)
    }
    names(plot.data) = lbls
    
    if (missing (zlim)) zlim = z.rng
    for (i in 1:length(lbls))
        plot.data[[i]]$zlim = zlim
    
    ### If plots requested, do plots with a common image scale
    if (plot.it) for (i in 1:length(form)) {
        dat = plot.data[[i]]
        if (!missing(hook))
            if (!is.null(hook$pre.plot)) hook$pre.plot(dat$labs)
        args = list(x=dat$x, y=dat$y, z=dat$z, col=img.col, zlim = dat$zlim, ...)
        args$xlab = dat$labs[1]
        args$ylab = dat$labs[2]
        if (image) {
          do.call("image", args)
          args$add = TRUE
          args$xlab = args$ylab = args$col = NULL
          do.call("contour", args)
            #image(dat$x, dat$y, dat$z, col=img.col, 
            #      xlab = dat$labs[1], ylab = dat$labs[2], zlim = dat$zlim, ...)
            #contour(dat$x, dat$y, dat$z, add=TRUE, ...)
        }
        else {
            args$col = NULL
            do.call("contour", args)
            #contour(dat$x, dat$y, dat$z, 
            #        xlab = dat$labs[1], ylab = dat$labs[2], zlim = dat$zlim, ...)
            if (atpos == 3)
                title(sub = labs[5])
        }
        if (!missing(hook))
            if (!is.null(hook$post.plot)) hook$post.plot(dat$labs)
    }
    
    invisible(plot.data)
}


# Image plot for a lm
image.lm = function(x, form, at, bounds, zlim, xlabs, hook, atpos=1, decode=TRUE, ...)  {
    plot.data = contour.lm(x, form, at, bounds, zlim, xlabs, atpos=atpos, decode=decode, plot.it=FALSE)
    for (i in 1:length(plot.data)) {
        dat = plot.data[[i]]
        if (!missing(hook))
            if (!is.null(hook$pre.plot)) hook$pre.plot(dat$labs)
        
        image(dat$x, dat$y, dat$z, 
              xlab = dat$labs[1], ylab = dat$labs[2], zlim = dat$zlim, ...)
        if (atpos == 3)
            title(sub = dat$labs[5])
        
        if (!missing(hook))
            if (!is.null(hook$post.plot)) hook$post.plot(dat$labs)
    }
    
    invisible(plot.data)
}


# Perspective plot(s) for 'lm' objects
# arg notes:
# col: facet colors; if null, default, else color palette based on z value
# contours: if TRUE, black contours.  Can also be a list with elements
#   z="bottom" (or "top" or value), col="black", lwd=1
persp.lm = function(x, form, at, bounds, zlim, zlab, 
                    xlabs, col = "white", contours=NULL, hook, atpos=3, decode = TRUE,
                    theta = -25, phi = 20, r = 4, border = NULL, box = TRUE,
                    ticktype = "detailed", ...) 
{
    draw.cont.line = function(line) {
        if (cont.varycol) {
            cont.col = col
            if (length(col) > 1) cont.col = col[cut(c(line$level, dat$zlim), length(col))][1]
        }
        lines(trans3d(line$x, line$y, cont.z, transf),
              col=cont.col, lwd=cont.lwd)
    }
    plot.data = contour.lm(x, form, at, bounds, zlim, xlabs, atpos=atpos, plot.it=FALSE)
    transf = list()
    if (missing(zlab)) zlab = ""
    
    facet.col = col
    
    cont = !is.null(contours)
    if (mode(contours) == "logical") cont = contours
    cont.first = cont
    cont.z = cz = plot.data[[1]]$zlim[1]
    cont.col = 1
    cont.varycol = FALSE
    cont.lwd = 1
    if (is.character(contours)) {
        idx = charmatch(contours, c("top","bottom", "colors"), 0)
        if (idx == 1) {
            cont.first = FALSE
            cont.z = plot.data[[1]]$zlim[2]
        }
        else if (idx == 2) {}
        else if (idx == 3) {
            cont.varycol = TRUE
            if (length(col) < 2) col = rainbow(40)
        }
        else
            cont.col = contours
    }
    else if (is.list(contours)) {
        if(!is.null(contours$z)) cz = contours$z
        if (is.numeric(cz)) cont.z = cz
        else if (cz=="top") {
            cont.first = FALSE
            cont.z = plot.data[[1]]$zlim[2]
        }
        if(!is.null(contours$col)) cont.col = contours$col
        if(!is.null(contours$lwd)) cont.lwd = contours$lwd
        if(charmatch(cont.col, "colors", 0) == 1) {
            cont.varycol = TRUE
            if (length(col) < 2) col = rainbow(40)
        }
    }
    
    # Loop through the plots
    for (i in 1:length(plot.data)) {
        dat = plot.data[[i]]
        cont.lines = NULL
        if (!missing(hook))
            if (!is.null(hook$pre.plot)) hook$pre.plot(dat$labs)
        if (cont) cont.lines = contourLines(dat$x, dat$y, dat$z)
        if (cont && cont.first) {
            transf = persp(dat$x, dat$y, dat$z, zlim=dat$zlim, theta=theta, phi=phi, r=r, col = NA, border=NA, box=FALSE, ...)
            lapply(cont.lines, draw.cont.line)
            par(new=TRUE)
        }
        if (length(col) > 1) {
            nrz = nrow(dat$z)
            ncz = ncol(dat$z)
            zfacet = dat$z[-1,-1] + dat$z[-1,-ncz] + dat$z[-nrz,-1] + dat$z[-nrz,-ncz]
            zfacet = c(zfacet/4, dat$zlim)
            facet.col = cut(zfacet, length(col))
            facet.col = col[facet.col]
        }
        transf = persp(dat$x, dat$y, dat$z,
                       xlab=dat$labs[1], ylab=dat$labs[2], zlab=zlab,
                       zlim=dat$zlim, col=facet.col, border=border, box=box, theta=theta, phi=phi, r=r, ticktype=ticktype, ...)
        if (atpos == 3)
            title(sub = dat$labs[5], ...)
        
        if (cont && !cont.first)
            lapply(cont.lines, draw.cont.line)
        if (!missing(hook))
            if (!is.null(hook$post.plot)) hook$post.plot(dat$labs)
        plot.data[[i]]$transf = transf
    }
    invisible(plot.data)
}
