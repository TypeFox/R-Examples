# Helper functions for plotting.

.GraphParamFromSpec <- function(property, spec, default)
{
    # Returns a vector of parameters. 

    # property - the property by which to decide parameters. If NULL, default 
    # is returned and the returned vector will be of length one. 

    # spec - a named vector mapping elements in property to parameter values. 
    # An error will be raised if any element in property is not in spec.
    # Spec can contain an unnamed default parameter, which matches 
    # is.na(property) | 0==property | ''==property. 
    # If NULL and propery is not NULL, the returned vector is 
    # as.integer(as.factor(property)).

    # default should be NULL or of length one. 

    stopifnot(1==length(default) || is.null(default))
    if(!is.null(property))
    {
        if(is.vector(spec))
        {
            # A named vector that is a mapping from property values to 
            # graphical parameter values.
            # Property values that are 0, empty string, .UnnamedString() or NA 
            # match unnamed values in spec.
            stopifnot(length(unique(names(spec)))==length(spec))
            param <- spec[property]

            # Fix empty 
            if(any(''==names(spec)) && 
                   any(is.na(property) | 0==property | ''==property | 
                       property==.UnnamedString()))
            {
                d <- which(''==names(spec))
                param[is.na(property) | 0==property | ''==property | 
                      property==.UnnamedString()] <- spec[d]
            }

            missing <- is.na(param)   # !property %in% names(spec)

            # What is the right thing to do here?
            if(TRUE && any(missing))
            {
                stop(paste('Properties missing from spec:', 
                            '[', paste(unique(property[missing]),collapse=','), 
                            ']', sep=''))
            }
            else
            {
                param[missing] <- default
            }
        }
        else if(FALSE && is.function(spec))
        {
            # Is this a good idea?
            # A function that will return a vector of length(property)
            param <- spec(property)
        }
        else if(is.null(spec))
        {
            # Map logicals to TRUE=1, FALSE=2
            if(is.logical(property))
            {
                param <- rep(1, length(property))
                param[!property] <- 2
            }
            else
            {
                # Get an integer value for each of the unique values in property
                param <- as.integer(as.factor(property))
            }
        }
        else
        {
            stop(paste("Can't use object of class [", class(spec), 
                       "] as spec for graphical parameters", sep=''))
        }
    }
    else
    {
        param <- default
    }

    return (param)
}

.NodeGraphParamFromSpec <- function(community, param.name, param, property, 
                                   spec, default=NULL)
{
    if(is.null(param))
    {
        if(is.function(property))
        {
            property <- property(community)
            if(length(property)!=NumberOfNodes(community))
            {
                stop(paste('The function given for', param.name, 'returned', 
                           length(property), 'values but should return', 
                           'NumberOfNodes() values'))
            }
        }
        else if(1==length(property))
        {
            # This breaks if the community has only 1 node
            property <- NPS(community, property)[,property]
        }

        return (.GraphParamFromSpec(property, spec, default))
    }
    else
    {
        if(!length(param) %in% c(1, NumberOfNodes(community)))
        {
            stop(paste(param.name, 'contains', length(param), 'values but', 
                       'should be either of length 1 or of length', 
                       'NumberOfNodes()'))
        }
        else
        {
            return (param)
        }
    }
}

.GraphParNodeLabels <- function(community, node.labels)
{
    # Returns a vector of node labels
    # node.labels should be a function, NULL, a node property name or a vector 
    # of length NumberOfNodes()
    if(is.function(node.labels))
    {
        node.labels <- node.labels(community)
        if(length(node.labels)!=NumberOfNodes(community))
        {
            stop(paste('The function given for node.labels returned', 
                       length(node.labels), 'values but should return', 
                       'NumberOfNodes() values'))
        }
    }
    else if(is.null(node.labels))
    {
        node.labels <- 1:NumberOfNodes(community)
    }
    else if(1==length(node.labels))
    {
        # This breaks if the community has only 1 node
        node.labels <- NPS(community, node.labels)[,node.labels]
    }
    else if(length(node.labels)!=NumberOfNodes(community))
    {
        stop(paste('node.labels contains', length(node.labels), 'values', 
                   'but should be either of length 1 or of length', 
                   'NumberOfNodes()'))
    }

    return (node.labels)
}

.TrophicLinkGraphParamFromSpec <- function(community, param.name, param, 
                                           property, spec, default=NULL)
{
    if(is.null(param))
    {
        if(is.function(property))
        {
            property <- property(community)
            if(length(property)!=NumberOfTrophicLinks(community))
            {
                stop(paste('The function given for', param.name, 'returned', 
                           length(property), 'values but should return', 
                           'NumberOfTrophicLinks() values'))
            }
        }
        else if(1==length(property))
        {
            # This breaks if the community has only 1 node
            property <- .ResolveTrophicLinkProperty(community, property)
        }

        return (.GraphParamFromSpec(property, spec, default))
    }
    else
    {
        if(!length(param) %in% c(1, NumberOfTrophicLinks(community)))
        {
            stop(paste(param.name, 'contains', length(param), 'values but', 
                       'should be either of length 1 or of length', 
                       'NumberOfTrophicLinks()'))
        }
        else
        {
            return (param)
        }
    }
}

.BlendColour <- function(a, b, weight=0.5)
{
    # Returns a colour between a and b.
    stopifnot(0<weight && weight<1)

    a <- col2rgb(a)
    b <- col2rgb(b)
    lower <- pmin(a,b)
    upper <- pmax(a,b)
    colour <- lower + weight*(upper-lower)
    to.rgb <- function(c) { return (rgb(c[1],c[2],c[3],maxColorValue=0xff)) }
    return (apply(colour, FUN=to.rgb, MARGIN=2))
}

.ColourAlpha <- function(col, alpha=127)
{
    # Returns col with the alpha channel set
    # Either 0<=alpha<1, in which case it is taken as a proportion or 
    # 0<=alpha<256.
    stopifnot(0<=alpha && alpha<256)
    if(alpha<1)
    {
        alpha <- alpha * 256
    }

    col <- lapply(col, col2rgb)
    return (sapply(col, function(c) rgb(c[1], c[2], c[3], alpha=alpha, 
                                        maxColorValue=0xff)))
}

DefaultLinkColour <- function()
{
    return ('#c7c7c7')
}

DefaultCategoryColours <- function()
{
    return (c('black', 
              producer='#33bd2b', 
              invertebrate='#4d82ff', 
              vert.ecto='#c175f0', 
              vert.endo='#ff3852'))
}

DefaultCategorySymbols <- function()
{
    return (c(19, 
              producer=21, 
              invertebrate=22, 
              vert.ecto=23, 
              vert.endo=24))
}

DefaultCategoryLabelColours <- function()
{
    # TODO Decide on these colours. They must work with a range of node colours
    # and on a white background.
    return (c('black', 
              producer='black', 
              invertebrate='black', 
              vert.ecto='black', 
              vert.endo='black'))
}

Log10MLabel <- function(community, name='italic(M)', 
                        units=CPS(community)$M.units)
{
    lab <- paste('log[10](', name, ') ~ (', units, ')', sep='')
    return (as.formula(lab))
}

Log10NLabel <- function(community, name='italic(N)', 
                        units=CPS(community)$N.units)
{
    lab <- paste('log[10](', name, ') ~ (', units, ')', sep='')
    return (as.formula(lab))
}

Log10BLabel <- function(community, name='italic(B)', units=with(CPS(community), 
                        paste(M.units, '~', N.units)))
{
    cp <- CPS(community)
    lab <- paste('log[10](', name, ') ~ (', units, ')', sep='')
    return (as.formula(lab))
}

LMabline <- function(model, ...)
{
    # Like abline(model) but restricts the x extent of the line to the x 
    # values that the model was fitted to.

    if(!is.null(model))
    {
        if('lm'!=class(model)) stop('Not an lm object')
        # TODO Error if model is not an ordinary linear regression

        x <- model$model[, 2]
        y <- fitted(model)
        y <- c(y[which(x==min(x))[1]], y[which(x==max(x))[1]])
        lines(range(x), y, ...)
    }
}

PlotLinearModels <- function(models, colour.spec, col, ...)
{
    # A helper that plots the list of models using LMabline(). List elements 
    # can be NULL. Returns a vector of colours used for each model.

    # If col is missing then colour.spec is used to provide colours using the 
    # names of the models list.
    # If col is missing and colour.spec is missing, DefaultCategoryColours()
    # if all of the model names are in the DefaultCategoryColours() names or 
    # are '', 'all' or .UnnamedString().

    if('lm'==class(models))
    {
        models <- list(all==models)
    }

    stopifnot('list'==class(models))
    models <- models[!sapply(models, is.null)]

    model.name <- names(models)

    if(missing(col))
    {
        if(missing(colour.spec))
        {
            # Use the default category colours if no spec provided
            colour.spec <- DefaultCategoryColours()
        }

        # A list of linear models fitted by LinearRegressionsByClass() 
        # will have an 'all' model that has been fitted to all data points.
        # TODO This feels like a bit of a hack - find a better solution.
        if('all' %in% model.name && !'all' %in% names(colour.spec) && 
           '' %in% names(colour.spec))
        {
            model.name['all'==model.name] <- ''
        }

        col <- .GraphParamFromSpec(model.name, colour.spec, 'black')
    }

    mapply(LMabline, model=models, col=col, ...)
    return (col)
}

PlaceMissingPoints <- function(x, xlim, y, ylim)
{
    # Returns a matrix of two columns that contain x and y values respectively. 
    # x and y must be of the same length. 
    # Points with either x or y of NA are placed at bottom of either the range 
    # given by xlim (ylim) or the range of x (y) if xlim (ylim) is NULL
    F <- function(values, lim)
    {
        # This catches two cases - an element in values is 0 and an element 
        # in values is NA. In both cases log10() returns NA.
        valid <- !is.na(values)
        if(any(!valid))
        {
            # Limits come the user-provided limits; if these are NULL, 
            # limits are taken from the extent of the data.
            if(is.null(lim))
            {
                lim <- range(values[valid])
            }

            # Place values that are NA at lower limit of the axis. 
            values[!valid] <- lim[1]
        }
        return (values)
    }

    processed.x <- F(x, xlim)
    processed.y <- F(y, ylim)

    # Jitter values that are NA in both x and y
    jitter <- which(is.na(x) & is.na(y))
    if(length(jitter)>1)
    {
        # Jitter in x
        x.amount <- diff(range(processed.x)) / 50

        # Jitter no more than 1/3 of the way across the x axis
        max.x <- diff(range(processed.x)) / 3

        x.from <- min(processed.x)
        x.to <- x.from + min(length(jitter)*x.amount, max.x)
        processed.x[jitter] <- 
                   seq(x.from, x.to,length.out=length(jitter))

        if(length(jitter)>2)
        {
            # Jitter every other node in y
            y.fraction <- 50
            jitter.y <- jitter[0==(1:length(jitter) %% 2)]
            processed.y[jitter.y] <- processed.y[jitter.y] + 
                                     diff(range(processed.y))/y.fraction
        }
    }

    return (cbind(processed.x, processed.y))
}

.HighlightColours <- function(n, col, highlight, lowlight)
{
    # Returns a vector of colours of length n or NULL 
    # col should be NULL, of length 1 or of length n.
    # highlight and lowlight should contain integer indices. 

    # Nodes to be highlighted are blended towards black.
    # Nodes to be lowlighted are blended towards white. 
    if(length(highlight)>0 || length(lowlight)>0)
    {
        if(is.null(col))         col <- rep(par('col'), n)
        else if(1==length(col))  col <- rep(col, n)
        stopifnot(length(col)==n)

        if(length(highlight)>0)
        {
            # Blend highlight to black
            col[highlight] <- .BlendColour(col[highlight], "black")
        }

        if(length(lowlight)>0)
        {
            # Blend lowlight to white
            col[lowlight] <- .BlendColour(col[lowlight], "white")
        }
    }

    return (col)
}

.HighlightSymbols <- function(n, pch, highlight, lowlight)
{
    if(length(highlight)>0)
    {
        if(is.null(pch))         pch <- rep(par('pch'), n)
        else if(1==length(pch))  pch <- rep(pch, n)
        stopifnot(length(pch)==n)

        # Show highlight as circle in col filled with bg
        pch[highlight] <- 21
    }

    return (pch)
}

.HighlightBG <- function(n, bg, highlight, lowlight)
{
    if(length(highlight)>0 || length(lowlight)>0)
    {
        if(is.null(bg))         bg <- rep(par('bg'), n)
        else if(1==length(bg))  bg <- rep(bg, n)
        stopifnot(length(bg)==n)

        if(length(highlight)>0)
        {
            # TODO Which of these is a better scheme?
            if(TRUE)
            {
                bg[highlight] <- .BlendColour(bg[highlight], "white")
            }
            else
            {
                bg[highlight] <- 'white'
            }
        }

        if(length(lowlight)>0)
        {
            bg[lowlight] <- .BlendColour(bg[lowlight], "white")
        }
    }

    return (bg)
}

.HighlightLabelColours <- function(n, label.col, highlight, lowlight)
{
    # TODO Decide what to do here
    if(FALSE && length(highlight)>0)
    {
        if(is.null(label.col))         label.col <- rep(par('col'), n)
        else if(1==length(label.col))  label.col <- rep(label.col, n)
        stopifnot(length(label.col)==n)

        if(length(highlight)>0)
        {
            label.col[highlight] <- 'black'
        }
    }

    return (label.col)
}

.HighlightCex <- function(n, cex, highlight, lowlight)
{
    # Don't do anything special with cex.
    return (cex)
}


.AddAxisTicks <- function(...)
{
    # A private helper that adds tick marks to the top and right of the plot, 
    # if appropriate. 
    if(getOption('cheddarTopAndRightTicks', TRUE))
    {
        dots <- list(...)

        # Use x %in% y form (rather than equality test) to check for xaxt in 
        # dots, as it handles NULL case
        if('n'!=par('xaxt') && !'n' %in% dots[['xaxt']])
        {
            axis(3, labels=FALSE)
        }

        if('n'!=par('yaxt') && !'n' %in% dots[['yaxt']])
        {
            axis(4, labels=FALSE)
        }
    }
}

.PlotHighlightedPoints <- function(x, y, col, pch, bg, cex, 
                                  labels, label.col, label.cex, 
                                  show.points, show.labels, 
                                  highlight, ...)
{
    # A helper that plots points / text not in 'highlight' and then 
    # points / text in 'highlight'.
    stopifnot(length(x)==length(y))
    stopifnot(is.null(highlight) || all(highlight>0 & highlight<=length(y)))

    F <- function(to.plot)
    {
        # A helper that plots the subset given by 'to.plot'
        # The relevant subset of any params (col, pch etc) is used
        r <- list()
        r$x <- x[to.plot]
        r$y <- y[to.plot]
        if(length(col)==length(x))  r$col <- col[to.plot]
        else                        r$col <- col
        if(length(pch)==length(x))  r$pch <- pch[to.plot]
        else                        r$pch <- pch
        if(length(bg)==length(x))   r$bg <-  bg[to.plot]
        else                        r$bg <- bg
        if(length(cex)==length(x))  r$cex <- cex[to.plot]
        else                        r$cex <- cex

        if(length(labels)==length(x))  r$labels <- labels[to.plot]
        else                           r$labels <- labels
        if(length(label.col)==length(x))  r$label.col <- label.col[to.plot]
        else                              r$label.col <- label.col
        if(length(label.cex)==length(x))  r$label.cex <- label.cex[to.plot]
        else                              r$label.cex <- label.cex

        if(show.points)
        {
           points(r$x, r$y, pch=r$pch, col=r$col, bg=r$bg, cex=r$cex, ...)
        }

        if(show.labels)
        {
           text(r$x, r$y, r$labels, col=r$label.col, cex=r$label.cex, ...)
        }
    }

    # Plot non-highlighted points first
    first <- setdiff(1:length(x), highlight)
    if(length(first))
    {
        F(first)
    }

    # Then plot highlighted points first
    if(length(highlight))
    {
        F(setdiff(highlight, first))
    }

    if(FALSE)
    {
    # An alternative
    # Reorder parameters such that highlighted points are plotted last
    if(length(highlight)>0)
    {
        x <- c(x[-highlight], x[highlight])
        y <- c(y[-highlight], y[highlight])
        if(length(col)==length(x))  col <- c(col[-highlight], col[highlight])
        if(length(pch)==length(x))  pch <- c(pch[-highlight], pch[highlight])
        if(length(bg)==length(x))    bg <- c(bg[-highlight], bg[highlight])
        if(length(cex)==length(x))  cex <- c(cex[-highlight], cex[highlight])

        if(length(labels)==length(x))  labels <- c(labels[-highlight], 
                                                   labels[highlight])
        if(length(label.col)==length(x))  label.col <- c(label.col[-highlight], 
                                                         label.col[highlight])
        if(length(label.cex)==length(x))  label.cex <- c(label.cex[-highlight], 
                                                         label.cex[highlight])
        }

    if(show.points)
    {
       points(x, y, pch=pch, col=col, bg=bg, cex=cex, ...)
    }

    if(show.labels)
    {
       text(x, y, labels, col=label.col, cex=label.cex, ...)
    }
    }
}

.PlotNetwork <- function(x, y, network, col, lty, lwd, highlight, ...)
{
    # A helper that plots lines connecting network using points in x and y. 
    # x and y should be named vectors containing the locations on nodes on the 
    # graph. network should be a matrix of two columns that are indices 
    # of x and y respectively. highlight should be a vector of row indices 
    # of network.

    F <- function(to.plot)
    {
        # A helper that plots the subset given by 'to.plot'
        # The relevant subset of any params (col, pch etc) is used
        r <- list()

        if(length(col)==nrow(network))    r$col <- col[to.plot]
        else                              r$col <- col
        if(length(lty)==nrow(network))    r$lty <- lty[to.plot]
        else                              r$lty <- lty
        if(length(lwd)==nrow(network))    r$lwd <- lwd[to.plot]
        else                              r$lwd <- lwd

        segments(x0=x[network[to.plot,1]], 
                 y0=y[network[to.plot,1]], 
                 x1=x[network[to.plot,2]],
                 y1=y[network[to.plot,2]],
                 col=r$col,
                 lty=r$lty,
                 lwd=r$lwd)
    }

    # Plot non-highlighted points first
    first <- setdiff(1:nrow(network), highlight)
    if(length(first))
    {
        F(first)
    }

    # Then plot highlighted points first
    if(length(highlight))
    {
        F(setdiff(highlight, first))
    }
}

.ResolveLevel <- function(community, level)
{
    # level should be either a function, the name of a property or a vector of 
    # levels.
    # A vector of levels of length(NumberOfNodes) is returned. 
    if(is.function(level))
    {
        level <- level(community)
        if(length(level)!=NumberOfNodes(community))
        {
            stop(paste('The function given for level returned', 
                       length(level), 'values but should return', 
                       'NumberOfNodes values'))
        }
    }
    else if(1==length(level))
    {
        # This breaks if the community has only 1 node
        level <- NPS(community, level)[,level]
    }

    if(length(level)!=NumberOfNodes(community))
    {
        stop(paste('level contains', length(level), 'values but', 
                   'should be of length NumberOfNodes'))
    }

    if(!all(level>0 & level<Inf & !is.na(level)))
    {
        stop('All values in level should be >0 and <Inf and not NA')
    }

    return (level)
}

.ResolveTrophicLinkProperty <- function(community, property)
{
    sub <- substr(property, '10', nchar(property))
    if(.StringStartsWith(property, 'resource.') || 
       .StringStartsWith(property, 'consumer.'))
    {
        return (TLPS(community, node.properties=sub)[,property])
    }
    else
    {
        return (TLPS(community, link.properties=property)[,property])
    }
}

.ResolveTrophicLinkAxisLabel <- function(property)
{
    sub <- substr(property, '10', nchar(property))
    if(.StringStartsWith(property, 'resource.') || 
       .StringStartsWith(property, 'consumer.'))
    {
        return(as.expression(substitute(property[end], 
                                        list(property=sub, 
                                             end=substr(property, 1, 8)))))
    }
    else
    {
        return (property)
    }
}

