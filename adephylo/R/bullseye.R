##
## PLOT A FAN TREE, WITH BULLSEYE LEGEND AND AXIS, AND OPTIONAL COLORS
## FOR TIPS
##
## Author: Thibaut Jombart, May 2013.
## t.jombart@imperial.ac.uk
##

############
## bullseye
############
bullseye <- function(phy, traits=NULL, col.tips.by=NULL, col.pal=spectral,
                     circ.n=6, circ.bg=transp("royalblue",.1), circ.unit=NULL,
                     legend=TRUE, leg.posi="bottomleft", leg.title="", leg.bg="white",
                     traits.inset=1.1, traits.space=0.05, traits.pch=19, traits.cex=1,
                     alpha=1, axis=TRUE, ...){
    ## CHECKS ##
    if(inherits(phy, c("phylo4","phylo4d"))) phy <- as(phy, "phylo")
    if(!is.list(col.pal)) col.pal <- c(col.pal)
    leg.info <- NULL

    ## REORDER DATA BY TIP LABEL ##
    ## make sure traits is a data.frame
    if(!is.null(traits)) traits <- as.data.frame(traits)
    if(!is.null(traits) && !is.null(row.names(traits))){
        if(!all(phy$tip.label %in% row.names(traits))){
            warning("tip labels and names of the traits matrix do not match")
        } else {
            traits <- traits[phy$tip.label,,drop=FALSE]
        }
    }

    ## col.tips.by
    if(!is.null(col.tips.by) && is.data.frame(col.tips.by)){
        old.names <- row.names(col.tips.by)
        col.tips.by <- unlist(col.tips.by)
        names(col.tips.by) <- old.names
    }
    if(!is.null(col.tips.by) && !is.null(names(col.tips.by))){
        col.tips.by <- col.tips.by[phy$tip.label]
    }

    ## recycle col.pal
    pal.length <- 0
    if(!is.null(traits)) pal.length <- pal.length + ncol(traits)
    if(!is.null(col.tips.by)) pal.length <- pal.length + 1
    col.pal <- rep(col.pal, length=pal.length)


    ## PLOT THE PHYLOGENY
    ## window setting
    oxpd <- par("xpd")
    par(xpd=TRUE)
    on.exit(par(oxpd))

    ## handle color info
    if(!is.null(col.tips.by)){
        tip.col.info <- any2col(col.tips.by, col.pal=col.pal[[1]])
        plot(phy, type="fan", tip.col=transp(tip.col.info$col,alpha), ...)
    } else{
        plot(phy, type="fan", ...)
    }


    ## HANDLE THE 'BULLSEYE' ##
    ## annot info
    if(is.null(circ.unit)){
        annot.max <- 0.5*diff(par("usr")[1:2])
        annot.dist <- seq(from=0, to=annot.max, length=circ.n)
    } else {
        annot.dist <- seq(from=0, by=circ.unit, length=circ.n)
        annot.max <- max(annot.dist)
    }

    ## trace the disks
    symbols(rep(0,circ.n), rep(0,circ.n), circles=annot.dist, inches=FALSE,
        bg=circ.bg, fg=NA, add=TRUE)

    ## axis annotation
    if(axis){
        segments(-annot.dist[2],0,-annot.dist[3],0)
        text(-mean(annot.dist[2:3]),-annot.dist[2]/5,
             label=format(annot.dist[2], scientific=TRUE, digits=3),cex=.7)
    }


    ## PLOT TRAITS ##
    if(!is.null(traits)){
        ## recycle pch and cex
        traits.pch <- rep(traits.pch, length=ncol(traits))
        traits.cex <- rep(traits.cex, length=ncol(traits))

        ## get tips coordinates
        tips.x <- get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[1:length(phy$tip.label)]
        tips.y <- get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[1:length(phy$tip.label)]

        ## use furthest tip from the root to define new base coords
        vec.length <- sqrt(tips.x^2 + tips.y^2)

        x.base <- (tips.x/vec.length) * max(vec.length) * traits.inset
        y.base <- (tips.y/vec.length) * max(vec.length) * traits.inset

        ## plot traits
        for(i in 1:ncol(traits)){
            col.info <- any2col(traits[,i], col.pal=col.pal[[i]])
            temp.x <- x.base * (traits.inset + i*traits.space)
            temp.y <- y.base * (traits.inset + i*traits.space)
            points(temp.x, temp.y, pch=traits.pch[i], col=transp(col.info$col,alpha), cex=traits.cex[i])

            ## save info for legend if needed
            if(is.null(col.tips.by) && i==1){
                leg.info <- list(col=transp(col.info$leg.col,alpha), txt=col.info$leg.txt)
            }
        }
    }


    ## ADD LEGEND ##
    ## legend info
    if(!is.null(legend)){
        ## legend for tip colors
        if(!is.null(col.tips.by)){
            leg.col <- transp(tip.col.info$leg.col,alpha)
            leg.txt <- tip.col.info$leg.txt
            leg.info <- list(col=transp(tip.col.info$leg.col,alpha), txt=tip.col.info$leg.txt)
        }

        ## plot legend
        if(!is.null(leg.info) && legend){
            leg.info$posi <- leg.posi
            legend(x=leg.info$posi, legend=leg.info$txt, fill=leg.info$col, title=leg.title, bg=leg.bg)
            return(invisible(leg.info))
        }
    }

    return(invisible())
} # end bullseye
