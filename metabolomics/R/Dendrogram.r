Dendrogram <- function(inputdata, distmethod="manhattan", aggmethod="ward", 
    main="Dendrogram", cex=1, ...)
{
    # Get groups information
    groups <- factor(inputdata[, 1], levels=unique(inputdata[, 1]))

    # Remove groups for data processing
    hcadata <- editcolnames(inputdata[, -1])

    # Prepare distance matrix for samples
    dist_sample <- dist(hcadata, method=distmethod, ...)
    hca_sample <- hclust(dist_sample, method=aggmethod, ...)

    smpl_labels<-paste(groups, ':', rownames(hcadata), sep='')
    
    par_defs<-par(mai=par()$mai,
        mar=par()$mar,
        mex=par()$mex,
        mgp=par()$mgp,
        oma=par()$oma,
        omd=par()$omd,
        omi=par()$omi,
        plt=par()$plt,
        usr=par()$usr,
        xaxp=par()$xaxp,
        xpd=par()$xpd,
        yaxp=par()$yaxp
    )
    on.exit(par(par_defs))
    
    lmat<-matrix(c(4,0,0,2,1,3), ncol=2)
    ### These need to be determined properly
    lwid<-c(0.6,4.25)
    lhei<-c(4.25, 0.25, 1)
    col_list<-ColList(length(levels(groups)))
    cols_used<-col_list[groups[hca_sample$order]]
    
    # Show where things are supposed to go
    #dev.new()
    #plotmat<-layout(lmat, lwid, lhei)
    #layout.show(plotmat)
    
    # Do plots
    #dev.new()
    layout(lmat, lwid, lhei)
    # 1) group colours
    # Use rect() for correct placement of group colours (image() width wrong)
    xpos<-seq(0, 1, length.out=length(smpl_labels))
    xwid<-xpos[2]/2
    par(mar=c(0,0,0,1))
    plot(matrix(c(0,0,1,1),nrow=2,byrow=TRUE), type="n", axes=FALSE, 
        xlab="", ylab=""
    )
    for (ii in 1:length(xpos)) {
        rect(xpos-xwid, 0, xpos+xwid, 1, density=NA, col=cols_used)
    }
    
    # 2) Dendro
    par(mar=c(0,0,0,1))
    plot(hca_sample, 
        labels=rep("", length(groups)),        # clear labels here
        hang=(-1),                             # even-length ends
        main=NULL,                             # no plot title 
                                               ########### fix (oma? mtext?)
        axes=FALSE,                            # suppress axes (see below)
        xlab="",                               # suppress all axis titles
        ylab="", 
        sub="",
        ...
    )
    
    axmap<-par()$usr
    dendyaxp<-par()$yaxp
    # 3) axis labels (x-axis)
    par(cex=cex, mar=c(0,0,0,1), usr=axmap)
    axis(1, 1:length(smpl_labels), labels=smpl_labels[hca_sample$order], 
        las=2, tick=0, 
        line=1.5 # if this is not displaced, it is placed _over_ the group colours
    )
    
    # 4) height bar (y-axis)
    par(cex=cex, mar=c(0,0,0,0), usr=axmap, yaxp=dendyaxp)
    ## interestingly, for very large values of height (i.e. 0 to 1.2e+09; from
    ## raw data), this hits a memory limit. 
    axis(2, pretty(dendyaxp[1]:dendyaxp[2], n=dendyaxp[3]+1), 
        line=0.5, las=2
    )
    
    par(par_defs)
}

