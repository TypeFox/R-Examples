#
#
#
idendro<-structure(function# Interactive Dendrogram
###
### 'idendro' is a plot enabling users to visualize a dendrogram and
### inspect it interactively: to select and color clusters anywhere in
### the dendrogram, to zoom and pan the dendrogram, and to visualize
### the clustered data not only in a built-in heat map, but also in any
### interactive plot implemented in GGobi (as available using the
### 'rggobi' package). The integration with GGobi (enabled using the
### 'ggobi' argument), but also with the user's code is implemented in
### terms of two callbacks (see the 'colorizeCallback' and
### 'fetchSelectedCallback' arguments).
### 'idendro' can be used to inspect quite large dendrograms (tens
### of thousands of observations, at least).
###
### The 'idendr0' package is a lightweight backport of the 'idendro'
### package. While the 'idendro' package depends on libraries not
### easily available on some platforms (e.g. Windows), the 'idendr0'
### package is based on platform-independent Tcl/Tk graphic widget
### toolkit, and thus made widely available. However, the 'idendro'
### package should be preferred, if available, for its better
### interactivity and performance.
###
##details<<
## 'idendro' displays an interactive dendrogram enriched, optionally,
## with a heat map and/or a brushed map.
##
## The dendrogram represents the result of a hierarchical cluster
## analysis performed on a set of observations (see e.g. 'hclust').
## There is an axis drawn below the dendrogram displaying the "height"
## of the clusters in the dendrogram.
##
## The heat map visualizes the observations living in k-dimensional
## feature space by mapping their features onto a color scale and
## displaying them as rows of 'k' colored rectangles. By default,
## normalization (scaling) of individual features to a common visual
## scale is enabled. Scaling of observations is also supported (see the
## 'doScaleHeatmapByRows' argument).
##
## The brushed map can indicate which observations are currently
## selected in some external plot/tool 'idendro' is integrated
## with (e.g. a GGobi scatter plot matrix). Technically speaking,
## the current selection must be determined explicitly by clicking the
## "fetCh selected" button (or pressing the 'Alt+C' shortcut), which
## results in calling the 'fetchSelectedCallback' function (see
## arguments).
##
## The dendrogram can be zoomed and panned. To zoom in a
## specific region, right click and drag in the dendrogram.
## Mouse wheel can also be used to zoom in and out. To pan a zoomed
## dendrogram, middle click and drag the mouse. Zooming and panning
## history is available (see 'GUI').
##
## User can select clusters manually one by one (by clicking
## at individual clusters in the dendrogram), or automatically by
## "cutting" the dendrogram at a specified height. To cut the
## dendrogram, navigate the mouse close to the dendrogram axis
## (a dashed line will appear across the dendrogram at a specified
## height), and left click. Clusters just beneath the cutting
## height will get selected, replacing the clusters currently
## selected. Selection history is available (see 'GUI').
##
##   \emph{Graphic User interface (GUI):}
##
## In the left part of the dendrogram window, there is a simple GUI.
## In the top part of the GUI come cluster-specific controls and info
## panels arranged in rows. (The number of rows is determined by the
## 'maxClusterCount' argument.)
## In each row, there is the current cluster selector (a radio button
## decorated with a cluster ID and a color code (determined by the
## 'clusterColors' argument)), and cluster-specific statistics: the
## total number (and the ratio) of the observations in that specific
## cluster out of the total number of observations, and the number
## (and the ratio) of the observations in that cluster out of the
## observations brushed.
## The current cluster determines which color and ID will be
## associated with a cluster selected in the dendrogram,
## At any time, exactly one cluster is selected as the current
## cluster. 
##
## At the bottom of the GUI window, there are buttons controling
## zooming, cluster selection, and heat map smoothing:
##
## "Undo zoom" - retrieves the previous zoom region from history
##
## "Full view" - zooms the dendrogram out maximally
##
## "Undo selection" - retrieves the previous cluster selection from
##     history
##
## "Unselect" - unselects the current cluster in the dendrogram
##
## "Unselect all" - unselects all clusters
##
## The "heat map smoothing" mode can be set to one of:
##
##   "none" - the heat map gets never smoothed, it displays the
##      features of all the individual observations
##
##   "cluster" - the heat map displays the average features for the 
##       currently selected clusters
##
##   "zoom" - the heat map displays the average feature for each
##      elementary (i.e. the finest) cluster seen in the dendrogram
##      currently. When the dendrogram is zoomed out maximally,
##      the features of all the elementary clusters (i.e. the
##      individual observations) are displayed. When the user zooms in
##      the dendrogram, such that some clusters get hidden, the
##      features of the observations forming the hidden clusters get
##      averaged.
##
##  "Quit"
##
##seealso<< idendro::idendro, stats::hclust, stats::plot.hclust
##
(   h, ##<< object of class 'stats::hclust' (or other class
    ## convertible to class 'hclust' by the 'as.hclust' function)
    ## describing a hierarchical clustering.
    ## If _inversions_ in heights (see 'hclust') is detected,
    ## the heights get fixed in a simple naive way by preserving
    ## non-negative relative differences in the heights, but changing
    ## negative differences to zero. Using clustering with monotone
    ## distance measure should be considered in that case.

    x=NULL, ##<< data frame holding observations tha were clustered
    ## giving rise to 'h'.
    ## The heat map will depict this data. (The heat map can be scaled
    ## - see the 'doScaleHeatmap' and 'doScaleHeatmapByRows' arguments.)
    ## Non-numeric types will get converted to numeric using 'as.numeric'.
    ## This parameter is optional.

    qx=NULL, ##<< (unused, appears for compatibility with
    ## idendro::idendro).

    clusters=NULL, ##<< the assignment of observations to clusters to start
    ## with, typically the value of a previous call to 'idendro'.
    ## A numeric vector of length of the number of observations is
    ## expected, in which 0s denote unselected observations, and values
    ## of i > 0 mark members of the cluster `i'.

    hscale=1.5, ##<< horizontal scaling factor of the dendrogram
    ## figure. As the dendrogram is implemented as a Tcl/Tk image, and
    ## rtcltk does not support image resizing (e.g. on window
    ## maximization), the dendrogram keeps its original size regardless
    ## of the size of its enclosing window. Thus specifying the
    ## hscale of more than 100% is preferred to make the dendrogram
    ## large enough.

    vscale=1.5, ##<< vertical scaling factor of the dendrogram
    ## figure. See 'hscale'.

    silent=FALSE, ##<< if TRUE, no informative messages will be shown

    zoomFactor=1/240, ##<<the amount of zoom in/out as controlled by the
    ## mouse wheel

    observationAnnotationEnabled=TRUE, ##<< shall the names of individual
    ## observations (rownames of 'x') be shown next to the
    ## dendrogram/heat map?

    clusterColors=c('red','green','blue','yellow','magenta','cyan','darkred','darkgreen','purple','darkcyan'),##<< colors
    ## of individual clusters

    unselectedClusterColor='black',##<< the color of unselected dendrogram
    ## branches

    maxClusterCount=max(length(clusterColors),ifelse(!is.null(clusters),max(clusters),0)), ##<< the
    ## maximum number of clusters user can select. If greater than the number
    ## of 'clusterColors', cluster colors will get recycled.
    ## This parameter affects the size of the GUI and the number of
    ## clusters that can be selected automatically by "cutting" the
    ## dendrogram.

    heatmapEnabled=TRUE, ##<< shall the heat map be drawn?

    heatmapSmoothing=c('none','cluster','zoom'),##<< heat map smoothing mode,
    ## one of
    ## 'none' - the heat map gets never smoothed, it displays the
    ##      features of all the individual observations
    ## 'cluster' - the heat map depicts the average features
    ##      for the currently selected clusters,
    ## 'zoom' - the heat map displays the average feature for each
    ##      elementary (i.e. the finest) cluster seen in the
    ##      dendrogram currently.

    heatmapColors=colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))(10), ##<< heat map
    ## color palette represented by a list of colors, e.g.
    ## a sequential palette generated by `brewer.pal', or
    ## `colorRampPalette(.)(.)', `gray.colors(.)', or `hsv(.)'.

    doScaleHeatmap=TRUE, ##<< scale each heat map column to the <0,1> range?
    ## (The default is TRUE.)

    doScaleHeatmapByRows=FALSE, ##<< scale heat map rows, not columns
    ## (The default is FALSE.)

    heatmapRelSize=.2, ##<< relative size of the heat map - the ratio
    ## of the heat map width to the width of the dendrogram, the heat
    ## map, and the brushed map. The default is 20%.

    colorizeCallback=NULL, ##<< callback function called when cluster
    ## selection changes; the argument is a vector assigning color
    ## indices (0=no color, >0 colors) to individual observations.

    fetchSelectedCallback=NULL, ##<< callback function used to fetch
    ## observation selection made externally. The callback must return
    ## a boolean vector of length of the number of observations in `x'.
    ## i-th element in the vector specifies whether given observation
    ## is selected.

    brushedmapEnabled=!is.null(fetchSelectedCallback), ##<< shall
    ## brushed map be drawn? If TRUE, a column vector is drawn next to
    ## dendrogram (and heatmap, if there is one) depicting observation
    ## that were fetched by fetchSelectedCallback. The color of the
    ## observations is the color of the cluster used to fetch
    ## observations into.

    brushedmapRelSize=ifelse(!is.null(x),heatmapRelSize/ncol(x),.05), ##<< relative
    ## size of the brushed map - the ratio of the brushed map width to
    ## the width of the dendrogram, the heat map, and the brushed map.
    ## The default is the size of a single column in the heat map, or
    ## 5% if there is no heatmap.

    geometry=NULL, ##<< window geometry ("width x height + xoffset +
    ## yoffset"). Almost useless as the dendrogram does not resize, see
    ## the 'hscale' and 'vscale' arguments instead.

    ggobi=FALSE, ##<< plot feature space projections of `x' in ggobi
    ## and bidirectionally integrate with the plot? (defaults to FALSE
    ## as some users may not have ggobi available)

    ggobiGlyphType=1, ##<< ggobi glyph type used to draw observations in
    ## ggobi (defaults to a single pixel; see rggobi::glyph_type)

    ggobiGlyphSize=1, ##<< size of ggobi glyphs (see rggobi::glyph_size)

    ggobiFetchingStyle='selected', ##<< how should we recognize
    ## ggobi-selected observations to be fetched to idendro?
    ## Use 'selected' to fetch observations selected by ggobi brush,
    ## or glyph type number 2-6 to fetch observations selected by ggobi
    ## persistent brushing with a specific glyph type.

    ggobiColorScheme='Paired 12', ##<< GGobi color scheme used to
    ## color observations in ggobi plots according to the clusters
    ## selected in the dendrogram

    dbg=0, ##<< debug level (0=none, 1=brief, 2=verbose)

    ... ##<< additional graphical parameters to be passed to the
    ## dendrogram plot.
    ) {

    #### debugs
    ####
    # numeric (0=none, 1=brief, 2=verbose)
    dbg.args<-1*dbg
    dbg.dendro.zoom<-1*dbg
    dbg.pan<-dbg*1
    dbg.heatmap<-1*dbg
    dbg.heatmap.smooth<-1*dbg
    dbg.tx<-0*dbg
    dbg.dendro<-1*dbg
    dbg.dendro.cut<-1*dbg
    dbg.dendro.select<-1*dbg
    dbg.clustersArg<-1*dbg
    dbg.colorizeCb<-1*dbg
    dbg.fetched<-1*dbg
    dbg.mouse<-0*dbg
    dbg.geometry<-0*dbg
    dbg.ggobi<-1*dbg

    #### arguments handling
    ####
    if (ggobi) {
        if (ggobiFetchingStyle!='selected') {
            selectionType<-suppressWarnings(as.numeric(ggobiFetchingStyle))
            if (dbg.ggobi) printVar(selectionType)
            if (is.na(selectionType) || selectionType!=round(selectionType) || selectionType<2 || selectionType>6) {
                stop(paste('invalid ggobi fetching type \'',ggobiFetchingStyle,'\', use one \'selected\' or glyph type 2 to 6',sep=''))
            }
        }
        ggobiColorizeCallback<-function(colors) {
            if (dbg.ggobi>1) cat('ggobi colorizeCallback called\n')
            if (dbg.ggobi>1) printVar(colors)
            rggobi::glyph_color(g[1])<-(1+colors)
            if (!is.null(colorizeCallback)) colorizeCallback(colors)
        }
        ggobiFetchSelectedCallback<-function() {
            if (dbg.ggobi>1) cat('ggobi fetchSelectedCallback called\n')
            if (dbg.ggobi>1) cat(paste('fetching type: ',ggobiFetchingStyle,'\n'))
            if (ggobiFetchingStyle=='selected') {
               selection<-rggobi::selected(g[1])
            } else {
                selection<-(rggobi::glyph_type(g[1])==ggobiFetchingStyle)
            }
            if (dbg.ggobi>1) printVar(selection)
            return(selection)
        }

        if (!requireNamespace("rggobi",quietly=TRUE)) {
            stop('The \'rggobi\' package is not installed, can\'t integrate with GGobi.')
        }
        # We need to attach "rggobi", otherwise thw following error gets thrown from ggobi:
        #  Error in .RGtkCall("R_setGObjectProps", obj, value, PACKAGE = "RGtk2") :
        #    Invalid property 1!
        # However, in order to keep the search path intact after the call to idendro(),
        # we attempt to unload "rggobi" later (see below).
        if (!"package:rggobi"%in%search()) {
          attachNamespace("rggobi")
          ggobiAttached<-TRUE
        } else {
          ggobiAttached<-FALSE
        }

        if (!silent) {
            message('Note: integrating with GGobi, ignoring the \'clusterColors\' argument: using colors from the \'',
                ggobiColorScheme,'\' GGobi color scheme specified using the \'ggobiColorScheme\' argument.')
        }
        g<-rggobi::ggobi(x)
        # set color scheme
        rggobi::colorscheme(g)<-ggobiColorScheme
        # read the ggobi color scheme colors in order to use the same colors in the dendrogram
        cols<-sapply(rggobi::colorscheme(g)$colors,function(x)rgb(x[1],x[2],x[3]))
        if (dbg.ggobi) printVar(cols)
        # close the default display(s)
        sapply(rggobi::displays(g),close)
        # set the "pixel" draw style
        rggobi::glyph_type(g[1])<-ggobiGlyphType
        rggobi::glyph_size(g[1])<-ggobiGlyphSize
        # draw scatterplot matrix of all parameters
        rggobi::display(g[1],"Scatterplot Matrix")

        rv<-idendro(
            h=h,
            x=x,
            qx=qx,
            clusters=clusters,
            hscale=hscale,
            vscale=vscale,
            silent=silent,
            zoomFactor=zoomFactor,
            observationAnnotationEnabled=observationAnnotationEnabled,
            clusterColors=cols[-1],
            unselectedClusterColor=unselectedClusterColor,
            maxClusterCount=maxClusterCount,
            heatmapEnabled=heatmapEnabled,
            heatmapSmoothing=heatmapSmoothing,
            heatmapColors=heatmapColors,
            doScaleHeatmap=doScaleHeatmap,
            doScaleHeatmapByRows=doScaleHeatmapByRows,
            heatmapRelSize=heatmapRelSize,
            colorizeCallback=ggobiColorizeCallback,
            fetchSelectedCallback=ggobiFetchSelectedCallback,
            # brushedmapEnabled not set - will be set accoriding to fetchSelectedCallback
            brushedmapRelSize=brushedmapRelSize,
            geometry=geometry,
            ggobi=FALSE, # disabled on the recursive call
            ggobiGlyphType=ggobiGlyphType,
            ggobiFetchingStyle=ggobiFetchingStyle,
            ggobiColorScheme=ggobiColorScheme,
            dbg=dbg)

        # wait until the user closes the dendrogram,
        # and close ggobi as well
        close(g)

        # unload "rggobi" not to alter the search path
        if (ggobiAttached) {
          unloadNamespace("rggobi")
        }

        return(invisible(rv))
    }

    if (dbg.args) cat('--- user supplied arguments: ---\n')
    if (dbg.args) printVar(!is.null(x))
    if (dbg.args) printVar(heatmapEnabled)
    if (dbg.args) printVar(doScaleHeatmap)
    #if (dbg.args) printVar(brushedmapEnabled)
    #if (dbg.args) printVar(observationAnnotationEnabled)

    if (!inherits(h,'hclust')) {
        h<-as.hclust(h)
    }

    if (inherits(x,'dist')) {
        stop('\'x\' argument of invalid type of class \'dist\'')
    }

    if (heatmapEnabled && is.null(x) || heatmapRelSize<=0) {
        # can't draw heat map if we have no data
        heatmapEnabled<-FALSE
        heatmapRelSize<-0
    }

    if (brushedmapRelSize<=0 || !brushedmapEnabled) {
        brushedmapEnabled<-FALSE
        brushedmapRelSize<-0
    }

    if (heatmapRelSize+brushedmapRelSize>.9) {
        stop(sprintf('heatmapRelSize %.2f + brushedmapRelSize %.2f too large',heatmapRelSize,brushedmapRelSize))
    }

    if (is.unsorted(h$height)) {
        if (!silent) {
            message('Note: non-monotone distance detected, applying a simple workaround. Consider using clustering with monotone distance.')
        }
        # 1  4  2  7  6  5  8  9  # h$height
        #    3 -2  5 -1 -1  3  1  # tmp<-diff(h$height),  min(tmp[tmp>0]) = 1
        #    2 -3  4 -2 -2  2  0  # tmp2<-tmp-min(tmp[tmp>0]
        #    0 -3  0 -2 -2  0  0  # tmp2*I(tmp<0)
        #    0 -3 -3 -5 -7 -7 -7  # cumsum(tmp2*I(tmp<0))
        #    0  3  3  5  7  7  7  # -cumsum(tmp2*I(tmp<0))
        # 1  4  2  7  6  5  8  9  # h$height
        # 1  4  5 10 11 12 15 16  # h$height + c(0,-cumsum(tmp*I(tmp<0)))
        tmp<-diff(h$height)
        tmp2<-tmp-min(tmp[tmp>0])
        h$height<-h$height+c(0,-cumsum(tmp2*I(tmp<0)))
    }

    n<-length(h$height)+1
    if (!is.null(x) && nrow(x)!=n) {
        stop(paste('Clustering (of ',n,' objects) does not fit data (',nrow(x),' rows).',sep=''))
    }

    if (is.function(heatmapColors)) {
        stop('\'heatmapColors\' argument of type \'function\' has been deprecated, please supply a color list instead')
    }

    # convert non-numeric data to numeric, if necessary
    if (heatmapEnabled) {
        xOrig<-x
        nonNumericColumnFound<-FALSE
        for (i in 1:ncol(x)) {
            if (!is.numeric(x[,i])) {
                x[,i]<-as.numeric(x[,i])
                nonNumericColumnFound<-TRUE
            }
        }
        if (nonNumericColumnFound) {
            if (!silent) {
                message('Note: non-numeric data found, converting to numeric (in order to enable heatmap drawing).')
            }
        }

        if (is.data.frame(x)) {
            tmp<-names(x)
            x<-as.matrix(x)
            names(x)<-tmp
        }
    } else {
        xOrig<-NULL
    }

    # scale heat map
    if (heatmapEnabled && doScaleHeatmap) {
        scaleVector<-function(x) {
            mn<-min(x,na.rm=TRUE)
            mx<-max(x,na.rm=TRUE)
            if (mn!=mx) {
                x<-(x-mn)/(mx-mn)
            } else {
                x<-x-mn+.5
            }
            x
        }
        if (doScaleHeatmapByRows) {
            for (i in 1:nrow(x)) {
                x[i,]<-scaleVector(x[i,])
            }
        } else {
            for (i in 1:ncol(x)) {
                x[,i]<-scaleVector(x[,i])
            }
        }
    }

    adaptHeatmapForDrawing<-function(x) {
        if (is.null(x)) {
            return(NULL)
        } else {
            return(t(x[,ncol(x):1]))
        }
    }

    df<-prepareDendro(h,x,xOrig,FALSE,dbg.dendro)
    df$heatmap<-adaptHeatmapForDrawing(df$xOrdered)
    df$fetchedBranches<-df$noBranches<-list(indices=c(),branches=data.frame(x1s=c(),x2s=c(),y1s=c(),y2s=c()))
    df$fetchedLeafCount<-0
    fetchedInfo<-NULL
    df$fetchedMap<-df$emptyFetchedMap<-matrix(rep(0,n),nrow=1)
    # initialize clusters from user-supplied clusters argument, if any
    if (!is.null(clusters)) {
        df<-createClustersFromLeafColors(df,clusters,maxClusterCount,dbg.dendro)
    }
    df$leafColorIdxs<-computeLeafColorIdxs(df)

    if (dbg.dendro>1) printVar(df)

    # observation annotations
    if (!is.null(x) && !is.null(rownames(x))) {
        df$observationLabelsOrdered<-rownames(x)
    } else {
        df$observationLabelsOrdered<-h$labels
    }
    if (!is.null(df$observationLabelsOrdered)) {
        df$observationLabelsOrdered<-df$observationLabelsOrdered[df$leafOrder]
    }
    if (observationAnnotationEnabled && is.null(df$observationLabelsOrdered)) {
        # if annotation not available, nothing to draw
        observationAnnotationEnabled<-FALSE
    }

    # dim annotations
    if (!is.null(x) && !is.null(colnames(x))) {
        df$dimLabels<-colnames(x)
    } else {
        df$dimLabels<-NULL
    }

    df$dendroZoomHistory<-list()
    lastDendroZoomHistorySaver<-'none'

    if (dbg.args) cat('--- consolidated arguments: ---\n')
    if (dbg.args) printVar(!is.null(x))
    if (dbg.args) printVar(heatmapEnabled)
    if (dbg.args) printVar(brushedmapEnabled)
    if (dbg.args) printVar(observationAnnotationEnabled)

    heatmapSmoothing<-match.arg(heatmapSmoothing)

    # Determine the color for cluster of given ID (starting at 1).
    clusterColor<-function(id) {
         clusterColors[((id-1)%%length(clusterColors))+1]
    }

    # Callbacks invoked when clusters change.
    updateBrushedClusterInfosOnChange<-function(df) {
        if (brushedmapEnabled) {
            # recompute number of leafs in each cluster
            if (length(df$clusters)==0 || length(df$brushed)==0) {
                txt<-printClusterInfo(0)
                for (i in 1:maxClusterCount) {
                    if (dbg) cat(sprintf('cluster %i: 0 brushed leafs\n',i))
                    eval(parse(text=paste('tclvalue(clusterBrushedInfo',as.character(i),')<-',txt,sep='')))
                }
            } else {
                for (i in 1:maxClusterCount) {
                    if (i>length(df$clusters) || is.null(df$clusters[[i]])) {
                        cnt<-0
                    } else {
                        cnt<-sum(df$leafColorIdxs==i & df$brushed)
                    }
                    if (dbg) cat(sprintf('cluster %i: %s brushed leafs\n',i,cnt))
                    txt<-printClusterInfo(cnt)
                    eval(parse(text=paste('tclvalue(clusterBrushedInfo',as.character(i),')<-',txt,sep='')))
                }
            }
        }
    }
    updateClusterInfosOnChange<-function(df) {
        # recompute number of leafs in each cluster
        if (length(df$clusters)==0) {
            txt<-printClusterInfo(0)
            for (i in 1:maxClusterCount) {
                if (dbg) cat(sprintf('cluster %i: 0 leafs\n',i))
                eval(parse(text=paste('tclvalue(clusterTotalInfo',as.character(i),')<-',txt,sep='')))
            }
        } else {
            for (i in 1:maxClusterCount) {
                if (i>length(df$clusters) || is.null(df$clusters[[i]])) {
                    cnt<-0
                } else if (length(df$clusters[[i]]$indices)==0) {
                    cnt<-0
                } else {
                    # note the '+1' as the most bottom observation does not
                    # form a (sub)cluster, so it is listed in the indices,
                    # but it must be counted
                    cnt<-length(df$clusters[[i]]$indices)+1
                }
                if (dbg) cat(sprintf('cluster %i: %s leafs\n',i,cnt))
                txt<-printClusterInfo(cnt)
                eval(parse(text=paste('tclvalue(clusterTotalInfo',as.character(i),')<-',txt,sep='')))
            }
        }
        updateBrushedClusterInfosOnChange(df)
    }

    # Callbacks invoked when clusters change.
    updateClustersOnChangePreDraw<-function(df) {
        updateClusterInfosOnChange(df)
        
        df$fetchedBranches<-df$noBranches

        if (tclvalue(heatmapSmoothingMode)==HEATMAP_SMOOTHING_CLUSTER) {
            df<-smoothHeatmapAccordingToClusters(df,dbg.heatmap.smooth)
            if (dbg.heatmap.smooth>1) printVar(df$xOrderedSmoothed)
            df$heatmap<-adaptHeatmapForDrawing(as.matrix(df$xOrderedSmoothed))
        }
        return(df)
    }
    updateClustersOnChangePostDraw<-function(df) {
        # colorize external plotmatrix, if any
        if (!is.null(colorizeCallback)) {
            if (dbg.colorizeCb) cat(' calling colorizeCallback\n')
            colorizeCallback(df$leafColorIdxs)
            if (dbg.colorizeCb) cat(' colorizeCallback called\n')
        }
    }

    updateHeatmap<-function(df) {
        if (dbg.heatmap) cat('updateHeatmap called\n')
        if (dbg.heatmap) printVar(tclvalue(heatmapSmoothingMode))
        if (tclvalue(heatmapSmoothingMode)==HEATMAP_SMOOTHING_CLUSTER) {
            if (dbg.heatmap.smooth) cat('smoothing according to clusters\n')
            df$elemClusterCount<-df$n
            df<-smoothHeatmapAccordingToClusters(df,dbg.heatmap.smooth)
            if (dbg.heatmap.smooth>1) printVar(df$xOrderedSmoothed)
            df$heatmap<-adaptHeatmapForDrawing(as.matrix(df$xOrderedSmoothed))
        } else if (tclvalue(heatmapSmoothingMode)==HEATMAP_SMOOTHING_ZOOM) {
            if (dbg.heatmap.smooth) cat('smoothing according to zoom\n')
            ch<-cutree(df$h,h=stillXlim[2])
            if (dbg.heatmap.smooth) printVar(max(ch))
            if (dbg.heatmap.smooth) printVar(df$elemClusterCount)
            if (max(ch)!=df$elemClusterCount) {
                if (dbg.heatmap) cat('smoothing heat map\n')
                df$xOrderedSmoothed<-smoothHeatmap(df$xOrdered,ch[df$leafOrder],dbg.heatmap.smooth)
                df$elemClusterCount<-max(ch)
                df$heatmap<-adaptHeatmapForDrawing(as.matrix(df$xOrderedSmoothed))
            }
        } else if (tclvalue(heatmapSmoothingMode)==HEATMAP_SMOOTHING_NONE) {
            if (dbg.heatmap.smooth) cat('reinitialing heatmap\n')
            df$elemClusterCount<-df$n
            df$heatmap<-adaptHeatmapForDrawing(df$xOrdered)
        } else {
            stop('invalid heatmap smoothing mode')
        }
        if (dbg.heatmap) printVar(dim(df$heatmap))
        if (dbg.heatmap.smooth>1) printVar(df$heatmap)

        return(df)
    }

    # Callback invoked when heat map smoothing mode gets changed (via
    # GUI buttons).
    heatmapSmoothingChanged<-function(mode) {
        #heatmapSmoothingModeChanged<<-TRUE##RMV
        df<<-updateHeatmap(df)
        scalingReplot(img)
    }

    ###################################################################
    ## start of graphics
    ###################################################################
    tt<-tcltk::tktoplevel()

    done<-tclVar(0)

    ACTION_TYPE_SELECT<-'select'
    ACTION_TYPE_ZOOM<-'zoom'
    ACTION_TYPE_PAN<-'pan'
    actionType<-tclVar()
    tclvalue(actionType)<-ACTION_TYPE_SELECT

    PAN_TYPE_VERTICAL<-'vertical'
    PAN_TYPE_VERTICAL_AND_HORIZONTAL<-'vertical+horizontal'

    currentCluster<-tclVar()
    tclvalue(currentCluster)<-1

    zoomMode<-tclVar()
    ZOOM_MODE_IN<-0
    ZOOM_MODE_OUT<-1
    tclvalue(zoomMode)<-ZOOM_MODE_IN

    heatmapSmoothingMode<-tclVar()
    HEATMAP_SMOOTHING_NONE<-'heatmapSmoothingNone'
    HEATMAP_SMOOTHING_CLUSTER<-'heatmapSmoothingCluster'
    HEATMAP_SMOOTHING_ZOOM<-'heatmapSmoothingZoom'
    switch (heatmapSmoothing,
      'none'={tclvalue(heatmapSmoothingMode)<-HEATMAP_SMOOTHING_NONE},
      'cluster'={tclvalue(heatmapSmoothingMode)<-HEATMAP_SMOOTHING_CLUSTER},
      'zoom'={tclvalue(heatmapSmoothingMode)<-HEATMAP_SMOOTHING_ZOOM}
    )

    ## full (maximal) limits of full zoom-out
    fullXlim<-c(max(h$height),0)
    fullYlim<-.5+c(0,n)
    # current drawing limits
    xlim<-stillXlim<-fullXlim
    ylim<-stillYlim<-fullYlim
    # The rectangle to zoom to gets defined by right clicking the top
    # left corner, dragging the mouse and releasing the right mouse
    # button at the bottom right corner. zoomTopLeftXY defines the top
    # left corner, zoomBottomRightXY defines the bottom right corner,
    # which is set to c(NA,NA) during the mouse dragging, and finally
    # set once the mouse is released.
    zoomTopLeftXY<-c(NA,NA)
    zoomBottomRightXY<-c(NA,NA)

    # has plot been zoomed or panned (if yes, we need to refresh the
    # physical<->user coordinate mapping)
    zoomed<-FALSE
    panned<-FALSE

    # during the "zoom rectangle selection" this is
    # the current zooming bottom right corner:
    zoomingBottomRightXY<-c(NA,NA)

    heatmapSmoothingModeChanged<-FALSE

    # panning
    panStartXY<-c(NA,NA)
    panEndXY<-c(NA,NA)
    panShiftXY<-c(NA,NA)

    plotRegionUsrCoords<-rep(NA,4)
    plotRegionDevRel<-rep(NA,4)
    # the number the replot() has been called (if more than 1 times, fixed to 2)
    callCount<-0

    devWidth<-NA
    devHeight<-NA

    clusterCuttingHeight<-NA

    #
    # graphics coordinates conversions
    #
    my_grconvertX<-function(x,from,to) {
        usrMinX<-plotRegionUsrCoords[1]
        usrMaxX<-plotRegionUsrCoords[2]
        usrRangeX<-usrMaxX-usrMinX
        devMinX<-plotRegionDevRel[1]*devWidth
        devMaxX<-plotRegionDevRel[2]*devWidth
        devRangeX<-devMaxX-devMinX
        if (from=='device' && to=='user') {
            return(usrMinX+(x-devMinX)*usrRangeX/devRangeX)
        } else if (from=='user' && to=='device') {
            return(devMinX+(x-usrMinX)*devRangeX/usrRangeX)
        } else stop(paste('unsupported request \'from\' =',from,', \'to\' =',to))
    }
    my_grconvertY<-function(y,from,to) {
        usrMinY<-plotRegionUsrCoords[3]
        usrMaxY<-plotRegionUsrCoords[4]
        usrRangeY<-usrMaxY-usrMinY
        devMinY<-plotRegionDevRel[3]*devHeight
        devMaxY<-plotRegionDevRel[4]*devHeight
        devRangeY<-devMaxY-devMinY
        if (from=='device' && to=='user') {
            y<-devHeight-y
            y<-usrMinY+(y-devMinY)*usrRangeY/devRangeY
        } else if (from=='user' && to=='device') {
            y<-devMinY+(y-usrMinY)*devRangeY/usrRangeY
            y<-devHeight-y
        } else stop(paste('unsupported request \'from\' =',from,', \'to\' =',to))
    }
    dev2usrX<-function(x) my_grconvertX(x,from='device',to='user')
    dev2usrY<-function(y) my_grconvertY(y,from='device',to='user')
    dev2usrXY<-function(x,y) c(dev2usrX(x),dev2usrY(y))
    usr2devX<-function(x) my_grconvertX(x,from='user',to='device')
    usr2devY<-function(y) my_grconvertY(y,from='user',to='device')

    printClusterInfo<-function(leafCount) {
        txt<-paste('"',as.character(leafCount),ifelse(leafCount>0,sprintf('(%d%%)',round(100*leafCount/n)),''),'"',sep='')
    }
    printClusterInfoTotalOverview<-function(df) {
        txt<-paste(df$n,' total',sep='')
    }
    printClusterInfoBrushedOverview<-function(df) {
        txt<-paste(df$fetchedLeafCount,' brushed',sep='')
    }

    panXlim<-function(xlim,shift,fullXlim) {
      xlim<-xlim-shift
      if (xlim[1]>fullXlim[1]) xlim<-xlim-(xlim[1]-fullXlim[1])
      if (xlim[2]<fullXlim[2]) xlim<-xlim-(xlim[2]-fullXlim[2])
      return(xlim)
    }
    panYlim<-function(ylim,shift,fullYlim) {
      ylim<-ylim-shift
      if (ylim[1]<fullYlim[1]) ylim<-ylim-(ylim[1]-fullYlim[1])
      if (ylim[2]>fullYlim[2]) ylim<-ylim-(ylim[2]-fullYlim[2])
      return(ylim)
    }
    zoomOut<-function(xy,zoomRatio=1.41) {
        if (dbg.dendro.zoom) cat(sprintf('zoom-out, zoomRatio %.3f\n',zoomRatio))

        if (zoomRatio>1) {
            tmp<-c(stillXlim[1],max(0,stillXlim[2])) # ignore the heatmap
            if (dbg.dendro.zoom) printVar(stillXlim)
            if (dbg.dendro.zoom) printVar(stillYlim)
            if (dbg.dendro.zoom) printVar(tmp)
            xlimHalfSize<-diff(range(tmp))/2*zoomRatio
            ylimHalfSize<-diff(range(stillYlim))/2*zoomRatio
            if (dbg.dendro.zoom) printVar(xlimHalfSize)
            if (dbg.dendro.zoom) printVar(ylimHalfSize)

            # adapt x limits
            xlim<-xy[1]+c(xlimHalfSize,-xlimHalfSize)
            if (dbg.dendro.zoom) printVar(xlim)
            # make the whole X view to double in size, i.e. do
            # not necessarily center around the click point if
            # it is too close to some border
            restX<-max(max(0,xlim[1]-max(h$height)),max(0,-xlim[2]))
            if (dbg.dendro.zoom) printVar(restX)
            xlim<-xlim+c(restX,-restX)
            if (dbg.dendro.zoom) printVar(xlim)
            xlim[1]<-min(xlim[1],max(h$height))
            xlim[2]<-max(xlim[2],0)
            if (dbg.dendro.zoom) printVar(xlim)

            # adapt y limits
            ylim<-xy[2]+c(-ylimHalfSize,ylimHalfSize)
            if (dbg.dendro.zoom) printVar(ylim)
            restY<-max(max(0,.5-ylim[1]),max(ylim[2]-(.5+n)))
            if (dbg.dendro.zoom) printVar(restY)
            ylim<-ylim+c(-restY,restY)
            if (dbg.dendro.zoom) printVar(ylim)
            ylim[1]<-max(ylim[1],.5)
            ylim[2]<-min(ylim[2],.5+n)
            if (dbg.dendro.zoom) printVar(ylim)

            # propagate to the environment of the idendro function
            stillXlim<<-xlim<<-xlim
            stillYlim<<-ylim<<-ylim
            df<<-updateHeatmap(df)
         }
    }
    zoomIn<-function(xy,zoomRatio=.71) {
        if (dbg.dendro.zoom) cat(sprintf('zoom-in, zoomRatio %.3f\n',zoomRatio))

        if (zoomRatio<1) {
            xlimHalfSize<-diff(range(stillXlim))*zoomRatio/2
            ylimHalfSize<-diff(range(stillYlim))*zoomRatio/2

            # adapt x limits
            xlim<-xy[1]+c(xlimHalfSize,-xlimHalfSize)
            if (dbg.dendro.zoom) printVar(xlim)
            restX<-max(max(0,xlim[1]-max(h$height)),max(0,-xlim[2]))
            if (dbg.dendro.zoom) printVar(restX)
            xlim<-xlim+c(restX,-restX)
            if (dbg.dendro.zoom) printVar(xlim)
            xlim[1]<-min(xlim[1],stillXlim[1])
            xlim[2]<-max(xlim[2],stillXlim[2])
            if (dbg.dendro.zoom) printVar(xlim)

            # adapt y limits
            ylim<-xy[2]+c(-ylimHalfSize,ylimHalfSize)
            if (dbg.dendro.zoom) printVar(ylim)
            restY<-max(max(0,.5-ylim[1]),max(ylim[2]-(.5+n)))
            if (dbg.dendro.zoom) printVar(restY)
            ylim<-ylim+c(-restY,restY)
            if (dbg.dendro.zoom) printVar(ylim)
            ylim[1]<-max(ylim[1],stillYlim[1])
            ylim[2]<-min(ylim[2],stillYlim[2])
            if (dbg.dendro.zoom) printVar(ylim)

            # propagate to the environment of the idendro function
            stillXlim<<-xlim<<-xlim
            stillYlim<<-ylim<<-ylim
            df<<-updateHeatmap(df)
        }
    }

    zoom<-function(zoomTopLeftXY,zoomBottomRightXY) {
        if (dbg.dendro.zoom) cat(' zooming to user-specified rectangle\n')

        xlim<-c(zoomTopLeftXY[1],zoomBottomRightXY[1])
        ylim<-c(zoomBottomRightXY[2],zoomTopLeftXY[2])
        ylim[1]<-max(ylim[1],stillYlim[1])
        ylim[2]<-min(ylim[2],stillYlim[2])

        # propagate to the environment of the idendro function
        stillXlim<<-xlim<<-xlim
        stillYlim<<-ylim<<-ylim
        df<<-updateHeatmap(df)
    }

    # w=width,h=height,i=window, see https://www.tcl.tk/man/tcl8.4/TkCmd/bind.htm
    handleResize<-function(w,h,i) {
        if (dbg.geometry) cat('handleResize called\n')
        w<-as.numeric(w)
        h<-as.numeric(h)
        if (dbg.geometry) printVar(w)
        if (dbg.geometry) printVar(h)
        if (dbg.geometry) printVar(i)
        if (i==as.character(tkwinfo('id',tt))) {
            # configuring the main window `tt'
            windowWidth<-as.numeric(tclvalue(tkwinfo("reqwidth",tt)))
            windowHeight<-as.numeric(tclvalue(tkwinfo("reqheight",tt)))
            imgWidth<-as.numeric(tclvalue(tkwinfo("reqwidth",img)))
            imgHeight<-as.numeric(tclvalue(tkwinfo("reqheight",img)))
            if (dbg.geometry) printVar(windowWidth)
            if (dbg.geometry) printVar(windowHeight)
            if (dbg.geometry) printVar(imgWidth)
            if (dbg.geometry) printVar(imgHeight)
            #if (dbg.geometry) printVar(guiWidth)
            # disabled, not working, TODO
            # the problem is that one an image is created by tkrplot,
            # it is of fixed size, see help for tkrplot:
            #
            #   The function 'tkrplot' creates and returns a Tk label
            #   widget containing a Tk image of type Rplot.  For now
            #   the size is hard-wired. 
            #if (abs(w-guiWidth-imgWidth)>10 || abs(h-imgHeight)>10) {
            #    hscale<<-vscale<<-min((w-guiWidth)/imgWidth,h/imgHeight)
            #    if (dbg.geometry) printVar(vscale)
            #    if (dbg.geometry) printVar(hscale)
            #    scalingReplot(img)
            #}
        }
    }

    ########################################
    # The core plotting function.
    #######################################
    scalingReplot<-function(img) {
        if (dbg.geometry) cat(sprintf('scalingReplot: hscale %.3f, vscale %.3f\n',hscale,vscale))
        tkrreplot(img,hscale=hscale,vscale=vscale)
    }

    replot<-function () {
        if (dbg) cat('-------- replot called --------\n')
        if (dbg.dendro.zoom) printVar(zoomTopLeftXY)
        if (dbg.dendro.zoom) printVar(zoomBottomRightXY)

        ##
        ## compute plotting limits
        ##
        if (dbg.dendro.zoom) printVar(xlim)
        if (dbg.dendro.zoom) printVar(ylim)
        # compute xlim and ylim to include heatmap and heatmap legend
        xlimIncludingMaps<-xlim
        if (heatmapEnabled) {
            # heatmapRelSize is relative to dendro+heatmap+brushedmap
            # heatmapRelSizeRelToDendro si relative to dendro only
            heatmapRelSizeRelToDendro<-heatmapRelSize/(1-heatmapRelSize-brushedmapRelSize)
            # heatmap will extend to the very right border (use `df$k-1')
            # (otherwise it might be separated by a tiny strip (use `df$k-.5'))
            xlimIncludingMaps[2]<-xlim[2]+diff(xlim)*heatmapRelSizeRelToDendro#*.99
            if (df$k>1) {
                heatmapXs<-heatmapXs4Labels<-((df$k:1)-.5)/df$k*diff(xlim)*heatmapRelSizeRelToDendro+xlim[2]
            } else {
                heatmapXs<-c(1,0)*diff(xlim)*heatmapRelSizeRelToDendro+xlim[2]
                heatmapXs4Labels<-.5*diff(xlim)*heatmapRelSizeRelToDendro+xlim[2]
            }
            if (dbg.dendro.zoom) printVar(heatmapXs)
        } else {
            heatmapRelSizeRelToDendro<-heatmapRelSize
        }
        if (brushedmapEnabled) {
            # brushedmapRelSize is relative to dendro+heatmap+brushedmap
            # brushedmapRelSizeRelToDendro si relative to dendro only
            brushedmapRelSizeRelToDendro<-brushedmapRelSize/(1-heatmapRelSize-brushedmapRelSize)
            if (heatmapEnabled) {
                fetchedMapXs<-c(1,0)*diff(xlim)*brushedmapRelSizeRelToDendro+xlimIncludingMaps[2]
                fetchedMapXs4Labels<-.5*diff(xlim)*brushedmapRelSizeRelToDendro+xlimIncludingMaps[2]
            } else {
                fetchedMapXs<-c(1,0)*diff(xlim)*brushedmapRelSizeRelToDendro+xlim[2]
                fetchedMapXs4Labels<-.5*diff(xlim)*brushedmapRelSizeRelToDendro+xlim[2]
            }
            xlimIncludingMaps[2]<-xlimIncludingMaps[2]+diff(xlim)*brushedmapRelSizeRelToDendro
        }
        if (dbg.dendro.zoom) printVar(xlimIncludingMaps)

        ##
        ## setup plot
        ##
        xlab<-paste('height',ifelse(!is.na(clusterCuttingHeight),paste(' (cutting at ',format(clusterCuttingHeight),')',sep=''),''),sep='')
        opar<-par(ask=FALSE,...)
        plot(c(max(h$height),-max(h$height)*heatmapRelSizeRelToDendro),c(0,n*1),type='n',
            frame.plot=T,fg='gray',ylab='',xlim=xlimIncludingMaps,xlab=xlab,ylim=ylim,yaxt="n",xaxs='i',yaxs='i',xaxt='n')
        # only positive x axis ticks make sense
        xl<-axTicks(1)
        axis(1,xl[xl>=0])

        if (observationAnnotationEnabled) {
            axis(4,at=1:n,labels=df$observationLabelsOrdered,tick=FALSE,las=1,cex.axis=.7)
        }

        if (callCount<2) {
            callCount<<-callCount+1
            if (callCount==1) {
                plotRegionDevRel<<-par("plt")
                if (dbg.dendro.zoom) printVar(plotRegionDevRel)
            }
        }

        if (zoomed || panned || callCount==2) {
            plotRegionUsrCoords<<-par("usr")
            if (dbg.dendro.zoom) printVar(plotRegionUsrCoords)
        }

        if (dbg.geometry>1) {
            print((tkwinfo("geometry",tt)))
            print((tkwinfo("reqwidth",tt)))
            print((tkwinfo("geometry",img)))
            print((tkwinfo("reqwidth",img)))
        }

        ##
        ## dendrogram
        ##
        #with(visibleUnselectedBranches$branches,segments(x1s,y1s,x2s,y2s))
        if (dbg) cat('plotting dendrogram\n')
        with(df$unselectedBranches$branches,segments(x1s,y1s,x2s,y2s))
        if (dbg) cat(' plotting subclusters\n')
        for (i in seq(along=df$clusters)) {
            if (!is.null(df$clusters[[i]]) && length(df$clusters[[i]]$branches)>0) {
                if (dbg) cat(sprintf('cluster %i: color %s\n',i,clusterColor(i)))
                with(df$clusters[[i]]$branches,segments(x1s,y1s,x2s,y2s,col=clusterColor(i)))
            }
        }

        if (length(df$fetchedBranches$indices)>0) {
            if (dbg) cat(' plotting externally selected branches\n')
            if (dbg.fetched>1) printVar(df$fetchedBranches)

            # in which color to draw the externally selected branches?
            if (0) {
                # get color of the first unused cluster
                col<-clusterColors[1]
                clusterFound<-0
                for (i in seq(along=df$clusters)) {
                    if (is.null(df$clusters[[i]]) || length(df$clusters[[i]]$indices)==0) {
                        if (dbg.fetched) cat(sprintf('cluster %d not used\n',i))
                        clusterFound<-i
                        col<-clusterColor(i)
                        break
                    }
                }
            } else {
                # use color of the current cluster
                clusterFound<-as.numeric(tclvalue(currentCluster))
                col<-clusterColor(clusterFound)
            }

            if (dbg.fetched) printVar(col)
            with(df$fetchedBranches$branches,segments(x1s,y1s,x2s,y2s,col=col))
            if (clusterFound>0) {
                fetchedInfo<<-sprintf("%d points selected (%d subclusters), drawn in the color of cluster %d",df$fetchedLeafCount,length(df$fetchedBranches$indices),clusterFound)
            } else {
                fetchedInfo<<-sprintf("%d points selected (%d subclusters), but not drawn: there is no free cluster",df$fetchedLeafCount,length(df$fetchedBranches$indices))
            }
        } else if (df$fetchedLeafCount>0) {
            fetchedInfo<<-sprintf("%d point(s) selected, which formed no subclusters - nothing to draw.",df$fetchedLeafCount)
        }

        if (!is.na(clusterCuttingHeight)) {
            lines(rep(clusterCuttingHeight,2),ylim,col='gray',lwd=2,lty=2)
        }

        if (all(!is.na(zoomingBottomRightXY))) {
            if (dbg) cat('drawing current zooming rectangle\n')
            lines(c(zoomTopLeftXY[1],zoomingBottomRightXY[1],zoomingBottomRightXY[1],zoomTopLeftXY[1],zoomTopLeftXY[1]),
                c(zoomTopLeftXY[2],zoomTopLeftXY[2],zoomingBottomRightXY[2],zoomingBottomRightXY[2],zoomTopLeftXY[2]),
                col='gray',lwd=2,lty=2)
        }

        ##
        ## heatmap
        ##
        if (heatmapEnabled) {
            if (dbg) cat(' plotting heatmap\n')
            if (dbg.heatmap) printVar(dim(heatmap))
            if (dbg.heatmap) printVar(dim(heatmapXs))
            if (dbg.heatmap) printVar(length(heatmapXs))
            if (dbg.heatmap) printVar(length(1:n))
            image(x=heatmapXs,y=1:n,df$heatmap,add=T,col=heatmapColors)
            dnames<-dimnames(df$heatmap)
            if (dbg.heatmap) printVar(dnames)
            if (!is.null(dnames) && !is.null(dnames[[1]])) {
                #text(heatmapXs,ylimIncludingHeatMap[2],dnames[[2]],srt=90,pos=3,cex=.7,offset=c(1,0))
                axis(3,at=heatmapXs4Labels,labels=dnames[[1]],las=2,cex.axis=.7)
            }
            if (dbg) cat(' heatmap plotted\n')
        }

        ##
        ## brushedmap
        ##
        if (brushedmapEnabled) {
            if (dbg) cat(' plotting ext. selected map\n')
            image(x=fetchedMapXs,y=1:n,df$fetchedMap,add=T,col=c('white','black'))#clusterColor(as.numeric(tclvalue(currentCluster)))))
            #axis(3,at=fetchedMapXs4Labels,labels=paste('fetched',df$fetchedLeafCount),las=2,cex.axis=.7)
            axis(3,at=fetchedMapXs4Labels,labels='fetched',las=2,cex.axis=.7)
            if (dbg) cat(' ext.selected map plotted\n')
        }
        par(opar)
    }

    ################
    ## zooming stuff
    undoZoom<-function() {
        if (dbg.dendro.zoom) cat('undoZoom called\n')
        rv<-popDendroZoomHistory(df,dbg.dendro.zoom)
        if (!is.null(rv$dendroZoom)) {
            if (dbg) cat('using zoom history\n')
            xlim<<-stillXlim<<-rv$dendroZoom$xlim
            ylim<<-stillYlim<<-rv$dendroZoom$ylim
            plotRegionUsrCoords<<-rv$dendroZoom$plotRegionUsrCoords
            df<<-updateHeatmap(rv$df)
            scalingReplot(img)
        } else {
            if (dbg) cat('no zoom lim history\n')
        }
    }
    fullView<-function() {
        if (dbg) cat('fullView called\n')
        df<<-pushDendroZoomHistory(df,list(xlim=stillXlim,ylim=stillYlim,
            plotRegionUsrCoords=plotRegionUsrCoords),dbg.dendro.zoom)
        #zoomTopLeftXY<<-zoomBottomRightXY<<-c(NA,NA)
        xlim<<-stillXlim<<-fullXlim
        ylim<<-stillYlim<<-fullYlim
        scalingReplot(img)
    }
    deselectCurrentCluster<-function() {
        if (dbg) cat('deselectCurrentCluster called\n')
        df$currentCluster<<-as.numeric(tclvalue(currentCluster))
        rv<-unselectCurrentCluster(df,dbg.dendro.select)
        if (rv$selectionChanged) {
            df<<-updateClustersOnChangePreDraw(rv$df)
            scalingReplot(img)
            updateClustersOnChangePostDraw(df)
        }
    }
    deselectAllClusters<-function() {
        if (dbg) cat('deselectAll called\n')
        rv<-unselectAllClusters(df,dbg.dendro.select)
        if (rv$selectionChanged) {
            df<<-updateClustersOnChangePreDraw(rv$df)
            scalingReplot(img)
            updateClustersOnChangePostDraw(df)
        }
    }
    undoSelection<-function() {
        if (dbg) cat('undoSelection called\n')
        rv<-popSelectionHistory(df,dbg.dendro.select)
        if (!is.null(rv)) {
            df<<-updateClustersOnChangePreDraw(rv)
            # instruct replot to redraw clusters and colorize external plotmatrix, if any
            scalingReplot(img)
            updateClustersOnChangePostDraw(df)
        }
    }
    fetchSelected<-function() {
        if (dbg) cat('fetchSelected called\n')
        if (!is.null(fetchSelectedCallback)) {
            if (TRUE) {
                # fetching in the brushed map only
                df$brushed<<-fetchedLeafs<-fetchSelectedCallback()
                if (dbg.fetched) printVar(fetchedLeafs)
                if (length(fetchedLeafs)!=df$n || !is.logical(fetchedLeafs)) {
                    cat('error: fetchSelectedCallback returned invalid value\n')
                } else {
                    df$brushed<<-fetchedLeafs
                    df$fetchedLeafCount<<-sum(fetchedLeafs)
                    tclvalue(clusterInfoBrushedOverviewVar)<-printClusterInfoBrushedOverview(df)
                    if (dbg.fetched) printVar(df$fetchedLeafCount)
                    df$fetchedMap<<-matrix(fetchedLeafs[df$leafOrder],nrow=1)
                    updateBrushedClusterInfosOnChange(df)
                    scalingReplot(img)
                }
            } else {
                # fetching into the current cluster - DISABLED

                # first check that the current clusters is empty
                i<-as.numeric(tclvalue(currentCluster))
                if (!is.null(df$clusters) && length(df$clusters)>=i && length(df$clusters[[i]]$indices)>0) {
                    if (!silent) tkmessageBox(message=sprintf("Can't fetch into the current cluster, it is not empty. Select another cluster."))
                } else {
                    fetchedLeafs<-fetchSelectedCallback()
                    df$fetchedLeafCount<<-sum(fetchedLeafs)
                    if (dbg.fetched) printVar(df$fetchedLeafCount)
                    if (df$fetchedLeafCount) {
                        # find clusters whose members are all brushed
                        fetchedClusters<-rep(FALSE,df$clusterCount)
                        for (i in 1:df$clusterCount) {
                            if ((h$merge[i,1]<0 && fetchedLeafs[-h$merge[i,1]] || h$merge[i,1]>0 && fetchedClusters[h$merge[i,1]]) &&
                                (h$merge[i,2]<0 && fetchedLeafs[-h$merge[i,2]] || h$merge[i,2]>0 && fetchedClusters[h$merge[i,2]])) {
                                fetchedClusters[i]<-TRUE
                            }
                        }
                        fetchedIdx<-which(fetchedClusters)
                        if (dbg.fetched>1) printVar(fetchedIdx)
                        idx<-clusterId2SegmentIds(fetchedIdx)
                        df$fetchedBranches<<-list(indices=fetchedIdx,branches=with(df$allBranches$branches,data.frame(x1s=x1s[idx],y1s=y1s[idx],x2s=x2s[idx],y2s=y2s[idx])))
                        if (dbg.fetched>1) printVar(df$fetchedBranches)
                        df$fetchedMap<<-matrix(fetchedLeafs[df$leafOrder],nrow=1)
                        scalingReplot(img)
                        if (!is.null(fetchedInfo)) if (!silent) tkmessageBox(message=fetchedInfo)
                        fetchedInfo<<-NULL
                    } else {
                        if (!silent) tkmessageBox(message=sprintf("No points selected."))
                    }
                }
            }
        }
    }
    doQuit<-function() {
        if (dbg) cat('doQuit called\n')
        #tkconfigure(tt,cursor="top_left_arrow")
        tclvalue(done)<-1
    }

    img<-tkrplot(tt,replot,vscale=vscale,hscale=hscale)
    # TODO: set geometry according to hscale,vscale?
    if (!is.null(geometry)) {
        tkwm.geometry(tt,sprintf('%dx%d+%d+%d',geometry[3],geometry[4],geometry[1],geometry[2]))
    }
    # to get a list of possible options, use:
    #printVar(tkwinfo(tt,"screenwidth"))

    devWidth<-as.numeric(tclvalue(tkwinfo("reqwidth",img)))
    devHeight<-as.numeric(tclvalue(tkwinfo("reqheight",img)))
    #tkconfigure(tt,cursor="ul_angle")
    clusterHeadingVar<-tclVar()
    tclvalue(clusterHeadingVar)<-"cluster"
    tkgrid(tklabel(tt,textvariable=clusterHeadingVar),row=1,column=1)#,sticky='n')

    clusterInfoTotalOverviewVar<-tclVar()
    tclvalue(clusterInfoTotalOverviewVar)<-printClusterInfoTotalOverview(df)
    tkgrid(tklabel(tt,textvariable=clusterInfoTotalOverviewVar),row=1,column=2)

    if (brushedmapEnabled) {
        clusterInfoBrushedOverviewVar<-tclVar()
        tclvalue(clusterInfoBrushedOverviewVar)<-printClusterInfoBrushedOverview(df)
        tkgrid(tklabel(tt,textvariable=clusterInfoBrushedOverviewVar),row=1,column=3)
    }

    for (i in 1:maxClusterCount) {
        tkgrid(tkradiobutton(tt,variable=currentCluster,value=i,text=paste(i),background=clusterColor(i)),row=i+1,column=1)
    }
    for (i in 1:maxClusterCount) {
        clusterTotalInfoTxt<-sprintf('clusterTotalInfo%d',i)
        clusterTotalInfoExpr<-parse(text=clusterTotalInfoTxt)
        eval(parse(text=paste(clusterTotalInfoTxt,'<-tclVar()')))
        leafCount<-ifelse(length(df$clusters)<i || is.null(df$clusters[[i]]),0,length(df$clusters[[i]]$indices))
        clusterTotalInfo<-printClusterInfo(leafCount)
        eval(parse(text=paste('tclvalue(',clusterTotalInfoTxt,')<-',clusterTotalInfo)))
        tkgrid(tklabel(tt,textvariable=eval(clusterTotalInfoExpr)),row=i+1,column=2)

        if (brushedmapEnabled) {
            clusterBrushedInfoTxt<-sprintf('clusterBrushedInfo%d',i)
            clusterBrushedInfoExpr<-parse(text=clusterBrushedInfoTxt)
            eval(parse(text=paste(clusterBrushedInfoTxt,'<-tclVar()')))
            leafCount<-0
            clusterBrushedInfo<-printClusterInfo(leafCount)
            eval(parse(text=paste('tclvalue(',clusterBrushedInfoTxt,')<-',clusterBrushedInfo)))
            tkgrid(tklabel(tt,textvariable=eval(clusterBrushedInfoExpr)),row=i+1,column=3)
        }
    }

    zoomModeHeading<-tclVar()
    rowCurr<-maxClusterCount+2
    tkgrid(tkbutton(tt,text="undo Zoom",command=undoZoom),row=rowCurr,column=1)
    tkgrid(tkbutton(tt,text="Full view",command=fullView),row=rowCurr,column=2)
    rowCurr<-rowCurr+1
    tkgrid(tkbutton(tt,text="undo Selection",command=undoSelection),row=rowCurr,column=1,columnspan=2)
    rowCurr<-rowCurr+1
    tkgrid(tkbutton(tt,text="Unselect",command=deselectCurrentCluster),row=rowCurr,column=1)
    tkgrid(tkbutton(tt,text="unselect All",command=deselectAllClusters),row=rowCurr,column=2)
    rowCurr<-rowCurr+1
    fetchColor<-ifelse(is.null(fetchSelectedCallback),'gray','black')
    tkgrid(tkbutton(tt,text="fetCh selected",command=fetchSelected,foreground=fetchColor,activeforeground=fetchColor),row=rowCurr,column=1,columnspan=2)
    rowCurr<-rowCurr+1
    heatmapSmoothingModeHeading<-tclVar()
    tclvalue(heatmapSmoothingModeHeading)<-"heatmap smoothing:"
    tkgrid(tklabel(tt,textvariable=heatmapSmoothingModeHeading),row=rowCurr,column=1,columnspan=3)
    rowCurr<-rowCurr+1
    heatmapSmoothingModeNoneButton<-tkradiobutton(tt,variable=heatmapSmoothingMode,value=HEATMAP_SMOOTHING_NONE,text="none",command=heatmapSmoothingChanged)
    tkgrid(heatmapSmoothingModeNoneButton,row=rowCurr,column=1)
    heatmapSmoothingModeClusterButton<-tkradiobutton(tt,variable=heatmapSmoothingMode,value=HEATMAP_SMOOTHING_CLUSTER,text="cluster",command=heatmapSmoothingChanged)
    tkgrid(heatmapSmoothingModeClusterButton,row=rowCurr,column=2)
    heatmapSmoothingModeZoomButton<-tkradiobutton(tt,variable=heatmapSmoothingMode,value=HEATMAP_SMOOTHING_ZOOM,text="zoom",command=heatmapSmoothingChanged)
    tkgrid(heatmapSmoothingModeZoomButton,row=rowCurr,column=3)
    rowCurr<-rowCurr+1
    tkgrid(tkbutton(tt,text="Quit",command=doQuit),row=rowCurr,column=1,columnspan=2)
    #rowCurr<-rowCurr+1
    #tkgrid(tkframe(tt),row=rowCurr,column=1,rowspan=1,columnspan=3)
    tkgrid(img,row=1,column=4,rowspan=rowCurr)

    mouse.down.impl <- function(x,y,actionType) {
        if (dbg) cat(paste('mouse.down.impl called, actionType',actionType,'\n'))
        if (dbg) printVar(c(x,y))
        x<-dev2usrX(as.numeric(x))
        y<-dev2usrY(as.numeric(y))
        if (dbg) printVar(c(x,y))
        if (actionType==ACTION_TYPE_SELECT) {
            if (dbg) cat('action:select\n')
            if (y>=stillYlim[1]) {
                if (dbg) cat(' selecting cluster\n')
                dendroZoom<-xy2gw(list(x=stillXlim,y=stillYlim))
                pos.xy<-list(x=x,y=y)
                origLeafColorIdxs<-df$leafColorIdxs
                df$currentCluster<<-as.numeric(tclvalue(currentCluster))
                df<<-selectCluster(pos.xy,df,dendroZoom,dbg.dendro.select)
                if (any(df$leafColorIdxs!=origLeafColorIdxs)) {
                    df<<-updateClustersOnChangePreDraw(df)
                    scalingReplot(img)
                    updateClustersOnChangePostDraw(df)
                }
            } else {
              if (dbg) cat(' cutting dendrogram\n')
              dendroZoom<-xy2gw(list(x=stillXlim,y=stillYlim))
              x.gw<-xy2gw(list(x=x,y=NA))
              rv<-cutDendro(df,x.gw$g,dendroZoom,dbg.dendro.cut)
              if (dbg.dendro.cut) printVar(length(rv$clusters))
              if (dbg.dendro.cut>1) printVar(rv)
              if (rv$selectedClusterCount>maxClusterCount) {
                if (!silent) tkmessageBox(message=
                  paste('You have selected ',rv$selectedClusterCount,
                    ' clusters, but the maximal configured number of clusters (maxClusterCount) is ',
                    maxClusterCount,'.',sep=''))
              } else {
                df<<-updateClustersOnChangePreDraw(rv$df)
                scalingReplot(img)
                updateClustersOnChangePostDraw(df)
              }
            }
        } else if (actionType==ACTION_TYPE_ZOOM) {
            if (dbg) cat('action:zoom\n')
            zoomTopLeftXY<<-c(x,y)
            zoomBottomRightXY<<-c(NA,NA)
            if (dbg.dendro.zoom) printVar(zoomTopLeftXY)
            #tkconfigure(tt,cursor="lr_angle")
            zoomingBottomRightXY<<-c(x,y)
            scalingReplot(img)
        } else if (actionType==ACTION_TYPE_PAN) {
            if (dbg) cat('action:pan\n')
            panStartXY<<-c(x,y)
            if (dbg.pan) printVar(panStartXY)
            tkconfigure(tt,cursor="fleur")
        } else {
            stop(paste('invalid action type',actionType))
        }
    }
    mouse1.down <- function(x,y) {
        if (dbg) cat('\nmouse1.down called\n')
        mouse.down.impl(x,y,ACTION_TYPE_SELECT)
    }
    mouse2.down <- function(x,y) {
        if (dbg) cat('\nmouse2.down called\n')
        mouse.down.impl(x,y,ACTION_TYPE_PAN)
    }
    mouse3.down <- function(x,y) {
        if (dbg) cat('\nmouse3.down called\n')
        mouse.down.impl(x,y,ACTION_TYPE_ZOOM)
    }
    mouse.up.impl <- function(x,y,actionType) {
        if (dbg) cat(paste('mouse.up.impl called, actionType',actionType,'\n'))
        if (dbg) printVar(c(x,y))
        x<-dev2usrX(as.numeric(x))
        y<-dev2usrY(as.numeric(y))
        if (dbg) printVar(c(x,y))
        if (actionType==ACTION_TYPE_SELECT) {
            if (dbg) cat('action:select\n')
        } else if (actionType==ACTION_TYPE_ZOOM) {
            if (dbg) cat('action:zoom\n')
            zoomBottomRightXY<<-c(x,y)
            if (dbg) printVar(zoomTopLeftXY)
            if (dbg) printVar(zoomBottomRightXY)
            tmp<-c(min(zoomTopLeftXY[1],zoomBottomRightXY[1]),min(zoomTopLeftXY[2],zoomBottomRightXY[2]))
            zoomTopLeftXY<<-c(max(zoomTopLeftXY[1],zoomBottomRightXY[1]),max(zoomTopLeftXY[2],zoomBottomRightXY[2]))
            zoomBottomRightXY<<-tmp

            if (dbg) printVar(stillXlim)
            if (dbg) printVar(stillYlim)
            if (dbg) printVar(zoomTopLeftXY)
            if (dbg) printVar(zoomBottomRightXY)

            if (zoomTopLeftXY[1]-zoomBottomRightXY[1] < -diff(stillXlim)/10 && zoomTopLeftXY[2]-zoomBottomRightXY[2]<diff(stillYlim)/10) {
                if (dbg) cat(' zooming rectangle too small, converting to half-zoom\n')
                zoomBottomRightXY<<-zoomTopLeftXY
            }
            # clip to dendrogram area
            if (zoomTopLeftXY[1]<=0) {
                if (dbg) cat(' can\'t zoom the heatmap')
                zoomTopLeftXY<<-zoomBottomRightXY<-c(NA,NA)
            } else {
                zoomBottomRightXY[1]<<-max(zoomBottomRightXY[1],0)
            }
            #tkconfigure(tt,cursor="ul_angle")
            # stop drawing the current zooming rectangle
            zoomingBottomRightXY<<-c(NA,NA)
            if (all(!is.na(zoomTopLeftXY))) {
                df<<-pushDendroZoomHistory(df,list(xlim=stillXlim,ylim=stillYlim,
                    plotRegionUsrCoords=plotRegionUsrCoords),dbg.dendro.zoom)
                if (any(zoomTopLeftXY!=zoomBottomRightXY)) {
                    zoom(zoomTopLeftXY,zoomBottomRightXY)
                } else {
                    if (as.numeric(tclvalue(zoomMode))==ZOOM_MODE_OUT) {
                        zoomOut(zoomTopLeftXY)
                    } else {
                        zoomIn(zoomTopLeftXY)
                    }
                }
                #rm(df) # get rid of local `df' copy, get access to the global one
                zoomed<<-TRUE
                scalingReplot(img)
            }
            zoomTopLeftXY<<-zoomBottomRightXY<<-zoomingBottomRightXY<<-c(NA,NA)
        } else if (actionType==ACTION_TYPE_PAN) {
            if (dbg) cat('action:pan\n')
            panShiftXY<<-c(x,y)-panStartXY
            if (dbg.pan) printVar(panShiftXY)
            tkconfigure(tt,cursor="top_left_arrow")
            if (any(panShiftXY!=0)) {
                if (dbg.dendro.zoom) cat('panning current limits\n')
                df<<-pushDendroZoomHistory(df,list(xlim=stillXlim,ylim=stillYlim,
                    plotRegionUsrCoords=plotRegionUsrCoords),dbg.dendro.zoom)
                xlim<<-stillXlim<<-panXlim(stillXlim,panShiftXY[1],fullXlim)
                ylim<<-stillYlim<<-panYlim(stillYlim,panShiftXY[2],fullYlim)
                if (dbg.pan) cat('clearing panning variables\n')
                panStartXY<<-panEndXY<<-panShiftXY<<-c(NA,NA)
                df<<-updateHeatmap(df)
                panned<-TRUE
                scalingReplot(img)
            } else {
                if (dbg.pan) cat('clearing panning variables\n')
                panStartXY<<-panEndXY<<-panShiftXY<<-c(NA,NA)
            }
        } else {
            stop(paste('invalid action type',actionType))
        }
    }
    mouse1.up <- function(x,y) {
        if (dbg) cat('\nmouse1.up called\n')
        mouse.up.impl(x,y,ACTION_TYPE_SELECT)
    }
    mouse2.up <- function(x,y) {
        if (dbg) cat('\nmouse2.up called\n')
        mouse.up.impl(x,y,ACTION_TYPE_PAN)
    }
    mouse3.up <- function(x,y) {
        if (dbg) cat('\nmouse3.up called\n')
        mouse.up.impl(x,y,ACTION_TYPE_ZOOM)
    }
    mouse.move.impl <- function(x,y,actionType) {
        if (dbg.mouse) cat('mouse.move.impl called, actionType',actionType,'\n')
        if (dbg.mouse) printVar(c(x,y))
        x<-dev2usrX(as.numeric(x))
        y<-dev2usrY(as.numeric(y))
        if (dbg.mouse) printVar(c(x,y))
        lastClusterCuttingHeight<-clusterCuttingHeight
        clusterCuttingHeight<<-NA
        if (actionType==ACTION_TYPE_SELECT) {
            if (dbg.mouse) cat('action:select\n')
            if (x>stillXlim[2] && y<stillYlim[1]) {
                clusterCuttingHeight<<-x
                if (dbg.dendro.cut) printVar(clusterCuttingHeight)
                if (is.na(clusterCuttingHeight)!=is.na(lastClusterCuttingHeight) ||
                    !is.na(clusterCuttingHeight) && clusterCuttingHeight!=lastClusterCuttingHeight) {
                    scalingReplot(img)
                }
            } else if (!is.na(lastClusterCuttingHeight)) {
                scalingReplot(img)
            }
        } else if (actionType==ACTION_TYPE_ZOOM) {
            if (dbg.mouse) cat('action:zoom\n')
            if (all(!is.na(zoomTopLeftXY))) {
                zoomingBottomRightXY<<-c(x,y)
                scalingReplot(img)
            }
        } else if (actionType==ACTION_TYPE_PAN) {
            if (dbg.mouse) cat('action:pan\n')
            if (!any(is.na(panStartXY))) {
                panEndXY<<-c(x,y)
                if (dbg.pan) printVar(panEndXY)
                if (any(panEndXY-panStartXY!=0)) {
                    # keep panning the dendrogram on the fly
                    if (dbg.dendro.zoom) cat(paste('temporarily padding current limits by',panEndXY[1]-panStartXY[1],'and',panEndXY[2]-panStartXY[2],'\n'))
                    xlim<<-panXlim(stillXlim,panEndXY[1]-panStartXY[1],fullXlim)
                    ylim<<-panYlim(stillYlim,panEndXY[2]-panStartXY[2],fullYlim)
                    ## TODO: set stillXlim as well and recompute
                    # hetmap, or it would be too slow?
                    scalingReplot(img)
                }
            }
        } else {
            stop(paste('invalid action type',actionType))
        }
    }
    mouse.move <- function(x,y) {
        if (dbg.mouse) cat('\nmouse.move called\n')
        mouse.move.impl(x,y,ACTION_TYPE_SELECT)
    }
    mouse2.move <- function(x,y) {
        if (dbg.mouse) cat('\nmouse2.move called\n')
        mouse.move.impl(x,y,ACTION_TYPE_PAN)
    }
    mouse3.move <- function(x,y) {
        if (dbg.mouse) cat('\nmouse3.move called\n')
        mouse.move.impl(x,y,ACTION_TYPE_ZOOM)
    }
    mouse.wheel.up <- function(x,y) {
        if (dbg.mouse) cat('\nmouse.wheel.up called\n')
        mouse.wheel.impl(x,y,1)
    }
    mouse.wheel.down <- function(x,y) {
        if (dbg.mouse) cat('\nmouse.wheel.down called\n')
        mouse.wheel.impl(x,y,-1)
    }
    mouse.wheel.impl <- function(x,y,D) {
        # TODO: handle the value of `D' (the delta value of a
        # MouseWheel event)
        # see https://www.tcl.tk/man/tcl8.4/TkCmd/bind.htm
        if (dbg.mouse) cat('\nmouse.wheel.impl called\n')
        x<-dev2usrX(as.numeric(x))
        y<-dev2usrY(as.numeric(y))
        df<<-pushDendroZoomHistory(df,list(xlim=stillXlim,ylim=stillYlim,
            plotRegionUsrCoords=plotRegionUsrCoords),dbg.dendro.zoom)
        # delta D value ignored for now (not sure what is the typical value,
        # https://www.tcl.tk/man/tcl8.4/TkCmd/bind.htm says it is 120 on Windows,
        # while we observe the value of 1 on linux using the tcltk package)
        if (dbg.dendro.zoom) printVar(D)
        if (dbg.dendro.zoom) printVar(zoomFactor*D)
        if (D>0) {
            zoomIn(c(x,y))
        } else {
            zoomOut(c(x,y))
        }
        zoomed<<-TRUE
        scalingReplot(img)
    }
    tkbind(img, "<ButtonPress-1>", mouse1.down)
    tkbind(img, "<ButtonRelease-1>", mouse1.up)
    tkbind(img, "<ButtonPress-2>", mouse2.down)
    tkbind(img, "<ButtonRelease-2>", mouse2.up)
    tkbind(img, "<ButtonPress-3>", mouse3.down)
    tkbind(img, "<ButtonRelease-3>", mouse3.up)
    tkbind(img, "<Motion>", mouse.move)
    tkbind(img, "<B2-Motion>", mouse2.move)
    tkbind(img, "<B3-Motion>", mouse3.move)
    tkbind(img, "<Button-4>", mouse.wheel.up)
    tkbind(img, "<Button-5>", mouse.wheel.down)
    #tkbind(tt, "<Configure>", handleResize) # unused, as the `img' can't be resized
    tkbind(tt, "f", fullView)
    tkbind(tt, "F", fullView)
    tkbind(tt, "z", undoZoom)
    tkbind(tt, "Z", undoZoom)
    tkbind(tt, "c", fetchSelected)
    tkbind(tt, "C", fetchSelected)
    tkbind(tt, "u", deselectCurrentCluster)
    tkbind(tt, "U", deselectCurrentCluster)
    tkbind(tt, "a", deselectAllClusters)
    tkbind(tt, "A", deselectAllClusters)
    tkbind(tt, "s", undoSelection)
    tkbind(tt, "S", undoSelection)
    tkbind(tt, "q", doQuit) ##TODO: send event instead of destroying in order to get rid of the error "Bad window path name"?
    tkbind(tt, "Q", doQuit)
    for (i in 1:max(9,maxClusterCount)) {
        tkbind(tt, as.character(i), eval(parse(text=paste("function() {tclvalue(currentCluster)<-",i,"}"))))
    }

    scalingReplot(img)
    tkwait.variable(done)
    tkdestroy(tt)

    if (dbg>1) printVar(h$height)
    return(invisible(df$leafColorIdxs))
    ### vector of colors assigned to observations. 0s denote unselected
    ### observations, while values of i > 0 denote the cluster `i'.
},ex=function() {
    if (interactive()) {
        data(iris, envir = environment())
        hc <- hclust(dist(iris[, 1:4]))
        idendro(hc, iris)
    }
    # see demos for more examples
})
