
# authors:  Norm Matloff and Yingkang Xie

# avoids the "black screen problem" in parallel coordinates, using
# various methods: 

#    method = "maxdens":  plots only the most frequent lines (uses
#    density estimation); by plotting only the most typical patterns, we
#    can discern relationships among the variables by avoiding filling
#    up the screen; density values are computed within groups, not overall;
#    one can also plot the LEAST typical patterns, to hunt for outliers

#    method = "locmax":  plots only those lines whose density values
#    are largest in their neighborhood, for the purpose of cluster
#    hunting

#    method = "randsamp":  the parallel coordinates graph will be
#    plotted for a random sampling of the data points, i.e. for a random
#    selection of the rows of x, but the lines will be color-coded
#    according to the density estimate, so that one sees which lines are
#    more "typical"

# arguments:

# x:  data, and number of variables/columns (more than 2)
# m:  number of lines to plot; action depends of method:
#     "maxdens": the most frequent rows of x (via density estimate from x) 
#                are plotted, or for m < 0, the least frequent ones 
#                (outlier hunting)
#     "locmax":  rows corresponding to local maxima of the density are
#                plotted (cluster hunting); m and grpvar are forced to 1
#                and NULL, respectively
#     "randsamp":  random rows from x are plotted, m lines from each group 
#                  (the nongrouping case is treated as 1 group)
# dispcols:  indices of columns of x to plot
# k:  number of nearest neighbors for density estimation
# klm:  number of nearest neighbors for finding local maxima points;
#       only makes sense if method = "locmax"
# grpvar:  grouping variable (use interaction() to combine several);
# method:  "maxdens", "locmax" or "randsamp"; see above comments
# faceting:  "vert" for vertical stacking, "horiz" for horizontal,
#            "none" to draw all groups together but code group by color; 
#            ignored if have no grouping
# keepidxs:  if not NULL, retain the indices of the data points that are
#            plotted, and the points themselves; the latter will be
#            ordered by data column number keepidxs
# plotidxs:  if TRUE, the case numbers, i.e. row numbers within x, will
#            be shown in the plot (subject to overplotting; useful only 
#            if m is small)
# cls:  snow cluster, if present, for parallel computation

# return value:

#    object of class "gg", with additional components idxs and xdisp if
#    keepidxs is not NULL; if the object is printed, the graph is displayed;
#    idxs, if included, will be the indices within x of the data points
#    plotted, and xdisp will be the corresponding rows of x, ordered by
#    x[,keepidxs]

freqparcoord <- function(x,m,dispcols=1:ncol(x),grpvar=NULL,
      method="maxdens",faceting="vert",k=50,klm=5*k,
      keepidxs=NULL,plotidxs=FALSE,cls=NULL) {
   if (!method %in% c("randsamp","maxdens","locmax")) 
      stop('method must be either "randsamp", "locmax" or "maxdens"')
   if (!is.numeric(keepidxs) && !is.null(keepidxs)) 
      stop('keepidxs must be either NULL or a column number')
   mneg <- (m < 0)
   if (mneg) m <- -m
   if (method == "locmax") {
      m <- 1
      grpvar <- NULL
   }
   rownames(x) <- NULL

   # to avoid the rescaling in ggparcoord(); this way the scale will be
   # reflect the entire data set, not just the few rows we choose below;
   # also, retain the original x if needed
   if (!is.null(keepidxs)) xorig <- x
   x[,dispcols] <- scale(x[,dispcols])

   nrx <- nrow(x)
   # determine group indices, i.e. which rows of x belong with which
   # groups; place in list grpidxs, 1 vector per group
   if (is.null(grpvar)) {
      # no grouping, so all rows are in one big "group"
      grpidxs <- list(1:nrx)
   } else {
      grpidxs <- split(1:nrx,x[,grpvar])
   }
   ngrps <- length(grpidxs)
   dens <- vector(length=nrx)
   # dens will contain the estimated density value for each row
   # of x; find it group by group (which means that densities will be
   # within-group)
   for (i in 1:ngrps) {
      ix <- grpidxs[[i]]  # row indices within x for this group
      dens[ix] <-
         smoothz(x[ix,dispcols],knndens,k,cls=cls,scalefirst=T)
   }
   x <- as.data.frame(x)
   x$dens <- dens
   # if (!is.null(grpvar)) x$grp <- x[,grpvar]
   if (!is.null(grpvar)) {
      grp <- x[,grpvar]
      x <- cbind(x,grp)
   }
   dispnms <- names(x[,dispcols])
   # determine which rows to plot
   rws <- NULL
   for (i in 1:ngrps) {
      ix <- grpidxs[[i]]
      if (method == "maxdens") { 
         dns <- if (mneg) -dens else dens
         tmp <- findtop(dns,ix,m)
      } else if (method == "locmax") {
         tmp <- findlocmax(x,dispcols,klm)
      } else {  # "randsamp" case
         tmp <- sample(ix,m,replace=F)
      }
      rws <- c(rws,tmp)
   }
   # drawing a single line in the entire plot can cause problems with
   # some graphic systems 
   xsub <- x[rws,]
   # call our altered ggparcoord() with scale = globalminmax so that
   # rescaler() is NOT called within ggparcoord()
   p <- ggparcoord(xsub,columns=dispcols,groupColumn=grpvar,
           scale="globalminmax")
   if (!is.null(keepidxs)) {
      p$idxs <- rws
      xdisp <- xorig[rws,]
      tmp <- order(xdisp[,keepidxs])
      p$xdisp <- xdisp[tmp,]  
   }
   p$data$dens <- xsub$dens
   if (!is.null(grpvar)) {
      p$data$grp <- xsub$grp
      p$labels$colour <- NULL
      if (faceting == "none")  {
         p$mapping$colour <- NULL
         p <- p + geom_line(aes(colour=grp))
      } else  {
         p$mapping$colour <- p$data$dens
         p <- p + geom_line(aes(colour=dens))
         if (faceting == "vert") {
            p <- p + facet_grid(grp ~ .)
         } else p <- p + facet_grid(. ~ grp)
      }
   } else p <- p + geom_line(aes(colour=dens))
   if (plotidxs) {
      # p$data$rownum<-rownames(xsub)
      # p$data$rownum<-rws
      rownum <- rws
      p$data <- cbind(p$data,rownum)
      p <- p + geom_text(aes(label=rownum))
   }
#    if (length(rws) <= 1) {
#       print("plot consists of a single line")
#       print("this may crash some graphics systems")
#       print("graph object returned in list for safety")
#       return(list(p))
#    }
   p
}

findlocmax <- function(x,dispcols,klm) {
   dta <- x[,dispcols]
   knninst <- get.knn(dta,k=klm)
   # row i of densmatrix gives the density values for dta[i,]
   densmatrix <- matrix(x$dens[knninst$nn.index],ncol=klm)
   # for each dta[i,], find max density among neighbors
   maxdens<-apply(densmatrix,1,max)
   # return those i for which density of dta[i,] is max in its
   # neighborhood
   tmp <- which(x$dens >= maxdens)
   tmp
}

# finds the indices within z[idxs] of the elements that have
# the m highest values; indices are respect to z overall, not respect to
# z[idxs]
findtop <- function(z,idxs,m) {
   n <- length(z)
   notidxs <- setdiff(1:n,idxs)
   tmp <- z
   mv <- min(z[idxs])
   tmp[notidxs] <- mv - 1
   order(tmp)[(n-m+1):n]
}

