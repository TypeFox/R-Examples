
#' Create plot of cluster centroids.
#' 
#' Visualize centroids produced by clustering function like k-means.
#' Plots available are line plot, bar plot, or heatmap. Parameter \code{format}
#' specifies which one to create.
#' 
#' @param km an object of class \code{"toakmeans"} returned by \code{\link{computeKmeans}}.
#' @param format type of plot to use: \code{"line"}, \code{"bar"}, \code{"bar_dodge"}, 
#'   \code{"bar_facet"} (same as \code{"bar"}) or \code{"heatmap"}.
#' @param groupByCluster logical: indicates if centroids are grouped by clusters or variables. \code{groupByCluster} 
#'   has no effect when \code{format="heatmap"}.
#' @param baseSize \code{\link{theme}} base font size.
#' @param baseFamily \code{\link{theme}} base font family.
#' @param title plot title.
#' @param xlab a label for the x axis, defaults to a description of x.
#' @param ylab a label for the y axis, defaults to a description of y.
#' @param legendPosition the position of legends. ("left", "right", "bottom", "top", or two-element numeric 
#'   vector). "none" is no legend.
#' @param coordFlip logical flipped cartesian coordinates so that horizontal becomes vertical, and vertical horizontal (see 
#'   \link{coord_flip}).
#' @param ticks \code{logical} Show axis ticks using default theme settings (see \code{defaultTheme})? 
#' @param defaultTheme plot theme settings with default value \code{\link[ggthemes]{theme_tufte}}. More themes
#'   are available here: \code{\link[ggplot2]{ggtheme}} (by \href{http://ggplot2.org/}{ggplot2}) 
#'   and \code{\link[ggthemes]{ggthemes}}.
#' @param themeExtra any additional \code{\link[ggplot2]{theme}} settings that override default theme.
#' 
#' @return ggplot object
#' @export
#' @examples 
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#'                          
#' km = computeKmeans(conn, "batting", centers=5, iterMax = 25,
#'                    aggregates = c("COUNT(*) cnt", "AVG(g) avg_g", "AVG(r) avg_r", "AVG(h) avg_h"),
#'                    id="playerid || '-' || stint || '-' || teamid || '-' || yearid", 
#'                    include=c('g','r','h'), scaledTableName='kmeans_test_scaled', 
#'                    centroidTableName='kmeans_test_centroids',
#'                    where="yearid > 2000")
#' createCentroidPlot(km)
#' createCentroidPlot(km, format="bar_dodge")
#' createCentroidPlot(km, format="heatmap", coordFlip=TRUE)
#' }
createCentroidPlot <- function(km, format='line', groupByCluster=TRUE, 
                               baseSize = 12, baseFamily = "serif",
                               title = paste("Cluster Centroids", format, "Plot"), 
                               xlab, ylab = ifelse(format=="heatmap", "cluster", ifelse(km$scale, "scaled value", "value")), 
                               legendPosition = ifelse(format=="bar", "none", "right"),
                               coordFlip = FALSE, ticks = FALSE,
                               defaultTheme=theme_tufte(base_size = baseSize, base_family = baseFamily, ticks=ticks),
                                     themeExtra = NULL) {
  
  # match argument values
  format = match.arg(format, c('line', 'bar', 'heatmap','bar_dodge'))
  
  if (missing(km) || !is.object(km) || !inherits(km, "toakmeans")) {
    stop("Kmeans object must be specified.")
  }
  
  if(is.null(km$centers))
    stop("Kmeans object is missing cluster centers.")
  
  clusterid = "clusterid"
  
  centroids = data.frame(km$centers, stringsAsFactors = TRUE)
  centroids[, clusterid] = factor(rownames(km$centers))
  data = melt(centroids,id.vars=clusterid)
  
  if(groupByCluster) {
    x = "variable"
    group = clusterid
    if(missing(xlab) && format!='bar_dodge') xlab = "variable" else xlab = "cluster"
  }else {
    x = clusterid
    group = "variable"
    if(missing(xlab) && format!='bar_dodge') xlab = "cluster" else xlab = "variable"
  }
  
  if (format=='line') {
    p = plotLineCentroids(data, x, group)
  }else if (format=='bar') {
    p = plotBarCentroids(data, x, group)
  }else if (format=='bar_dodge') {
    p = plotBarDodgeCentroids(data, x, group)
  }else {
    p = plotHeatmapCentroids(data, clusterid)
  }
  
  border_element = if(format=='bar') element_rect(fill=NA) else element_blank()
  
  p = p +
    labs(title=title, x=xlab, y=ylab) +
    defaultTheme + 
    theme(legend.position=legendPosition,
          panel.border = border_element) +
    themeExtra
  
  if (format!='heatmap')
    p = p + scale_y_continuous(labels=scales::comma)
  
  if (coordFlip)
    p = p + coord_flip()
  
  return(p)
}


# Barplot with facets
plotBarCentroids <- function(data, x, group) {
  
  facet_formula = stats::as.formula(paste("~", group))
  
  ggplot(data) +
    geom_bar(aes_string(x, "value", fill=group), stat="identity", position="dodge") +
    if (group == "clusterid")
      facet_wrap(facet_formula, scales="fixed", dir="h", labeller=labeller(.default=cluster_labeller))
    else
      facet_wrap(facet_formula, scales="fixed", dir="h", labeller=labeller(.default=agg_labeller))
}

cluster_labeller <- function(value) {
  paste("Cluster", value)
}

# Barplot dodged
plotBarDodgeCentroids <- function(data, x, group) {
  
  ggplot(data) +
    geom_bar(aes_string(group, "value", fill=x), 
             stat="identity", position="dodge", color="black") 
}

# Lineplot
plotLineCentroids <- function(data, x, group) {
  
  ggplot(data=data, aes_string(x, "value", color=group)) +
    geom_line(aes_string(group=group)) +
    geom_point(size=3)
                
}

# Heatmp
plotHeatmapCentroids <- function(data, id) {
  
  ggplot(data) +
    geom_tile(aes_string("variable", id, fill="value")) +
    scale_fill_gradient2()
  
}


#' Create clusters' properties plot.
#' 
#' @param km an object of class \code{"toakmeans"} returned by \code{\link{computeKmeans}}.
#' @param baseSize \code{\link{theme}} base font size.
#' @param baseFamily \code{\link{theme}} base font family.
#' @param title plot title.
#' @param xlab a label for the x axis, defaults to a description of x.
#' @param ylab a label for the y axis, defaults to a description of y.
#' @param border boolean indicates to use border around plotting area. In case of facets border is around each facet.
#' @param colorByCluster logical: color corresponds to clusters or properties.
#' @param ticks \code{logical} Show axis ticks using default theme settings (see \code{defaultTheme})?
#' @param defaultTheme plot theme settings with default value \code{\link[ggthemes]{theme_tufte}}. More themes
#'   are available here: \code{\link[ggplot2]{ggtheme}} (by \href{http://ggplot2.org/}{ggplot2}) 
#'   and \code{\link[ggthemes]{ggthemes}}.
#' @param themeExtra any additional \code{\link[ggplot2]{theme}} settings that override default theme.
#' 
#' @return ggplot object
#' @seealso \code{\link{computeKmeans}}
#' @export
#' @examples 
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#'                          
#' km = computeKmeans(conn, "batting", centers=5, iterMax = 25,
#'                    aggregates = c("COUNT(*) cnt", "AVG(g) avg_g", "AVG(r) avg_r", "AVG(h) avg_h"),
#'                    id="playerid || '-' || stint || '-' || teamid || '-' || yearid", 
#'                    include=c('g','r','h'), scaledTableName='kmeans_test_scaled', 
#'                    centroidTableName='kmeans_test_centroids',
#'                    where="yearid > 2000")
#' createClusterPlot(km)
#' }
createClusterPlot <- function(km, baseSize = 12, baseFamily = "serif",
                              title = paste("Cluster Properties Plot"), xlab = "cluster", ylab = "value", 
                              border=TRUE, colorByCluster=TRUE, ticks=FALSE,
                              defaultTheme=theme_tufte(base_size = baseSize, base_family = baseFamily, ticks=ticks),
                              themeExtra = NULL) {
  
  if (missing(km) || !is.object(km) || !inherits(km, "toakmeans")) {
    stop("Kmeans object must be specified.")
  }
  
  if(is.null(km$aggregates))
    stop("Kmeans object is missing cluster aggregates.")
  
  clusterid="clusterid"
  aggregates = km$aggregates
  
  if (!is.factor(aggregates$clusterid)) 
    aggregates$clusterid = factor(aggregates$clusterid)
  
  if (all(c("cnt","withinss") %in% names(aggregates))) {
    aggregates$unit_withinss = aggregates$withinss / aggregates$cnt
  }
  
  data = melt(aggregates,id.vars=clusterid)
  
  facet_formula = stats::as.formula(paste("~", "variable"))
  border_element = if(border) element_rect(fill=NA) else element_blank()
  fill = ifelse(colorByCluster, clusterid, "variable")
  
  p = ggplot(data) +
    geom_bar(aes_string(clusterid, "value", fill=fill), stat="identity", position="dodge") +
    facet_wrap(facet_formula, scales="free", dir="h", labeller=labeller(.default=agg_labeller)) +
    labs(title=title, x=xlab, y=ylab) +
    defaultTheme + 
    theme(legend.position="none",
          panel.border = border_element) +
    themeExtra
  
  return(p)
}

agg_labeller <- function(value) {
  paste("Property", value)
}


#' Create cluster variable plot.
#' 
#' @param km an object of class \code{"toakmeans"} returned by \code{\link{computeKmeans}}.
#' @param baseSize \code{\link{theme}} base font size.
#' @param baseFamily \code{\link{theme}} base font family.
#' @param title plot title.
#' @param ticks \code{logical} Show axis ticks using default theme settings (see \code{defaultTheme})?
#' @param defaultTheme plot theme settings with default value \code{\link[ggthemes]{theme_tufte}}. More themes
#'   are available here: \code{\link[ggplot2]{ggtheme}} (by \href{http://ggplot2.org/}{ggplot2}) 
#'   and \code{\link[ggthemes]{ggthemes}}.
#' @param themeExtra any additional \code{\link[ggplot2]{theme}} settings that override default theme.
#' @param ... other parameters being suplied to geom's \code{aes}.
#' 
#' @return ggplot object
#' @export
#' @examples 
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#'                          
#' km = computeKmeans(conn, "batting", centers=5, iterMax = 25,
#'                    aggregates = c("COUNT(*) cnt", "AVG(g) avg_g", "AVG(r) avg_r", "AVG(h) avg_h"),
#'                    id="playerid || '-' || stint || '-' || teamid || '-' || yearid", 
#'                    include=c('g','r','h'), scaledTableName='kmeans_test_scaled', 
#'                    centroidTableName='kmeans_test_centroids',
#'                    where="yearid > 2000")
#' km = computeClusterSample(conn, km, 0.01)
#' createClusterPairsPlot(km, title="Batters Clustered by G, H, R", ticks=FALSE)
#' }
createClusterPairsPlot <- function(km, baseSize = 12, baseFamily = "serif",
                                   title="Cluster Variable Pairs", ticks=FALSE,
                                   defaultTheme=theme_tufte(base_size = baseSize, base_family = baseFamily, ticks = ticks),
                                   themeExtra = theme(), ...) {
  
  if (missing(km) || !is.object(km) || !inherits(km, "toakmeans")) {
    stop("Kmeans object must be specified.")
  }
  
  if(is.null(km$data))
    stop("Kmeans object is missing sample data.")
  
  kms = km$data
  
  if (!is.factor(kms$clusterid)) 
    kms$clusterid = factor(kms$clusterid)
  
  p = GGally::ggpairs(kms, color='clusterid', title=title, ...) +
    defaultTheme +
    themeExtra
  
  return(p)
}


#' Create cluster silhouette profile plot.
#' 
#' @param km an object of class \code{"toakmeans"} returned by \code{\link{computeKmeans}}.
#' @param baseSize \code{\link{theme}} base font size.
#' @param baseFamily \code{\link{theme}} base font family.
#' @param title plot title.
#' @param xlab a label for the x axis, defaults to a description of x.
#' @param ylab a label for the y axis, defaults to a description of y.
#' @param coordFlip logical flipped cartesian coordinates so that horizontal becomes vertical, and vertical horizontal (see 
#'   \link{coord_flip}).
#' @param ticks \code{logical} Show axis ticks using default theme settings (see \code{defaultTheme})?
#' @param defaultTheme plot theme settings with default value \code{\link[ggthemes]{theme_tufte}}. More themes
#'   are available here: \code{\link[ggplot2]{ggtheme}} (by \href{http://ggplot2.org/}{ggplot2}) 
#'   and \code{\link[ggthemes]{ggthemes}}.
#' @param themeExtra any additional \code{\link[ggplot2]{theme}} settings that override default theme.
#'  
#' @return ggplot object
#' @export
#' @examples 
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#'                          
#' km = computeKmeans(conn, "batting", centers=5, iterMax = 25,
#'                    aggregates = c("COUNT(*) cnt", "AVG(g) avg_g", "AVG(r) avg_r", "AVG(h) avg_h"),
#'                    id="playerid || '-' || stint || '-' || teamid || '-' || yearid", 
#'                    include=c('g','r','h'), scaledTableName='kmeans_test_scaled', 
#'                    centroidTableName='kmeans_test_centroids',
#'                    where="yearid > 2000")
#' km = computeSilhouette(conn, km)
#' createSilhouetteProfile(km, title="Cluster Silhouette Histograms (Profiles)")
#' }
createSilhouetteProfile <- function(km, baseSize = 12, baseFamily = "serif",
                                   title="Cluster Silhouette Profile (Histogram)", xlab="Silhouette Value", ylab="Count",
                                   coordFlip = TRUE, ticks=FALSE,
                                   defaultTheme=theme_tufte(base_size = baseSize, base_family = baseFamily, ticks = ticks),
                                   themeExtra = NULL) {
  
  if (missing(km) || !is.object(km) || !inherits(km, "toakmeans")) {
    stop("Kmeans object must be specified.")
  }
  
  if(is.null(km$sil) || is.null(km$sil$profile))
    stop("Kmeans object is missing silhouette data.")
  
  silprofile = km$sil$profile
  
  if (!is.factor(silprofile$clusterid)) 
    silprofile$clusterid = factor(silprofile$clusterid)
  
  silprofile = silprofile[silprofile$bin_percent != 0, ]
  p = ggplot(silprofile) +
    geom_bar(aes_string("bin_start","bin_percent",group="clusterid",fill="clusterid"), stat="identity", position="dodge") +
    facet_wrap(~clusterid, ncol=1) +
    labs(title=title, x=xlab, y=ylab) +
    defaultTheme +
    themeExtra +
    theme(legend.position="none")
  
  if (coordFlip)
    p = p + coord_flip()
  
  return(p)
}