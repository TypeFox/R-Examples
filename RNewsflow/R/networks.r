#' Create a document similarity network
#'
#' Combines document similarity data (d) with document meta data (meta) into an \link[igraph]{igraph} network/graph.
#' 
#' This function is mainly offered to mimic the output of the \link[RNewsflow]{newsflow.compare} function when using imported document similarity data.
#' This way the functions for transforming, aggregating and visualizing the document similarity data can be used.
#'
#' @param d A data.frame with three columns, that represents an edgelist with weight values. 
#' The first and second column represent the names/ids of the 'from' and 'to' documents/vertices. 
#' The third column represents the similarity score. Column names are ignored
#' @param meta A data.frame where rows are documents and columns are document meta information. 
#' Should at least contain 2 columns: the document name/id and date. 
#' The name/id column should match the document names/ids of the edgelist, and its label is specified in the `id.var` argument. 
#' The date column should be intepretable with \link[base]{as.POSIXct}, and its label is specified in the `date.var` argument.            
#' @param id.var The label for the document name/id column in the `meta` data.frame. Default is "document_id"
#' @param date.var The label for the document date column in the `meta` data.frame . default is "date"
#' @param min.similarity For convenience, ignore all edges where the weight is below `min.similarity`. Default is zero.
#'
#' @return A network/graph in the \link[igraph]{igraph} class
#' @export
#'
#' @examples
#' d = data.frame(x = c(1,1,1,2,2,3),
#'                y = c(2,3,5,4,5,6),
#'                v = c(0.3,0.4,0.7,0.5,0.2,0.9))
#' 
#' meta = data.frame(document_id = 1:8,
#'                   date = seq.POSIXt(from = as.POSIXct('2010-01-01 12:00:00'), 
#'                          by='hour', length.out = 8),
#'                   medium = c(rep('Newspapers', 4), rep('Blog', 4)))
#' 
#' g = document.network(d, meta)
#' 
#' igraph::get.data.frame(g, 'both')
#' igraph::plot.igraph(g)
document.network <- function(d, meta, id.var='document_id', date.var='date', min.similarity=0){
  confirm.dtm.meta(meta, id.var, date.var)
  
  colnames(d) = c('x','y','similarity')
  d = d[d$similarity >= min.similarity, c('x','y','similarity')]
  d = d[order(-d$similarity),]
  
  g = igraph::graph.data.frame(d[,c('x','y')])
  igraph::E(g)$weight = d$similarity
  
  if(mean(igraph::V(g)$name %in% meta[,id.var]) < 1) stop("Not all documents in d match with an 'id' in the meta information")
  
  ## add documents in meta data.frame that do not appear in the edgelist (in other words, isolates)
  missingmeta = as.character(meta[!meta[,id.var] %in% igraph::V(g)$name, id.var])
  g = igraph::add.vertices(g, nv = length(missingmeta), attr = list(name=missingmeta))

  meta = meta[match(igraph::V(g)$name, meta[,id.var]),]
  attribs = colnames(meta)[!colnames(meta) == id.var]
  for(attrib in attribs){
    g = igraph::set.vertex.attribute(g, attrib, value=as.character(meta[,attrib]))
  }
  
  edgelist = igraph::get.edges(g, igraph::E(g))
  dates = as.POSIXct(igraph::V(g)$date) 
  igraph::E(g)$hourdiff = round(difftime(dates[edgelist[,2]], dates[edgelist[,1]], units = 'hours'),3)
  
  g
}

#' Visualize (a subcomponent) of the document similarity network 
#'
#' @param g A document similarity network, as created with \link[RNewsflow]{newsflow.compare} or \link[RNewsflow]{document.network}
#' @param date.attribute The label of the vertex/document date attribute. Default is "date"
#' @param source.attribute The label of the vertex/document source attribute. Default is "source"
#' @param subcomp_i Optional. If an integer is given, the network is decomposed into subcomponents (i.e. unconnected components) and only the i-th component is visualized.
#' @param dtm Optional. If a document-term matrix that contains the documents in g is given, a wordcloud with the most common words of the network is added. 
#' @param sources Optional. Use a character vector to select only certain sources
#' @param only.outer.date If TRUE, only the labels for the first and last date are reported on the x-axis
#' @param date.format The date format of the date labels (see \link[base]{format.POSIXct})
#' @param margins The margins of the network plot. The four values represent bottom, left, top and right margin. 
#' @param isolate.color Optional. Set a custom color for isolates
#' @param source.loops If set to FALSE, all edges between vertices/documents of the same source are ignored.
#' @param ... Additional arguments for the network plotting function \link[igraph]{plot.igraph} 
#'
#' @return Nothing. 
#' @export
#'
#' @examples
#' data(docnet)
#' data(dtm)
#' 
#' docnet_comps = igraph::decompose.graph(docnet) # get subcomponents
#' 
#' # subcomponent 1
#' document.network.plot(docnet_comps[[1]]) 
#' 
#' # subcomponent 2 with wordcloud
#' document.network.plot(docnet_comps[[2]], dtm=dtm) 
#' 
#' # subcomponent 3 with additional arguments passed to plot.igraph 
#' document.network.plot(docnet_comps[[3]], dtm=dtm, vertex.color='red') 
document.network.plot <- function(g, date.attribute='date', source.attribute='source', subcomp_i=NULL, dtm=NULL, sources=NULL, only.outer.date=FALSE, date.format='%Y-%m-%d %H:%M', margins=c(5,8,1,13), isolate.color=NULL, source.loops=TRUE, ...){
  g = igraph::set.vertex.attribute(g, 'date', value= igraph::get.vertex.attribute(g, date.attribute))
  g = igraph::set.vertex.attribute(g, 'source', value= igraph::get.vertex.attribute(g, source.attribute))
  
  if(is.null(subcomp_i)){
    cluster = g
  } else {
    clusters = igraph::decompose.graph(g)
    cluster = clusters[[subcomp_i]]
  }
  
  if(is.null(sources)) sources = sort(unique(igraph::V(cluster)$source))
  cluster = igraph::delete.vertices(cluster, which(!igraph::V(cluster)$source %in% sources))
  
  if(!is.null(dtm)){
    graphics::layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE), widths = c(1.5, 2.5), heights = c(1, 2))  
    topwords = slam::col_sums(dtm[rownames(dtm) %in% igraph::V(cluster)$name,])
    topwords = topwords[topwords > 0, drop=FALSE]
    topwords = head(topwords[order(-topwords)], 10)
    graphics::par(mar = c(0,0,0,0))
    suppressWarnings(
    wordcloud::wordcloud(names(topwords), sqrt(as.numeric(topwords)), scale=c(2,0.7))
    )
  } 
  
  plotDocNet(cluster, sources, only.outer.date, date.format, margins, isolate.color=NULL, source.loops=TRUE, ...)
  if(!is.null(dtm)) graphics::layout(matrix(c(1, 1, 1, 1), 2, 2, byrow = TRUE))
}

plotDocNet <- function(cluster, sources=NULL, only.outer.date=FALSE, date.format='%Y-%m-%d %H:%M', margins=c(5,8,1,13), isolate.color=NULL, source.loops=TRUE, ...){
  graphics::par(mar = margins)
  if(is.null(sources)) sources = sort(unique(igraph::V(cluster)$source))
  igraph::E(cluster)$width = scales::rescale(igraph::E(cluster)$weight, c(2,4), from = c(0,1))
  igraph::E(cluster)$arrow.size = 0.45
  igraph::E(cluster)$color = grDevices::grey(scales::rescale(igraph::E(cluster)$weight, c(0.8,0.2), from = c(0,1)))
  
  if(!source.loops){
    edgelist = igraph::get.edgelist(cluster)
    xmatch = match(edgelist[,1], igraph::V(cluster)$name)
    ymatch = match(edgelist[,2], igraph::V(cluster)$name)
    sourceloop = which(igraph::V(cluster)$source[xmatch] == igraph::V(cluster)$source[ymatch])
    cluster = igraph::delete.edges(cluster, sourceloop)
  }
  
  igraph::V(cluster)$x = as.POSIXct(igraph::V(cluster)$date)
  igraph::V(cluster)$y = length(sources) + 1 - match(igraph::V(cluster)$source, sources)
  
  vertex.atts = names(igraph::vertex.attributes(cluster))
  if(!'size' %in% vertex.atts) igraph::V(cluster)$size = 7.5
  if(!'color' %in% vertex.atts) igraph::V(cluster)$color = 'lightgrey'   
  if(!'label' %in% vertex.atts) igraph::V(cluster)$label = ''
  if(!'label.color' %in% vertex.atts) igraph::V(cluster)$label.color = 'black'
  
  if(!is.null(isolate.color)) igraph::V(cluster)$color[igraph::degree(cluster) == 0] = isolate.color
  
  ymin = 1 - max(igraph::V(cluster)$y)*(1/max(igraph::V(cluster)$y))
  cluster$layout = igraph::layout.norm(cbind(as.numeric(igraph::V(cluster)$x), igraph::V(cluster)$y), -1, 1, ymin, 1)
  
  if(length(unique(igraph::V(cluster)$x)) == 1) cluster$layout[,1] = 0
  igraph::plot.igraph(cluster, rescale=FALSE, asp=0, ...)

  for(source in sources){
    ycoord = cluster$layout[igraph::V(cluster)$source == source,2][1]
    graphics::text(1.1, ycoord, source, adj=0)
  }
  
  ## determine inner vertical lines (to indicate date)
  inner_date = seq.POSIXt(as.POSIXct(min(igraph::V(cluster)$date)), as.POSIXct(max(igraph::V(cluster)$date)), by='day')
  inner_date = as.POSIXct(as.Date(inner_date))
  inner_date = inner_date[as.numeric(inner_date) > min(igraph::V(cluster)$x) & as.numeric(inner_date) < max(igraph::V(cluster)$x)]
  inner_date_i = scales::rescale(as.numeric(inner_date), to=c(-1,1), from=as.numeric(c(min(igraph::V(cluster)$x), max(igraph::V(cluster)$x))))
  
  graphics::clip(-1.05,1,ymin,1.05)
  graphics::abline(h=unique(cluster$layout[,2]), col="grey10", lty='dotted')
  
  ## determine outer vertical lines
  outer_date_i = c(-1,1)
  outer_date = c(min(igraph::V(cluster)$date), max(igraph::V(cluster)$date))
  if(length(unique(outer_date)) == 1){
    outer_date_i = 0
    outer_date = outer_date[1]
  }
  ## draw lines, add text
  graphics::abline(v=outer_date_i, col="grey10", lty='dotted')
  graphics::abline(v=inner_date_i, col="grey10", lty='dotted')
  
  graphics::clip(-2,2,-100,2)
  outer_date_print = format(as.POSIXct(outer_date), date.format)
  if(!only.outer.date & length(inner_date) > 0){
    graphics::text(x=outer_date_i, y=ymin-0.2, label=outer_date_print, srt=30, adj=1) 
    inner_date_print = format(as.POSIXct(inner_date), date.format)
    graphics::text(x=inner_date_i, y=ymin-0.2, label=inner_date_print, srt=30, adj=1) 
  } else graphics::text(x=outer_date_i, y=ymin-0.2, label=outer_date_print, srt=30, adj=1) 
}

#' Show time window of document pairs 
#'
#' This function aggregates the edges for all combinations of attributes specified in `from.attribute` and `to.attribute`, and shows the minimum and maximum hour difference for each combination.
#'
#' The \link[RNewsflow]{filter.window} function can be used to filter edges that fall outside of the intended time window. 
#'
#' @param g A document similarity network, as created with \link[RNewsflow]{newsflow.compare} or \link[RNewsflow]{document.network}  
#' @param to.attribute The vertex attribute to aggregate the `to` group of the edges
#' @param from.attribute The vertex attribute to aggregate the `from` group of the edges
#'
#' @return A data.frame showing the left and right edges of the window for each unique group.
#' @export
#'
#' @examples
#' data(docnet)
#' show.window(docnet, to.attribute = 'source')
#' show.window(docnet, to.attribute = 'sourcetype')
#' show.window(docnet, to.attribute = 'sourcetype', from.attribute = 'sourcetype')
show.window <- function(g, to.attribute=NULL, from.attribute=NULL){
  from.vertices = if(is.null(from.attribute)) rep('any document', igraph::vcount(g)) else igraph::get.vertex.attribute(g, from.attribute)
  to.vertices = if(is.null(to.attribute)) rep('any document', igraph::vcount(g)) else igraph::get.vertex.attribute(g, to.attribute)
  
  e = data.frame(igraph::get.edges(g, igraph::E(g)))
  medwindow = stats::aggregate(igraph::E(g)$hourdiff, by=list(from=from.vertices[e$X1], to=to.vertices[e$X2]), FUN = function(x) cbind(min(x), max(x)))
  d = data.frame(from = medwindow$from, to=medwindow$to, window.left=medwindow$x[,1], window.right=medwindow$x[,2])
  d
}


#' Filter edges from the document similarity network based on hour difference
#'
#' The `filter.window` function can be used to filter the document pairs (i.e. edges) using the `hour.window` parameter, which works identical to the `hour.window` parameter in the `newsflow.compare` function. 
#' In addition, the `from.vertices` and `to.vertices` parameters can be used to select the vertices (i.e. documents) for which this filter is applied.
#'
#' It is recommended to use the \link[RNewsflow]{show.window} function to verify whether the hour windows are correct according to the assumptions and focus of the study.
#'
#' @param g A document similarity network, as created with \link[RNewsflow]{newsflow.compare} or \link[RNewsflow]{document.network}  
#' @param hour.window A vector of length 2, in which the first and second value determine the left and right side of the window, respectively. For example, c(-10, 36) will compare each document to all documents between the previous 10 and the next 36 hours.
#' @param to.vertices A filter to select the vertices `to` which an edge is filtered. 
#' For example, if `V(g)$sourcetype == "newspaper"` is used, then the hour.window filter is only applied for edges `to` newspaper documents (specifically, where the sourcetype attribute is "newspaper"). 
#' @param from.vertices A filter to select the vertices `from` which an edge is filtered. 
#' Works identical to `to.vertices`.  
#'
#' @return A network/graph in the \link[igraph]{igraph} class
#' @export
#'
#' @examples
#' data(docnet)
#' show.window(docnet, to.attribute = 'source') # before filtering
#' 
#' docnet = filter.window(docnet, hour.window = c(0.1,24))
#' 
#' docnet = filter.window(docnet, hour.window = c(6,36), 
#'                        to.vertices = V(docnet)$sourcetype == 'Print NP')
#' 
#' show.window(docnet, to.attribute = 'sourcetype') # after filtering per sourcetype
#' show.window(docnet, to.attribute = 'source') # after filtering per source
filter.window <- function(g, hour.window, to.vertices=NULL, from.vertices=NULL){
  if(is.null(from.vertices)) from.vertices = rep(TRUE, igraph::vcount(g))
  if(is.null(to.vertices)) to.vertices = rep(TRUE, igraph::vcount(g))
  
  # get vector indicating which edges are selected given the vertex.filter
  from.edge.i = igraph::get.edges(g, igraph::E(g))[,1]
  to.edge.i = igraph::get.edges(g, igraph::E(g))[,2]
  selected_edge = from.vertices[from.edge.i] & to.vertices[to.edge.i]
  
  # filter the selected edges
  delete.i = selected_edge & (igraph::E(g)$hourdiff < hour.window[1] | igraph::E(g)$hourdiff > hour.window[2])
  igraph::delete.edges(g, which(delete.i))
}

########

#' Transform document network so that each document only matches the earliest dated matching document
#'
#' Transforms the network so that a document only has an edge to the earliest dated document it matches within the specified time window[^duplicate].
#' 
#' If there are multiple earliest dated documents (that is, having the same publication date) then edges to all earliest dated documents are kept.
#' 
#' @param g A document similarity network, as created with \link[RNewsflow]{newsflow.compare} or \link[RNewsflow]{document.network}  
#'
#' @return A network/graph in the \link[igraph]{igraph} class
#' @export
#'
#' @examples
#' data(docnet)
#' 
#' subcomp1 = igraph::decompose.graph(docnet)[[2]]
#' subcomp2 = only.first.match(subcomp1)
#' 
#' igraph::get.data.frame(subcomp1)
#' igraph::get.data.frame(subcomp2)
#' 
#' graphics::par(mfrow=c(2,1))
#' document.network.plot(subcomp1, main='All matches')
#' document.network.plot(subcomp2, main='Only first match')
#' graphics::par(mfrow=c(1,1))
only.first.match <- function(g){
  e = data.frame(igraph::get.edges(g, igraph::E(g)))
  
  # first make from and to columns of same temporal order (from before to) by switching columns with negative hourdiff.
  switch_from.to = igraph::get.data.frame(g, 'edges')$hourdiff < 0
  e = data.frame(from = ifelse(switch_from.to, e$X2, e$X1),
                 to   = ifelse(switch_from.to, e$X1, e$X2))
  
  e$i = 1:nrow(e)
  e$from.date = igraph::V(g)$date[e$from]
  e = e[order(e$from.date),]
  
  first = !duplicated(e$to)
  first = paste(e$to, e$from.date, sep='.') %in% paste(e$to[first], e$from.date[first], sep='.') # also first if same date as first
  
  igraph::delete.edges(g, e$i[!first])
}




#' Aggregate the edges of a network by vertex attributes 
#' 
#' This function offers a versatile way to aggregate the edges of a network based on the vertex attributes.
#' Although it was designed specifically for document similarity networks, it can be used for any network in the \link[igraph]{igraph} class.
#' 
#' The first argument is the network (in the `igraph` class). 
#' The second argument, for the `by` parameter, is a character vector to indicate one or more vertex attributes based on which the edges are aggregated.
#' Optionally, the `by` parameter can also be specified separately for `by.from` and `by.to`. 
#' 
#' By default, the function returns the aggregated network as an \link[igraph]{igraph} class.
#' The edges in the aggregated network have five standard attributes. 
#' The `edges` attribute counts the number of edges from the `from` group to the `to` group. 
#' The `from.V` attribute shows the number of vertices in the `from` group that matched with a vertex in the `to` group.
#' The `from.Vprop attribute shows this as the proportion of all vertices in the `from` group.
#' The `to.V` and `to.Vprop` attributes show the same for the `to` group.
#' 
#' In addition, one of the edge attributes of the original network can be aggregated with a given function.
#' These are specified in the `edge.attribute` and `agg.FUN` parameters.
#'
#' @param g A network/graph in the \link[igraph]{igraph} class
#' @param by A character string indicating the vertex attributes by which the edges will be aggregated. 
#' @param by.from Optionally, specify different vertex attributes to aggregate the `from` side of edges 
#' @param by.to Optionally, specify different vertex attributes to aggregate the `to` side of edges 
#' @param edge.attribute Select an edge attribute to aggregate using the function specified in `agg.FUN`. Defaults to 'weight'
#' @param agg.FUN The function used to aggregate the edge attribute
#' @param return.df Optional. If TRUE, the results are returned as a data.frame. This can in particular be convenient if by.from and by.to are used.
#'
#' @return A network/graph in the \link[igraph]{igraph} class, or a data.frame if return.df is TRUE.
#' @export
#'
#' @import data.table
#'
#' @examples
#' data(docnet)
#' aggdocnet = network.aggregate(docnet, by='sourcetype')
#' igraph::get.data.frame(aggdocnet, 'both')
#' 
#' aggdocdf = network.aggregate(docnet, by.from='sourcetype', by.to='source', return.df = TRUE)
#' head(aggdocdf)
network.aggregate <- function(g, by=NULL, by.from=by, by.to=by, edge.attribute='weight', agg.FUN=mean, return.df=FALSE){
  igraph::V(g)$any_vertex = '---'
  if(is.null(by.from)) by.from = 'any_vertex'
  if(is.null(by.to)) by.to = 'any_vertex'
  
  e = data.frame(igraph::get.edges(g, igraph::E(g)))
  v = igraph::get.data.frame(g, 'vertices')  
  
  e = cbind(e, v[e$X1, by.from, drop=FALSE])
  e = cbind(e, v[e$X2, by.to, drop=FALSE])
  colnames(e) = c('from','to', paste('from', by.from, sep='.'), paste('to', by.to, sep='.'))
  aggvars = colnames(e)[3:ncol(e)]
  
  e$from_unique_match = as.numeric(!duplicated(e[,colnames(e)[-2]]))
  e$to_unique_match = as.numeric(!duplicated(e[,colnames(e)[-1]]))
  e$var = igraph::get.edge.attribute(g, edge.attribute)
  
  e = data.table::setDT(e)
  e = e[, list(edges=.N,
               edge.attribute=agg.FUN(get('var')), # note that get('var') is used because the CRAN check would otherwise wrongfully interpret 'var' as a global variable (yielding a NOTE)  
               from.V=sum(get('from_unique_match')),
               to.V=sum(get('to_unique_match'))),
               by= eval(paste(aggvars, collapse=','))]
  e = as.data.frame(e)
  

  
  colnames(e)[colnames(e) == 'edge.attribute'] = paste('agg',edge.attribute, sep='.')
  
  # match total vertices
  v = data.table::setDT(v)
  from.totals = v[, list(from.N=.N), by= eval(paste(by.from, collapse=','))]
  to.totals = v[, list(to.N=.N), by= eval(paste(by.to, collapse=','))]
  
  e = merge(e, from.totals, by.x=paste('from', by.from, sep='.'), by.y=by.from, all.x=TRUE)
  e = merge(e, to.totals, by.x=paste('to', by.to, sep='.'), by.y=by.to, all.x=TRUE)
  e$from.Vprop = e$from.V / e$from.N
  e$to.Vprop = e$to.V / e$to.N
  
  if(!return.df) {
    e = return.network.aggregate(e)
    if('any_vertex' %in% names(igraph::vertex.attributes(e))) e = igraph::delete_vertex_attr(e, 'any_vertex')
    return(e)
  }
  if(return.df) {
    e = e[,!colnames(e) %in% c('from.N', 'to.N', 'from.any_vertex', 'to.any_vertex')]
    
    edge_col = which(colnames(e) == 'edges') # silly but necessary way to order the from and to attribute columns ()
    ord = c(order(colnames(e)[1:(edge_col-1)]), edge_col:ncol(e))
    return(e[,ord]) 
  }
}

return.network.aggregate <- function(g.df){
  vnames = colnames(g.df)[1:(grep('^edges$', colnames(g.df))-1)]
  
  # create network
  from.v = g.df[,grep('^from\\.', vnames), drop=FALSE]
  to.v = g.df[,grep('^to\\.', vnames), drop=FALSE]
  from.name = if(ncol(from.v) > 1) apply(from.v, MARGIN = 1, paste, collapse='; ') else from.v[,1]
  to.name = if(ncol(to.v) > 1) apply(to.v, MARGIN = 1, paste, collapse='; ') else to.v[,1]
  g.agg = igraph::graph.data.frame(cbind(from.name,to.name))
  
  aggvar = colnames(g.df)[grep('agg\\.', colnames(g.df))]
  for(edge.attrib in c('edges',aggvar,'from.V','from.Vprop','to.V','to.Vprop')){
    g.agg = igraph::set.edge.attribute(g.agg, edge.attrib, value=g.df[,edge.attrib])
  }
  
  # create meta table
  from.v = cbind(from.name, from.v, from.N=g.df$from.N)
  to.v = cbind(to.name, to.v, to.N=g.df$to.N)
  
  colnames(from.v) = gsub('^from\\.', '', colnames(from.v))
  colnames(to.v) = gsub('^to\\.', '', colnames(to.v))
  metavars = unique(c(colnames(from.v), colnames(to.v)))
  for(metavar in metavars){
    if(!metavar %in% colnames(from.v)) from.v[,metavar] = NA
    if(!metavar %in% colnames(to.v)) to.v[,metavar] = NA
  }

  meta = unique(rbind(from.v[,metavars], to.v[,metavars]))
  
  # add meta attributes
  meta_i = match(igraph::V(g.agg)$name, meta$name)
  for(metavar in metavars[!metavars == 'name']){
    g.agg = igraph::set.vertex.attribute(g.agg, metavar, value=meta[meta_i,metavar])
  }
  g.agg
}


graph.plot.presets <- function(g){
  igraph::E(g)$curved = 0.2
  igraph::E(g)$width = igraph::E(g)$weight * 5
  igraph::E(g)$arrow.size = 0.7
  igraph::E(g)$label = round(igraph::E(g)$weight, 2)
  
  igraph::V(g)$label.color = 'black'
  igraph::V(g)$color = 'white'
  igraph::V(g)$size = 10
  g
}


#' A wrapper for \link[igraph]{plot.igraph} for visualizing directed networks.
#'
#' This is a convenience function for visualizing directed networks with edge labels using \link[igraph]{plot.igraph}. 
#' It was designed specifically for visualizing aggregated document similarity networks in the RNewsflow package, but works with any network in the \link[igraph]{igraph} class.  
#'
#' @param g A network/graph in the \link[igraph]{igraph} class 
#' @param weight.var The edge attribute that is used to specify the edges
#' @param weight.thres A threshold for weight. Edges below the threshold are ignored 
#' @param delete.isolates If TRUE, isolates (i.e. vertices without edges) are ignored.
#' @param vertex.size The size of the verticex/nodes. Defaults to 30. Can be a vector with values per vertex.
#' @param vertex.color Color of vertices/nodes. Default is lightblue. Can be a vector with values per vertex.
#' @param vertex.label.color Color of labels for vertices/nodes. Defaults to black. Can be a vector with values per vertex. 
#' @param vertex.label.cex Size of the labels for vertices/nodes. Defaults to 0.7. Can be a vector with values per vertex. 
#' @param edge.color Color of the edges. Defaults to grey. Can be a vector with values per edge.
#' @param show.edge.labels Logical. Should edge labels be displayed? Default is TRUE.
#' @param edge.label.color Color of the edge labels. Defaults to black. Can be a vector with values per edge. 
#' @param edge.label.cex Size of the edge labels. Defaults to 0.6. Can be a vector with values per edge. 
#' @param edge.arrow.size Size of the edge arrows. Defaults to 1. Can only be set globally (igraph might update this at some point)
#' @param layout The igraph layout used to plot the network. Defaults to \link[igraph]{layout.davidson.harel}
#' @param ... Arguments to be passed to the \link[igraph]{plot.igraph} function. 
#'
#' @return Nothing
#' @export
#'
#' @examples
#' data(docnet)
#' aggdocnet = network.aggregate(docnet, by='source')
#' directed.network.plot(aggdocnet, weight.var = 'to.Vprop', weight.thres = 0.2)
directed.network.plot <- function(g, weight.var='from.Vprop', weight.thres=NULL, delete.isolates=FALSE, 
                                   vertex.size=30, vertex.color='lightblue', vertex.label.color='black', vertex.label.cex=0.7, 
                                   edge.color = 'grey', show.edge.labels=TRUE, edge.label.color = 'black', edge.label.cex = 0.6, edge.arrow.size=1, 
                                   layout=igraph::layout.davidson.harel, ...){
  igraph::E(g)$weight = igraph::get.edge.attribute(g, weight.var)
  if(!is.null(weight.thres)) g = igraph::delete.edges(g, which(igraph::E(g)$weight < weight.thres))
  if(delete.isolates) g = igraph::delete.vertices(g, igraph::degree(g) == 0)
  
  argnames = names(list(...))  
  igraph::E(g)$color = edge.color
  igraph::E(g)$arrow.size = edge.arrow.size
  igraph::E(g)$label.cex = edge.label.cex
  igraph::E(g)$label.color = edge.label.color
  if(!'edge.width' %in% argnames) igraph::E(g)$width = scales::rescale(igraph::E(g)$weight, to = c(1,5))
  if(!'edge.label' %in% argnames & show.edge.labels) igraph::E(g)$label = round(igraph::E(g)$weight, 2)
  if(!'edge.label.font' %in% argnames) igraph::E(g)$label.font = 3
  if(!'edge.curved' %in% argnames) {
    e = igraph::get.edges(g, igraph::E(g))
    twoway = paste(e[,1], e[,2], sep='.') %in% paste(e[,2], e[,1], sep='.')
    igraph::E(g)$curved = ifelse(twoway, 0.35, 0.15)
  } 
  
  igraph::V(g)$size = vertex.size
  igraph::V(g)$label.cex = vertex.label.cex
  igraph::V(g)$color = vertex.color
  igraph::V(g)$label.color = vertex.label.color
  if(!'vertex.label' %in% argnames) igraph::V(g)$label = gsub(' ', '\n', igraph::V(g)$name)  
  if(!'vertex.frame.color' %in% argnames) igraph::V(g)$frame.color = igraph::V(g)$color
  
  g$layout = layout
  igraph::plot.igraph(g, ...)
}

