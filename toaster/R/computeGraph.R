#' Define an object corresponding to a graph in Aster database.
#' 
#' In Aster Database, to process graphs using SQL-GR, it is recommended to represent
#' a graph using two tables:
#' - Vertices table
#' - Edges table
#' Vertices table must contain a unique key so that each row represents a vertex.
#' Edges table must contain a pair of source and target keys (from vertices table)
#' so that each row represents an edge.
#' Both vertices and edges tables may contain additional columns representing 
#' optional attributes. For example if edges table has column 'weight' it can
#' correspond a graph with edge weights.
#' @param vertices A table, view, or query of a collection of vertices in the graph.
#' @param edges A table, view, or query of a collection of edges of the graph. 
#'   The collection must contain at least two lists of columns, one list that represents 
#'   the source vertex key and another list that represents the target vertex key.
#' @param directed logical: should edges be interpreted as directed?
#' @param key name of the column with vertex unique id (in the table \code{vertices}).
#' @param source name of the column with the from vertex (in the table \code{edges}).
#' @param target name of the column with the to vertex (in the tbale \code{edges}).
#' @param vertexAttrnames optionally, a list of columns containing vertex 
#'   attribute names.
#' @param edgeAttrnames optionally, a list of columns containing edge 
#'   attribute names.
#' @param vertexWhere optionally, a \code{SQL WHERE} clause to subset vertex table (use SQL 
#'   as if in \code{WHERE} clause but omit the keyword \code{WHERE}).
#' @param edgeWhere optionally, a \code{SQL WHERE} clause to subset edge table (use SQL 
#'   as if in \code{WHERE} clause but omit the keyword \code{WHERE}).
#' 
#' @export
#' @examples 
#' # undirected graph
#' policeGraphUn = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_un", 
#'      directed = FALSE, key = "officer", source = "officer1", target = "officer2", 
#'      vertexAttrnames = c("offense_count"), edgeAttrnames = c("weight"))
#'                          
#' # directed graph with the vertex filter
#' policeGraphDi = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_di", 
#'      directed = TRUE, key = "officer", source = "officer1", target = "officer2", 
#'      vertexAttrnames = c("offense_count"), edgeAttrnames = c("weight"),
#'      vertexWhere = "officer ~ '[A-Z].*'")
#'      
toaGraph <- function(vertices, edges, directed=FALSE, 
                     key='id', source='source', target='target', 
                     vertexAttrnames=NULL, edgeAttrnames=NULL, 
                     vertexWhere = NULL, edgeWhere = NULL) {
  
  if(is.null(vertices) || is.null(edges))
    stop("Both vertices and edges must be defined.")
  
  z <- structure(list(vertices = vertices,
                      edges = edges,
                      directed = directed,
                      key = key,
                      source = source,
                      target = target,
                      vertexAttrnames = vertexAttrnames,
                      edgeAttrnames = edgeAttrnames,
                      vertexWhere = vertexWhere,
                      edgeWhere = edgeWhere
  ),
  class = "toagraph")
 
  z 
}

#' Materialize Aster graph as network object in R.
#'
#' Results in \code{\link{network}} object representation of the graph 
#' stored in Aster tables. Usually in Aster database a graph is represented
#' using a pair of vertice and edge tables.
#' 
#' Use caution when computing network objects stored in Aster with this function 
#' as data may include considerable amount of vertices and edges which are too large to
#' load into a memory.  
#'
#'
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param graph an object of class \code{'toagraph'} referencing graph 
#'   tables in Aster database.
#' @param v a SQL \code{SELECT} that returns key values or a list of key values (corresponding 
#'   to the \code{vertex.names} attribute) of the vertices to include in the graph. 
#'   When not \code{NULL} this guarentees that no other vertices or edges between other vertices 
#'   are included in the resulting network.
#' @param vertexWhere optionally, a \code{SQL WHERE} clause to subset vertex table. When not \code{NULL}
#'   it overrides \code{vertexWhere} condition from the \code{graph}.
#' @param edgeWhere optionally, a \code{SQL WHERE} clause to subset edge table. When not \code{NULL}
#'   it overrides \code{edgeWhere} condition from the \code{graph}.
#' @param allTables pre-built information about existing tables.
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \pkg{RODBC} 
#'   functions: \link{sqlQuery} and \link{sqlSave}).
#'   
#' @export
#' @examples 
#' if(interactive()) {
#' library(GGally)
#'
#' policeGraphUn = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_un", 
#'                          directed = FALSE, key = "officer", 
#'                          source = "officer1", target = "officer2", 
#'                          vertexAttrnames = c("offense_count"), edgeAttrnames = c("weight"))
#'                
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#'                              
#' # create network object and visualize with ggplot2
#' net1 = computeGraph(conn, policeGraphUn)
#' ggnet2(net1, node.label="vertex.names", node.size="offense_count", 
#'        legend.position="none")
#'        
#' # network object with filters and color attribute
#' net2 = computeGraph(conn, policeGraphUn, vertexWhere = "officer ~ '[A-Z].*'", 
#'                     edgeWhere = "weight > 0.36")
#' net2 %v% "color" = substr(get.vertex.attribute(net2, "vertex.names"), 1, 1)
#' ggnet2(net2, node.label="vertex.names", node.size="offense_count", 
#'        size.cut=TRUE, node.color="color", legend.position="none", 
#'        palette = "Set2")
#' 
#' # networ object for subgraph of top degree vertices
#' topDegree = computeGraphMetric(conn, policeGraphUn, type="degree", top=50)
#' net3 = computeGraph(conn, policeGraphUn, v=as.list(as.character(topDegree$key)))
#' net3 %v% "degree" = topDegree[match(get.vertex.attribute(net3, "vertex.names"), 
#'                                             topDegree$key), "degree"]
#' ggnet2(net3, node.label="vertex.names", node.size="degree", 
#'        legend.position="none")
#'                          
#' }                          
computeGraph <- function(channel, graph, v=NULL,
                         vertexWhere=graph$vertexWhere, 
                         edgeWhere=graph$edgeWhere, 
                         allTables=NULL, test=FALSE) {
  
  if (missing(graph) || !is.object(graph) || !inherits(graph, "toagraph"))
    stop("Graph object must be specified.")
  
  if (test && is.null(allTables))
    stop("Must provide allTables when test==TRUE.")
  
  isValidConnection(channel, test)
  
  isTableFlag = isTable(channel, c(vertices=graph$vertices, edges=graph$edges), allTables=allTables)
  
  if(!all(isTableFlag | is.na(isTableFlag)))
    stop("Both vertices and edges must exist as tables or views.")
  
  # Handle vertex list if defined
  vertexWhere = addVerticesInVertexWhere(graph, v, vertexWhere)
  
  if(test) {
    emptyLine = "--"
    sqlText = ""
  }
  
  # Edges select
  sqlComment = "-- Edges Select"
  edgesSql = makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere)

  if(test)
    sqlText = paste(sqlComment, edgesSql, sep='\n')
  else
    e = toaSqlQuery(channel, edgesSql, stringsAsFactors=FALSE)
  
  # Vertices select
  if ((!is.null(graph$vertexAttrnames) && length(graph$vertexAttrnames) > 0) ||
      !is.null(vertexWhere)) {
    sqlComment = "-- Vertices Select"
    verticesSql = makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE)
      
    if(test)
      sqlText = paste(sqlText, paste(emptyLine, sqlComment, verticesSql, sep='\n'), sep=';\n')
    else
      vx = toaSqlQuery(channel, verticesSql, stringsAsFactors=FALSE)
  }else
    vx = NULL
  
  # result
  if (test) {
    return(sqlText)
  }else {
    net = makeNetworkResult(graph, vx, e)

    return(net)
  }
  
}


#' Find the vertices not farther than a given limit from another fixed vertex, 
#' and create egographs (subgraphs) with the given order parameter.
#'
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param graph an object of class \code{'toagraph'} referencing graph 
#'   tables in Aster database.
#' @param ego list of vertices for which the calculation of corresponding ego graphs is performed.
#' @param order	integer giving the order of the ego graph neighborhood.
#' @param mode character constant, it specifies how to use the direction of the edges if a directed graph is analyzed. 
#'   For \code{'out'} only the outgoing edges are followed, so all vertices reachable from the source vertex in at 
#'   most order steps are counted. For \code{'in'} all vertices from which the source vertex is reachable in at most 
#'   \code{order} steps are counted. \code{'all'} ignores the direction of the edges. This argument is ignored 
#'   for undirected graphs.
#' @param createDistanceAttr logical: indicates if vertices should receive attribute with the distance
#'   to ego graph centeral vertex.
#' @param distanceAttrname name of the vertex distance attribute.
#' @param vertexWhere SQL WHERE clause limiting data from the vertex table. This value when not null
#'   overrides corresponding value \code{vertexWhere} from \code{graph} (use SQL as if in WHERE clause but 
#'   omit keyword WHERE).
#' @param edgeWhere SQL WHERE clause limiting data from the edge table. This value when not null
#'   overrides corresponding value \code{edgeWhere} from \code{graph} (use SQL as if in WHERE clause but 
#'   omit keyword WHERE).
#' @param allTables pre-built information about existing tables.
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \pkg{RODBC} 
#'   functions: \link{sqlQuery} and \link{sqlSave}).
#' @export
#' @examples 
#' if(interactive()) {
#' library(GGally)
#' 
#' policeGraphDi = toaGraph(vertices = "dallaspolice_officer_vertices", 
#'                          edges = "dallaspolice_officer_edges_di", 
#'                          directed = TRUE,
#'                          key = "officer", source = "officer1", target = "officer2", 
#'                          vertexAttrnames = c("offense_count"),
#'                          edgeAttrnames = c("weight"))
#'                
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#'                          
#' setVertexColor <- function(graph, vertex, color="red", default="grey") {
#'   graph %v% "color" = 
#'   ifelse(get.vertex.attribute(graph, "vertex.names") == as.character(vertex),
#'          color, default)
#'   
#'   return(graph)
#' }
#' 
#' topPagerankPolice = computeGraphMetric(conn, policeGraphDi, type='pagerank', top=3)
#' egoCenters = as.list(as.character(topPagerankPolice$key))
#' 
#' egoGraphsTopPagerank = computeEgoGraph(conn, policeGraphDi, order = 1, ego = egoCenters)
#' 
#' egoGraph = setVertexColor(egoGraphsTopPagerank[[1]], egoCenters[[1]])
#' ggnet2(egoGraph, node.label="vertex.names",  node.size="offense_count",
#'        legend.position="none", color="color")
#'
#' egoGraph = setVertexColor(egoGraphsTopPagerank[[2]], egoCenters[[2]])
#' ggnet2(egoGraph, node.label="vertex.names",  node.size="offense_count",
#'        legend.position="none", color="color")
#'
#' egoGraph = setVertexColor(egoGraphsTopPagerank[[3]], egoCenters[[3]])
#' ggnet2(egoGraph, node.label="vertex.names",  node.size="offense_count",
#'        legend.position="none", color="color")

#' }
computeEgoGraph <- function(channel, graph, ego, order=1, mode="all",
                            createDistanceAttr=TRUE, distanceAttrname="ego.distance",
                            vertexWhere=graph$vertexWhere, edgeWhere=graph$edgeWhere,
                            allTables=NULL, test=FALSE) {
  
  # match argument values
  mode = match.arg(mode, c('all','both','in','out'))
  
  if (missing(graph) || !is.object(graph) || !inherits(graph, "toagraph"))
    stop("Graph object must be specified.")
  
  if (test && is.null(allTables))
    stop("Must provide allTables when test==TRUE.")
  
  if (is.null(ego) || length(ego) == 0)
    stop("Must have at least one ego vertex defined.")
  
  if (!graph$directed && mode %in% c('in','out'))
    stop("Must be a directed graph when mode is 'in' or 'out'.")
  
  isValidConnection(channel, test)
  
  isTableFlag = isTable(channel, c(vertices=graph$vertices, edges=graph$edges), allTables=allTables)
  
  if(!all(isTableFlag | is.na(isTableFlag)))
    stop("Both vertices and edges must exist as tables or views.")
  
  if(test) {
    emptyLine = "--"
    sqlText = ""
  }
  
  egoVertexWhere = addVerticesInVertexWhere(graph, ego, vertexWhere)
  selfTargetOrSource = ifelse(mode == 'in', 'target', 'source')
  otherTargetOrSource = ifelse(mode == 'in', 'source', 'target')
  sqlComment = "-- Create temp table of the shortest paths from ego vertices"
  sqlBeginTran = "BEGIN"
  shortestPathSql = paste0(
      "CREATE TEMP FACT TABLE egographtemp 
       DISTRIBUTE BY HASH(source) 
       AS
       SELECT source, target, distance FROM AllPairsShortestPath(
         ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key , "
         ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
         ON (", makeVerticesSql(graph, isTableFlag, egoVertexWhere, FALSE), ") AS ", selfTargetOrSource,"s PARTITION BY ", graph$key, "
         TARGETKEY('",graph$target,"')
         DIRECTED('", ifelse(graph$directed, 'true', 'false'), "')
         MAXDISTANCE('",order,"')
       )"
  )
  
  if(test) {
    sqlText = paste(sqlBeginTran, paste(emptyLine, sqlComment, shortestPathSql, sep='\n'), sep=';\n')
  }else {
    odbcSetAutoCommit(channel, autoCommit = FALSE)
    toaSqlQuery(channel, sqlBeginTran)
    toaSqlQuery(channel, shortestPathSql)
  }
  
  if (createDistanceAttr) {
    distanceColumnSql = paste0(", eg.distance ", "__distance_attr__")
    distance0ColumnSql = paste0(", 0 ", "__distance_attr__")
  }else {
    distanceColumnSql = ""
    distance0ColumnSql = ""
  }

  egoGraphs = list()
  for(i in 1:length(ego)) {
    
    key = ego[[i]]
    ego.graph = graph
    
    sqlComment = "-- Edges Select"
    edgesSql = paste0(
      "SELECT e.*
         FROM (", makeEdgesSql(ego.graph, isTableFlag, vertexWhere, edgeWhere), ") e 
        WHERE ", graph$source, " IN (SELECT ", otherTargetOrSource, " FROM egographtemp WHERE ", selfTargetOrSource, " = '",key,"')
          AND ", graph$target, " IN (SELECT ", otherTargetOrSource, " FROM egographtemp WHERE ", selfTargetOrSource, " = '",key,"')
       UNION
       SELECT e.*
         FROM (", makeEdgesSql(ego.graph, isTableFlag, vertexWhere, edgeWhere), ") e
        WHERE ", makeEgoSelfEdgeWhereSql(graph, key, mode)
    )
    
    if(test)
      sqlText = paste(sqlText, paste(emptyLine, sqlComment, edgesSql, sep='\n'), sep=';\n')
    else
      e = toaSqlQuery(channel, edgesSql, stringsAsFactors=FALSE)
    
    if ((!is.null(ego.graph$vertexAttrnames) && length(ego.graph$vertexAttrnames) > 0) ||
        createDistanceAttr) {
      sqlComment = "-- Vertices Select"
      verticesSql = paste0(
        "SELECT v.*", distanceColumnSql, " 
           FROM egographtemp eg JOIN
                (", makeVerticesSql(ego.graph, isTableFlag, vertexWhere, FALSE), ") v ON (eg.",otherTargetOrSource," = v.",graph$key,")
          WHERE eg.",selfTargetOrSource," = '",key,"'
         UNION
         SELECT v.*", distance0ColumnSql, " 
           FROM egographtemp eg JOIN
                (", makeVerticesSql(ego.graph, isTableFlag, vertexWhere, FALSE), ") v ON (eg.",selfTargetOrSource," = v.",graph$key,")
          WHERE eg.",selfTargetOrSource," = '",key,"'"
      )
      
      if(test)
        sqlText = paste(sqlText, paste(emptyLine, sqlComment, verticesSql, sep='\n'), sep=';\n')
      else {
        v = toaSqlQuery(channel, verticesSql, stringsAsFactors=FALSE)
        if (createDistanceAttr) {
          ego.graph$vertexAttrnames = c(ego.graph$vertexAttrnames, distanceAttrname)
          names(v)[[length(v)]] = distanceAttrname
        }
      }
    }else
      v = NULL
    
    if(!test)
      egoGraphs[[i]] = makeNetworkResult(ego.graph, v, e)

  }
  
  sqlEndTran = "END"
  if(test) {
    sqlText = paste(sqlText, paste(emptyLine, sqlEndTran, sep='\n'), sep=';\n')
    
    return(sqlText)
  }else {
    toaSqlQuery(channel, sqlEndTran)
    odbcSetAutoCommit(channel, autoCommit = TRUE)
  
    return(egoGraphs)
  }
  
}


#' Compute various statistic distributions on graph edges and vertices.
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param graph an object of class \code{'toagraph'} referencing graph 
#'   tables in Aster database.
#' @param type choose between graph measures to compute histogram distribution for: 
#'   \code{'degree', 'clustering', 'shortestpath', 'pagerank', 'betweenness', 'eigenvector'}.
#' @param weight logical or character: if logical then \code{TRUE} indicates using \code{'weight'} edge
#'   attribute, otherwise no weight used. If character then use as a name for the edge weight attribute. 
#'   The edge weight may apply with types \code{'clustering', 'shortestpath'} and centrality measures.
#' @param binMethod one of several methods to determine number and size of bins: \code{'manual'} indicates to use 
#'   paramters below, both \code{'Sturges'} or \code{'Scott'} will use corresponding methods of computing number
#'   of bins and width (see \url{http://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width}).
#' @param binsize size (width) of discrete intervals defining histogram (all bins are equal).
#' @param startvalue lower end (bound) of values to include in histogram.
#' @param endvalue upper end (bound) of values to include in histogram.
#' @param numbins number of bins to use in histogram.
#' @param vertexWhere SQL WHERE clause limiting data from the vertex table. This value when not null
#'   overrides corresponding value \code{vertexWhere} from \code{graph} (use SQL as if in WHERE clause but 
#'   omit keyword WHERE).
#' @param edgeWhere SQL WHERE clause limiting data from the edge table. This value when not null
#'   overrides corresponding value \code{edgeWhere} from \code{graph} (use SQL as if in WHERE clause but 
#'   omit keyword WHERE).
#' @param allTables pre-built information about existing tables.
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \pkg{RODBC} 
#'   functions: \link{sqlQuery} and \link{sqlSave}).
#' @param ... other arguments passed on to Aster graph functions except for \code{EDGEWEIGHT} argument -
#'   use argument \code{weight} instead. Aster function areguments are not casesensetive 
#' 
#' @export
#' @examples 
#' if(interactive()) {
#' 
#' policeGraphUn = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_un", 
#'                          directed = FALSE, key = "officer", 
#'                          source = "officer1", target = "officer2", 
#'                          vertexAttrnames = c("offense_count"), edgeAttrnames = c("weight"))
#'                
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#'                          
#' hdegreePolice = computeGraphHistogram(conn, policeGraphUn, type='degree', numbins=36) 
#' createHistogram(hdegreePolice, 
#'                 title = "Dallas Police Graph Degree Distribution", 
#'                 xlab='Degree', ylab='Count')
#'
#' hshortestpathPolice = computeGraphHistogram(conn, policeGraphUn, type='shortestpath',
#'                               numbins = 10)
#' createHistogram(hshortestpathPolice, 
#'                 title = "Dallas Police Shortest Path Distribution", 
#'                 xlab = "Distance", ylab = "Count")
#' }
computeGraphHistogram <- function(channel, graph, type='degree', weight=FALSE,
                                  binMethod='manual', numbins = NULL, binsize = NULL, 
                                  startvalue = NULL, endvalue = NULL,
                                  vertexWhere=graph$vertexWhere, edgeWhere=graph$edgeWhere,
                                  allTables=NULL, test=FALSE, ...) {
  
  # match argument values
  type = match.arg(type, c('degree', 'clustering', 'shortestpath', 'pagerank', 'betweenness',
                           'eigenvector', 'closeness', 'avg-closeness', 'k-degree', 'alt-closeness'))
  binMethod = match.arg(binMethod, c('manual', 'Sturges', 'Scott'))
  
  if (missing(graph) || !is.object(graph) || !inherits(graph, "toagraph"))
    stop("Graph object must be specified.")
  
  if (test && is.null(allTables))
    stop("Must provide allTables when test==TRUE.")
  
  isValidConnection(channel, test)
  
  isTableFlag = isTable(channel, c(vertices=graph$vertices, edges=graph$edges), allTables=allTables)
  
  if(!all(isTableFlag | is.na(isTableFlag)))
    stop("Both vertices and edges must exist as tables or views.")
  
  # weight attribute if present
  weight = parseWeightArgument(graph, weight)
  
  argsSql = makeGraphFunctionArgumentsSql(...)
  
  # Graph sql, names
  if (type=='degree') {
    histogramSelectTempTableSql = makeDegreeSelectSql(graph, isTableFlag, NULL, vertexWhere, edgeWhere)
    histogramValue = 'degree'
    histogramGroup = 'degree_type'
  }else if (type=='clustering') {
    histogramSelectTempTableSql = makeClusteringSelectSql(graph, isTableFlag, weight, NULL, vertexWhere, edgeWhere, argsSql)
    histogramValue = 'cc'
    histogramGroup = 'cc_type'
  }else if (type=='closeness') {
    histogramSelectTempTableSql = makeClosenessSelectSql(graph, isTableFlag, weight, NULL, vertexWhere, edgeWhere, argsSql)
    histogramValue = 'closeness'
    histogramGroup = NULL
  }else if (type=='avg-closeness') {
    histogramSelectTempTableSql = makeAvgClosenessSelectSql(graph, isTableFlag, weight, NULL, vertexWhere, edgeWhere, argsSql)
    histogramValue = 'closeness'
    histogramGroup = NULL
  }else if (type=='alt-closeness') {
    histogramSelectTempTableSql = makeAltClosenessSelectSql(graph, isTableFlag, weight, NULL, vertexWhere, edgeWhere, argsSql)
    histogramValue = 'closeness'
    histogramGroup = NULL
  }else if (type=='k-degree') {
    histogramSelectTempTableSql = makeKdegreeSelectSql(graph, isTableFlag, weight, NULL, vertexWhere, edgeWhere, argsSql)
    histogramValue = 'kdegree'
    histogramGroup = NULL
  }else if (type=='shortestpath') {
    histogramSelectTempTableSql = makeShortestPathSelectSql(graph, isTableFlag, weight, vertexWhere, edgeWhere, argsSql)
    histogramValue = 'distance'
    histogramGroup = NULL
  }else if (type=='pagerank') {
    histogramSelectTempTableSql = makePageRankSelectSql(graph, isTableFlag, weight, NULL, vertexWhere, edgeWhere, argsSql)
    histogramValue = 'pagerank'
    histogramGroup = NULL
  }else if (type == 'betweenness') {
    histogramSelectTempTableSql = makeBetweennessSelectSql(graph, isTableFlag, weight, NULL, vertexWhere, edgeWhere, argsSql)
    histogramValue = 'betweenness'
    histogramGroup = NULL
  }else if (type == 'eigenvector') {
    histogramSelectTempTableSql = makeEigenVectorSelectSql(graph, isTableFlag, weight, NULL, vertexWhere, edgeWhere, argsSql)
    histogramValue = 'centrality'
    histogramGroup = NULL
  }

  if (binMethod=='manual') {
    if (is.null(startvalue) && is.null(endvalue) && is.null(numbins))
      stop("Number of bins or/and at least startvalue and endvalue must be defined when method is 'manual'.")
    
    # set number of bins if manual
    if( is.null(startvalue))
      startvalue = 0
    
    if( is.null(endvalue)) 
      endvalue = startvalue + numbins * binsize
    
    if(is.null(startvalue) || is.null(endvalue) || length(startvalue)==0 || length(endvalue)==0) {
      histPrep = paste0("ON hist_prep(
                           ON graphdataforhisttemp VALUE_COLUMN('",histogramValue,"')) as data_stat DIMENSION
                           BIN_SELECT('",numbins,"')")
    }else {
      if (startvalue >= endvalue)
        stop("End value should be greater than start value.")
      
      if(is.null(binsize))
        binsize = (endvalue - startvalue) / numbins
      
      histPrep = paste0("binsize('", binsize, "')
                         startvalue('", startvalue, "')
                         endvalue('", endvalue, "')")
    }
  }else {
    # compute histogram parameters if missing
    histPrep = paste0("ON hist_prep(
                         ON graphdataforhisttemp VALUE_COLUMN('",histogramValue,"')) as data_stat DIMENSION
                         BIN_SELECT('",binMethod,"')")
  }
  
  if(test) {
    emptyLine = "--"
    sqlText = ""
  }
  
  # Create temp table with the node metric
  sqlComment = paste("-- Compute", type, "into temp table with all vertices")
  sqlBeginTran = "BEGIN"
  graphDegreeSql = paste0(
    "CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     ", histogramSelectTempTableSql
  )
  
  if(test) {
    sqlText = paste(sqlBeginTran, paste(emptyLine, sqlComment, graphDegreeSql, sep='\n'), sep=';\n')
  }else {
    odbcSetAutoCommit(channel, autoCommit = FALSE)
    toaSqlQuery(channel, sqlBeginTran)
    toaSqlQuery(channel, graphDegreeSql)
  }
  
  # Compute metric histogram
  sqlComment = "-- Compute metric histogram"
  histDegreeSql = paste0(
    "SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          ", histPrep, "
             VALUE_COLUMN('",histogramValue,"')
            ",ifelse(is.null(histogramGroup),"",paste0("GROUP_COLUMNS('",histogramGroup,"')")),"
         ) PARTITION BY ",ifelse(is.null(histogramGroup),"1",histogramGroup),"
    )"
  )
  
  if(test)
    sqlText = paste(sqlText, paste(emptyLine, sqlComment, histDegreeSql, sep='\n'), sep=';\n')
  else
    histograms = toaSqlQuery(channel, histDegreeSql, stringsAsFactors=FALSE)
  
  sqlEndTran = "END"
  if(test) {
    sqlText = paste(sqlText, paste(emptyLine, sqlEndTran, sep='\n'), sep=';\n')
    
    return(sqlText)
  }else {
    toaSqlQuery(channel, sqlEndTran)
    odbcSetAutoCommit(channel, autoCommit = TRUE)
  
    return(histograms)
  }
  
}


#' Compute top vertices by the metric values on a graph.
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param graph an object of class \code{'toagraph'} referencing graph 
#'   tables in Aster database.
#' @param type choose between graph metrics to compute: \code{'degree', 'in-degree', 
#'   'out-degree', 'clustering', 'shortestpath', 'pagerank', 'betweenness', 
#'   'eigenvector', 'closeness', 'avg-closeness', 'k-degree', 'alt-closeness'}.
#' @param top the number of vertices to return. If \code{top >= 0} then \code{top} vertices 
#'   sorted by the metric value are returned, otherwise all vertices are returned. 
#'   Returned vertices are ordered by the computed graph metric only when \code{top >= 0}. 
#' @param rankFunction one of \code{rownumber, rank, denserank, percentrank}. Rank computed and
#'   returned for each vertex and each metric type. \code{rankFunction} determines which SQL window 
#'   function computes vertex rank value (default \code{rank} corresponds to SQL \code{RANK()} window function). 
#'   When threshold \code{top} is greater than 0 ranking function used to limit number of 
#'   vertices returned (see details).
#' @param weight logical or character: if logical then \code{TRUE} indicates using \code{'weight'} edge
#'   attribute, otherwise no weight used. If character then use as a name for the edge weight attribute. 
#'   The edge weight may apply with types \code{'clustering', 'shortestpath'} and centrality measures.
#' @param vertexWhere SQL WHERE clause limiting data from the vertex table. This value when not null
#'   overrides corresponding value \code{vertexWhere} from \code{graph} (use SQL as if in WHERE clause but 
#'   omit keyword WHERE).
#' @param edgeWhere SQL WHERE clause limiting data from the edge table. This value when not null
#'   overrides corresponding value \code{edgeWhere} from \code{graph} (use SQL as if in WHERE clause but 
#'   omit keyword WHERE).
#' @param keyAsFactor logical: should key column be converted to factor? If \code{TRUE} then conversion
#'   will always take place for any of integer, numeric, or character data types.
#' @param allTables pre-built information about existing tables.
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \pkg{RODBC} 
#'   functions: \link{sqlQuery} and \link{sqlSave}).
#' @param ... other arguments passed on to Aster graph functions except for \code{EDGEWEIGHT} argument -
#'   use argument \code{weight} instead. Aster function areguments are not case-sensetive.
#'   
#' @return dataframe containing one vertice per row with key value, computed metric value, and its rank 
#'   using \code{rankFunction}.
#' 
#' @export
#' @examples 
#' if(interactive()) {
#' library(ggplot2)
#' 
#' policeGraphUn = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_un", 
#'                          directed = FALSE, key = "officer", 
#'                          source = "officer1", target = "officer2", 
#'                          vertexAttrnames = c("offense_count"), edgeAttrnames = c("weight"))
#'                
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#'                          
#' createTopMetricPlot <- function(data, metric, xlab='Officer', ylab='Degree', title) {
#'    p = ggplot(data) +
#'        geom_bar(aes_string("key", metric, fill="key"), stat='identity') +
#'        labs(x=xlab,y=ylab,title=title) +
#'        ggthemes::theme_tufte() + 
#'        theme(legend.position='none', 
#'              axis.text.x = element_text(size=16, angle = 315, vjust = 1), 
#'              plot.title = element_text(size=20),
#'              axis.ticks = element_blank())
#'    
#'    return(p)
#' }
#' 
#' # top degree officers
#' topDegree = computeGraphMetric(conn, policeGraphUn, type="degree", top=30)
#' createTopMetricPlot(topDegree, 'degree', ylab='Degree', title='Top 30 Officers by Degree') 
#' 
#' # top betweenness officers
#' topbetweenness = computeGraphMetric(conn, policeGraphUn, type='betweenness', top=25)
#' createTopMetricPlot(topbetweenness, 'betweenness', ylab='Betweenness', 
#'                     title='Top 25 Officers (Betweenness)')                     
#'                          
#' }
computeGraphMetric <- function(channel, graph, type='degree', top=10, rankFunction='rank',
                               weight=FALSE, vertexWhere=graph$vertexWhere, edgeWhere=graph$edgeWhere,
                               keyAsFactor=TRUE, allTables=NULL, test=FALSE, ...) {
  
  # match argument values
  type = match.arg(type, c('degree', 'in-degree', 'out-degree', 'clustering', 'shortestpath', 'pagerank',
                           'betweenness', 'eigenvector', 'closeness', 'avg-closeness', 'k-degree', 'alt-closeness'))
  rankFunction = match.arg(rankFunction, c('rank', 'rownumber', 'row', 'denserank', 'percentrank'))
  
  if (missing(graph) || !is.object(graph) || !inherits(graph, "toagraph"))
    stop("Graph object must be specified.")
  
  if (test && is.null(allTables))
    stop("Must provide allTables when test==TRUE.")
  
  isValidConnection(channel, test)
  
  isTableFlag = isTable(channel, c(vertices=graph$vertices, edges=graph$edges), allTables=allTables)
  
  if(!all(isTableFlag | is.na(isTableFlag)))
    stop("Both vertices and edges must exist as tables or views.")
  
  # weight attribute if present
  weight = parseWeightArgument(graph, weight)
  
  argsSql = makeGraphFunctionArgumentsSql(...)
  windowFunction = getWindowFunction(rankFunction)
  
  if (type=='degree') {
    selectSql = makeDegreeSelectSql(graph, isTableFlag, windowFunction, vertexWhere, edgeWhere)
    value = 'degree'
    if (graph$directed) {
      group = 'degree_type'
      groupValue = 'degree'
    }else
      group = NULL
  }else if(type=='in-degree') {
    selectSql = makeDegreeSelectSql(graph, isTableFlag, windowFunction, vertexWhere, edgeWhere)
    value = 'degree'
    group = 'degree_type'
    groupValue = 'indegree'
  }else if(type=='out-degree') {
    selectSql = makeDegreeSelectSql(graph, isTableFlag, windowFunction, vertexWhere, edgeWhere)
    value = 'degree'
    group = 'degree_type'
    groupValue = 'outdegree'
  }else if (type=='clustering') {
    selectSql = makeClusteringSelectSql(graph, isTableFlag, weight, windowFunction, vertexWhere, edgeWhere, argsSql)
    value = 'cc'
    group = 'cc_type'
    groupValue = ifelse(graph$directed, 'avg_cc', 'cc')
  }else if (type=='pagerank') {
    selectSql = makePageRankSelectSql(graph, isTableFlag, weight, windowFunction, vertexWhere, edgeWhere, argsSql)
    value = 'pagerank'
    group = NULL
  }else if (type == 'betweenness') {
    selectSql = makeBetweennessSelectSql(graph, isTableFlag, weight, windowFunction, vertexWhere, edgeWhere, argsSql)
    value = 'betweenness'
    group = NULL
  }else if (type == 'eigenvector') {
    selectSql = makeEigenVectorSelectSql(graph, isTableFlag, weight, windowFunction, vertexWhere, edgeWhere, argsSql)
    value = 'centrality'
    group = NULL
  }else if (type == 'closeness') {
    selectSql = makeClosenessSelectSql(graph, isTableFlag, weight, windowFunction, vertexWhere, edgeWhere, argsSql)
    value = 'closeness'
    group = NULL
  }else if (type == 'avg-closeness') {
    selectSql = makeAvgClosenessSelectSql(graph, isTableFlag, weight, windowFunction, vertexWhere, edgeWhere, argsSql)
    value = 'closeness'
    group = NULL
  }else if (type == 'alt-closeness') {
    selectSql = makeAltClosenessSelectSql(graph, isTableFlag, weight, windowFunction, vertexWhere, edgeWhere, argsSql)
    value = 'closeness'
    group = NULL
  }else if (type == 'k-degree') {
    selectSql = makeKdegreeSelectSql(graph, isTableFlag, weight, windowFunction, vertexWhere, edgeWhere, argsSql)
    value = 'kdegree'
    group = NULL
  }
  
  selectSql = paste0(
      selectSql, "
     ",ifelse(is.null(group), "", paste0(" WHERE ", group, " = '", groupValue, "'")), "
     ",ifelse(!is.null(top) && is.numeric(top) && top >= 0, 
              paste0("ORDER BY ", value, " DESC LIMIT ", as.character(top)), "")
  )

  if (test)
    return(selectSql)
  else
    result = toaSqlQuery(channel, selectSql, stringsAsFactors=FALSE)
  
  # make a factor ordered by rank
  if(keyAsFactor) {
    result$key = factor(result$key)
    result$key = factor(result$key, levels=result$key[order(result$rank)], ordered = TRUE)
  }
  
  return(result)
}


parseWeightArgument <- function(graph, weight) {
  if (is.null(weight) || weight==FALSE) weight = NULL
  if (is.logical(weight) && weight) weight = 'weight'
  if (!is.null(weight) && !weight %in% graph$edgeAttrnames)
    stop(paste0("No edge attribute '", weight, "' found in graph."))
  
  return(weight)
}


makeGraphFunctionArgumentsSql <- function(...) {
  
  args = list(...)
  
  if (is.null(args) || length(args)==0) return("")
  
  argNames = names(args)
  result = ""
  for(i in 1:length(args)) {
      result = paste(result, "
           ", argNames[[i]],"('",args[[i]],"')")
  }
  
  return(result)
}


makeDegreeSelectSql <- function(graph, isTableFlag, rankFunction, vertexWhere, edgeWhere) {
  
  if(graph$directed) {
    sql = paste0(
    "SELECT key, degree_type, degree_long degree", 
            ifelse(!is.null(rankFunction), paste0(", ", rankFunction, " OVER (PARTITION BY degree_type ORDER BY degree_long DESC) rank"), ""), "
       FROM unpivot(
         ON (SELECT COALESCE(s.key, t.key) key, 
                    COALESCE(s.cnt_source,0) outdegree, 
                    COALESCE(t.cnt_target,0) indegree,  
                    COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) degree,
                    (COALESCE(t.cnt_target,0) + 1)/(COALESCE(s.cnt_source,0) + 1) inbyoutdegree
               FROM (SELECT ",graph$source," key, COUNT(*) cnt_source 
                       FROM (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") e 
                      GROUP BY 1) s FULL JOIN
                    (SELECT ",graph$target," key, COUNT(*) cnt_target 
                       FROM (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") e
                      GROUP BY 1) t ON (s.key = t.key)
         )
         colsToUnpivot('outdegree','indegree','degree','inbyoutdegree')
         colsToAccumulate('key')
         keepInputColumnTypes('true')
         ATTRIBUTECOLUMNNAME('degree_type')
         VALUECOLUMNNAME('degree')
     )"
    )
  }else {
    sql = paste0(
    "SELECT COALESCE(s.key, t.key) key,
            'degree'::varchar degree_type,
            COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) degree",
            ifelse(!is.null(rankFunction), 
                   paste0(", ", rankFunction, " OVER (PARTITION BY 1 ORDER BY COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) DESC) rank"), 
                   ""), "
      FROM (SELECT ",graph$source," key, COUNT(*) cnt_source 
              FROM (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") e 
             GROUP BY 1) s FULL JOIN
           (SELECT ",graph$target," key, COUNT(*) cnt_target 
              FROM (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") e
             GROUP BY 1) t ON (s.key = t.key)"
    )
  }
  
  return(sql)
}


makeClusteringSelectSql <- function(graph, isTableFlag, weight, rankFunction, vertexWhere, edgeWhere, argsSql) {
  
  sql = paste0(
    "SELECT ", graph$key, " key, cc_type, coalesce(cc_double, cc_str::double) cc",
               ifelse(!is.null(rankFunction),
                      paste0(", ", rankFunction, " OVER (PARTITION BY cc_type ORDER BY coalesce(cc_double, cc_str::double) DESC) rank"), ""), "
       FROM unpivot(
         ON (SELECT * FROM LocalClusteringCoefficient(
               ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key, "
               ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
               targetKey('",graph$target,"')
             ",ifelse(is.null(weight), "", paste0("edgeweight('",weight,"')")),"
               directed('",ifelse(graph$directed,"true","false"),"')
               accumulate('",graph$key,"')",
               argsSql,"
         ))
         colsToUnpivot(",ifelse(graph$directed,"'cyc_cc','mid_cc','in_cc','out_cc','avg_cc'","'cc'"),")
         colsToAccumulate('",graph$key,"')
         keepInputColumnTypes('true')
         ATTRIBUTECOLUMNNAME('cc_type')
         VALUECOLUMNNAME('cc')
     )"
  )
}


makeClosenessSelectSql <- function(graph, isTableFlag, weight, rankFunction, vertexWhere, edgeWhere, argsSql) {
  
  sql = paste0(
    "SELECT ", graph$key, " key, inv_sum_dist closeness",
            ifelse(!is.null(rankFunction),
                      paste0(", ", rankFunction, " OVER (PARTITION BY 1 ORDER BY inv_sum_dist DESC) rank"), ""), "
               FROM Closeness(
               ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key, "
               ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
               targetKey('",graph$target,"')
             ",ifelse(is.null(weight), "", paste0("edgeweight('",weight,"')")),"
               directed('",ifelse(graph$directed,"true","false"),"')
               accumulate('",graph$key,"')",
               argsSql,"
         )"
  )
}


makeAvgClosenessSelectSql <- function(graph, isTableFlag, weight, rankFunction, vertexWhere, edgeWhere, argsSql) {
  
  sql = paste0(
    "SELECT ", graph$key, " key, inv_avg_dist closeness",
            ifelse(!is.null(rankFunction),
                      paste0(", ", rankFunction, " OVER (PARTITION BY 1 ORDER BY inv_avg_dist DESC) rank"), ""), "
               FROM Closeness(
               ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key, "
               ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
               targetKey('",graph$target,"')
             ",ifelse(is.null(weight), "", paste0("edgeweight('",weight,"')")),"
               directed('",ifelse(graph$directed,"true","false"),"')
               accumulate('",graph$key,"')",
               argsSql,"
         )"
  )
}


makeAltClosenessSelectSql <- function(graph, isTableFlag, weight, rankFunction, vertexWhere, edgeWhere, argsSql) {
  
  sql = paste0(
    "SELECT ", graph$key, " key, sum_inv_dist closeness",
            ifelse(!is.null(rankFunction),
                      paste0(", ", rankFunction, " OVER (PARTITION BY 1 ORDER BY sum_inv_dist DESC) rank"), ""), "
               FROM Closeness(
               ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key, "
               ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
               targetKey('",graph$target,"')
             ",ifelse(is.null(weight), "", paste0("edgeweight('",weight,"')")),"
               directed('",ifelse(graph$directed,"true","false"),"')
               accumulate('",graph$key,"')",
               argsSql,"
         )"
  )
}


makeKdegreeSelectSql <- function(graph, isTableFlag, weight, rankFunction, vertexWhere, edgeWhere, argsSql) {
  
  sql = paste0(
    "SELECT ", graph$key, " key, kdegree",
            ifelse(!is.null(rankFunction),
                      paste0(", ", rankFunction, " OVER (PARTITION BY 1 ORDER BY kdegree DESC) rank"), ""), "
               FROM Closeness(
               ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key, "
               ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
               targetKey('",graph$target,"')
             ",ifelse(is.null(weight), "", paste0("edgeweight('",weight,"')")),"
               directed('",ifelse(graph$directed,"true","false"),"')
               accumulate('",graph$key,"')",
               argsSql,"
         )"
  )
}


makeShortestPathSelectSql <- function(graph, isTableFlag, weight, vertexWhere, edgeWhere, argsSql) {
  
  sql = paste0(
    "SELECT source key, target, distance FROM AllPairsShortestPath(
       ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key, "
       ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
       targetKey('",graph$target,"')
       directed('",ifelse(graph$directed,"true","false"),"')
     ",ifelse(is.null(weight), "", paste0("edgeweight('",weight,"')")),
       argsSql,"
)"
  )
}


makePageRankSelectSql <- function(graph, isTableFlag, weight, rankFunction, vertexWhere, edgeWhere, argsSql) {
  
  sql = paste0(
    "SELECT ", graph$key, " key, pagerank",
       ifelse(!is.null(rankFunction),
                      paste0(", ", rankFunction, " OVER (PARTITION BY 1 ORDER BY pagerank DESC) rank"), ""), "
       FROM PageRank(
       ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key, "
       ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
       targetKey('",graph$target,"')
       accumulate('",graph$key,"')
     ",ifelse(is.null(weight), "", paste0("edgeweight('",weight,"')")),
       argsSql,"
)"
  )
}


makeBetweennessSelectSql <- function(graph, isTableFlag, weight, rankFunction, vertexWhere, edgeWhere, argsSql) {
  
  sql = paste0(
    "SELECT ", graph$key, " key, betweenness", 
       ifelse(!is.null(rankFunction),
                      paste0(", ", rankFunction, " OVER (PARTITION BY 1 ORDER BY betweenness DESC) rank"), ""), "
       FROM Betweenness(
       ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key, "
       ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
       targetKey('",graph$target,"')
       accumulate('",graph$key,"')
     ",ifelse(is.null(weight), "", paste0("edgeweight('",weight,"')")),
       argsSql,"
)"
  )
}


makeEigenVectorSelectSql <- function(graph, isTableFlag, weight, rankFunction, vertexWhere, edgeWhere, argsSql) {
  
  sql = paste0(
    "SELECT ", graph$key, " key, centrality",
       ifelse(!is.null(rankFunction),
                      paste0(", ", rankFunction, " OVER (PARTITION BY 1 ORDER BY centrality DESC) rank"), ""), "
       FROM EigenVectorCentrality(
       ON (", makeVerticesSql(graph, isTableFlag, vertexWhere, FALSE), ") AS vertices PARTITION BY ", graph$key, "
       ON (", makeEdgesSql(graph, isTableFlag, vertexWhere, edgeWhere), ") AS edges PARTITION BY ", graph$source, "
       targetKey('",graph$target,"')
       accumulate('",graph$key,"')
       directed('",ifelse(graph$directed,"true","false"),"')
     ",ifelse(is.null(weight), "", paste0("edgeweight('",weight,"')")),
       argsSql,"
     )"
  )
}

addVerticesInVertexWhere <- function(graph, v, vertexWhere) {
  
  if(!is.null(v) && length(v) > 0) {
    if(is.list(v))
      vertexValueWhere = paste0(graph$key, " IN (", makeSqlValueList(unlist(v)) , ")")
    else
      vertexValueWhere = paste0(graph$key, " IN (", v[[1]], ")")
    
    if(is.null(vertexWhere))
      vertexWhere = vertexValueWhere
    else
      vertexWhere = paste0("(", vertexWhere, ") AND ", vertexValueWhere)
  }
  
  return(vertexWhere)
}


makeVerticesSql <- function(graph, isTableFlag, vertexWhere, keyOnlyFlag) {
  
  if(keyOnlyFlag) 
    selectList = graph$key
  else
    selectList = makeSqlColumnList(c(graph$key, graph$vertexAttrnames))
  
  paste0(
    "SELECT ", selectList, " 
       FROM ", makeFromClause(graph$vertices, isTableFlag[['vertices']], "t"),
    makeWhereClause(vertexWhere)
  )
}
 

makeEdgesSql <- function(graph, isTableFlag, vertexWhere, edgeWhere) {
  
  if (!is.null(vertexWhere)) {
    verticesSql = makeVerticesSql(graph, isTableFlag, vertexWhere, TRUE)
    
    if(!is.null(edgeWhere)) {
      edgeWhere = paste0(
        "(", edgeWhere, ") AND 
         ", graph$source, " IN (", verticesSql, ") AND 
         ", graph$target, " IN (", verticesSql, ")"
      )
    }else {
      edgeWhere = paste0(
        graph$source, " IN (", verticesSql, ") AND 
     ", graph$target, " IN (", verticesSql, ")"
      )
    }
  }
  
  paste0(
      "SELECT ", makeSqlColumnList(c(graph$source, graph$target, graph$edgeAttrnames)), " 
         FROM ", makeFromClause(graph$edges, isTableFlag[['edges']], "t"),
      makeWhereClause(edgeWhere)
  )
}


makeEgoSelfEdgeWhereSql <- function(graph, key, mode) {
  
  if (mode %in% c('all','both'))
    whereSql = paste0(graph$source," = '",key,"'
                 OR ",graph$target," = '",key,"'")
  else if (mode == 'in')
    whereSql = paste0(graph$target, " = '",key,"'")
  else
    whereSql = paste0(graph$source, " = '",key,"'")
  
  return(whereSql)
}


makeNetworkResult <- function(graph, v, e){

  # this step is necessary to eliminate integer values which are processed 
  # not like character values by network constructor
  if (is.numeric(e[,graph$source]))
    e[,graph$source] = as.character(e[,graph$source])
  if (is.numeric(e[,graph$target]))
    e[,graph$target] = as.character(e[,graph$target])
  if (!is.null(v) && is.numeric(v[,graph$key]))
    v[,graph$key] = as.character(v[,graph$key])
  
  # create network object using edge list
  net = network(e, directed=graph$directed, matrix.type="edgelist", ignore.eval=FALSE)
  
  # handle the case of nodes without edges
  if (!is.null(v)) {
    vnames = v[, graph$key]
    evnames = get.vertex.attribute(net, "vertex.names")
    if(!all(vnames %in% evnames)) {
      add.vertices(net, length(v[!vnames %in% evnames, graph$key]), 
                   lapply(v[!vnames %in% evnames, graph$key], 
                          FUN = function(x) {list(na=FALSE, vertex.names=x)}))
    }
  }
  
  # vertex data is optional; if it's not null then it contains vertex attributes
  if(!is.null(v)) {
    net.v = data.frame(id=1:length(net$val), vertex.name=matrix(unlist(net$val), ncol=2, byrow=TRUE)[,2],
                       stringsAsFactors=FALSE)
    net.v = merge(net.v, v, by.x="vertex.name", by.y=graph$key, all=FALSE, sort=FALSE)
      
    for(attrname in graph$vertexAttrnames) {
      set.vertex.attribute(net, attrname, net.v[, attrname], net.v[, "id"])
    }
  }
  
  return(net)
}