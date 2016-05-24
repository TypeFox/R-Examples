#' Perform k-means clustering on the table.
#' 
#' K-means clustering algorithm runs in-database, returns object compatible with \code{\link{kmeans}} and 
#' includes arbitrary aggregate metrics computed on resulting clusters.
#' 
#' The function fist scales not-null data (if \code{scale=TRUE}) or just eliminate nulls without scaling. After 
#' that the data given (table \code{tableName} with option of filering with \code{where}) are clustered by the 
#' k-means in Aster. Next, all standard metrics of k-means clusters plus additional aggregates provided with
#' \code{aggregates} are calculated again in-database.
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}.
#' @param tableName Aster table name.
#' @param tableInfo pre-built summary of data to use (require when \code{test=TRUE}). See \code{\link{getTableSummary}}.
#' @param id column name or SQL expression containing unique table key.
#' @param idAlias SQL alias for table id. This is required when SQL expression is given for \code{id}.
#' @param include a vector of column names with variables (must be numeric). Model never contains variables other than in the list.
#' @param except a vector of column names to exclude from variables. Model never contains variables from the list.
#' @param centers either the number of clusters, say \code{k}, or a matrix of initial (distinct) cluster centres. 
#'   If a number, a random set of (distinct) rows in x is chosen as the initial centres. If a matrix then number 
#'   of rows determines the number of clusters as each row determines initial center.
#' @param threshold the convergence threshold. When the centroids move by less than this amount, 
#'   the algorithm has converged.
#' @param iterMax the maximum number of iterations the algorithm will run before quitting if the convergence 
#'   threshold has not been met.
#' @param aggregates vector with SQL aggregates that define arbitrary aggreate metrics to be computed on each cluster 
#'   after running k-means. Aggregates may have optional aliases like in \code{"AVG(era) avg_era"}. 
#'   Subsequently, used in \code{\link{createClusterPlot}} as cluster properties.
#' @param scale logical if TRUE then scale each variable in-database before clustering. Scaling performed results in 0 mean and unit
#'   standard deviation for each of input variables.
#' @param where specifies criteria to satisfy by the table rows before applying
#'   computation. The creteria are expressed in the form of SQL predicates (inside
#'   \code{WHERE} clause).
#' @param scaledTableName name of Aster table with results of scaling
#' @param centroidTableName name of Aster table with centroids found by kmeans  
#' @param schema name of Aster schema tables \code{scaledTableName} and \code{centroidTableName} belong.
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \pkg{RODBC} 
#'   functions: \link{sqlQuery} and \link{sqlSave}).
#' @return \code{computeKmeans} returns an object of class \code{"toakmeans"} (compatible with class \code{"kmeans"}).
#' It is a list with at least the following components:
#' \describe{
#'   \item{\code{cluster}}{A vector of integers (from 0:K-1) indicating the cluster to which each point is allocated. 
#'     \code{\link{computeKmeans}} leaves this component empty. Use function \code{\link{computeClusterSample}} to set this compoenent.}
#'   \item{\code{centers}}{A matrix of cluster centres.}
#'   \item{\code{totss}}{The total sum of squares.}
#'   \item{\code{withinss}}{Vector of within-cluster sum of squares, one component per cluster.}
#'   \item{\code{tot.withinss}}{Total within-cluster sum of squares, i.e. \code{sum(withinss)}.}
#'   \item{\code{betweenss}}{The between-cluster sum of squares, i.e. \code{totss-tot.withinss}.}
#'   \item{\code{size}}{The number of points in each cluster. These includes all points in the Aster table specified that 
#'     satisfy optional \code{where} condition.}
#'   \item{\code{iter}}{The number of (outer) iterations.}
#'   \item{\code{ifault}}{integer: indicator of a possible algorithm problem (always 0).}
#'   \item{\code{scale}}{logical: indicates if variable scaling was performed before clustering.}
#'   \item{\code{aggregates}}{Vectors (dataframe) of aggregates computed on each cluster.}
#'   \item{\code{tableName}}{Aster table name containing data for clustering.}
#'   \item{\code{columns}}{Vector of column names with variables used for clustering.}
#'   \item{\code{scaledTableName}}{Aster table name containing scaled data for clustering.}
#'   \item{\code{centroidTableName}}{Aster table name containing cluster centroids.}
#'   \item{\code{id}}{Column name or SQL expression containing unique table key.}
#'   \item{\code{idAlias}}{SQL alias for table id.}
#'   \item{\code{whereClause}}{SQL \code{WHERE} clause expression used (if any).}
#'   \item{\code{time}}{An object of class \code{proc_time} with user, system, and total elapsed times
#'     for the \code{computeKmeans} function call.} 
#' }
#' 
#' @export
#' @seealso \code{\link{computeClusterSample}}, \code{\link{computeSilhouette}}
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
#' km
#' createCentroidPlot(km)
#' createClusterPlot(km)
#' }
computeKmeans <- function(channel, tableName, centers, threshold=0.0395, iterMax=10, 
                          tableInfo, id, include=NULL, except=NULL, 
                          aggregates="COUNT(*) cnt", scale=TRUE, idAlias=gsub("[^0-9a-zA-Z]+", "_", id), 
                          where=NULL, scaledTableName=NULL, centroidTableName=NULL, schema=NULL,
                          test=FALSE) {
  
  ptm = proc.time()
  
  if (test & missing(tableInfo)) {
    stop("Must provide tableInfo when test==TRUE")
  }
  
  isValidConnection(channel, test)
  
  # validate centers (initial clusters)
  if (is.matrix(centers)) 
    K = nrow(centers)
  else if (is.numeric(centers)) 
    K = as.integer(centers)
  else 
    stop("Parameter centers must be numeric.")
  
  if (K < 1) 
    stop("Number of clusters must be greater or equal to 1.")
  
  tableName = normalizeTableName(tableName)
  
  if (missing(tableInfo)) {
    tableInfo = sqlColumns(channel, tableName)
  }
  
  columns = getNumericColumns(tableInfo, names.only=TRUE, include=include, except=except)
  columns = sort(setdiff(columns, id))
  
  if (is.null(columns) || length(columns) < 1) {
    stop("Kmeans operates on one or more numeric variables.")
  }
  
  # check if id alias is not one of independent variables
  if(idAlias %in% columns)
    stop(paste0("Id alias '", idAlias, "' can't be one of variable names."))
  
  # adjust id alias if it's exactly one of the table columns
  if(idAlias %in% tableInfo$COLUMN_NAME)
    idAlias = paste("_", idAlias, "_", sep="_")
  
  if (is.matrix(centers))
    if (length(columns) != ncol(centers))
      stop(paste0("Kmeans received incompatible parameters: dimension of initial cluster centers doesn't match variables: '", 
                  paste0(columns, collapse = "', '"), "'"))
  
  aggregates = makeAggregatesAlwaysContainCount(aggregates)
  
  if (is.null(scaledTableName))
    scaledTableName = makeTempTableName('scaled', 30, schema)
  else if (!is.null(schema))
    scaledTableName = paste0(schema, ".", scaledTableName)
  
  if(is.null(centroidTableName))
    centroidTableName = makeTempTableName('centroids', 30, schema)
  else if (!is.null(schema))
    centroidTableName = paste0(schema, ".", centroidTableName)
  
  where_clause = makeWhereClause(where)
  
  emptyLine = "--"
  
  if(test)
    sqlText = ""
  
  # scale data or just eliminate incomplete observations (if not scaling)
  sqlComment = "-- Data Prep: scale"
  sqlDrop = paste("DROP TABLE IF EXISTS", scaledTableName)
  sql = getDataPrepSql(scale, tableName, scaledTableName, columns, id, idAlias, where_clause)
  if(test) {
    sqlText = paste(sqlComment, sqlDrop, sep='\n')
    sqlText = paste(sqlText, sql, sep=';\n')
  }else {
    toaSqlQuery(channel, sqlDrop)
    toaSqlQuery(channel, sql)
  }
    
  
  # run kmeans
  sqlComment = "-- Run k-means"
  sqlDrop = paste("DROP TABLE IF EXISTS", centroidTableName)
  sql = getKmeansSql(scaledTableName, centroidTableName, centers, threshold, iterMax)
  if(test) {
    sqlText = paste(sqlText, emptyLine, sqlComment, sqlDrop, sql, sep=';\n')
  }else {
    toaSqlQuery(channel, sqlDrop)
    kmeansResultStr = toaSqlQuery(channel, sql, stringsAsFactors=FALSE)
    if (kmeansResultStr[2,'message'] == "Successful!" &&
        kmeansResultStr[3,'message'] == "Algorithm converged.") {
      iter = as.integer(gsub("[^0-9]", "", kmeansResultStr[[4,'message']]))
    }else {
      msg = paste(kmeansResultStr[,'message'], collapse="\n")
      stop(msg)
    }
  }
  
  
  # compute cluster stats
  sqlComment = "-- Run cluster assignment, cluster stats, and within-cluster sum of squares"
  sql = getKmeansStatsSql(tableName, scaledTableName, centroidTableName, columns, 
                         K, id, idAlias, aggregates, where_clause)
  if(test)
    sqlText = paste(sqlText, emptyLine, sqlComment, sql, sep=';\n')
  else
    kmeansstats = toaSqlQuery(channel, sql, stringsAsFactors=FALSE)
  
  
  # compute total sum of squares
  sqlComment = "-- Compute Total Sum of Squares"
  sql = getTotalSumOfSquaresSql(scaledTableName, columns, idAlias, scale)
  
  if(test)
    sqlText = paste(sqlText, emptyLine, sqlComment, sql, sep=';\n')
  else {
    rs = toaSqlQuery(channel, sql)
    totss = rs$totss[[1]]
  }

  
  # return sql
  if(test) {
    sqlText = paste0(sqlText, ';')
    return(sqlText)
  }
  
  result = makeKmeansResult(kmeansstats, K, totss, iter, tableName, columns, scale,
                            scaledTableName, centroidTableName, id, idAlias, 
                            where_clause, ptm)
  
  return(result)
}


# Phase 1: Data Prep
getDataPrepSql <- function(scale, tableName, tempTableName, columns, id, idAlias, whereClause) {
  
  dataPrepSql = ifelse(scale, 
                      getDataScaledSql(tableName, columns, id, idAlias, whereClause),
                      getDataNoNullsSql(tableName, columns, id, idAlias, whereClause))
  
  tempTableSql = paste0(
    "CREATE FACT TABLE ", tempTableName, " DISTRIBUTE BY HASH(", idAlias, ") AS 
       ", dataPrepSql
  )
}

getDataScaledSql <- function(tableName, columns, id, idAlias, whereClause) {
  
  sqlmr_column_list = makeSqlMrColumnList(columns)
  query_as_table = getDataSql(tableName, columns, id, idAlias, whereClause)
  
  scaleMapSql = paste0(
    "SELECT * FROM ScaleMap (
       ON (", query_as_table, ")
    InputColumns (", sqlmr_column_list, ")
    -- MissValue ('OMIT')
    )"
  )
  
  scaleSql = paste0(
    "SELECT * FROM Scale(
       ON (", query_as_table, ") AS input PARTITION BY ANY
       ON (", scaleMapSql, ") AS STATISTIC DIMENSION
       Method ('STD')
       Accumulate('", idAlias, "')
       GlobalScale ('false')
       InputColumns (", sqlmr_column_list, ")
     )"
  )
  
  return(scaleSql)
}

getDataNoNullsSql <- function(tableName, columns, id, idAlias, whereClause) {
  
  not_null_clause = paste0(c(idAlias, columns), " IS NOT NULL", collapse=" AND ")
  query_as_table = getDataSql(tableName, columns, id, idAlias, whereClause)
  
  noNullsSql = paste0(
    "SELECT * FROM (", query_as_table, ") d
      WHERE ", not_null_clause)
  
}

# Phase: kmeans 
getKmeansSql <- function(scaledTableName, centroidTableName, centers, threshold, maxiternum) {
  
  if (is.matrix(centers)) {
    initCenters = paste0("MEANS(", paste0("'", paste0(apply(centers, 1, paste0, collapse='_'), collapse="', '"), "'"), ")")
  }else
    initCenters = paste0("NUMBERK('", centers, "')")
  
  kmeansSql = paste0(
    "SELECT * FROM kmeans(
      ON (SELECT 1)
      PARTITION BY 1
      INPUTTABLE('", scaledTableName, "')
      OUTPUTTABLE('", centroidTableName, "')
   ", initCenters, "
      THRESHOLD('", threshold, "')
      MAXITERNUM('", maxiternum, "')
    )")
}

getKmeansStatsSql <- function(tableName, scaledTableName, centroidTableName, columns, K, id, idAlias, aggregates, whereClause) {
  
  clustersWithValuesSql = paste0(
    "SELECT c1.*, c2.withinss  
       FROM (SELECT clusterid, means, ", paste(aggregates, collapse=", "), 
    "          FROM (", getClusteredDataSql(tableName, scaledTableName, centroidTableName, columns, id, idAlias, whereClause), "
                    ) clustered_data
              GROUP BY clusterid, means
            ) c1 JOIN ( 
            ", paste(sapply(1:K, FUN=getClusterSumOfSquaresSql, scaledTableName, centroidTableName, columns, idAlias), 
                     collapse="\nUNION ALL\n"),
    "
            ) c2 ON (c1.clusterid = c2.clusterid)
      ORDER BY clusterid"
  )
}

getDataSql <- function(tableName, columns, id, idAlias, whereClause) {
  
  paste0("SELECT ", id, " ", idAlias, ", ", makeSqlColumnList(columns), " FROM ", tableName, whereClause)
}


getTableDataSql <- function(tableName, id, idAlias, whereClause) {
  
  paste0("SELECT ", id, " ", idAlias, ", * FROM ", tableName, whereClause)
}


getClusteredDataSql <- function(tableName, scaledTableName, centroidTableName, columns, id, idAlias, whereClause) {
  
  query_as_table = getTableDataSql(tableName, id, idAlias, whereClause)
  
  paste0(
    "SELECT c.clusterid, c.means, d.* 
      FROM ", centroidTableName, " c JOIN 
    kmeansplot (
      ON ", scaledTableName, " PARTITION BY ANY
      ON ", centroidTableName, " DIMENSION
      centroidsTable('",centroidTableName,"')
    ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
    (", query_as_table, ") d on (kmp.", idAlias, " = d.", idAlias, ")"
  )
}


getClusterSumOfSquaresSql <- function(clusterid, scaledTableName, centroidTableName, columns, idAlias) {
  
  clusterid = as.character(clusterid - 1)
  
  sql = paste0(
    "SELECT ", clusterid, " clusterid, SUM(distance::double ^ 2) withinss FROM ", 
    getUnpivotedClusterSql(clusterid, scaledTableName, centroidTableName, columns, idAlias)
  )
}


getUnpivotedClusterSql <- function(clusterid, scaledTableName, centroidTableName, columns, idAlias) {
  
  sql_column_list = makeSqlColumnList(columns)
  sqlmr_column_list = makeSqlMrColumnList(columns)
  
  sql = paste0(
    
    "VectorDistance(
       ON (
         SELECT clusterid, ", idAlias, ", variable, coalesce(value_double, value_long, value_str::double) value
           FROM unpivot(
                  ON (SELECT d.* 
                        FROM kmeansplot (
                               ON ", scaledTableName, " PARTITION BY ANY
                               ON ", centroidTableName, " DIMENSION
                               centroidsTable('",centroidTableName,"')
                             ) d 
                       WHERE clusterid = ", clusterid, "
                  )
                  COLSTOUNPIVOT(", sqlmr_column_list, ")
                  COLSTOACCUMULATE('",idAlias,"','clusterid')
                  ATTRIBUTECOLUMNNAME('variable')
                  VALUECOLUMNNAME('value')
                  KEEPINPUTCOLUMNTYPES('true')
                )
       ) AS target PARTITION BY ", idAlias, "
       ON (
         ", getCentroidTableSql(centroidTableName, sql_column_list, clusterid), "
       ) AS ref DIMENSION
       TARGETIDCOLUMNS('",idAlias,"')
       TARGETFEATURECOLUMN('variable')
       TARGETVALUECOLUMN('value')
       REFIDCOLUMNS('clusterid')
       REFFEATURECOLUMN('variable')
       REFVALUECOLUMN('value')
       MEASURE('Euclidean')
     )"
  )
  
}


getCentroidTableSql <- function(centroidTableName, sql_column_list, clusterid = NULL) {
  
  whereClause = ifelse(is.null(clusterid), 
                       "", 
                       paste0(" WHERE clusterid = ", clusterid))

  paste0(
    "SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('", sql_column_list, "', ', ') variable 
           FROM ", centroidTableName, whereClause
    
  )
}


getTotalSumOfSquaresSql <- function(scaledTableName, columns, idAlias, scale) {
  
  sqlmr_column_list = makeSqlMrColumnList(columns)
  
  if (scale) {
    sqlagg_column_list = paste0("0.0::double ", columns, collapse = ", ")
    globalCenterSql = paste0("SELECT 1 id, ", sqlagg_column_list)
  }else {
    sqlagg_column_list = makeSqlAggregateColumnList(columns, "avg", FALSE, cast="::double")
    globalCenterSql = paste0("SELECT 1 id, ", sqlagg_column_list, " FROM ", scaledTableName)
  }

  sql = paste0(
    "SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
       ON (SELECT ", idAlias, ", variable, coalesce(value_double, value_long, value_str::double) value
             FROM unpivot(
               ON ", scaledTableName , "
               COLSTOUNPIVOT(", sqlmr_column_list, ")
               COLSTOACCUMULATE('",idAlias,"')
               ATTRIBUTECOLUMNNAME('variable')
               VALUECOLUMNNAME('value')
               KEEPINPUTCOLUMNTYPES('true')
             ) 
       ) AS target PARTITION BY ", idAlias, "
       ON (SELECT id, variable, value_double
             FROM unpivot(
               ON (", globalCenterSql, ")
               COLSTOUNPIVOT(", sqlmr_column_list, ")
               COLSTOACCUMULATE('id')
               ATTRIBUTECOLUMNNAME('variable')
               VALUECOLUMNNAME('value')
               KEEPINPUTCOLUMNTYPES('true')
             )
       ) AS ref DIMENSION
       TARGETIDCOLUMNS('",idAlias,"')
       TARGETFEATURECOLUMN('variable')
       TARGETVALUECOLUMN('value')
       REFIDCOLUMNS('id')
       REFFEATURECOLUMN('variable')
       REFVALUECOLUMN('value_double')
       MEASURE('Euclidean')
     )"
  )
}


makeKmeansResult <- function(data, K, totss, iter, tableName, columns, scale,
                             scaledTableName, centroidTableName, id, idAlias, 
                             whereClause, ptm) {
  
  # parse data (kmeansplot) and form kmeans object
  centers = matrix(as.numeric(unlist(strsplit(as.vector(data$means), split = " "))), 
                   ncol=length(columns), nrow=length(data$means), byrow=TRUE)
  colnames(centers) = columns
  rownames(centers) = data$clusterid
  
  aggregates = data.frame(clusterid=data$clusterid, data[, c(-1,-2)])
  
  tot_withinss = sum(data$withinss)
  
  z <- structure(list(cluster=integer(0),
                      centers=centers,
                      totss=totss,
                      withinss = data$withinss,
                      tot.withinss = tot_withinss,
                      betweenss = totss - tot_withinss,
                      size = aggregates$cnt,
                      iter=iter,
                      ifault = 0,
                      
                      scale=scale,
                      aggregates=aggregates,
                      tableName=tableName,
                      columns=columns,
                      scaledTableName=scaledTableName,
                      centroidTableName=centroidTableName,
                      id=id,
                      idAlias=idAlias,
                      whereClause=whereClause,
                      time=proc.time() - ptm
  ),
  class = c("toakmeans", "kmeans"))
  
  return (z)
}


makeAggregatesAlwaysContainCount <- function(aggregates){
  
  # empty or NULL list 
  if (length(aggregates) == 0)
    return("COUNT(*) cnt")
 
  # parse aggregates into tuples of function and alias
  aggFun = unlist(sapply(strsplit(aggregates, '[[:space:]]'), 
                 FUN=function(x) {
                   s = paste0(x[1:length(x)-1], collapse = ' ')
                   if (nchar(s)==0) NULL else s
                 }))
 
  aggAlias = unlist(sapply(strsplit(aggregates, '[[:space:]]'), 
                   FUN=function(x) {
                     s = x[length(x)]
                     if (nchar(s)==0) NULL else s
                  }))
 
  # detect missing alias
  if (length(aggFun) != length(aggAlias))
   stop("Check aggregates: at least one missing alias found.")
 
  # detect if 'COUNT(*) cnt' present
  missingCount = TRUE
  for(i in 1:length(aggregates)) {
    fun = aggFun[[i]]
    alias = aggAlias[[i]]
   
    if(tolower(fun) == 'count(*)' && alias == 'cnt')
      missingCount = FALSE
 }
 
  # form final aggregates
  aggregates = paste(aggFun, aggAlias)
  if (missingCount)
    aggregates = c(aggregates, "COUNT(*) cnt")
 
  return(aggregates)
}


#' Random sample of clustered data
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}.
#' @param km an object of class \code{"toakmeans"} obtained with \code{computeKmeans}.
#' @param sampleFraction one or more sample fractions to use in the sampling of data. (multipe 
#'   sampling fractions are not yet supported.)
#' @param sampleSize total sample size (applies only when \code{sampleFraction} is missing).
#' @param scaled logical: indicates if original (default) or scaled data returned.
#' @param includeId logical indicates if sample should include the key uniquely identifying
#'   each data row.
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \pkg{RODBC} 
#'   functions: \link{sqlQuery} and \link{sqlSave}).
#' @return \code{computeClusterSample} returns an object of class \code{"toakmeans"} (compatible with class \code{"kmeans"}).
#' @seealso \code{\link{computeKmeans}}
#' 
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
#' km
#' createClusterPairsPlot(km, title="Batters Clustered by G, H, R", ticks=FALSE)
#' }
computeClusterSample <- function(channel, km, sampleFraction, sampleSize, scaled=FALSE, includeId=FALSE, test=FALSE) {
  
  isValidConnection(channel, test)
  
  if (missing(km) || !is.object(km) || !inherits(km, "toakmeans")) {
    stop("Kmeans object must be specified.")
  }
  
  if ((missing(sampleFraction) || is.null(sampleFraction)) && 
      (missing(sampleSize) || is.null(sampleSize))) {
    stop("Sample fraction or sample size must be specified.")
  }
  
  table_name = km$tableName
  columns = km$columns
  scaled_table_name = km$scaledTableName
  centroid_table_name = km$centroidTableName
  id = km$id
  idAlias = km$idAlias
  where_clause = km$whereClause
  centers = nrow(km$centers)
  
  conditionOnSql = paste0("'",paste0(1:centers-1, collapse="','"),"'")
  query_as_table = getDataSql(table_name, columns, id, idAlias, where_clause)
  
  
  if (!missing(sampleFraction) && !is.null(sampleFraction)) {
    
    # using sample fraction
    stopifnot(sampleFraction >= 0, sampleFraction <= 1)

    sql = paste0(
      "SELECT * FROM sample(
             ON (", getKmeansplotDataSql(scaled_table_name, centroid_table_name, scaled, query_as_table, idAlias), "
             )
             CONDITIONONCOLUMN('clusterid')
             CONDITIONON(",conditionOnSql,")
             SAMPLEFRACTION('", as.character(sampleFraction), "')
       )")
  }else {
    # using sample size
    sql = paste0(
      "WITH stratum_counts AS (
         SELECT clusterid stratum, count(*) stratum_count 
           FROM kmeansplot(
             ON ", scaled_table_name, " PARTITION BY ANY
             ON ", centroid_table_name, " DIMENSION
             centroidsTable('baseball.kmeans_test_centroids')
           ) 
          WHERE clusterid != -1
         GROUP BY 1
       )
       SELECT * FROM sample (
         ON (", getKmeansplotDataSql(scaled_table_name, centroid_table_name, scaled, query_as_table, idAlias), "
            ) AS data PARTITION BY ANY
         ON stratum_counts AS summary DIMENSION
         CONDITIONONCOLUMN('clusterid')
         CONDITIONON(",conditionOnSql,")
         ApproximateSampleSize('", as.character(sampleSize), "')
      )"
    )
  }
  
  if (!includeId) {
    sql = paste0(
      "SELECT * FROM antiselect(
         ON 
           (",sql,"
           )
         EXCLUDE('",idAlias,"')
       )")
  }
  
  if(test) {
    return(sql)
  }else {
    data = toaSqlQuery(channel, sql)
    km$cluster = data$clusterid
    km$data = data
    return(km)
  }
}


getKmeansplotDataSql <- function(scaled_table_name, centroid_table_name, scaled, query_as_table, idAlias) {
  
  sql = paste0(
                "SELECT ", ifelse(scaled,  " d.* ", " clusterid, d.* "), "
                   FROM kmeansplot(
                     ON ", scaled_table_name, " PARTITION BY ANY
                     ON ", centroid_table_name, " DIMENSION
                     centroidsTable('",centroid_table_name,"')
                   ) ", ifelse(scaled, " d ", 
                               paste0( " kmp JOIN (", query_as_table, ") d ON (kmp.", idAlias, " = d.", idAlias, ")")), "
                  WHERE clusterid != -1"
  )
}


#' Compute Silhouette (k-means clustering).
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}.
#' @param km an object of class \code{"toakmeans"} obtained with \code{computeKmeans}.
#' @param scaled logical: indicates if computation performed on original (default) or scaled values.
#' @param silhouetteTableName name of the Aster table to hold silhouette scores. The table persists silhoutte scores 
#'   for all clustered elements. Set parameter \code{drop=F} to keep the table.
#' @param drop logical: indicates if the table \code{silhouetteTableName}
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \pkg{RODBC} 
#'   functions: \link{sqlQuery} and \link{sqlSave}).
#' @return \code{computeSilhouette} returns an object of class \code{"toakmeans"} (compatible with class \code{"kmeans"}).
#'   It adds a named list \code{sil} the \code{km} containing couple of elements: average value of silhouette \code{value} and silhouette profile  
#'   (distribution of silhouette values on each cluster) \code{profile}
#' @seealso \code{\link{computeKmeans}}
#'
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
#' km$sil
#' createSilhouetteProfile(km, title="Cluster Silhouette Histograms (Profiles)")
#' }
computeSilhouette <- function(channel, km, scaled=TRUE, silhouetteTableName=NULL, drop=TRUE, test=FALSE) {
  
  ptm = proc.time()
  
  isValidConnection(channel, test)
  
  if(test && is.null(silhouetteTableName)){
    stop("Silhouette table name is required when test=TRUE.")
  }
    
  if (missing(km) || !is.object(km) || !inherits(km, "toakmeans")) {
    stop("Kmeans object must be specified.")
  }
  
  if (nrow(km$centers) == 1)
    stop("Silhouette values are trivial in case of single cluster model.")
  
  if (is.null(silhouetteTableName))
    silhouetteTableName = makeTempTableName('silhouette', 30)
  
  table_name = km$tableName
  columns = km$columns
  scaled_table_name = km$scaledTableName
  centroid_table_name = km$centroidTableName
  id = km$id
  idAlias = km$idAlias
  where_clause = km$whereClause
  
  emptyLine = "--"
  
  if(test)
    sqlText = ""
  
  # make silhouette data
  sqlComment = "-- Create Analytical Table with Silhouette Data"
  sqlDrop = paste("DROP TABLE IF EXISTS", silhouetteTableName)
  sql = makeSilhouetteDataSql(table_name, silhouetteTableName, columns, id, idAlias, where_clause, scaled_table_name, centroid_table_name, scaled)
  if(test) {
    sqlText = paste(sqlComment, sqlDrop, sep='\n')
    sqlText = paste(sqlText, sql, sep=';\n')
  }else {
    toaSqlQuery(channel, sqlDrop)
    toaSqlQuery(channel, sql)
  }
  
  # Compute overall silhouette value
  sqlComment = "-- Compute overall silhouette value"
  sql = paste0("SELECT AVG((b-a)/greatest(a,b)) silhouette_value FROM ", silhouetteTableName)
  if(test) {
    sqlText = paste(sqlText, emptyLine, sqlComment, sql, sep=';\n')
  }else {
    data = toaSqlQuery(channel, sql)
    sil_value = data[[1,'silhouette_value']]
  }
  
  # Compute silhouette profiles (histograms by cluster)
  sqlComment = "-- Compute silhouette cluster profiles"
  sql = paste0(
    "SELECT * FROM Hist_Reduce(
       ON Hist_Map(
         ON (SELECT clusterid::varchar clusterid, (b-a)/greatest(a,b) silhouette_value FROM ", silhouetteTableName, "
         )
         STARTVALUE('-1')
         BINSIZE('0.05')
         ENDVALUE('1')
         VALUE_COLUMN('silhouette_value')
         GROUP_COLUMNS('clusterid')
       ) PARTITION BY clusterid
     )"
  )
  if(test){
    sqlText = paste(sqlText, emptyLine, sqlComment, sql, sep=';\n')
  }else {
    data = toaSqlQuery(channel, sql)
    sil_profile = toaSqlQuery(channel, sql, stringsAsFactors=TRUE)
  }
  
  # Drop Analytical Table with Silhouette Data 
  if(drop) {
    sqlComment = "-- Drop Analytical Table with Silhouette Data"
    sql = paste0("DROP TABLE IF EXISTS ", silhouetteTableName)
    if(test){
      sqlText = paste(sqlText, emptyLine, sqlComment, sql, sep=';\n')
    }else {
      toaSqlQuery(channel, sql)
    }
  }
  
  # get result back
  if(test){
    sqlText = paste0(sqlText, ';')
    return(sqlText)
  }
  
  sil = list(value=sil_value, profile=sil_profile)
  if(!drop) {
    sil = c(sil, tableName=silhouetteTableName)
  }
  sil$time = proc.time() - ptm
  km$sil = sil
  
  return(km)
}


makeSilhouetteDataSql <- function(table_name, temp_table_name, columns, id, idAlias, 
                                  where_clause, scaled_table_name, centroid_table_name, scaled) {
  
  sqlmr_column_list = makeSqlMrColumnList(columns)
  query_as_table = getDataSql(table_name, columns, id, idAlias, where_clause)
  
  sql = paste0(
    "CREATE ANALYTIC TABLE ", temp_table_name, "
     DISTRIBUTE BY HASH(clusterid)
     AS
     WITH kmeansplotresult AS (
         SELECT clusterid, ", idAlias, ", variable, coalesce(value_double, value_long, value_str::double) value
           FROM unpivot(
                  ON (", getKmeansplotDataSql(scaled_table_name, centroid_table_name, scaled, query_as_table, idAlias), "
                  )
                  COLSTOUNPIVOT(", sqlmr_column_list, ")
                  COLSTOACCUMULATE('",idAlias,"','clusterid')
                  ATTRIBUTECOLUMNNAME('variable')
                  VALUECOLUMNNAME('value')
                  KEEPINPUTCOLUMNTYPES('true')
                ) 
     )
     SELECT target_clusterid clusterid, target_", idAlias, " ", idAlias, ", a, b 
       FROM (
         SELECT target_clusterid, target_", idAlias, ",
                MAX(CASE WHEN target_clusterid = ref_clusterid THEN dissimilarity ELSE 0 END) a,
                MIN(CASE WHEN target_clusterid = ref_clusterid THEN 'Infinity' ELSE dissimilarity END) b 
           FROM
             (SELECT target_clusterid, target_", idAlias, ", ref_clusterid, avg(distance) dissimilarity 
                FROM VectorDistance(
                  ON kmeansplotresult AS target PARTITION BY ", idAlias, "
                  ON kmeansplotresult AS ref DIMENSION
                  TARGETIDCOLUMNS('clusterid','", idAlias, "')
                  TARGETFEATURECOLUMN('variable')
                  TARGETVALUECOLUMN('value')
                  REFIDCOLUMNS('clusterid','", idAlias, "')
                  REFFEATURECOLUMN('variable')
                  REFVALUECOLUMN('value')
                  MEASURE('Euclidean')
                ) 
          WHERE target_", idAlias," != ref_", idAlias, " 
          GROUP BY 1,2,3
         ) agg  
       GROUP BY 1,2
     ) sil"
  )
}