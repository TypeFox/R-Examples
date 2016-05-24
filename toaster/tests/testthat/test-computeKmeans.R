context("computeKmeans") 

batting_info = dget("_battingInfo.dat")

kmeans_obj <- structure(list(
    cluster=integer(0),
    centers=matrix(c(1,2,3,11,22,22,111,222,333,1111,2222,3333,11111,22222,33333,111111,222222,333333), nrow=6, byrow = TRUE),
    totss=100000,
    withinss = c(10,20,30,40,50,60),
    tot.withinss = sum(c(10,20,30,40,50,60)),
    betweenss = 99999,
    size = c(1,2,3,4,5,6),
    iter=2,
    ifault = 0,
    scale=TRUE,
    tableName="batting",
    columns=c('g', 'h', 'r'),
    scaledTableName="baseball.kmeans_test_scaled",
    centroidTableName="baseball.kmeans_test_centroids",
    id="playerid || '-' || stint || '-' || teamid || '-' || yearid",
    idAlias="id",
    whereClause=" WHERE yearid > 2010  "),
    class = c("toakmeans", "kmeans"))


test_that("computeKmeans throws errors", {
  
  expect_error(computeKmeans(NULL),
               "Connection is not valid RODBC object.")
  
  expect_error(computeKmeans(NULL, tableName="fielding", test=TRUE),
               "Must provide tableInfo when test==TRUE")
  
  expect_error(computeKmeans(NULL, tableName='XXXXX', centers="a", tableInfo=batting_info, test=TRUE),
               "Parameter centers must be numeric.")
  
  expect_error(computeKmeans(NULL, tableName='XXXXX', centers=0.1, tableInfo=batting_info, test=TRUE),
               "Number of clusters must be greater or equal to 1.")
  
  expect_error(computeKmeans(NULL, tableName="batting", centers=4, tableInfo=batting_info, id="id", include=c('lgid','playerid'),
                             test=TRUE),
               "Kmeans operates on one or more numeric variables.")
  
  expect_error(computeKmeans(NULL, tableName="batting", centers=4, tableInfo=batting_info, id="id", idAlias="g",
                             test=TRUE),
               "Id alias 'g' can't be one of variable names")
  
  expect_error(computeKmeans(NULL, tableName="batting", centers=4, tableInfo=batting_info, id="id", include=c('g','h','r'),
                             idAlias="g", test=TRUE), "Id alias 'g' can't be one of variable names")
  
  expect_error(computeKmeans(NULL, tableName="batting", centers=matrix(c(10,20,30,40,50,60,70,80), nrow=2, byrow=TRUE),
                             tableInfo = batting_info, id="id", include=c('g','h','r'), test=TRUE),
               "Kmeans received incompatible parameters: dimension of initial cluster centers doesn't match variables: 'g', 'h', 'r'")
  
  expect_error(computeKmeans(NULL, tableName="batting", centers=4, tableInfo=batting_info, id="id", include=c('g','h','r'),
                             aggregates=c("AVG(a) a", "a"), test=TRUE),
               "Check aggregates: at least one missing alias found.")
  
})


test_that("computeClusterSample throws errors", {
  
  expect_error(computeClusterSample(NULL),
               "Connection is not valid RODBC object.")
  
  expect_error(computeClusterSample(NULL, test=TRUE),
               "Kmeans object must be specified.")
  
  expect_error(computeClusterSample(NULL, NULL, test=TRUE),
               "Kmeans object must be specified.")
  
  expect_error(computeClusterSample(NULL, character(1), test=TRUE),
               "Kmeans object must be specified.")
  
  expect_error(computeClusterSample(NULL, data.frame(1:5), test=TRUE),
               "Kmeans object must be specified.")
  
  expect_error(computeClusterSample(NULL, structure(list(cluster=integer(0)), class = c("toakmeans", "kmeans")), test=TRUE),
               "Sample fraction or sample size must be specified.")
})



test_that("computeSilhouette throws errors", {
  
  expect_error(computeSilhouette(NULL),
               "Connection is not valid RODBC object.")
  
  expect_error(computeSilhouette(NULL, test=TRUE),
               "Silhouette table name is required when test=TRUE.")
  
  expect_error(computeSilhouette(NULL, silhouetteTableName='name', test=TRUE),
               "Kmeans object must be specified.")
  
  expect_error(computeSilhouette(NULL, NULL, silhouetteTableName='name', test=TRUE),
               "Kmeans object must be specified.")
  
  expect_error(computeSilhouette(NULL, character(1), silhouetteTableName='name', test=TRUE),
               "Kmeans object must be specified.")
  
  expect_error(computeSilhouette(NULL, data.frame(1:5), silhouetteTableName='name', test=TRUE),
               "Kmeans object must be specified.")
  
  expect_error(computeSilhouette(NULL, km=structure(list(centers=matrix(c(1), nrow=1, byrow = TRUE)),
                                                    class = c("toakmeans", "kmeans")), 
                                 silhouetteTableName = 'name', test=TRUE),
               "Silhouette values are trivial in case of single cluster model.")
})


test_that("computeKmeans SQL is correct", {
  
  expect_equal_normalized(computeKmeans(NULL, "batting", centers=3, iterMax = 25,
                                        tableInfo=batting_info, include=c('g','ab','r','h'),
                                        id="playerid || '-' || stint || '-' || teamid || '-' || yearid", idAlias="id",
                                        aggregates = c("COUNT(*) cnt", "AVG(g) avg_g", "AVG(ab) avg_ab", "AVG(r) avg_r", "AVG(h) avg_h"),
                                        scaledTableName='kmeans_test_scaled', centroidTableName='kmeans_test_centroids', 
                                        schema='baseball', test=TRUE),
                         "-- Data Prep: scale
                          DROP TABLE IF EXISTS baseball.kmeans_test_scaled;
                          CREATE FACT TABLE baseball.kmeans_test_scaled DISTRIBUTE BY HASH(id) AS 
                          SELECT * FROM Scale(
                            ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, ab, g, h, r FROM batting ) AS input PARTITION BY ANY
                            ON (SELECT * FROM ScaleMap (
                              ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, ab, g, h, r FROM batting )
                              InputColumns ('ab', 'g', 'h', 'r')
                              -- MissValue ('OMIT')
                            )) AS STATISTIC DIMENSION
                            Method ('STD')
                            Accumulate('id')
                            GlobalScale ('false')
                            InputColumns ('ab', 'g', 'h', 'r')
                          );
                          --;
                          -- Run k-means;
                          DROP TABLE IF EXISTS baseball.kmeans_test_centroids;
                          SELECT * FROM kmeans(
                            ON (SELECT 1)
                            PARTITION BY 1
                            INPUTTABLE('baseball.kmeans_test_scaled')
                            OUTPUTTABLE('baseball.kmeans_test_centroids')
                            NUMBERK('3')
                            THRESHOLD('0.0395')
                            MAXITERNUM('25')
                          );
                          --;
                          -- Run cluster assignment, cluster stats, and within-cluster sum of squares;
                          SELECT c1.*, c2.withinss  
                          FROM (SELECT clusterid, means, COUNT(*) cnt, AVG(g) avg_g, AVG(ab) avg_ab, AVG(r) avg_r, AVG(h) avg_h          
                                  FROM (SELECT c.clusterid, c.means, d.* 
                                          FROM baseball.kmeans_test_centroids c JOIN 
                                               kmeansplot (
                                                 ON baseball.kmeans_test_scaled PARTITION BY ANY
                                                 ON baseball.kmeans_test_centroids DIMENSION
                                                 centroidsTable('baseball.kmeans_test_centroids')
                                               ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
                                               (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, * 
                                                  FROM batting ) d on (kmp.id = d.id)
                                       ) clustered_data
                                 GROUP BY clusterid, means
                                ) c1 JOIN 
                               ( SELECT 0 clusterid, SUM(distance::double ^ 2) withinss 
                                   FROM VectorDistance(
                                     ON ( SELECT clusterid, id, variable, coalesce(value_double, value_long, value_str::double) value 
                                            FROM unpivot(
                                              ON (SELECT d.* 
                                                    FROM kmeansplot (
                                                      ON baseball.kmeans_test_scaled PARTITION BY ANY
                                                      ON baseball.kmeans_test_centroids DIMENSION
                                                      centroidsTable('baseball.kmeans_test_centroids')
                                                    ) d 
                                                   WHERE clusterid = 0
                                              )
                                              COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                              COLSTOACCUMULATE('id','clusterid')
                                              ATTRIBUTECOLUMNNAME('variable')
                                              VALUECOLUMNNAME('value')
                                              KEEPINPUTCOLUMNTYPES('true')
                                            )
                                         ) AS target PARTITION BY id
                                      ON ( SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                FROM baseball.kmeans_test_centroids
                                WHERE clusterid = 0
                              ) AS ref DIMENSION
                              TARGETIDCOLUMNS('id')
                              TARGETFEATURECOLUMN('variable')
                              TARGETVALUECOLUMN('value')
                              REFIDCOLUMNS('clusterid')
                              REFFEATURECOLUMN('variable')
                              REFVALUECOLUMN('value')
                              MEASURE('Euclidean')
                            )
                            UNION ALL
                            SELECT 1 clusterid, SUM(distance::double ^ 2) withinss FROM VectorDistance(
                              ON (
                                SELECT clusterid, id, variable, coalesce(value_double, value_long, value_str::double) value
                                FROM unpivot(
                                  ON (SELECT d.* 
                                        FROM kmeansplot (
                                          ON baseball.kmeans_test_scaled PARTITION BY ANY
                                          ON baseball.kmeans_test_centroids DIMENSION
                                          centroidsTable('baseball.kmeans_test_centroids')
                                        ) d 
                                      WHERE clusterid = 1
                                  )
                                  COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                  COLSTOACCUMULATE('id','clusterid')
                                  ATTRIBUTECOLUMNNAME('variable')
                                  VALUECOLUMNNAME('value')
                                  KEEPINPUTCOLUMNTYPES('true')
                                )
                              ) AS target PARTITION BY id
                              ON (
                                SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                FROM baseball.kmeans_test_centroids
                                WHERE clusterid = 1
                              ) AS ref DIMENSION
                              TARGETIDCOLUMNS('id')
                              TARGETFEATURECOLUMN('variable')
                              TARGETVALUECOLUMN('value')
                              REFIDCOLUMNS('clusterid')
                              REFFEATURECOLUMN('variable')
                              REFVALUECOLUMN('value')
                              MEASURE('Euclidean')
                            )
                            UNION ALL
                            SELECT 2 clusterid, SUM(distance::double ^ 2) withinss FROM VectorDistance(
                              ON (
                                SELECT clusterid, id, variable, coalesce(value_double, value_long, value_str::double) value
                                FROM unpivot(
                                  ON (SELECT d.* 
                                        FROM kmeansplot (
                                          ON baseball.kmeans_test_scaled PARTITION BY ANY
                                          ON baseball.kmeans_test_centroids DIMENSION
                                          centroidsTable('baseball.kmeans_test_centroids')
                                        ) d 
                                      WHERE clusterid = 2
                                  )
                                  COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                  COLSTOACCUMULATE('id','clusterid')
                                  ATTRIBUTECOLUMNNAME('variable')
                                  VALUECOLUMNNAME('value')
                                  KEEPINPUTCOLUMNTYPES('true')
                                )
                              ) AS target PARTITION BY id
                              ON (
                                SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                FROM baseball.kmeans_test_centroids
                                WHERE clusterid = 2
                              ) AS ref DIMENSION
                              TARGETIDCOLUMNS('id')
                              TARGETFEATURECOLUMN('variable')
                              TARGETVALUECOLUMN('value')
                              REFIDCOLUMNS('clusterid')
                              REFFEATURECOLUMN('variable')
                              REFVALUECOLUMN('value')
                              MEASURE('Euclidean')
                            )
                          ) c2 ON (c1.clusterid = c2.clusterid)
                          ORDER BY clusterid;
                          --;
                          -- Compute Total Sum of Squares;
                          SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
                            ON (SELECT id, variable, coalesce(value_double, value_long, value_str::double) value
                                FROM unpivot(
                                  ON baseball.kmeans_test_scaled
                                  COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                  COLSTOACCUMULATE('id')
                                  ATTRIBUTECOLUMNNAME('variable')
                                  VALUECOLUMNNAME('value')
                                  KEEPINPUTCOLUMNTYPES('true')
                                ) 
                            ) AS target PARTITION BY id
                            ON (SELECT id, variable, value_double
                                FROM unpivot(
                                  ON (SELECT 1 id, 0.0::double ab, 0.0::double g, 0.0::double h, 0.0::double r)
                                  COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                  COLSTOACCUMULATE('id')
                                  ATTRIBUTECOLUMNNAME('variable')
                                  VALUECOLUMNNAME('value')
                                  KEEPINPUTCOLUMNTYPES('true')
                                )
                            ) AS ref DIMENSION
                            TARGETIDCOLUMNS('id')
                            TARGETFEATURECOLUMN('variable')
                            TARGETVALUECOLUMN('value')
                            REFIDCOLUMNS('id')
                            REFFEATURECOLUMN('variable')
                            REFVALUECOLUMN('value_double')
                            MEASURE('Euclidean')
                          );",
                         "kmeans for 3 clusters with aggregates without WHERE clause")
  
  expect_equal_normalized(computeKmeans(NULL, "batting", centers=1, tableInfo=batting_info, 
                                        include=c('g','ab','r','h'),
                                        id="playerid || '-' || stint || '-' || teamid || '-' || yearid", idAlias="id",
                                        scaledTableName='kmeans_test_scaled', centroidTableName='kmeans_test_centroids', 
                                        where="yearid > 2000", test=TRUE),
                         "-- Data Prep: scale
                          DROP TABLE IF EXISTS kmeans_test_scaled;
                          CREATE FACT TABLE kmeans_test_scaled DISTRIBUTE BY HASH(id) AS 
                          SELECT * FROM Scale(
                            ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, ab, g, h, r FROM batting WHERE yearid > 2000  ) AS input PARTITION BY ANY
                            ON (SELECT * FROM ScaleMap (
                              ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, ab, g, h, r FROM batting WHERE yearid > 2000  )
                              InputColumns ('ab', 'g', 'h', 'r')
                              -- MissValue ('OMIT')
                            )) AS STATISTIC DIMENSION
                            Method ('STD')
                            Accumulate('id')
                            GlobalScale ('false')
                            InputColumns ('ab', 'g', 'h', 'r')
                          );
                          --;
                          -- Run k-means;
                          DROP TABLE IF EXISTS kmeans_test_centroids;
                          SELECT * FROM kmeans(
                            ON (SELECT 1)
                            PARTITION BY 1
                            INPUTTABLE('kmeans_test_scaled')
                            OUTPUTTABLE('kmeans_test_centroids')
                            NUMBERK('1')
                            THRESHOLD('0.0395')
                            MAXITERNUM('10')
                          );
                          --;
                          -- Run cluster assignment, cluster stats, and within-cluster sum of squares;
                          SELECT c1.*, c2.withinss  
                            FROM (SELECT clusterid, means, COUNT(*) cnt          
                                    FROM (SELECT c.clusterid, c.means, d.* 
                                      FROM kmeans_test_centroids c JOIN 
                                           kmeansplot (
                                             ON kmeans_test_scaled PARTITION BY ANY
                                             ON kmeans_test_centroids DIMENSION
                                             centroidsTable('kmeans_test_centroids')
                                           ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
                                           (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, * 
                                              FROM batting WHERE yearid > 2000 ) d on (kmp.id = d.id)
                                          ) clustered_data
                                    GROUP BY clusterid, means
                                  ) c1 JOIN ( 
                                  SELECT 0 clusterid, SUM(distance::double ^ 2) withinss 
                                    FROM VectorDistance(
                                      ON ( SELECT clusterid, id, variable, coalesce(value_double, value_long, value_str::double) value
                                             FROM unpivot(
                                               ON (SELECT d.* 
                                                     FROM kmeansplot (
                                                       ON kmeans_test_scaled PARTITION BY ANY
                                                       ON kmeans_test_centroids DIMENSION
                                                       centroidsTable('kmeans_test_centroids')
                                                   ) d 
                                            WHERE clusterid = 0
                                               )
                                               COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                               COLSTOACCUMULATE('id','clusterid')
                                               ATTRIBUTECOLUMNNAME('variable')
                                               VALUECOLUMNNAME('value')
                                               KEEPINPUTCOLUMNTYPES('true')
                                            )
                                      ) AS target PARTITION BY id
                                      ON (
                                          SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                            FROM kmeans_test_centroids
                                           WHERE clusterid = 0
                                      ) AS ref DIMENSION
                                      TARGETIDCOLUMNS('id')
                                      TARGETFEATURECOLUMN('variable')
                                      TARGETVALUECOLUMN('value')
                                      REFIDCOLUMNS('clusterid')
                                      REFFEATURECOLUMN('variable')
                                      REFVALUECOLUMN('value')
                                      MEASURE('Euclidean')
                                )
                              ) c2 ON (c1.clusterid = c2.clusterid)
                              ORDER BY clusterid;
                          --;
                          -- Compute Total Sum of Squares;
                          SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
                            ON (SELECT id, variable, coalesce(value_double, value_long, value_str::double) value
                                  FROM unpivot(
                                    ON kmeans_test_scaled
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('id')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                  ) 
                            ) AS target PARTITION BY id
                            ON (SELECT id, variable, value_double
                                  FROM unpivot(
                                    ON (SELECT 1 id, 0.0::double ab, 0.0::double g, 0.0::double h, 0.0::double r)
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('id')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                )
                            ) AS ref DIMENSION
                            TARGETIDCOLUMNS('id')
                            TARGETFEATURECOLUMN('variable')
                            TARGETVALUECOLUMN('value')
                            REFIDCOLUMNS('id')
                            REFFEATURECOLUMN('variable')
                            REFVALUECOLUMN('value_double')
                            MEASURE('Euclidean')
                          );",
                          "kmeans for 1 cluster with no aggregates with WHERE clause")
  
  
  expect_equal_normalized(computeKmeans(NULL, "batting", centers=1, tableInfo=batting_info, 
                                        include=c('g','ab','r','h'),
                                        id="playerid || '-' || stint || '-' || teamid || '-' || yearid", 
                                        scaledTableName='kmeans_test_scaled', centroidTableName='kmeans_test_centroids', 
                                        where="yearid > 2000", test=TRUE),
                         "-- Data Prep: scale
                          DROP TABLE IF EXISTS kmeans_test_scaled;
                          CREATE FACT TABLE kmeans_test_scaled DISTRIBUTE BY HASH(playerid_stint_teamid_yearid) AS 
                          SELECT * FROM Scale(
                            ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, 
                                       ab, g, h, r FROM batting WHERE yearid > 2000  
                            ) AS input PARTITION BY ANY
                            ON (SELECT * FROM ScaleMap (
                              ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, 
                                         ab, g, h, r FROM batting WHERE yearid > 2000  )
                              InputColumns ('ab', 'g', 'h', 'r')
                              -- MissValue ('OMIT')
                            )) AS STATISTIC DIMENSION
                            Method ('STD')
                            Accumulate('playerid_stint_teamid_yearid')
                            GlobalScale ('false')
                            InputColumns ('ab', 'g', 'h', 'r')
                          );
                          --;
                          -- Run k-means;
                          DROP TABLE IF EXISTS kmeans_test_centroids;
                          SELECT * FROM kmeans(
                            ON (SELECT 1)
                            PARTITION BY 1
                            INPUTTABLE('kmeans_test_scaled')
                            OUTPUTTABLE('kmeans_test_centroids')
                            NUMBERK('1')
                            THRESHOLD('0.0395')
                            MAXITERNUM('10')
                          );
                          --;
                          -- Run cluster assignment, cluster stats, and within-cluster sum of squares;
                          SELECT c1.*, c2.withinss  
                            FROM (SELECT clusterid, means, COUNT(*) cnt          
                                    FROM (SELECT c.clusterid, c.means, d.* 
                                      FROM kmeans_test_centroids c JOIN 
                                           kmeansplot (
                                             ON kmeans_test_scaled PARTITION BY ANY
                                             ON kmeans_test_centroids DIMENSION
                                             centroidsTable('kmeans_test_centroids')
                                           ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
                                           (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, *
                                              FROM batting WHERE yearid > 2000 ) d on (kmp.playerid_stint_teamid_yearid = d.playerid_stint_teamid_yearid)
                                          ) clustered_data
                                    GROUP BY clusterid, means
                                  ) c1 JOIN ( 
                                  SELECT 0 clusterid, SUM(distance::double ^ 2) withinss 
                                    FROM VectorDistance(
                                      ON ( SELECT clusterid, playerid_stint_teamid_yearid, variable, coalesce(value_double, value_long, value_str::double) value
                                             FROM unpivot(
                                               ON (SELECT d.* 
                                                     FROM kmeansplot (
                                                       ON kmeans_test_scaled PARTITION BY ANY
                                                       ON kmeans_test_centroids DIMENSION
                                                       centroidsTable('kmeans_test_centroids')
                                                   ) d 
                                            WHERE clusterid = 0
                                               )
                                               COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                               COLSTOACCUMULATE('playerid_stint_teamid_yearid','clusterid')
                                               ATTRIBUTECOLUMNNAME('variable')
                                               VALUECOLUMNNAME('value')
                                               KEEPINPUTCOLUMNTYPES('true')
                                            )
                                      ) AS target PARTITION BY playerid_stint_teamid_yearid
                                      ON (
                                          SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                            FROM kmeans_test_centroids
                                           WHERE clusterid = 0
                                      ) AS ref DIMENSION
                                      TARGETIDCOLUMNS('playerid_stint_teamid_yearid')
                                      TARGETFEATURECOLUMN('variable')
                                      TARGETVALUECOLUMN('value')
                                      REFIDCOLUMNS('clusterid')
                                      REFFEATURECOLUMN('variable')
                                      REFVALUECOLUMN('value')
                                      MEASURE('Euclidean')
                                )
                              ) c2 ON (c1.clusterid = c2.clusterid)
                              ORDER BY clusterid;
                          --;
                          -- Compute Total Sum of Squares;
                          SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
                            ON (SELECT playerid_stint_teamid_yearid, variable, coalesce(value_double, value_long, value_str::double) value
                                  FROM unpivot(
                                    ON kmeans_test_scaled
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('playerid_stint_teamid_yearid')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                  ) 
                            ) AS target PARTITION BY playerid_stint_teamid_yearid
                            ON (SELECT id, variable, value_double
                                  FROM unpivot(
                                    ON (SELECT 1 id, 0.0::double ab, 0.0::double g, 0.0::double h, 0.0::double r)
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('id')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                )
                            ) AS ref DIMENSION
                            TARGETIDCOLUMNS('playerid_stint_teamid_yearid')
                            TARGETFEATURECOLUMN('variable')
                            TARGETVALUECOLUMN('value')
                            REFIDCOLUMNS('id')
                            REFFEATURECOLUMN('variable')
                            REFVALUECOLUMN('value_double')
                            MEASURE('Euclidean')
                          );",
                          "kmeans for 1 cluster with no aggregates with default id and WHERE clause")
  
  expect_equal_normalized(computeKmeans(NULL, "batting", centers=1, tableInfo=batting_info, 
                                        include=c('g','ab','r','h'),
                                        id="playerid || '-' || stint || '-' || teamid || '-' || yearid", 
                                        aggregates = NULL,
                                        scaledTableName='kmeans_test_scaled', centroidTableName='kmeans_test_centroids', 
                                        where="yearid > 2000", test=TRUE),
                         "-- Data Prep: scale
                          DROP TABLE IF EXISTS kmeans_test_scaled;
                          CREATE FACT TABLE kmeans_test_scaled DISTRIBUTE BY HASH(playerid_stint_teamid_yearid) AS 
                          SELECT * FROM Scale(
                            ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, 
                                       ab, g, h, r FROM batting WHERE yearid > 2000  
                            ) AS input PARTITION BY ANY
                            ON (SELECT * FROM ScaleMap (
                              ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, 
                                         ab, g, h, r FROM batting WHERE yearid > 2000  )
                              InputColumns ('ab', 'g', 'h', 'r')
                              -- MissValue ('OMIT')
                            )) AS STATISTIC DIMENSION
                            Method ('STD')
                            Accumulate('playerid_stint_teamid_yearid')
                            GlobalScale ('false')
                            InputColumns ('ab', 'g', 'h', 'r')
                          );
                          --;
                          -- Run k-means;
                          DROP TABLE IF EXISTS kmeans_test_centroids;
                          SELECT * FROM kmeans(
                            ON (SELECT 1)
                            PARTITION BY 1
                            INPUTTABLE('kmeans_test_scaled')
                            OUTPUTTABLE('kmeans_test_centroids')
                            NUMBERK('1')
                            THRESHOLD('0.0395')
                            MAXITERNUM('10')
                          );
                          --;
                          -- Run cluster assignment, cluster stats, and within-cluster sum of squares;
                          SELECT c1.*, c2.withinss  
                            FROM (SELECT clusterid, means, COUNT(*) cnt          
                                    FROM (SELECT c.clusterid, c.means, d.* 
                                      FROM kmeans_test_centroids c JOIN 
                                           kmeansplot (
                                             ON kmeans_test_scaled PARTITION BY ANY
                                             ON kmeans_test_centroids DIMENSION
                                             centroidsTable('kmeans_test_centroids')
                                           ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
                                           (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, *
                                              FROM batting WHERE yearid > 2000 ) d on (kmp.playerid_stint_teamid_yearid = d.playerid_stint_teamid_yearid)
                                          ) clustered_data
                                    GROUP BY clusterid, means
                                  ) c1 JOIN ( 
                                  SELECT 0 clusterid, SUM(distance::double ^ 2) withinss 
                                    FROM VectorDistance(
                                      ON ( SELECT clusterid, playerid_stint_teamid_yearid, variable, coalesce(value_double, value_long, value_str::double) value
                                             FROM unpivot(
                                               ON (SELECT d.* 
                                                     FROM kmeansplot (
                                                       ON kmeans_test_scaled PARTITION BY ANY
                                                       ON kmeans_test_centroids DIMENSION
                                                       centroidsTable('kmeans_test_centroids')
                                                   ) d 
                                            WHERE clusterid = 0
                                               )
                                               COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                               COLSTOACCUMULATE('playerid_stint_teamid_yearid','clusterid')
                                               ATTRIBUTECOLUMNNAME('variable')
                                               VALUECOLUMNNAME('value')
                                               KEEPINPUTCOLUMNTYPES('true')
                                            )
                                      ) AS target PARTITION BY playerid_stint_teamid_yearid
                                      ON (
                                          SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                            FROM kmeans_test_centroids
                                           WHERE clusterid = 0
                                      ) AS ref DIMENSION
                                      TARGETIDCOLUMNS('playerid_stint_teamid_yearid')
                                      TARGETFEATURECOLUMN('variable')
                                      TARGETVALUECOLUMN('value')
                                      REFIDCOLUMNS('clusterid')
                                      REFFEATURECOLUMN('variable')
                                      REFVALUECOLUMN('value')
                                      MEASURE('Euclidean')
                                )
                              ) c2 ON (c1.clusterid = c2.clusterid)
                              ORDER BY clusterid;
                          --;
                          -- Compute Total Sum of Squares;
                          SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
                            ON (SELECT playerid_stint_teamid_yearid, variable, coalesce(value_double, value_long, value_str::double) value
                                  FROM unpivot(
                                    ON kmeans_test_scaled
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('playerid_stint_teamid_yearid')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                  ) 
                            ) AS target PARTITION BY playerid_stint_teamid_yearid
                            ON (SELECT id, variable, value_double
                                  FROM unpivot(
                                    ON (SELECT 1 id, 0.0::double ab, 0.0::double g, 0.0::double h, 0.0::double r)
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('id')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                )
                            ) AS ref DIMENSION
                            TARGETIDCOLUMNS('playerid_stint_teamid_yearid')
                            TARGETFEATURECOLUMN('variable')
                            TARGETVALUECOLUMN('value')
                            REFIDCOLUMNS('id')
                            REFFEATURECOLUMN('variable')
                            REFVALUECOLUMN('value_double')
                            MEASURE('Euclidean')
                          );",
                          "kmeans for 1 cluster with NULL aggregates with default id and WHERE clause")
    
  expect_equal_normalized(computeKmeans(NULL, "batting", centers=1, tableInfo=batting_info, 
                                        include=c('g','ab','r','h'),
                                        id="playerid || '-' || stint || '-' || teamid || '-' || yearid", 
                                        aggregates = character(0),
                                        scaledTableName='kmeans_test_scaled', centroidTableName='kmeans_test_centroids', 
                                        where="yearid > 2000", test=TRUE),
                         "-- Data Prep: scale
                          DROP TABLE IF EXISTS kmeans_test_scaled;
                          CREATE FACT TABLE kmeans_test_scaled DISTRIBUTE BY HASH(playerid_stint_teamid_yearid) AS 
                          SELECT * FROM Scale(
                            ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, 
                                       ab, g, h, r FROM batting WHERE yearid > 2000  
                            ) AS input PARTITION BY ANY
                            ON (SELECT * FROM ScaleMap (
                              ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, 
                                         ab, g, h, r FROM batting WHERE yearid > 2000  )
                              InputColumns ('ab', 'g', 'h', 'r')
                              -- MissValue ('OMIT')
                            )) AS STATISTIC DIMENSION
                            Method ('STD')
                            Accumulate('playerid_stint_teamid_yearid')
                            GlobalScale ('false')
                            InputColumns ('ab', 'g', 'h', 'r')
                          );
                          --;
                          -- Run k-means;
                          DROP TABLE IF EXISTS kmeans_test_centroids;
                          SELECT * FROM kmeans(
                            ON (SELECT 1)
                            PARTITION BY 1
                            INPUTTABLE('kmeans_test_scaled')
                            OUTPUTTABLE('kmeans_test_centroids')
                            NUMBERK('1')
                            THRESHOLD('0.0395')
                            MAXITERNUM('10')
                          );
                          --;
                          -- Run cluster assignment, cluster stats, and within-cluster sum of squares;
                          SELECT c1.*, c2.withinss  
                            FROM (SELECT clusterid, means, COUNT(*) cnt          
                                    FROM (SELECT c.clusterid, c.means, d.* 
                                      FROM kmeans_test_centroids c JOIN 
                                           kmeansplot (
                                             ON kmeans_test_scaled PARTITION BY ANY
                                             ON kmeans_test_centroids DIMENSION
                                             centroidsTable('kmeans_test_centroids')
                                           ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
                                           (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, *
                                              FROM batting WHERE yearid > 2000 ) d on (kmp.playerid_stint_teamid_yearid = d.playerid_stint_teamid_yearid)
                                          ) clustered_data
                                    GROUP BY clusterid, means
                                  ) c1 JOIN ( 
                                  SELECT 0 clusterid, SUM(distance::double ^ 2) withinss 
                                    FROM VectorDistance(
                                      ON ( SELECT clusterid, playerid_stint_teamid_yearid, variable, coalesce(value_double, value_long, value_str::double) value
                                             FROM unpivot(
                                               ON (SELECT d.* 
                                                     FROM kmeansplot (
                                                       ON kmeans_test_scaled PARTITION BY ANY
                                                       ON kmeans_test_centroids DIMENSION
                                                       centroidsTable('kmeans_test_centroids')
                                                   ) d 
                                            WHERE clusterid = 0
                                               )
                                               COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                               COLSTOACCUMULATE('playerid_stint_teamid_yearid','clusterid')
                                               ATTRIBUTECOLUMNNAME('variable')
                                               VALUECOLUMNNAME('value')
                                               KEEPINPUTCOLUMNTYPES('true')
                                            )
                                      ) AS target PARTITION BY playerid_stint_teamid_yearid
                                      ON (
                                          SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                            FROM kmeans_test_centroids
                                           WHERE clusterid = 0
                                      ) AS ref DIMENSION
                                      TARGETIDCOLUMNS('playerid_stint_teamid_yearid')
                                      TARGETFEATURECOLUMN('variable')
                                      TARGETVALUECOLUMN('value')
                                      REFIDCOLUMNS('clusterid')
                                      REFFEATURECOLUMN('variable')
                                      REFVALUECOLUMN('value')
                                      MEASURE('Euclidean')
                                )
                              ) c2 ON (c1.clusterid = c2.clusterid)
                              ORDER BY clusterid;
                          --;
                          -- Compute Total Sum of Squares;
                          SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
                            ON (SELECT playerid_stint_teamid_yearid, variable, coalesce(value_double, value_long, value_str::double) value
                                  FROM unpivot(
                                    ON kmeans_test_scaled
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('playerid_stint_teamid_yearid')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                  ) 
                            ) AS target PARTITION BY playerid_stint_teamid_yearid
                            ON (SELECT id, variable, value_double
                                  FROM unpivot(
                                    ON (SELECT 1 id, 0.0::double ab, 0.0::double g, 0.0::double h, 0.0::double r)
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('id')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                )
                            ) AS ref DIMENSION
                            TARGETIDCOLUMNS('playerid_stint_teamid_yearid')
                            TARGETFEATURECOLUMN('variable')
                            TARGETVALUECOLUMN('value')
                            REFIDCOLUMNS('id')
                            REFFEATURECOLUMN('variable')
                            REFVALUECOLUMN('value_double')
                            MEASURE('Euclidean')
                          );",
                          "kmeans for 1 cluster with 0-length aggregates with default id and WHERE clause")
  
  expect_equal_normalized(computeKmeans(NULL, "batting", centers=1, tableInfo=batting_info, 
                                        include=c('g','ab','r','h'),
                                        id="playerid || '-' || stint || '-' || teamid || '-' || yearid", 
                                        aggregates = c("AVG( a ) a"),
                                        scaledTableName='kmeans_test_scaled', centroidTableName='kmeans_test_centroids', 
                                        where="yearid > 2000", test=TRUE),
                         "-- Data Prep: scale
                          DROP TABLE IF EXISTS kmeans_test_scaled;
                          CREATE FACT TABLE kmeans_test_scaled DISTRIBUTE BY HASH(playerid_stint_teamid_yearid) AS 
                          SELECT * FROM Scale(
                            ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, 
                                       ab, g, h, r FROM batting WHERE yearid > 2000  
                            ) AS input PARTITION BY ANY
                            ON (SELECT * FROM ScaleMap (
                              ON (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, 
                                         ab, g, h, r FROM batting WHERE yearid > 2000  )
                              InputColumns ('ab', 'g', 'h', 'r')
                              -- MissValue ('OMIT')
                            )) AS STATISTIC DIMENSION
                            Method ('STD')
                            Accumulate('playerid_stint_teamid_yearid')
                            GlobalScale ('false')
                            InputColumns ('ab', 'g', 'h', 'r')
                          );
                          --;
                          -- Run k-means;
                          DROP TABLE IF EXISTS kmeans_test_centroids;
                          SELECT * FROM kmeans(
                            ON (SELECT 1)
                            PARTITION BY 1
                            INPUTTABLE('kmeans_test_scaled')
                            OUTPUTTABLE('kmeans_test_centroids')
                            NUMBERK('1')
                            THRESHOLD('0.0395')
                            MAXITERNUM('10')
                          );
                          --;
                          -- Run cluster assignment, cluster stats, and within-cluster sum of squares;
                          SELECT c1.*, c2.withinss  
                            FROM (SELECT clusterid, means, AVG( a ) a, COUNT(*) cnt          
                                    FROM (SELECT c.clusterid, c.means, d.* 
                                      FROM kmeans_test_centroids c JOIN 
                                           kmeansplot (
                                             ON kmeans_test_scaled PARTITION BY ANY
                                             ON kmeans_test_centroids DIMENSION
                                             centroidsTable('kmeans_test_centroids')
                                           ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
                                           (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid playerid_stint_teamid_yearid, *
                                              FROM batting WHERE yearid > 2000 ) d on (kmp.playerid_stint_teamid_yearid = d.playerid_stint_teamid_yearid)
                                          ) clustered_data
                                    GROUP BY clusterid, means
                                  ) c1 JOIN ( 
                                  SELECT 0 clusterid, SUM(distance::double ^ 2) withinss 
                                    FROM VectorDistance(
                                      ON ( SELECT clusterid, playerid_stint_teamid_yearid, variable, coalesce(value_double, value_long, value_str::double) value
                                             FROM unpivot(
                                               ON (SELECT d.* 
                                                     FROM kmeansplot (
                                                       ON kmeans_test_scaled PARTITION BY ANY
                                                       ON kmeans_test_centroids DIMENSION
                                                       centroidsTable('kmeans_test_centroids')
                                                   ) d 
                                            WHERE clusterid = 0
                                               )
                                               COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                               COLSTOACCUMULATE('playerid_stint_teamid_yearid','clusterid')
                                               ATTRIBUTECOLUMNNAME('variable')
                                               VALUECOLUMNNAME('value')
                                               KEEPINPUTCOLUMNTYPES('true')
                                            )
                                      ) AS target PARTITION BY playerid_stint_teamid_yearid
                                      ON (
                                          SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                            FROM kmeans_test_centroids
                                           WHERE clusterid = 0
                                      ) AS ref DIMENSION
                                      TARGETIDCOLUMNS('playerid_stint_teamid_yearid')
                                      TARGETFEATURECOLUMN('variable')
                                      TARGETVALUECOLUMN('value')
                                      REFIDCOLUMNS('clusterid')
                                      REFFEATURECOLUMN('variable')
                                      REFVALUECOLUMN('value')
                                      MEASURE('Euclidean')
                                )
                              ) c2 ON (c1.clusterid = c2.clusterid)
                              ORDER BY clusterid;
                          --;
                          -- Compute Total Sum of Squares;
                          SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
                            ON (SELECT playerid_stint_teamid_yearid, variable, coalesce(value_double, value_long, value_str::double) value
                                  FROM unpivot(
                                    ON kmeans_test_scaled
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('playerid_stint_teamid_yearid')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                  ) 
                            ) AS target PARTITION BY playerid_stint_teamid_yearid
                            ON (SELECT id, variable, value_double
                                  FROM unpivot(
                                    ON (SELECT 1 id, 0.0::double ab, 0.0::double g, 0.0::double h, 0.0::double r)
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('id')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                )
                            ) AS ref DIMENSION
                            TARGETIDCOLUMNS('playerid_stint_teamid_yearid')
                            TARGETFEATURECOLUMN('variable')
                            TARGETVALUECOLUMN('value')
                            REFIDCOLUMNS('id')
                            REFFEATURECOLUMN('variable')
                            REFVALUECOLUMN('value_double')
                            MEASURE('Euclidean')
                          );",
                          "kmeans for 1 cluster with aggregates without COUNT with default id and WHERE clause")
  
  expect_equal_normalized(computeKmeans(NULL, "batting", centers=1, tableInfo=batting_info, 
                                        include=c('g','ab','r','h'), scale = FALSE,
                                        id="playerid || '-' || stint || '-' || teamid || '-' || yearid", idAlias="id",
                                        scaledTableName='kmeans_test_scaled', centroidTableName='kmeans_test_centroids', 
                                        where="yearid > 2000", test=TRUE),
                         "-- Data Prep: scale
DROP TABLE IF EXISTS kmeans_test_scaled;
CREATE FACT TABLE kmeans_test_scaled DISTRIBUTE BY HASH(id) AS 
       SELECT * FROM (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, ab, g, h, r FROM batting WHERE yearid > 2000  ) d
      WHERE id IS NOT NULL AND ab IS NOT NULL AND g IS NOT NULL AND h IS NOT NULL AND r IS NOT NULL;
--;
-- Run k-means;
DROP TABLE IF EXISTS kmeans_test_centroids;
SELECT * FROM kmeans(
      ON (SELECT 1)
      PARTITION BY 1
      INPUTTABLE('kmeans_test_scaled')
      OUTPUTTABLE('kmeans_test_centroids')
   NUMBERK('1')
      THRESHOLD('0.0395')
      MAXITERNUM('10')
    );
--;
-- Run cluster assignment, cluster stats, and within-cluster sum of squares;
SELECT c1.*, c2.withinss  
       FROM (SELECT clusterid, means, COUNT(*) cnt          FROM (SELECT c.clusterid, c.means, d.* 
      FROM kmeans_test_centroids c JOIN 
    kmeansplot (
      ON kmeans_test_scaled PARTITION BY ANY
      ON kmeans_test_centroids DIMENSION
      centroidsTable('kmeans_test_centroids')
    ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
    (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, * FROM batting WHERE yearid > 2000  ) d on (kmp.id = d.id)
                    ) clustered_data
              GROUP BY clusterid, means
            ) c1 JOIN ( 
            SELECT 0 clusterid, SUM(distance::double ^ 2) withinss FROM VectorDistance(
       ON (
         SELECT clusterid, id, variable, coalesce(value_double, value_long, value_str::double) value
           FROM unpivot(
                  ON (SELECT d.* 
                        FROM kmeansplot (
                               ON kmeans_test_scaled PARTITION BY ANY
                               ON kmeans_test_centroids DIMENSION
                               centroidsTable('kmeans_test_centroids')
                             ) d 
                       WHERE clusterid = 0
                  )
                  COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                  COLSTOACCUMULATE('id','clusterid')
                  ATTRIBUTECOLUMNNAME('variable')
                  VALUECOLUMNNAME('value')
                  KEEPINPUTCOLUMNTYPES('true')
                )
       ) AS target PARTITION BY id
       ON (
         SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
           FROM kmeans_test_centroids WHERE clusterid = 0
       ) AS ref DIMENSION
       TARGETIDCOLUMNS('id')
       TARGETFEATURECOLUMN('variable')
       TARGETVALUECOLUMN('value')
       REFIDCOLUMNS('clusterid')
       REFFEATURECOLUMN('variable')
       REFVALUECOLUMN('value')
       MEASURE('Euclidean')
     )
            ) c2 ON (c1.clusterid = c2.clusterid)
      ORDER BY clusterid;
--;
-- Compute Total Sum of Squares;
SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
       ON (SELECT id, variable, coalesce(value_double, value_long, value_str::double) value
             FROM unpivot(
               ON kmeans_test_scaled
               COLSTOUNPIVOT('ab', 'g', 'h', 'r')
               COLSTOACCUMULATE('id')
               ATTRIBUTECOLUMNNAME('variable')
               VALUECOLUMNNAME('value')
               KEEPINPUTCOLUMNTYPES('true')
             ) 
       ) AS target PARTITION BY id
       ON (SELECT id, variable, value_double
             FROM unpivot(
               ON (SELECT 1 id, AVG(ab)::double ab, AVG(g)::double g, AVG(h)::double h, AVG(r)::double r FROM kmeans_test_scaled)
               COLSTOUNPIVOT('ab', 'g', 'h', 'r')
               COLSTOACCUMULATE('id')
               ATTRIBUTECOLUMNNAME('variable')
               VALUECOLUMNNAME('value')
               KEEPINPUTCOLUMNTYPES('true')
             )
       ) AS ref DIMENSION
       TARGETIDCOLUMNS('id')
       TARGETFEATURECOLUMN('variable')
       TARGETVALUECOLUMN('value')
       REFIDCOLUMNS('id')
       REFFEATURECOLUMN('variable')
       REFVALUECOLUMN('value_double')
       MEASURE('Euclidean')
     );",
                          "kmeans for 1 cluster not scaled with no aggregates with WHERE clause")
  
  
  expect_equal_normalized(computeKmeans(NULL, "batting", centers=1, tableInfo=batting_info, 
                                        include=c('g','ab','r','h'),
                                        id="playerid",
                                        scaledTableName='kmeans_test_scaled', centroidTableName='kmeans_test_centroids', 
                                        where="yearid > 2000", test=TRUE),
                         "-- Data Prep: scale
                          DROP TABLE IF EXISTS kmeans_test_scaled;
                          CREATE FACT TABLE kmeans_test_scaled DISTRIBUTE BY HASH(__playerid__) AS 
                          SELECT * FROM Scale(
                            ON (SELECT playerid __playerid__, 
                                       ab, g, h, r FROM batting WHERE yearid > 2000  
                            ) AS input PARTITION BY ANY
                            ON (SELECT * FROM ScaleMap (
                              ON (SELECT playerid __playerid__, 
                                         ab, g, h, r FROM batting WHERE yearid > 2000  )
                              InputColumns ('ab', 'g', 'h', 'r')
                              -- MissValue ('OMIT')
                            )) AS STATISTIC DIMENSION
                            Method ('STD')
                            Accumulate('__playerid__')
                            GlobalScale ('false')
                            InputColumns ('ab', 'g', 'h', 'r')
                          );
                          --;
                          -- Run k-means;
                          DROP TABLE IF EXISTS kmeans_test_centroids;
                          SELECT * FROM kmeans(
                            ON (SELECT 1)
                            PARTITION BY 1
                            INPUTTABLE('kmeans_test_scaled')
                            OUTPUTTABLE('kmeans_test_centroids')
                            NUMBERK('1')
                            THRESHOLD('0.0395')
                            MAXITERNUM('10')
                          );
                          --;
                          -- Run cluster assignment, cluster stats, and within-cluster sum of squares;
                          SELECT c1.*, c2.withinss  
                            FROM (SELECT clusterid, means, COUNT(*) cnt          
                                    FROM (SELECT c.clusterid, c.means, d.* 
                                      FROM kmeans_test_centroids c JOIN 
                                           kmeansplot (
                                             ON kmeans_test_scaled PARTITION BY ANY
                                             ON kmeans_test_centroids DIMENSION
                                             centroidsTable('kmeans_test_centroids')
                                           ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
                                           (SELECT playerid __playerid__, *
                                              FROM batting WHERE yearid > 2000 ) d on (kmp.__playerid__ = d.__playerid__)
                                          ) clustered_data
                                    GROUP BY clusterid, means
                                  ) c1 JOIN ( 
                                  SELECT 0 clusterid, SUM(distance::double ^ 2) withinss 
                                    FROM VectorDistance(
                                      ON ( SELECT clusterid, __playerid__, variable, coalesce(value_double, value_long, value_str::double) value
                                             FROM unpivot(
                                               ON (SELECT d.* 
                                                     FROM kmeansplot (
                                                       ON kmeans_test_scaled PARTITION BY ANY
                                                       ON kmeans_test_centroids DIMENSION
                                                       centroidsTable('kmeans_test_centroids')
                                                   ) d 
                                            WHERE clusterid = 0
                                               )
                                               COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                               COLSTOACCUMULATE('__playerid__','clusterid')
                                               ATTRIBUTECOLUMNNAME('variable')
                                               VALUECOLUMNNAME('value')
                                               KEEPINPUTCOLUMNTYPES('true')
                                            )
                                      ) AS target PARTITION BY __playerid__
                                      ON (
                                          SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                            FROM kmeans_test_centroids
                                           WHERE clusterid = 0
                                      ) AS ref DIMENSION
                                      TARGETIDCOLUMNS('__playerid__')
                                      TARGETFEATURECOLUMN('variable')
                                      TARGETVALUECOLUMN('value')
                                      REFIDCOLUMNS('clusterid')
                                      REFFEATURECOLUMN('variable')
                                      REFVALUECOLUMN('value')
                                      MEASURE('Euclidean')
                                )
                              ) c2 ON (c1.clusterid = c2.clusterid)
                              ORDER BY clusterid;
                          --;
                          -- Compute Total Sum of Squares;
                          SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
                            ON (SELECT __playerid__, variable, coalesce(value_double, value_long, value_str::double) value
                                  FROM unpivot(
                                    ON kmeans_test_scaled
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('__playerid__')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                  ) 
                            ) AS target PARTITION BY __playerid__
                            ON (SELECT id, variable, value_double
                                  FROM unpivot(
                                    ON (SELECT 1 id, 0.0::double ab, 0.0::double g, 0.0::double h, 0.0::double r)
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('id')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                )
                            ) AS ref DIMENSION
                            TARGETIDCOLUMNS('__playerid__')
                            TARGETFEATURECOLUMN('variable')
                            TARGETVALUECOLUMN('value')
                            REFIDCOLUMNS('id')
                            REFFEATURECOLUMN('variable')
                            REFVALUECOLUMN('value_double')
                            MEASURE('Euclidean')
                          );",
                          "kmeans for 1 cluster with no aggregates with id one of the columns and WHERE clause")
  
  expect_equal_normalized(computeKmeans(NULL, "batting", centers=1, tableInfo=batting_info, 
                                        include=c('g','ab','r','h'),
                                        id="playerid", idAlias="playerid",
                                        scaledTableName='kmeans_test_scaled', centroidTableName='kmeans_test_centroids', 
                                        where="yearid > 2000", test=TRUE),
                         "-- Data Prep: scale
                          DROP TABLE IF EXISTS kmeans_test_scaled;
                          CREATE FACT TABLE kmeans_test_scaled DISTRIBUTE BY HASH(__playerid__) AS 
                          SELECT * FROM Scale(
                            ON (SELECT playerid __playerid__, 
                                       ab, g, h, r FROM batting WHERE yearid > 2000  
                            ) AS input PARTITION BY ANY
                            ON (SELECT * FROM ScaleMap (
                              ON (SELECT playerid __playerid__, 
                                         ab, g, h, r FROM batting WHERE yearid > 2000  )
                              InputColumns ('ab', 'g', 'h', 'r')
                              -- MissValue ('OMIT')
                            )) AS STATISTIC DIMENSION
                            Method ('STD')
                            Accumulate('__playerid__')
                            GlobalScale ('false')
                            InputColumns ('ab', 'g', 'h', 'r')
                          );
                          --;
                          -- Run k-means;
                          DROP TABLE IF EXISTS kmeans_test_centroids;
                          SELECT * FROM kmeans(
                            ON (SELECT 1)
                            PARTITION BY 1
                            INPUTTABLE('kmeans_test_scaled')
                            OUTPUTTABLE('kmeans_test_centroids')
                            NUMBERK('1')
                            THRESHOLD('0.0395')
                            MAXITERNUM('10')
                          );
                          --;
                          -- Run cluster assignment, cluster stats, and within-cluster sum of squares;
                          SELECT c1.*, c2.withinss  
                            FROM (SELECT clusterid, means, COUNT(*) cnt          
                                    FROM (SELECT c.clusterid, c.means, d.* 
                                      FROM kmeans_test_centroids c JOIN 
                                           kmeansplot (
                                             ON kmeans_test_scaled PARTITION BY ANY
                                             ON kmeans_test_centroids DIMENSION
                                             centroidsTable('kmeans_test_centroids')
                                           ) kmp ON (c.clusterid = kmp.clusterid) JOIN 
                                           (SELECT playerid __playerid__, *
                                              FROM batting WHERE yearid > 2000 ) d on (kmp.__playerid__ = d.__playerid__)
                                          ) clustered_data
                                    GROUP BY clusterid, means
                                  ) c1 JOIN ( 
                                  SELECT 0 clusterid, SUM(distance::double ^ 2) withinss 
                                    FROM VectorDistance(
                                      ON ( SELECT clusterid, __playerid__, variable, coalesce(value_double, value_long, value_str::double) value
                                             FROM unpivot(
                                               ON (SELECT d.* 
                                                     FROM kmeansplot (
                                                       ON kmeans_test_scaled PARTITION BY ANY
                                                       ON kmeans_test_centroids DIMENSION
                                                       centroidsTable('kmeans_test_centroids')
                                                   ) d 
                                            WHERE clusterid = 0
                                               )
                                               COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                               COLSTOACCUMULATE('__playerid__','clusterid')
                                               ATTRIBUTECOLUMNNAME('variable')
                                               VALUECOLUMNNAME('value')
                                               KEEPINPUTCOLUMNTYPES('true')
                                            )
                                      ) AS target PARTITION BY __playerid__
                                      ON (
                                          SELECT *, regexp_split_to_table(means, ' ')::numeric value, regexp_split_to_table('ab, g, h, r', ', ') variable 
                                            FROM kmeans_test_centroids
                                           WHERE clusterid = 0
                                      ) AS ref DIMENSION
                                      TARGETIDCOLUMNS('__playerid__')
                                      TARGETFEATURECOLUMN('variable')
                                      TARGETVALUECOLUMN('value')
                                      REFIDCOLUMNS('clusterid')
                                      REFFEATURECOLUMN('variable')
                                      REFVALUECOLUMN('value')
                                      MEASURE('Euclidean')
                                )
                              ) c2 ON (c1.clusterid = c2.clusterid)
                              ORDER BY clusterid;
                          --;
                          -- Compute Total Sum of Squares;
                          SELECT SUM(distance::double ^ 2) totss FROM VectorDistance(
                            ON (SELECT __playerid__, variable, coalesce(value_double, value_long, value_str::double) value
                                  FROM unpivot(
                                    ON kmeans_test_scaled
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('__playerid__')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                  ) 
                            ) AS target PARTITION BY __playerid__
                            ON (SELECT id, variable, value_double
                                  FROM unpivot(
                                    ON (SELECT 1 id, 0.0::double ab, 0.0::double g, 0.0::double h, 0.0::double r)
                                    COLSTOUNPIVOT('ab', 'g', 'h', 'r')
                                    COLSTOACCUMULATE('id')
                                    ATTRIBUTECOLUMNNAME('variable')
                                    VALUECOLUMNNAME('value')
                                    KEEPINPUTCOLUMNTYPES('true')
                                )
                            ) AS ref DIMENSION
                            TARGETIDCOLUMNS('__playerid__')
                            TARGETFEATURECOLUMN('variable')
                            TARGETVALUECOLUMN('value')
                            REFIDCOLUMNS('id')
                            REFFEATURECOLUMN('variable')
                            REFVALUECOLUMN('value_double')
                            MEASURE('Euclidean')
                          );",
                          "kmeans for 1 cluster with no aggregates with id one of the columns, explicit id alias, and WHERE clause")

})


test_that("computeClusterSample SQL is correct", {
  
  expect_equal_normalized(computeClusterSample(NULL, kmeans_obj, 0.01, test=TRUE), 
                          "SELECT * FROM antiselect(
         ON (SELECT * FROM sample(
             ON (SELECT clusterid, d.* FROM kmeansplot(
                   ON baseball.kmeans_test_scaled PARTITION BY ANY
                   ON baseball.kmeans_test_centroids DIMENSION
                   centroidsTable('baseball.kmeans_test_centroids')
                 ) kmp JOIN (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, g, h, r 
                               FROM batting 
                              WHERE yearid > 2010  ) d 
                       ON (kmp.id = d.id)
                 WHERE clusterid != -1
             )
             CONDITIONONCOLUMN('clusterid')
             CONDITIONON('0','1','2','3','4','5')
             SAMPLEFRACTION('0.01')
       )
         )
         EXCLUDE('id')
       )",
                          info="Sampling fraction unscaled data without row id.")
  
  expect_equal_normalized(computeClusterSample(NULL, kmeans_obj, 0.01, includeId=TRUE, test=TRUE), 
                          "SELECT * FROM sample(
             ON (SELECT clusterid, d.* FROM kmeansplot(
                   ON baseball.kmeans_test_scaled PARTITION BY ANY
                   ON baseball.kmeans_test_centroids DIMENSION
                   centroidsTable('baseball.kmeans_test_centroids')
                 ) kmp JOIN (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, g, h, r 
                               FROM batting 
                              WHERE yearid > 2010  ) d 
                       ON (kmp.id = d.id)
                 WHERE clusterid != -1
             )
             CONDITIONONCOLUMN('clusterid')
             CONDITIONON('0','1','2','3','4','5')
             SAMPLEFRACTION('0.01')
       )",
                          info="Sampling fraction unscaled data with row id.")
  
  expect_equal_normalized(computeClusterSample(NULL, kmeans_obj, 0.01, scaled=TRUE, test=TRUE), 
      "SELECT * FROM antiselect(
         ON (SELECT * FROM sample(
             ON (SELECT d.* FROM kmeansplot(
                   ON baseball.kmeans_test_scaled PARTITION BY ANY
                   ON baseball.kmeans_test_centroids DIMENSION
                   centroidsTable('baseball.kmeans_test_centroids')
                 ) d
                 WHERE clusterid != -1
             )
             CONDITIONONCOLUMN('clusterid')
             CONDITIONON('0','1','2','3','4','5')
             SAMPLEFRACTION('0.01')
       )
         )
         EXCLUDE('id')
       )",
      info="Sampling fraction scaled data without row id.")
  
  expect_equal_normalized(computeClusterSample(NULL, kmeans_obj, 0.01, scaled=TRUE, includeId=TRUE, test=TRUE), 
          "SELECT * FROM sample(
             ON (SELECT d.* FROM kmeansplot(
                   ON baseball.kmeans_test_scaled PARTITION BY ANY
                   ON baseball.kmeans_test_centroids DIMENSION
                   centroidsTable('baseball.kmeans_test_centroids')
                 ) d
                 WHERE clusterid != -1
             )
             CONDITIONONCOLUMN('clusterid')
             CONDITIONON('0','1','2','3','4','5')
             SAMPLEFRACTION('0.01')
       )",
          info="Sampling fraction scaled with row id.")
  
  expect_equal_normalized(computeClusterSample(NULL, kmeans_obj, sampleSize=1000, test=TRUE),
"SELECT * FROM antiselect(
          ON 
            (WITH stratum_counts AS (
          SELECT clusterid stratum, count(*) stratum_count 
            FROM kmeansplot(
              ON baseball.kmeans_test_scaled PARTITION BY ANY
              ON baseball.kmeans_test_centroids DIMENSION
              centroidsTable('baseball.kmeans_test_centroids')
            ) 
           WHERE clusterid != -1
           GROUP BY 1
         )
         SELECT * FROM sample (
           ON (SELECT clusterid, d.* FROM kmeansplot(
             ON baseball.kmeans_test_scaled PARTITION BY ANY
             ON baseball.kmeans_test_centroids DIMENSION
             centroidsTable('baseball.kmeans_test_centroids')
                ) kmp JOIN (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, g, h, r 
                                 FROM batting 
                                WHERE yearid > 2010  ) d 
                         ON (kmp.id = d.id)
                   WHERE clusterid != -1
               ) AS data PARTITION BY ANY
               ON stratum_counts AS summary DIMENSION
               CONDITIONONCOLUMN('clusterid')
               CONDITIONON('0','1','2','3','4','5')
               ApproximateSampleSize('1000')
             )
             )
             EXCLUDE('id')
          )",
          info="Sampling size unscaled data without row id.")
  
  expect_equal_normalized(computeClusterSample(NULL, kmeans_obj, sampleSize = 1000, includeId=TRUE, test=TRUE),
          "WITH stratum_counts AS (
            SELECT clusterid stratum, count(*) stratum_count 
              FROM kmeansplot(
                ON baseball.kmeans_test_scaled PARTITION BY ANY
                ON baseball.kmeans_test_centroids DIMENSION
                centroidsTable('baseball.kmeans_test_centroids')
              ) 
             WHERE clusterid != -1
            GROUP BY 1
           )
           SELECT * FROM sample (
             ON (SELECT clusterid, d.* FROM kmeansplot(
               ON baseball.kmeans_test_scaled PARTITION BY ANY
               ON baseball.kmeans_test_centroids DIMENSION
               centroidsTable('baseball.kmeans_test_centroids')
                 ) kmp JOIN (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, g, h, r 
                                 FROM batting 
                                WHERE yearid > 2010  ) d 
                         ON (kmp.id = d.id) 
                 WHERE clusterid != -1
             ) AS data PARTITION BY ANY
             ON stratum_counts AS summary DIMENSION
             CONDITIONONCOLUMN('clusterid')
             CONDITIONON('0','1','2','3','4','5')
             ApproximateSampleSize('1000')
           )",
          info="Sampling size unscaled data with id.")
  
  expect_equal_normalized(computeClusterSample(NULL, kmeans_obj, sampleSize=1000, scaled=TRUE, test=TRUE),
          "SELECT * FROM antiselect(
             ON 
             (WITH stratum_counts AS (
               SELECT clusterid stratum, count(*) stratum_count 
                 FROM kmeansplot(
                   ON baseball.kmeans_test_scaled PARTITION BY ANY
                   ON baseball.kmeans_test_centroids DIMENSION
                   centroidsTable('baseball.kmeans_test_centroids')
                 ) 
                WHERE clusterid != -1
               GROUP BY 1
             )
             SELECT * FROM sample (
               ON (SELECT d.* FROM kmeansplot(
                 ON baseball.kmeans_test_scaled PARTITION BY ANY
                 ON baseball.kmeans_test_centroids DIMENSION
                 centroidsTable('baseball.kmeans_test_centroids')
                  ) d
                  WHERE clusterid != -1
               ) AS data PARTITION BY ANY
               ON stratum_counts AS summary DIMENSION
               CONDITIONONCOLUMN('clusterid')
               CONDITIONON('0','1','2','3','4','5')
               ApproximateSampleSize('1000')
             )
             )
             EXCLUDE('id')
           )",
          info="Sampling size scaled data without id.")
  
  expect_equal_normalized(computeClusterSample(NULL, kmeans_obj, sampleSize = 1000, scaled=TRUE, includeId=TRUE, test=TRUE),
          "WITH stratum_counts AS (
            SELECT clusterid stratum, count(*) stratum_count 
              FROM kmeansplot(
                ON baseball.kmeans_test_scaled PARTITION BY ANY
                ON baseball.kmeans_test_centroids DIMENSION
                centroidsTable('baseball.kmeans_test_centroids')
              ) 
             WHERE clusterid != -1
            GROUP BY 1
           )
           SELECT * FROM sample (
             ON (SELECT d.* FROM kmeansplot(
               ON baseball.kmeans_test_scaled PARTITION BY ANY
               ON baseball.kmeans_test_centroids DIMENSION
               centroidsTable('baseball.kmeans_test_centroids')
                 ) d
                 WHERE clusterid != -1
             ) AS data PARTITION BY ANY
             ON stratum_counts AS summary DIMENSION
             CONDITIONONCOLUMN('clusterid')
             CONDITIONON('0','1','2','3','4','5')
             ApproximateSampleSize('1000')
           )",
          info="Sampling size scaled data with id.")
})


test_that("computeSilhouette SQL is correct", {
  
  kmeans_obj$scaledTableName="public.kmeans_test_scaled"
  kmeans_obj$centroidTableName="public.kmeans_test_centroids"
  
  expect_equal_normalized(computeSilhouette(NULL, kmeans_obj, scaled=TRUE, silhouetteTableName='public.kmeans_test_sil', test=TRUE),
    "-- Create Analytical Table with Silhouette Data
DROP TABLE IF EXISTS public.kmeans_test_sil;
CREATE ANALYTIC TABLE public.kmeans_test_sil
     DISTRIBUTE BY HASH(clusterid)
     AS
     WITH kmeansplotresult AS (
         SELECT clusterid, id, variable, coalesce(value_double, value_long, value_str::double) value 
           FROM unpivot(
                  ON (SELECT  d.* 
                   FROM kmeansplot(
                     ON public.kmeans_test_scaled PARTITION BY ANY
                     ON public.kmeans_test_centroids DIMENSION
                     centroidsTable('public.kmeans_test_centroids')
                   )  d 
                  WHERE clusterid != -1
                  )
                  COLSTOUNPIVOT('g', 'h', 'r')
                  COLSTOACCUMULATE('id','clusterid')
                  ATTRIBUTECOLUMNNAME('variable')
                  VALUECOLUMNNAME('value')
                  KEEPINPUTCOLUMNTYPES('true')
                ) 
     )
     SELECT target_clusterid clusterid, target_id id, a, b 
       FROM (
         SELECT target_clusterid, target_id,
                MAX(CASE WHEN target_clusterid = ref_clusterid THEN dissimilarity ELSE 0 END) a,
                MIN(CASE WHEN target_clusterid = ref_clusterid THEN 'Infinity' ELSE dissimilarity END) b 
           FROM
             (SELECT target_clusterid, target_id, ref_clusterid, avg(distance) dissimilarity 
                FROM VectorDistance(
                  ON kmeansplotresult AS target PARTITION BY id
                  ON kmeansplotresult AS ref DIMENSION
                  TARGETIDCOLUMNS('clusterid','id')
                  TARGETFEATURECOLUMN('variable')
                  TARGETVALUECOLUMN('value')
                  REFIDCOLUMNS('clusterid','id')
                  REFFEATURECOLUMN('variable')
                  REFVALUECOLUMN('value')
                  MEASURE('Euclidean')
                ) 
          WHERE target_id != ref_id 
          GROUP BY 1,2,3
         ) agg  
       GROUP BY 1,2
     ) sil;
--;
-- Compute overall silhouette value;
SELECT AVG((b-a)/greatest(a,b)) silhouette_value FROM public.kmeans_test_sil;
--;
-- Compute silhouette cluster profiles;
SELECT * FROM Hist_Reduce(
       ON Hist_Map(
         ON (SELECT clusterid::varchar clusterid, (b-a)/greatest(a,b) silhouette_value FROM public.kmeans_test_sil
         )
         STARTVALUE('-1')
         BINSIZE('0.05')
         ENDVALUE('1')
         VALUE_COLUMN('silhouette_value')
         GROUP_COLUMNS('clusterid')
       ) PARTITION BY clusterid
     );
--;
-- Drop Analytical Table with Silhouette Data;
DROP TABLE IF EXISTS public.kmeans_test_sil;",
          "compute Silhouette on scaled data")
  
  
  expect_equal_normalized(computeSilhouette(NULL, kmeans_obj, scaled=TRUE, silhouetteTableName='public.kmeans_test_sil', drop=FALSE, 
                                            test=TRUE),
    "-- Create Analytical Table with Silhouette Data
DROP TABLE IF EXISTS public.kmeans_test_sil;
CREATE ANALYTIC TABLE public.kmeans_test_sil
     DISTRIBUTE BY HASH(clusterid)
     AS
     WITH kmeansplotresult AS (
         SELECT clusterid, id, variable, coalesce(value_double, value_long, value_str::double) value 
           FROM unpivot(
                  ON (SELECT  d.* 
                   FROM kmeansplot(
                     ON public.kmeans_test_scaled PARTITION BY ANY
                     ON public.kmeans_test_centroids DIMENSION
                     centroidsTable('public.kmeans_test_centroids')
                   )  d 
                  WHERE clusterid != -1
                  )
                  COLSTOUNPIVOT('g', 'h', 'r')
                  COLSTOACCUMULATE('id','clusterid')
                  ATTRIBUTECOLUMNNAME('variable')
                  VALUECOLUMNNAME('value')
                  KEEPINPUTCOLUMNTYPES('true')
                ) 
     )
     SELECT target_clusterid clusterid, target_id id, a, b 
       FROM (
         SELECT target_clusterid, target_id,
                MAX(CASE WHEN target_clusterid = ref_clusterid THEN dissimilarity ELSE 0 END) a,
                MIN(CASE WHEN target_clusterid = ref_clusterid THEN 'Infinity' ELSE dissimilarity END) b 
           FROM
             (SELECT target_clusterid, target_id, ref_clusterid, avg(distance) dissimilarity 
                FROM VectorDistance(
                  ON kmeansplotresult AS target PARTITION BY id
                  ON kmeansplotresult AS ref DIMENSION
                  TARGETIDCOLUMNS('clusterid','id')
                  TARGETFEATURECOLUMN('variable')
                  TARGETVALUECOLUMN('value')
                  REFIDCOLUMNS('clusterid','id')
                  REFFEATURECOLUMN('variable')
                  REFVALUECOLUMN('value')
                  MEASURE('Euclidean')
                ) 
          WHERE target_id != ref_id 
          GROUP BY 1,2,3
         ) agg  
       GROUP BY 1,2
     ) sil;
--;
-- Compute overall silhouette value;
SELECT AVG((b-a)/greatest(a,b)) silhouette_value FROM public.kmeans_test_sil;
--;
-- Compute silhouette cluster profiles;
SELECT * FROM Hist_Reduce(
       ON Hist_Map(
         ON (SELECT clusterid::varchar clusterid, (b-a)/greatest(a,b) silhouette_value FROM public.kmeans_test_sil
         )
         STARTVALUE('-1')
         BINSIZE('0.05')
         ENDVALUE('1')
         VALUE_COLUMN('silhouette_value')
         GROUP_COLUMNS('clusterid')
       ) PARTITION BY clusterid
     );",
          "compute Silhouette on scaled data without drop of silhouette table")
  
  
  expect_equal_normalized(computeSilhouette(NULL, kmeans_obj, scaled=FALSE, silhouetteTableName='public.kmeans_test_sil', test=TRUE),
    "-- Create Analytical Table with Silhouette Data
DROP TABLE IF EXISTS public.kmeans_test_sil;
CREATE ANALYTIC TABLE public.kmeans_test_sil
     DISTRIBUTE BY HASH(clusterid)
     AS
     WITH kmeansplotresult AS (
         SELECT clusterid, id, variable, coalesce(value_double, value_long, value_str::double) value 
           FROM unpivot(
                  ON (SELECT  clusterid, d.* 
                   FROM kmeansplot(
                     ON public.kmeans_test_scaled PARTITION BY ANY
                     ON public.kmeans_test_centroids DIMENSION
                     centroidsTable('public.kmeans_test_centroids')
                   )  kmp JOIN (SELECT playerid || '-' || stint || '-' || teamid || '-' || yearid id, g, h, r FROM batting WHERE yearid > 2010  ) d ON (kmp.id = d.id)
                  WHERE clusterid != -1
                  )
                  COLSTOUNPIVOT('g', 'h', 'r')
                  COLSTOACCUMULATE('id','clusterid')
                  ATTRIBUTECOLUMNNAME('variable')
                  VALUECOLUMNNAME('value')
                  KEEPINPUTCOLUMNTYPES('true')
                ) 
     )
     SELECT target_clusterid clusterid, target_id id, a, b 
       FROM (
         SELECT target_clusterid, target_id,
                MAX(CASE WHEN target_clusterid = ref_clusterid THEN dissimilarity ELSE 0 END) a,
                MIN(CASE WHEN target_clusterid = ref_clusterid THEN 'Infinity' ELSE dissimilarity END) b 
           FROM
             (SELECT target_clusterid, target_id, ref_clusterid, avg(distance) dissimilarity 
                FROM VectorDistance(
                  ON kmeansplotresult AS target PARTITION BY id
                  ON kmeansplotresult AS ref DIMENSION
                  TARGETIDCOLUMNS('clusterid','id')
                  TARGETFEATURECOLUMN('variable')
                  TARGETVALUECOLUMN('value')
                  REFIDCOLUMNS('clusterid','id')
                  REFFEATURECOLUMN('variable')
                  REFVALUECOLUMN('value')
                  MEASURE('Euclidean')
                ) 
          WHERE target_id != ref_id 
          GROUP BY 1,2,3
         ) agg  
       GROUP BY 1,2
     ) sil;
--;
-- Compute overall silhouette value;
SELECT AVG((b-a)/greatest(a,b)) silhouette_value FROM public.kmeans_test_sil;
--;
-- Compute silhouette cluster profiles;
SELECT * FROM Hist_Reduce(
       ON Hist_Map(
         ON (SELECT clusterid::varchar clusterid, (b-a)/greatest(a,b) silhouette_value FROM public.kmeans_test_sil
         )
         STARTVALUE('-1')
         BINSIZE('0.05')
         ENDVALUE('1')
         VALUE_COLUMN('silhouette_value')
         GROUP_COLUMNS('clusterid')
       ) PARTITION BY clusterid
     );
--;
-- Drop Analytical Table with Silhouette Data;
DROP TABLE IF EXISTS public.kmeans_test_sil;",
          "compute Silhouette on non-scaled data")
})