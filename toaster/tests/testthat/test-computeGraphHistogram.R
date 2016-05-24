context("computeGraphHistogram")

test_that("computeGraphHistogram throws errors", {
  
  expect_error(computeGraphHistogram(type = 'not a type'),
               ".*'arg' should be one of \"degree\", \"clustering\", \"shortestpath\", \"pagerank\", \"betweenness\", \"eigenvector\"")
  
  expect_error(computeGraphHistogram(binMethod = 'not a method'),
               ".*'arg' should be one of \"manual\", \"Sturges\", \"Scott\"")
  
  expect_error(computeGraphHistogram(NULL, "not graph"),
               "Graph object must be specified.")
  
  expect_error(computeGraphHistogram(NULL, toaGraph("vs", "es"), test=TRUE),
               "Must provide allTables when test==TRUE.")
  
  expect_error(computeGraphHistogram(NULL, toaGraph("vs", "es"), type='shortestpath', numbins=30,
                          allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_di")),
                          test=TRUE),
               "Both vertices and edges must exist as tables or views.")
  
  expect_error(computeGraphHistogram(NULL, toaGraph("vs", "es"), weight=TRUE, 
                                     allTables = data.frame(TABLE_NAME=c("vs","es"), stringsAsFactors = FALSE), 
                                     test=TRUE),
               "No edge attribute 'weight' found in graph.")
  
  expect_error(computeGraphHistogram(NULL, toaGraph("vs", "es"), weight="notweight",
                                     allTables = data.frame(TABLE_NAME=c("vs","es"), stringsAsFactors = FALSE), 
                                     test=TRUE),
               "No edge attribute 'notweight' found in graph.")
  
  expect_error(computeGraphHistogram(NULL, toaGraph("vs", "es"), 
                                      binMethod = 'manual', numbins = NULL,
                                      allTables = data.frame(TABLE_NAME=c("vs","es"), stringsAsFactors = FALSE), 
                                      test=TRUE),
               "Number of bins or/and at least startvalue and endvalue must be defined when method is 'manual'.")
  
  expect_error(computeGraphHistogram(NULL, toaGraph("vs", "es"), 
                                      binMethod = 'manual', startvalue = 10, endvalue = 9,
                                      allTables = data.frame(TABLE_NAME=c("vs","es"), stringsAsFactors = FALSE), 
                                      test=TRUE),
               "End value should be greater than start value.")
})

policeGraphUn = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_un", FALSE,
                         "officer", "officer1", "officer2", vertexAttrnames = c("offense_count"),
                         edgeAttrnames = c("weight"))
policeGraphDi = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_di", TRUE,
                         "officer", "officer1", "officer2", vertexAttrnames = c("offense_count"),
                         edgeAttrnames = c("weight"))

test_that("computeGraphHistogram for degree works properly", {

  expect_equal_normalized(computeGraphHistogram(NULL, policeGraphUn,
             binMethod = 'manual', numbins = 30,
                      allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un")),
                      test=TRUE),
"BEGIN;
--
-- Compute degree into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT COALESCE(s.key, t.key) key,
            'degree'::varchar degree_type,
            COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) degree
      FROM (SELECT officer1 key, COUNT(*) cnt_source 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e 
             GROUP BY 1) s FULL JOIN
           (SELECT officer2 key, COUNT(*) cnt_target 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e
             GROUP BY 1) t ON (s.key = t.key);
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          ON hist_prep(
                           ON graphdataforhisttemp VALUE_COLUMN('degree')) as data_stat DIMENSION
                           BIN_SELECT('30')
             VALUE_COLUMN('degree')
            GROUP_COLUMNS('degree_type')
         ) PARTITION BY degree_type
    );
--
END", label="Degree histogram with numbins")
  
  expect_equal_normalized(computeGraphHistogram(NULL, policeGraphUn,
             binMethod = 'Scott', 
             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un")),
                      test=TRUE),
"BEGIN;
--
-- Compute degree into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT COALESCE(s.key, t.key) key,
            'degree'::varchar degree_type,
            COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) degree
      FROM (SELECT officer1 key, COUNT(*) cnt_source 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e 
             GROUP BY 1) s FULL JOIN
           (SELECT officer2 key, COUNT(*) cnt_target 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e
             GROUP BY 1) t ON (s.key = t.key);
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          ON hist_prep(
                         ON graphdataforhisttemp VALUE_COLUMN('degree')) as data_stat DIMENSION
                         BIN_SELECT('Scott')
             VALUE_COLUMN('degree')
            GROUP_COLUMNS('degree_type')
         ) PARTITION BY degree_type
    );
--
END", label="Degree histogram with Scott method")
  
  expect_equal_normalized(computeGraphHistogram(NULL, policeGraphUn,
             binMethod = 'manual', startvalue = 0, endvalue = 30, numbins = 25, 
             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un")),
             test=TRUE),
"BEGIN;
--
-- Compute degree into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT COALESCE(s.key, t.key) key,
            'degree'::varchar degree_type,
            COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) degree
      FROM (SELECT officer1 key, COUNT(*) cnt_source 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e 
             GROUP BY 1) s FULL JOIN
           (SELECT officer2 key, COUNT(*) cnt_target 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e
             GROUP BY 1) t ON (s.key = t.key);
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          binsize('1.2')
                         startvalue('0')
                         endvalue('30')
             VALUE_COLUMN('degree')
            GROUP_COLUMNS('degree_type')
         ) PARTITION BY degree_type
    );
--
END", label="Degree histogram with start value, end value, and numbins")
  
})


airlineGraphGlobal = toaGraph("karthik.global_vertex_view", "karthik.global_edges_view", TRUE, 
                              "city", "vertex1", "vertex2",
                              edgeAttrnames = c("score"))

test_that("computeGraphHistogram for clustering works properly", {
  
  expect_equal_normalized(computeGraphHistogram(conn, airlineGraphGlobal, type='clustering',
                           binMethod = 'manual', binsize = 0.01, endvalue = 1,
                           allTables = data.frame(TABLE_SCHEM=c("karthik","karthik"),
                                                  TABLE_NAME=c("global_vertex_view","global_edges_view")),
                           test=TRUE),
"BEGIN;
--
-- Compute clustering into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT city key, cc_type, coalesce(cc_double, cc_str::double) cc FROM unpivot(
       ON (SELECT * FROM LocalClusteringCoefficient(
             ON (SELECT city 
       FROM karthik.global_vertex_view ) AS vertices PARTITION BY city
             ON (SELECT vertex1, vertex2, score 
         FROM karthik.global_edges_view ) AS edges PARTITION BY vertex1
             targetKey('vertex2')
           
             directed('true')
             accumulate('city')
       ))
       colsToUnpivot('cyc_cc','mid_cc','in_cc','out_cc','avg_cc')
       colsToAccumulate('city')
       keepInputColumnTypes('true')
       ATTRIBUTECOLUMNNAME('cc_type')
       VALUECOLUMNNAME('cc')
     );
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          binsize('0.01')
                         startvalue('0')
                         endvalue('1')
             VALUE_COLUMN('cc')
            GROUP_COLUMNS('cc_type')
         ) PARTITION BY cc_type
    );
--
END", label="Clustering histogram for directed graph with manual binmethod, binsize, endvalue")
  
  expect_equal_normalized(computeGraphHistogram(NULL, policeGraphUn, type='clustering', numbins = 100,
                          allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un")),
                          test=TRUE),
"BEGIN;
--
-- Compute clustering into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT officer key, cc_type, coalesce(cc_double, cc_str::double) cc FROM unpivot(
       ON (SELECT * FROM LocalClusteringCoefficient(
             ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
             ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) AS edges PARTITION BY officer1
             targetKey('officer2')
           
             directed('false')
             accumulate('officer')
       ))
       colsToUnpivot('cc')
       colsToAccumulate('officer')
       keepInputColumnTypes('true')
       ATTRIBUTECOLUMNNAME('cc_type')
       VALUECOLUMNNAME('cc')
     );
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          ON hist_prep(
                           ON graphdataforhisttemp VALUE_COLUMN('cc')) as data_stat DIMENSION
                           BIN_SELECT('100')
             VALUE_COLUMN('cc')
            GROUP_COLUMNS('cc_type')
         ) PARTITION BY cc_type
    );
--
END", label="Clustering histogram with numbins (manual by default)")
  
})

test_that("computeGraphHistogram for shortestpath works properly", {
  
  expect_equal_normalized(computeGraphHistogram(conn, policeGraphUn, type='shortestpath',
                           binMethod = 'manual', binsize = 0.01, endvalue = 1,
                           allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un")),
                           test=TRUE),
"BEGIN;
--
-- Compute shortestpath into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT source key, target, distance FROM AllPairsShortestPath(
       ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
       ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) AS edges PARTITION BY officer1
       targetKey('officer2')
       directed('false')
     
);
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          binsize('0.01')
                         startvalue('0')
                         endvalue('1')
             VALUE_COLUMN('distance')
            
         ) PARTITION BY 1
    );
--
END", label="Shortestpath histogram with binsize and endvalue")
  
  expect_equal_normalized(computeGraphHistogram(NULL, policeGraphDi, type='shortestpath', numbins=15,
                          allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_di")),
                          test=TRUE),
"BEGIN;
--
-- Compute shortestpath into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT source key, target, distance FROM AllPairsShortestPath(
       ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
       ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_di ) AS edges PARTITION BY officer1
       targetKey('officer2')
       directed('true')
     
);
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          ON hist_prep(
                           ON graphdataforhisttemp VALUE_COLUMN('distance')) as data_stat DIMENSION
                           BIN_SELECT('15')
             VALUE_COLUMN('distance')
            
         ) PARTITION BY 1
    );
--
END", label="Shortestpath histggram for directed graph with numbins")
  
})


test_that("computeGraphHistogram for pagerank works properly", {
  
  expect_equal_normalized(computeGraphHistogram(NULL, policeGraphDi, type="pagerank", weight=TRUE, binMethod = "Scott",
                          allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_di")),
                          test=TRUE),
"BEGIN;
--
-- Compute pagerank into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT officer key, pagerank
       FROM PageRank(
       ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
       ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_di ) AS edges PARTITION BY officer1
       targetKey('officer2')
       accumulate('officer')
     edgeweight('weight')
);
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          ON hist_prep(
                         ON graphdataforhisttemp VALUE_COLUMN('pagerank')) as data_stat DIMENSION
                         BIN_SELECT('Scott')
             VALUE_COLUMN('pagerank')
            
         ) PARTITION BY 1
    );
--
END", label="Pagerank histogram for undirected graph with weight with binmethod 'Scott'")
  
})


test_that("computeGraphHistogram for betweenness works properly", {
  
  expect_equal_normalized(computeGraphHistogram(NULL, policeGraphUn, type="betweenness", weight=TRUE,
                          numbins = 100, startvalue = 1, endvalue = 5001,
                          allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un")),
                          test=TRUE),
"BEGIN;
--
-- Compute betweenness into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT officer key, betweenness
       FROM Betweenness(
       ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
       ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) AS edges PARTITION BY officer1
       targetKey('officer2')
       accumulate('officer')
     edgeweight('weight')
);
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          binsize('50')
                         startvalue('1')
                         endvalue('5001')
             VALUE_COLUMN('betweenness')
            
         ) PARTITION BY 1
    );
--
END", label="Betweenness histogram for undirected graph with weight and numbins, startvalue, endvalue")
  
})


test_that("computeGraphHistogram for eigenvector works properly", {
  
  expect_equal_normalized(computeGraphHistogram(NULL, policeGraphUn, type="eigenvector", weight=TRUE,
                          numbins = 100, startvalue = 0, endvalue = 5001,
                          allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un")),
                          test=TRUE),
"BEGIN;
--
-- Compute eigenvector into temp table with all vertices
CREATE TEMP FACT TABLE graphdataforhisttemp 
     DISTRIBUTE BY HASH(key) 
     AS
     SELECT officer key, centrality
       FROM EigenVectorCentrality(
       ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
       ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) AS edges PARTITION BY officer1
       targetKey('officer2')
       accumulate('officer')
       directed('false')
     edgeweight('weight')
     );
--
-- Compute metric histogram
SELECT * FROM Hist_Reduce(
           ON Hist_Map(
             ON graphdataforhisttemp as data_input PARTITION BY ANY 
          binsize('50.01')
                         startvalue('0')
                         endvalue('5001')
             VALUE_COLUMN('centrality')
            
         ) PARTITION BY 1
    );
--
END", label="Eigenvector histogram for undirected graph with weight and numbins, startvalue, endvalue")

})