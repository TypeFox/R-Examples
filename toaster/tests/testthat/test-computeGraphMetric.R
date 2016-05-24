context("computeGraphMetric")

test_that("computeGraphMetric throws errors", {
  
  expect_error(computeGraphMetric(type = 'not a type'),
               ".*'arg' should be one of \"degree\", \"in-degree\", \"out-degree\", \"clustering\", \"shortestpath\", \"pagerank\", \"betweenness\", \"eigenvector\"")
  
  expect_error(computeGraphMetric(type='degree', rank='not rank'),
               ".*'arg' should be one of \"rank\", \"rownumber\", \"row\", \"denserank\", \"percentrank\"")
  
  expect_error(computeGraphMetric(NULL, "not graph"),
               "Graph object must be specified.")
  
  expect_error(computeGraphMetric(NULL, toaGraph("vs", "es"), test=TRUE),
               "Must provide allTables when test==TRUE.")
  
  expect_error(computeGraphMetric(NULL, toaGraph("vs", "es"), 
                                  allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                                  test=TRUE),
               "Both vertices and edges must exist as tables or views.")

})


policeGraphUn = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_un", FALSE,
                         "officer", "officer1", "officer2", vertexAttrnames = c("offense_count"),
                         edgeAttrnames = c("weight"))
policeGraphDi = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_di", TRUE,
                         "officer", "officer1", "officer2", vertexAttrnames = c("offense_count"),
                         edgeAttrnames = c("weight"))

test_that("computeGraphMetric for degree works properly", {
  
  expect_equal_normalized(computeGraphMetric(NULL, policeGraphUn, type="degree", top=20,
                                             rankFunction="denserank", 
                                             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un",
                                                                                 "dallaspolice_officer_edges_di")),
                                             test=TRUE),
"SELECT COALESCE(s.key, t.key) key,
            'degree'::varchar degree_type,
            COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) degree, DENSE_RANK() OVER (PARTITION BY 1 ORDER BY COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) DESC) rank
      FROM (SELECT officer1 key, COUNT(*) cnt_source 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e 
             GROUP BY 1) s FULL JOIN
           (SELECT officer2 key, COUNT(*) cnt_target 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e
             GROUP BY 1) t ON (s.key = t.key)
     
      ORDER BY degree DESC LIMIT 20", label="degree metric for undirected graph")
  
  expect_equal_normalized(computeGraphMetric(NULL, policeGraphDi, type="out-degree", top=20,
                                             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un",
                                                                                 "dallaspolice_officer_edges_di")),
                                             test=TRUE),
"SELECT key, degree_type, degree_long degree, RANK() OVER (PARTITION BY degree_type ORDER BY degree_long DESC) rank
       FROM unpivot(
         ON (SELECT COALESCE(s.key, t.key) key, 
                    COALESCE(s.cnt_source,0) outdegree, 
                    COALESCE(t.cnt_target,0) indegree,  
                    COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) degree,
                    (COALESCE(t.cnt_target,0) + 1)/(COALESCE(s.cnt_source,0) + 1) inbyoutdegree
               FROM (SELECT officer1 key, COUNT(*) cnt_source 
                       FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_di ) e 
                      GROUP BY 1) s FULL JOIN
                    (SELECT officer2 key, COUNT(*) cnt_target 
                       FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_di ) e
                      GROUP BY 1) t ON (s.key = t.key)
         )
         colsToUnpivot('outdegree','indegree','degree','inbyoutdegree')
         colsToAccumulate('key')
         keepInputColumnTypes('true')
         ATTRIBUTECOLUMNNAME('degree_type')
         VALUECOLUMNNAME('degree')
     )
      WHERE degree_type = 'outdegree'
      ORDER BY degree DESC LIMIT 20", label="out-degree metric for directed graph")
  
  expect_equal_normalized(computeGraphMetric(NULL, policeGraphDi, type="in-degree", top=20,
                                             rankFunction = "percentrank",
                                             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un",
                                                                                 "dallaspolice_officer_edges_di")),
                                             test=TRUE),
"SELECT key, degree_type, degree_long degree, PERCENT_RANK() OVER (PARTITION BY degree_type ORDER BY degree_long DESC) rank
       FROM unpivot(
         ON (SELECT COALESCE(s.key, t.key) key, 
                    COALESCE(s.cnt_source,0) outdegree, 
                    COALESCE(t.cnt_target,0) indegree,  
                    COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) degree,
                    (COALESCE(t.cnt_target,0) + 1)/(COALESCE(s.cnt_source,0) + 1) inbyoutdegree
               FROM (SELECT officer1 key, COUNT(*) cnt_source 
                       FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_di ) e 
                      GROUP BY 1) s FULL JOIN
                    (SELECT officer2 key, COUNT(*) cnt_target 
                       FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_di ) e
                      GROUP BY 1) t ON (s.key = t.key)
         )
         colsToUnpivot('outdegree','indegree','degree','inbyoutdegree')
         colsToAccumulate('key')
         keepInputColumnTypes('true')
         ATTRIBUTECOLUMNNAME('degree_type')
         VALUECOLUMNNAME('degree')
     )
      WHERE degree_type = 'indegree'
      ORDER BY degree DESC LIMIT 20", label="in-degree metric for directed graph")
  
  expect_equal_normalized(computeGraphMetric(NULL, policeGraphUn, type="degree", top=-1,
                                             rankFunction="denserank", 
                                             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un",
                                                                                 "dallaspolice_officer_edges_di")),
                                             test=TRUE),
"SELECT COALESCE(s.key, t.key) key,
            'degree'::varchar degree_type,
            COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) degree, DENSE_RANK() OVER (PARTITION BY 1 ORDER BY COALESCE(s.cnt_source,0) + COALESCE(t.cnt_target,0) DESC) rank
      FROM (SELECT officer1 key, COUNT(*) cnt_source 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e 
             GROUP BY 1) s FULL JOIN
           (SELECT officer2 key, COUNT(*) cnt_target 
              FROM (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) e
             GROUP BY 1) t ON (s.key = t.key) ", label="degree metric for undirected graph")
  
})


test_that("computeGraphMetric for clustering works properly", {
  
  expect_equal_normalized(computeGraphMetric(NULL, policeGraphUn, type="clustering", top=40,
                                             rankFunction = "denserank",
                                             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un",
                                                                                 "dallaspolice_officer_edges_di")),
                                             test=TRUE),
"SELECT officer key, cc_type, coalesce(cc_double, cc_str::double) cc, DENSE_RANK() OVER (PARTITION BY cc_type ORDER BY coalesce(cc_double, cc_str::double) DESC) rank
       FROM unpivot(
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
     )
      WHERE cc_type = 'cc'
      ORDER BY cc DESC LIMIT 40", label="clustering metric for undirected graph")
  
  expect_equal_normalized(computeGraphMetric(NULL, policeGraphDi, type="clustering", top=40,
                                             rankFunction = "row",
                                             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un",
                                                                                 "dallaspolice_officer_edges_di")),
                                             test=TRUE),
"SELECT officer key, cc_type, coalesce(cc_double, cc_str::double) cc, ROW_NUMBER() OVER (PARTITION BY cc_type ORDER BY coalesce(cc_double, cc_str::double) DESC) rank
       FROM unpivot(
         ON (SELECT * FROM LocalClusteringCoefficient(
               ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
               ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_di ) AS edges PARTITION BY officer1
               targetKey('officer2')
             
               directed('true')
               accumulate('officer')
         ))
         colsToUnpivot('cyc_cc','mid_cc','in_cc','out_cc','avg_cc')
         colsToAccumulate('officer')
         keepInputColumnTypes('true')
         ATTRIBUTECOLUMNNAME('cc_type')
         VALUECOLUMNNAME('cc')
     )
      WHERE cc_type = 'avg_cc'
      ORDER BY cc DESC LIMIT 40", label="clustering metric for directed graph")
  
})

test_that("computeGraphMetric for pagerank works properly", {
  
  expect_equal_normalized(computeGraphMetric(NULL, policeGraphDi, type="pagerank", top=50,
                                             rankFunction = "rownumber",
                                             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un",
                                                                                 "dallaspolice_officer_edges_di")),
                                             test=TRUE),
"SELECT officer key, pagerank, ROW_NUMBER() OVER (PARTITION BY 1 ORDER BY pagerank DESC) rank
       FROM PageRank(
       ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
       ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_di ) AS edges PARTITION BY officer1
       targetKey('officer2')
       accumulate('officer')
     ) ORDER BY pagerank DESC LIMIT 50", label="pagerank metric for undirected graph")
  
})

test_that("computeGraphMetric for betweenness works properly", {
  
  expect_equal_normalized(computeGraphMetric(NULL, policeGraphUn, type="betweenness", top=25,
                                             rankFunction = "rank",
                                             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un",
                                                                                 "dallaspolice_officer_edges_di")),
                                             test=TRUE),
"SELECT officer key, betweenness, RANK() OVER (PARTITION BY 1 ORDER BY betweenness DESC) rank
       FROM Betweenness(
       ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
       ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) AS edges PARTITION BY officer1
       targetKey('officer2')
       accumulate('officer')
     ) ORDER BY betweenness DESC LIMIT 25", label="computeGraphMetric for betweenness works properly")
})

test_that("computeGraphMetric for eigenvector centrality works properly", {
  
  expect_equal_normalized(computeGraphMetric(NULL, policeGraphUn, type="eigenvector", top=10,
                                             rankFunction = "denserank",
                                             allTables = data.frame(TABLE_NAME=c("dallaspolice_officer_vertices","dallaspolice_officer_edges_un",
                                                                                 "dallaspolice_officer_edges_di")),
                                             test=TRUE),
"SELECT officer key, centrality, DENSE_RANK() OVER (PARTITION BY 1 ORDER BY centrality DESC) rank
       FROM EigenVectorCentrality(
       ON (SELECT officer, offense_count 
       FROM dallaspolice_officer_vertices ) AS vertices PARTITION BY officer
       ON (SELECT officer1, officer2, weight 
         FROM dallaspolice_officer_edges_un ) AS edges PARTITION BY officer1
       targetKey('officer2')
       accumulate('officer')
       directed('false')
     ) ORDER BY centrality DESC LIMIT 10", label="computeGraphMetric for eigenvector works properly")
})