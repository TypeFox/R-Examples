context("computeEgoGraph")

test_that("computeEgoGraph throws errors", {
  
  expect_error(computeEgoGraph(NULL, "not graph", mode="____"))
  
  expect_error(computeEgoGraph(NULL, "not graph"),
               "Graph object must be specified.")
  
  expect_error(computeEgoGraph(NULL, toaGraph("vs", "es"), test=TRUE),
               "Must provide allTables when test==TRUE.")
  
  expect_error(computeEgoGraph(NULL, toaGraph("vs","es"), ego=list("Tom Cruz"), 
                               allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                               test=TRUE),
               "Both vertices and edges must exist as tables or views.")
  
  expect_error(computeEgoGraph(NULL, toaGraph("vs", "es"), ego=list()),
               "Must have at least one ego vertex defined.")
  
  expect_error(computeEgoGraph(NULL, toaGraph("vs","es",FALSE), ego=list("Tom Cruz"), mode="in"),
               "Must be a directed graph when mode is 'in' or 'out'.")
})



filmGraph = toaGraph("graph.films_vertices", "graph.films_edges", FALSE,
                     "name", "name1", "name2")

test_that("computeEgoGraph works properly", {
  
  expect_equal_normalized(computeEgoGraph(NULL, filmGraph, 
                                          ego=list('Keanu Reeves'), 
                                          createDistanceAttr = FALSE,
                                          allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                                          test=TRUE),
                          computeEgoGraph(NULL, filmGraph, 
                                          ego=list('Keanu Reeves'), mode='both',
                                          createDistanceAttr = FALSE,
                                          allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                                          test=TRUE))
  
  expect_equal_normalized(computeEgoGraph(NULL, filmGraph, 
                                          ego=list('Keanu Reeves'), createDistanceAttr = FALSE,
                                          allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                                          test=TRUE),
                          
"BEGIN;
--
-- Create temp table of the shortest paths from ego vertices
CREATE TEMP FACT TABLE egographtemp 
       DISTRIBUTE BY HASH(source) 
       AS
       SELECT source, target, distance FROM AllPairsShortestPath(
         ON (SELECT name 
       FROM graph.films_vertices ) AS vertices PARTITION BY name
         ON (SELECT name1, name2 
         FROM graph.films_edges ) AS edges PARTITION BY name1
         ON (SELECT name 
       FROM graph.films_vertices WHERE name IN ('Keanu Reeves')  ) AS sources PARTITION BY name
         TARGETKEY('name2')
         DIRECTED('false')
         MAXDISTANCE('1')
       );
--
-- Edges Select
SELECT e.*
         FROM (SELECT name1, name2 
         FROM graph.films_edges ) e 
        WHERE name1 IN (SELECT target FROM egographtemp WHERE source = 'Keanu Reeves')
          AND name2 IN (SELECT target FROM egographtemp WHERE source = 'Keanu Reeves')
       UNION
       SELECT e.*
         FROM (SELECT name1, name2 
         FROM graph.films_edges ) e
        WHERE name1 = 'Keanu Reeves'
           OR name2 = 'Keanu Reeves';
--
END",
                          label="With Keanu Reeves without distance attr"
)
  
  expect_equal_normalized(computeEgoGraph(NULL, filmGraph, 
                                          ego=list('Keanu Reeves'),
                                          allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                                          test=TRUE),
"BEGIN;
--
-- Create temp table of the shortest paths from ego vertices
CREATE TEMP FACT TABLE egographtemp 
       DISTRIBUTE BY HASH(source) 
       AS
       SELECT source, target, distance FROM AllPairsShortestPath(
         ON (SELECT name 
       FROM graph.films_vertices ) AS vertices PARTITION BY name
         ON (SELECT name1, name2 
         FROM graph.films_edges ) AS edges PARTITION BY name1
         ON (SELECT name 
       FROM graph.films_vertices WHERE name IN ('Keanu Reeves')  ) AS sources PARTITION BY name
         TARGETKEY('name2')
         DIRECTED('false')
         MAXDISTANCE('1')
       );
--
-- Edges Select
SELECT e.*
         FROM (SELECT name1, name2 
         FROM graph.films_edges ) e 
        WHERE name1 IN (SELECT target FROM egographtemp WHERE source = 'Keanu Reeves')
          AND name2 IN (SELECT target FROM egographtemp WHERE source = 'Keanu Reeves')
       UNION
       SELECT e.*
         FROM (SELECT name1, name2 
         FROM graph.films_edges ) e
        WHERE name1 = 'Keanu Reeves'
           OR name2 = 'Keanu Reeves';
--
-- Vertices Select
SELECT v.*, eg.distance __distance_attr__ 
           FROM egographtemp eg JOIN
                (SELECT name 
       FROM graph.films_vertices ) v ON (eg.target = v.name)
          WHERE eg.source = 'Keanu Reeves'
         UNION
         SELECT v.*, 0 __distance_attr__ 
           FROM egographtemp eg JOIN
                (SELECT name 
       FROM graph.films_vertices ) v ON (eg.source = v.name)
          WHERE eg.source = 'Keanu Reeves';
--
END",
                          label="With Keanu Reeves with distance attr")
  
  expect_equal_normalized(computeEgoGraph(NULL, 
                                          toaGraph("graph.films_vertices", "graph.films_edges", FALSE,
                                                   "name", "name1", "name2",
                                                   vertexAttrnames = "role", edgeAttrnames = c("weight","years")),
                ego = list('Keanu Reeves','Takeshi Kitano'), order = 2,
                test=TRUE, allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                               TABLE_NAME=c("films_vertices","films_edges"))),
"BEGIN;
--
-- Create temp table of the shortest paths from ego vertices
CREATE TEMP FACT TABLE egographtemp 
       DISTRIBUTE BY HASH(source) 
       AS
       SELECT source, target, distance FROM AllPairsShortestPath(
         ON (SELECT name, role 
       FROM graph.films_vertices ) AS vertices PARTITION BY name
         ON (SELECT name1, name2, weight, years 
         FROM graph.films_edges ) AS edges PARTITION BY name1
         ON (SELECT name, role 
       FROM graph.films_vertices WHERE name IN ('Keanu Reeves', 'Takeshi Kitano')  ) AS sources PARTITION BY name
         TARGETKEY('name2')
         DIRECTED('false')
         MAXDISTANCE('2')
       );
--
-- Edges Select
SELECT e.*
         FROM (SELECT name1, name2, weight, years 
         FROM graph.films_edges ) e 
        WHERE name1 IN (SELECT target FROM egographtemp WHERE source = 'Keanu Reeves')
          AND name2 IN (SELECT target FROM egographtemp WHERE source = 'Keanu Reeves')
       UNION
       SELECT e.*
         FROM (SELECT name1, name2, weight, years 
         FROM graph.films_edges ) e
        WHERE name1 = 'Keanu Reeves'
           OR name2 = 'Keanu Reeves';
--
-- Vertices Select
SELECT v.*, eg.distance __distance_attr__ 
           FROM egographtemp eg JOIN
                (SELECT name, role 
       FROM graph.films_vertices ) v ON (eg.target = v.name)
          WHERE eg.source = 'Keanu Reeves'
         UNION
         SELECT v.*, 0 __distance_attr__ 
           FROM egographtemp eg JOIN
                (SELECT name, role 
       FROM graph.films_vertices ) v ON (eg.source = v.name)
          WHERE eg.source = 'Keanu Reeves';
--
-- Edges Select
SELECT e.*
         FROM (SELECT name1, name2, weight, years 
         FROM graph.films_edges ) e 
        WHERE name1 IN (SELECT target FROM egographtemp WHERE source = 'Takeshi Kitano')
          AND name2 IN (SELECT target FROM egographtemp WHERE source = 'Takeshi Kitano')
       UNION
       SELECT e.*
         FROM (SELECT name1, name2, weight, years 
         FROM graph.films_edges ) e
        WHERE name1 = 'Takeshi Kitano'
           OR name2 = 'Takeshi Kitano';
--
-- Vertices Select
SELECT v.*, eg.distance __distance_attr__ 
           FROM egographtemp eg JOIN
                (SELECT name, role 
       FROM graph.films_vertices ) v ON (eg.target = v.name)
          WHERE eg.source = 'Takeshi Kitano'
         UNION
         SELECT v.*, 0 __distance_attr__ 
           FROM egographtemp eg JOIN
                (SELECT name, role 
       FROM graph.films_vertices ) v ON (eg.source = v.name)
          WHERE eg.source = 'Takeshi Kitano';
--
END")
  
})