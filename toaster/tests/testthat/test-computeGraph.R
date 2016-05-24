context("computeGraph")

test_that("toaGraph throws errors", {
  
  expect_error(toaGraph(NULL, NULL),
               "Both vertices and edges must be defined.")
  
  expect_error(toaGraph("vertices", NULL),
               "Both vertices and edges must be defined.")
  
  expect_error(toaGraph(NULL, "edges"),
               "Both vertices and edges must be defined.")
})

test_that("computeGraph throws errors", {
  
  expect_error(computeGraph(NULL, "not graph"),
               "Graph object must be specified.")
  
  expect_error(computeGraph(NULL, toaGraph("vs", "es"), test=TRUE),
               "Must provide allTables when test==TRUE.")
  
  expect_error(computeGraph(NULL, toaGraph("vs", "es"), v=list(logical(1)), 
                            allTables = data.frame(TABLE_NAME=c("vs","_es_"), stringsAsFactors = FALSE), test=TRUE),
               "Both vertices and edges must exist as tables or views.")
  
  expect_error(computeGraph(NULL, toaGraph("vs", "es"), v=list(logical(1)), 
                            allTables = data.frame(TABLE_NAME=c("vs","es"), stringsAsFactors = FALSE), test=TRUE),
               ".*Values must be either numeric or character only.")
})


simplestGraph = toaGraph("vertices", "edges")
simplestGraphWithEdgeAttrs = toaGraph("vertices", "edges", edgeAttrnames = c("weight","cost"))
vertexWhereGraph = toaGraph("vertices", "edges", vertexWhere = "state='TX'")

test_that("computeGraph works properly", {
  
  expect_equal_normalized(computeGraph(NULL, simplestGraph, 
                                       allTables = data.frame(TABLE_NAME=c("vertices","edges"), stringsAsFactors = FALSE),
                                       test=TRUE),
                          "-- Edges Select
                           SELECT source, target FROM edges ")
  
  expect_equal_normalized(computeGraph(NULL, simplestGraph, v=character(0),
                                       allTables = data.frame(TABLE_NAME=c("vertices","edges"), stringsAsFactors = FALSE),
                                       test=TRUE),
                          "-- Edges Select
                           SELECT source, target FROM edges ")
  
  expect_equal_normalized(computeGraph(NULL, simplestGraph, v=list(),
                                       allTables = data.frame(TABLE_NAME=c("vertices","edges"), stringsAsFactors = FALSE),
                                       test=TRUE),
                          "-- Edges Select
                           SELECT source, target FROM edges ")
  
  expect_equal_normalized(computeGraph(NULL, simplestGraphWithEdgeAttrs,
                                       allTables = data.frame(TABLE_NAME=c("vertices","edges"), stringsAsFactors = FALSE),
                                       test=TRUE),
                          "-- Edges Select
                           SELECT source, target, weight, cost FROM edges ")
  
  expect_equal_normalized(computeGraph(NULL, vertexWhereGraph,
                                       allTables = data.frame(TABLE_NAME=c("vertices","edges"), stringsAsFactors = FALSE),
                                       test=TRUE),
                          "-- Edges Select
                            SELECT source, target FROM edges 
                            WHERE source IN (SELECT id FROM vertices WHERE state='TX' )
                              AND target IN (SELECT id FROM vertices WHERE state='TX' ) ;
                           -- -- Vertices Select 
                           SELECT id FROM vertices WHERE state='TX' ")
  
  expect_equal_normalized(computeGraph(NULL, simplestGraph, v=list(1,2,3),
                                       allTables = data.frame(TABLE_NAME=c("vertices","edges"), stringsAsFactors = FALSE),
                                       test=TRUE),
                          "-- Edges Select
                            SELECT source, target FROM edges 
                            WHERE source IN (SELECT id FROM vertices WHERE id IN (1, 2, 3) )
                              AND target IN (SELECT id FROM vertices WHERE id IN (1, 2, 3) ) ;
                           -- -- Vertices Select 
                           SELECT id FROM vertices WHERE id IN (1, 2, 3) ")
  
  expect_equal_normalized(computeGraph(NULL, simplestGraph, v=list('1','2','3'),
                                       allTables = data.frame(TABLE_NAME=c("vertices","edges"), stringsAsFactors = FALSE),
                                       test=TRUE),
                          "-- Edges Select
                            SELECT source, target FROM edges 
                            WHERE source IN (SELECT id FROM vertices WHERE id IN ('1', '2', '3') )
                              AND target IN (SELECT id FROM vertices WHERE id IN ('1', '2', '3') ) ;
                           -- -- Vertices Select 
                           SELECT id FROM vertices WHERE id IN ('1', '2', '3') ")
  
  expect_equal_normalized(computeGraph(NULL, simplestGraph, v=list('a','b','c'),
                                       allTables = data.frame(TABLE_NAME=c("vertices","edges"), stringsAsFactors = FALSE),
                                       test=TRUE),
                          "-- Edges Select
                            SELECT source, target FROM edges 
                            WHERE source IN (SELECT id FROM vertices WHERE id IN ('a', 'b', 'c') )
                              AND target IN (SELECT id FROM vertices WHERE id IN ('a', 'b', 'c') ) ;
                          -- -- Vertices Select 
                           SELECT id FROM vertices WHERE id IN ('a', 'b', 'c') ")
  
  expect_equal_normalized(computeGraph(NULL, 
             toaGraph("graph.films_vertices", "graph.films_edges", FALSE, "name", "name1", "name2",
                      vertexAttrnames = "role", 
                      edgeAttrnames = c("weight","years")),
                      allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                             TABLE_NAME=c("films_vertices","films_edges")),
                      test=TRUE),
             "-- Edges Select
             SELECT name1, name2, weight, years 
               FROM graph.films_edges ;
             --
             -- Vertices Select
             SELECT name, role FROM graph.films_vertices ")
  
  expect_equal_normalized(computeGraph(NULL,
                                       toaGraph("graph.films_vertices", "graph.films_edges", FALSE,
                                                "name", "name1", "name2", 
                                                vertexAttrnames = "role",
                                                vertexWhere = "role = 'Actor'"),
                                                allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                             TABLE_NAME=c("films_vertices","films_edges")),
                      test=TRUE),
                      
                      "-- Edges Select
                      SELECT name1, name2 
                        FROM graph.films_edges
                       WHERE name1 IN (SELECT name FROM graph.films_vertices WHERE role = 'Actor' )
                         AND name2 IN (SELECT name FROM graph.films_vertices WHERE role = 'Actor' ) ;
                      --
                      -- Vertices Select
                      SELECT name, role FROM graph.films_vertices
                       WHERE role = 'Actor' ")
  
  expect_equal_normalized(computeGraph(NULL,
                                       toaGraph("graph.films_vertices", "graph.films_edges", FALSE,
                                                "name", "name1", "name2",
                                                vertexAttrnames = "role",
                                                vertexWhere = "role = 'Actor'"),
                                       v = list("Cruz", "Pacino", "Beluci", "Portman"),
                                       allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                      test=TRUE),
                      
                      "-- Edges Select
                      SELECT name1, name2 
                        FROM graph.films_edges
                       WHERE name1 IN (SELECT name FROM graph.films_vertices WHERE (role = 'Actor') AND 
                                                                                    name IN ('Cruz', 'Pacino', 'Beluci', 'Portman') )
                         AND name2 IN (SELECT name FROM graph.films_vertices WHERE (role = 'Actor') AND
                                                                                    name IN ('Cruz', 'Pacino', 'Beluci', 'Portman') ) ;
                      --
                      -- Vertices Select
                      SELECT name, role FROM graph.films_vertices
                       WHERE (role = 'Actor') AND name IN ('Cruz', 'Pacino', 'Beluci', 'Portman') ",
                      label="With vertex attr, vertex where, v list")
  
  expect_equal_normalized(computeGraph(NULL,
                                       toaGraph("graph.films_vertices", "graph.films_edges", FALSE,
                                                "name", "name1", "name2",
                                                vertexAttrnames = "role",
                                                vertexWhere = "role = 'Actor'"),
                                       edgeWhere = "years > 0",
                                       v = list("Cruz", "Pacino", "Beluci", "Portman"),
                                       allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                      test=TRUE),
                      
                      "-- Edges Select
                      SELECT name1, name2 
                        FROM graph.films_edges
                       WHERE (years > 0) 
                         AND name1 IN (SELECT name FROM graph.films_vertices WHERE (role = 'Actor') AND 
                                                                                    name IN ('Cruz', 'Pacino', 'Beluci', 'Portman') )
                         AND name2 IN (SELECT name FROM graph.films_vertices WHERE (role = 'Actor') AND
                                                                                    name IN ('Cruz', 'Pacino', 'Beluci', 'Portman') ) ;
                      --
                      -- Vertices Select
                      SELECT name, role FROM graph.films_vertices
                       WHERE (role = 'Actor') AND name IN ('Cruz', 'Pacino', 'Beluci', 'Portman') ",
                      label="With vertex attr, vertex where, edge where, v list")
  
  expect_equal_normalized(computeGraph(NULL,
                                       toaGraph("graph.films_vertices", "graph.films_edges", FALSE,
                                                "name", "name1", "name2",
                                                vertexAttrnames = "role"),
                                       v = "SELECT name FROM graph.films_vertices WHERE name like '%Bill%' or name like '%Tom%'",
                                       allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                      test=TRUE),
                      
                      "-- Edges Select
                      SELECT name1, name2 
                        FROM graph.films_edges
                       WHERE name1 IN (SELECT name FROM graph.films_vertices
                                        WHERE name IN (SELECT name FROM graph.films_vertices WHERE name like '%Bill%' or name like '%Tom%') )
                         AND name2 IN (SELECT name FROM graph.films_vertices 
                                        WHERE name IN (SELECT name FROM graph.films_vertices WHERE name like '%Bill%' or name like '%Tom%') ) ;
                      --
                      -- Vertices Select
                      SELECT name, role FROM graph.films_vertices
                       WHERE name IN (SELECT name FROM graph.films_vertices WHERE name like '%Bill%' or name like '%Tom%') ",
                      label="With vertex attr, v SELECT")
  
  expect_equal_normalized(computeGraph(NULL,
                                       toaGraph("graph.films_vertices", "graph.films_edges", FALSE,
                                                "name", "name1", "name2",
                                                vertexAttrnames = "role",
                                                vertexWhere = "role = 'Actor'"),
                                       v = "SELECT name FROM graph.films_vertices WHERE name like '%Bill%' or name like '%Tom%'",
                                       allTables = data.frame(TABLE_SCHEM=c("graph","graph"),
                                                              TABLE_NAME=c("films_vertices","films_edges")),
                      test=TRUE),
                      
                      "-- Edges Select
                      SELECT name1, name2 
                        FROM graph.films_edges
                       WHERE name1 IN (SELECT name FROM graph.films_vertices 
                                        WHERE (role = 'Actor') 
                                          AND name IN (SELECT name FROM graph.films_vertices WHERE name like '%Bill%' or name like '%Tom%') )
                         AND name2 IN (SELECT name FROM graph.films_vertices 
                                        WHERE (role = 'Actor') 
                                          AND name IN (SELECT name FROM graph.films_vertices WHERE name like '%Bill%' or name like '%Tom%') ) ;
                      --
                      -- Vertices Select
                      SELECT name, role FROM graph.films_vertices
                       WHERE (role = 'Actor') 
                         AND name IN (SELECT name FROM graph.films_vertices WHERE name like '%Bill%' or name like '%Tom%') ",
                      label="With vertex attr, vertex where, v SELECT")
  
})