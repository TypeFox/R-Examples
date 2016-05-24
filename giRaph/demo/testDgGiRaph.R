
library(giRaph)
#library(dynamicGraph)

G.1 <- new("incidenceList", 
           E = list(u(1, 2), d(1, 3), u(3), 
                    d(2, 5), d(2, 5), d(3, c(1, 4), 5),
                    u(2, 4, 5), d(c(3, 4), c(2, 1)), r(1, 5)
                    ), 
           V = 1:5)
g.1 <- new("anyGraph", incidenceList = G.1)

G.2 <- new("incidenceList", E = list(), V = 5:10)
g.2 <- new("anyGraph", incidenceList = G.2)

incidenceList(as(g.1, "anyGraph"))
incidenceList(as(g.1, "generalGraph"))
incidenceList(as(g.1, "multiGraph"))
incidenceList(as(g.1, "simpleGraph"))

#str(as(as(g.1, "anyGraph"),     "dg.graph")) # 
#str(as(as(g.1, "generalGraph"), "dg.graph")) # hyper edges
#str(as(as(g.1, "multiGraph"),   "dg.simple.graph")) # multiple edges, no hyper
#str(as(as(g.1, "multiGraph"),   "dg.graph")) # multiple edges, no hyper
#str(as(as(g.1, "simpleGraph"),  "dg.simple.graph")) # simple, directed
#str(as(as(g.1, "simpleGraph"),  "dg.graph")) # simple, directed

#as(as(g.1, "anyGraph"),     "dg.graph")@edgeList
#as(as(g.1, "generalGraph"), "dg.graph")@edgeList
#as(as(g.1, "multiGraph"),   "dg.graph")@edgeList
#as(as(g.1, "simpleGraph"),  "dg.graph")@edgeList

#dynamic.Graph(as(g.1, "anyGraph"),     control = dg.control(label = "anyGraph, debug.edges = TRUE"))
#dynamic.Graph(as(g.1, "generalGraph"), control = dg.control(label = "generalGraph"))
#dynamic.Graph(as(g.1, "multiGraph"),   control = dg.control(label = "multiGraph"))
#dynamic.Graph(as(g.1, "simpleGraph"),  control = dg.control(label = "simpleGraph"))

#dg(as(as(g.1, "simpleGraph"), "dg.graph"))

#dynamic.Graph(as(as(g.1, "simpleGraph"), "anyGraph"))

#dynamic.Graph(as(as(g.1, "simpleGraph"), "anyGraph"), modelObject = new("dg.Model", name = "AnModelObject"))

#dynamic.Graph(as(as(g.1, "simpleGraph"), "anyGraph"), modelObject = new("dg.Model"))
