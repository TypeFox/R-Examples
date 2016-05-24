## f1-classDefinitions.R --- 
## Author          : Jens Henrik Badsberg, Claus Dethlefsen, Luca La Rocca
## Created On      : Tue Nov 30 13:35:00 2004
## Last Modified By: Luca La Rocca
## Last Modified On: Sat Feb 25 16:30:00 2006
## Update Count    : 4
## Status          : Unknown, Use with caution!
######################################################

## vertices and edges

# class for a vertex set
setClass("vertexSet",contains="character")
# a 'vertexSet' object consists of a vector
# of unique character identifiers
# that are syntactically valid names

# virtual class for all edges
setClass("edge",contains="VIRTUAL")
# no slots, just a container

# class for an undirected (hyper)edge
setClass("undirectedEdge", contains = c("edge","integer"))
# an 'undirectedEdge' consists of a vector
# of unique numeric identifiers
# that are strictly positive integers
# referring to some 'vertexSet'

# class for a directed (hyper)edge
setClass("directedEdge", contains = c("edge","list"))
# a 'directedEdge' consists of a list of vectors
# that indicate non-empty disjoint undirected edges
# by referring to some 'vertexSet'

# class for a multi-set of edges
setClass("edgeList",contains="list")
# an 'edgeList' is a list of 'edge' objects

## representations

# class for an incidence list
setClass("incidenceList", representation(V = "vertexSet",E= "edgeList"))
# "G=(V,E)" with integers in 'E' not exceeding the cardinality of 'V'

# class for an incidence matrix
setClass("incidenceMatrix", contains="matrix")
# Entries of "I" are positive integers.
# One row per edge and one column per vertex.
# Multiple edges are allowed.
# We also allow loops and edges involving more than two vertices.
# Note that if a row of "I" has other than 0's and 1's,
# its non-zero elements are interpreted as a partial ordering
# of the vertices forming the edge.

# class for an adjacency list
setClass("adjacencyList", contains="list")
# 'A@.Data[[i]]$ne' is a vector of strictly positive integers indicating the neighbours of 'i'
# 'A@.Data[[i]]$pa' is a vector of strictly positive integers indicating the parents of 'i'
# 'A@.Data[[i]]$ch' is a vector of strictly positive integers indicating the children of 'i'

# class for an adjacency matrix
setClass("adjacencyMatrix", contains="matrix")
# Entries of "X" are 0's and 1's. One row and one column per vertex.
# Only edges involving two vertices are allowed. No loops, nor multiple edges.

## graphs

# class for any graph (a vertex set together with a multi-set of edges)
setClass("anyGraph",representation(incidenceList="incidenceList"))
# one representation available

# class for the most general kind of graphs we consider
setClass("generalGraph",contains="anyGraph",representation(incidenceMatrix="incidenceMatrix"))
# two representations available

# class for graphs with multiple edges but without hyper-edges
setClass("multiGraph", contains="generalGraph",representation(adjacencyList="adjacencyList"))
# three representations available

# class for graphs that are "simple" according to the standard definition
setClass("simpleGraph", contains="multiGraph",representation(adjacencyMatrix = "adjacencyMatrix"))
# four representations available
