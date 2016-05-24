## it looks like <<- assignments here should actually be to env.

# Represent the tree as a flat collection of nodes
# but allocate the list ahead of time and grow it
# by doubling the space. This makes things a lot faster
# for large trees.

xmlFlatListTree =
function(nodes = list(),
         parents = character(), children = list(),
         env = new.env(),
         n = 200)
{
    # To make things reasonably fast, we store the nodes in a pre-allocated list

  env = structure(env, class = c("XMLFlatListTree", "XMLFlatTree"))
  
  assign("nodeSet", vector("list", n), env)
  assign("idx", 1, env)
  assign("parentCount", 0, env)

  assign("nodeNames", character(n), env)  
  assign("parents", character(n), env)


  #XXX Deal with this if parents is specified.  
  
  # Assign the parents and children values and fill in any orphans, etc.
  # after calling addNode for the different nodes.
  
  if(!exists(".nodes", env))
    env$.nodes <- env #?

    # function to generate a new node identifier.  Can be given the
    # proposed name and will then make one up if that conflicts with another
    # identifier.
  f = function(suggestion = "") {
     if(suggestion == "" || suggestion %in% nodeNames)
        as.character(idx + 1)  
     else
        suggestion
  }
  environment(f) = env
  
  assign( ".nodeIdGenerator", f, env)

  
  g = addParentNode  
  environment(g) = env
  assign(".addParentNode", g, env)

  
  assign(".this", env, env)
  assign("n", n, env)

  
  addNode = function(node, parentId) {
    node = asXMLTreeNode(node, .this)
    id = node$id

       # Put it in the nodeSet by position.
    nodeSet[[ idx ]] <<- node
    nodeNames[idx] <<- id
    
    idx <<- idx + 1

    if(inherits(parentId, "XMLTreeNode"))
      parentId = parentId$id
    
    if(length(parentId)) {
      parentCount <<- parentCount + 1
      .parents[ parentCount ] <<- parentId
      names(.parents)[parentCount] <<- id

      .children [[ parentId ]] <<- c(.children[[ parentId ]], id )
    }    
    if(idx == n) {
      n <<- 2*n
      length(nodeSet) <<- n
    }
      
    return(node)
  }
  environment(addNode)  = env
  env$.addNode <- addNode

  # Populate the tree with any initial nodes.
  # XXX putting these in .nodes and not nodeSet!
  ids = names(nodes)
  nodes = lapply(seq(along = nodes),
                  function(i) {
                         x = nodes[[ i ]]
                         if(!("id" %in% names(unclass(x))))
                            x$id = f( ifelse(ids[ i ] == "", xmlName(x), ids[i]) )

                         if(!inherits(x, "XMLTreeNode")) {
			    ## no 'e' is visible here
                            x$env = e
                            class(x) = c("XMLTreeNode", class(x))
                         }
                         x
                       })

  names(nodes) = sapply(nodes, function(x) x$id)
  env$.nodes <- nodes

  env$.parents = parents
  env$.children = children

  .tidy =
     # to be run when adding to the tree is complete.
     # This shrinks the vectors to their actual size
     # rather than their preallocated sizes.
   function() {
      idx <- idx - 1
      length(nodeSet) <- idx
      length(nodeNames) <- idx 
      names(nodeSet) <- nodeNames
      .nodes <<- nodeSet
      idx
   }
  .tidy
  environment(.tidy) <- env
  env$.tidy = .tidy

  env
}


xmlRoot.xmlFlatListTree =
function(x, skip = TRUE, ...)
{
  #XXX
   stop("not implemented")
}  



# Represent the tree as a flat collection of nodes
# combined with

# See tests/tree.R

# Use an environment within the node so that we can lookup the children and parent information
#  directly from within

#
#  provide tools to set parent and children relationship.
#
#  Validate entries for parents and children to ensure nodes exist.
#
#  as(, "XMLTreeNode") function to make certain environment, id and class are present.
#
#  Suppose we are given an empty xmlTree() object when parsing an XML document.
# Then when we are converting the results back to R, we need to add nodes as we traverse the tree.
#  Need to make no
#   see convertNode() called in createXMLNode()
#  Given out an id within this tree for each node
#



xmlFlatTree =
  #
  # This version just concatenates each node to an existing list and so suffers
  # horrifically from garbage collection.
  # We leave it here in case it is useful either directly to someone for use on
  # small documents, or for performance comparisons with other approaches.
  #
function(nodes = list(), parents = character(), children = list(), env = new.env())
{
  # Assign the parents and children values and fill in any orphans, etc.
  # after calling addNode for the different nodes.
  
  if(!exists(".nodes", env))
    env$.nodes <- env

    # function to generate a new node identifier.  Can be given the
    # proposed name and will then make one up if that conflicts with another
    # identifier.
  f = function(suggestion = "") {
     if(suggestion == "" || suggestion %in% names(.nodes))
        as.character(length(.nodes) + 1)
     else
        suggestion
  }
  environment(f) = env
  
  assign( ".nodeIdGenerator", f, env)

  g = addParentNode  
  environment(g) = env
  assign(".addParentNode", g, env)

  assign(".this", env, env)

  addNode = function(node, parentId) {
    node = asXMLTreeNode(node, .this)
    id = node$id
    
    if(length(parentId)) {
      .parents[ id ] <<- parentId
      .children [[ parentId ]] <<- c(.children[[ parentId ]], id )
    }
    .nodes[[ id ]] <<- node

    id
  }
  environment(addNode)  = env
  env$.addNode <- addNode
  
  ids = names(nodes)
  nodes = lapply(seq(along = nodes),
                  function(i) {
                         x = nodes[[ i ]]
                         if(!("id" %in% names(unclass(x))))
                            x$id = f( ifelse(ids[ i ] == "", xmlName(x), ids[i]) )

                         if(!inherits(x, "XMLTreeNode")) {
				## FIXME: there is no visible 'e' here
                            x$env = e
                            class(x) = c("XMLTreeNode", class(x))
                         }
                         x
                       })

  names(nodes) = sapply(nodes, function(x) x$id)
  env$.nodes <- nodes

  env$.parents = parents
  env$.children = children
  
  structure(env, class = c("XMLSimpleFlatTree", "XMLFlatTree"))
}





