#
#  This is an experiment to see if a simple hash of the values 
#  is faster.
#
#  Basically, we keep the parents, children and nodes
#  each as hash tables not a list.  Otherwise, this resembles
#  a flat tree
#
#  The notion is that we have a collection of nodes
#  and a collection of .parents and .children
#
#   Each element in .parents is assigned to the name of the node
#   whose parent is being stored. The value is the identifier for the
#   parent node. The top node has no entry in this collection.
#
#   The .children environment maintains a collection of entries
#   indexed by the identifier of the relevant node.
#   The value is a character vector containing the identifiers of the
#   nodes which  are children of this node.
#
xmlHashTree =
  #
  # Currently ignore the nodes, parents and children.
  #
function(nodes = list(), parents = character(), children = list(),
         env = new.env(TRUE, parent = emptyenv()))
{
    # function to generate a new node identifier.  Can be given the
    # proposed name and will then make one up if that conflicts with another
    # identifier.

  .count = 0

  # ability to be able to refer to this tree itself.
  # Not used here since the functions don't make use of the tt environment implicitly.
  # assign(".this", env, env)
  
    # We will store the children and parents as regular entries in these hash tables.
  env$.children = .children = new.env(TRUE)
  env$.parents = .parents = new.env(TRUE)

   #XXX we can do without this and make it a regular function
   # but we need to deal with XMLFlatListTree slightly different.
  f = function(suggestion = "") {
       # the check to see if suggestion is a name in env is very expensive.
     if(suggestion == "" || exists(suggestion, env, inherits = FALSE)) 
        as.character(.count + 1)  # can use length(tt)
     else
        suggestion
  }
  assign(".nodeIdGenerator", f, env)


  addNode =
     # This adds each new node.  
  function(node, parent = character(), ..., attrs = NULL, namespace = NULL,
           namespaceDefinitions = character(),
           .children = list(...),
           cdata = FALSE,
           suppressNamespaceWarning = getOption('suppressXMLNamespaceWarning', FALSE))
   {

    if(is.character(node))
       node = xmlNode(node, attrs = attrs, namespace = namespace, namespaceDefinitions = namespaceDefinitions)

     .kids = .children
     .children = .this$.children
     
     node = asXMLTreeNode(node, .this, className = "XMLHashTreeNode")
     id = node$id

     assign(id, node, env)
     .count <<- .count + 1

       # if no parent, 
     if(!inherits(parent, "XMLNode") && (!is.environment(parent) && length(parent) == 0) || parent == "")
       return(node)
    
     if(inherits(parent, "XMLHashTreeNode"))
        parent = parent$id     

    if(length(parent)) {
        assign(id, parent, envir = .parents)
        if(exists(parent, .children, inherits = FALSE))
           tmp = c(get(parent, .children), id)
        else
           tmp = id
        assign(parent, tmp, .children)
     }    
      
    return(node)
   }
   env$.addNode <- addNode

  
  # Create a .nodes vector with the names
  # of the node. And then makes a
  #

  .tidy = function() {
      idx <- idx - 1
      length(nodeSet) <- idx
      length(nodeNames) <- idx 
      names(nodeSet) <- nodeNames
      .nodes <<- nodeSet
      idx
  }
  #  environment(env$.tidy) <- env
  
  .this = structure(env, class = oldClass("XMLHashTree"))

  .this
}


# Example of looping over all elements
# table(unlist(eapply(tree, xmlName)))

getDescendants =
  #
  # This is trying to avoid recursion  and use iteration.
  #
function(id, tree, kids = tree$.children)
{
    # if no kids, then return empty list
  if(!exists(id, kids))
     return(character())

  ans = character()
  tmp = get(id, kids)
  hasKids = objects(kids)
  while( length( tmp  ) > 0) {
    ans = c(ans, tmp)
    tmp = tmp[ tmp %in% hasKids ]
    k = get(tmp[1], kids)
    tmp = c(tmp[-1], k)
  }
  ans
}

  


getDescendants =
  # Simple mechanism, for xmlHashTree trees.
  # This is recursive.
function(id, tree, kids = tree$.children)
{
   if(inherits(id, "XMLHashTreeNode")) {
      if(missing(tree))
         tree = id$env
     id = id$id
   }

  if(!exists(id, kids))
     return(character())

   ans = get(id, kids)
   c(ans, unlist(lapply(ans, getDescendants, tree, kids)))
   #Debugging
   # names(ans) = sapply(ans, function(i) xmlName(get(i, tree)))
}


subtree = copyXMLHashSubTree =
function(node)
{
    # find all the nodes below this node, i.e. in the subtree
  tree = node$env  
  ids = getDescendants(node$id, tree, tree$.children)

  newTree = xmlHashTree()

    # Now copy the parent & children information to the new tree
    # and also the modified nodes. The only modification necessary
    # is to set the env field of the original node to the new tree.
  sapply(c(node$id, ids), function(id) {
                 n = get(id, tree)
                 n$env = newTree
                 assign(id, n, newTree)
                 if(exists(id, tree$.children))
                    assign(id, get(id, tree$.children), newTree$.children)
                 if(exists(id, tree$.parents))                 
                    assign(id, get(id, tree$.parents), newTree$.parents)
              })
  
  remove(list = node$id, envir = newTree$.parents)
  newTree
}



xmlNamespaceDefinitions.XMLAbstractDocument =
function(x, addNames = TRUE, recursive = FALSE, simplify = FALSE, ...)  
{
   xmlNamespaces(as(x, "XMLAbstractNode"))
}

setAs("XMLAbstractDocument", "XMLAbstractNode",
        function(from)
           xmlRoot(from))

setAs("XMLHashTreeNode", "XMLHashTree",
       function(from)
         from$env
     )

"$.XMLHashTree" =
function(x, name)
  get(name, x, inherits = FALSE)


setMethod("xmlParent", "XMLHashTreeNode",
  # To get the parent of the node 'obj', we have to look in the .parents object
  # for the variable with obj's node identifier and then get the corresponding
  # value which is the identifier of the parent.
function(x, ...)
{
  p = get(".parents", x$env)
  idx = exists(x$id, p, inherits = FALSE)
  if(!idx)
      return(NULL)

  get(get(x$id, p), x$env)
}  )


xmlChildren.XMLHashTreeNode =
  #
  # For a given node 'obj', we have to use its id to find the entry
  # in the .children hash table and then the resulting entry is a character
  # vector giving the ids of the child nodes of obj. So we have to resolve those
  # children id's back in the hash table for the actual nodes.
function(x, addNames = TRUE, ...)
{
  e = x$env
  kids = get(".children", e)

  if(exists(x$id, kids, inherits = FALSE)) {
    ids = get(x$id,  kids, inherits = FALSE)
    nodes = lapply(ids, get, e, inherits = FALSE)
    names(nodes) = sapply(nodes, xmlName)
    nodes
  } else
    list()
}

if(useS4)
  setMethod("xmlChildren", "XMLHashTreeNode", xmlChildren.XMLHashTreeNode)

xmlSize.XMLHashTreeNode =
function(obj)
{
  length(xmlChildren(obj))
}


xmlSize.XMLHashTree =
function(obj)
{
  # 3 is the number of entries with a . prefix that we put there
  # for our own implementation purposes
  # We could calculate this as
  # length(grep("^.", objects(obj, all = TRUE))
  length(obj) - 3  
}  

#??? Currently overridden below
xmlRoot.XMLHashTree =
function(x, skip = TRUE, ...)
{  
  id = setdiff(objects(x), objects(x[[".parents"]]))
  get(id, x)
}


"[[.XMLHashTreeNode" =
function(x, ..., copy = FALSE, exact = TRUE)  
{
#   ans = NextMethod("[[")
   ans = xmlChildren(x)[[...]]  
   if(copy)
     xmlRoot(subtree(ans))
   else
     ans
}


addNode =
function(node, parent, to, ...)
{
  UseMethod("addNode", to)
}

addNode.XMLHashTree =
function(node, parent = character(), to, ...)
{
  to[[".addNode"]](node, parent, ...)
}

xmlRoot.XMLHashTree =
  #
  # This can return multiple roots
  #
  # Find all the identities of the nodes for which there is no 
  # corresponding entry n the .parents
  #
  #
  # If skip is TRUE, discard comment nodes. Leave  PI nodes, etc.
  #
  # If all is TRUE, return a list() with all the top-level nodes.
  # 
function(x, skip = TRUE, all = FALSE, ...)
{
  parents = get(".parents", x, inherits = FALSE)
  tops = objects(x)[ is.na(match(objects(x), objects(parents)))]

  if(length(tops) == 0)
    return(NULL)
  
  ans = mget(tops, x)

  if(skip)
    ans = ans[!sapply(ans, inherits, c("XMLCommentNode"))] #XXX names of XML hash tree nodes for comment, processing instruction, text node, etc.
  
  if(all)
    return(ans)

  ans[[1]]
}

getSibling =
  # Access the next field in the xmlNodePtr object.
  # not exported.
function(node, after = TRUE, ...)
  UseMethod("getSibling")

getSibling.XMLHashTreeNode =
function(node, after = TRUE, ...)
{
  .this = node$env
  parent = xmlParent(node)

  if(!is.null(parent)) {
      kids = xmlChildren(parent)
  } else
      kids = xmlRoot(.this, skip = FALSE, all = TRUE)

  i = match(node$id, sapply(kids, function(x) x$id))
  if(is.na(i))
    stop("shouldn't happen")

  if(after) {
    if(i < length(kids))
       kids[[i+1]]
    else
       NULL
  } else {
    if(i > 1)
      kids[[i-1]]
    else
      NULL
  }
}

print.XMLHashTree =
function(x, ...)
{
  print(xmlRoot(x), ...)
}  


xmlElementsByTagName.XMLHashTree =
  #
  # non-recursive version only at present
  #
function(el, name, recursive = FALSE)
{
  kids = xmlChildren(el)
  if(!recursive) 
    return(kids [ sapply(kids, xmlName) == name ])
}  


convertToHashTree =
function(from)
{
  xx = xmlHashTree()
  ans = .Call("R_convertDOMToHashTree", from, xx, xx$.children, xx$.parents, PACKAGE = "XML")
  docName(xx) = docName(from)
  xx
}  

setAs("XMLInternalDocument", "XMLHashTree",
       function(from) {
           convertToHashTree(xmlRoot(from, skip = FALSE))
       })

setAs("XMLInternalNode", "XMLHashTree",
       function(from) {
           ans = convertToHashTree(from)
           docName(ans) <- docName(from)
           ans           
       })

docName.XMLHashTree =
function(doc)
{
  if(exists(".doc", doc))
     doc$.doc
  else
     as.character(NA)
}

setMethod("docName", "XMLHashTree", docName.XMLHashTree)



