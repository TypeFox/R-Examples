
#XXX Deal with line numbers in the original document.

  #
  # We could also do things this way
  #  source(textConnection(saveXML(xsltApplyStyleSheet("cityTemperatures.xml", "~/Projects/org/omegahat/XML/Literate/segment.xsl"))))
  #
  #
  # Allow the user to specify a subset of nodes in which to find the code, etc. nodes.
  # Or an XPath query to restrict the search.
  # For example, suppose we have a document with two sections and we want to run the code
  # in only one of those sections.
  #  getNodeSet(section[@id='second one'])
  #

setOldClass("XMLNodeSet")

DefaultXPathNamespaces = 
                 c(r = "http://www.r-project.org",
                   s = "http://cm.bell-labs.com/stat/S4",
                   omg = "http://www.omegahat.net",
                 
                   mlb = "http://www.mathworks.com",  # matlab

                   sh="http://www.shell.org",
                   perl = "http://www.perl.org",
                   py = "http://www.python.org",  

                   fo="http://www.w3.org/1999/XSL/Format",
                   xsl="http://www.w3.org/1999/XSL/Transform",
                   xi="http://www.w3.org/2001/XInclude"                   
                  )

DefaultXMLSourceXPath =
   sprintf("%s[not(@eval='false') and not(ancestor::ignore) and not(ancestor::section[@r:eval = 'false'])]",
            c("//r:init", "//r:function", "//r:init", "//r:code", "//r:plot", "//r:expr"))

setGeneric("xmlSource",
function(url, ...,
          envir = globalenv(),
          xpath = character(),
          ids = character(),
          omit = character(),
          ask = FALSE,
          example = NA,
          fatal = TRUE, verbose = TRUE, echo = verbose, print = echo,
          xnodes = DefaultXMLSourceXPath,         
          namespaces = DefaultXPathNamespaces, section = character(), eval = TRUE, init = TRUE,
          setNodeNames = FALSE, parse = TRUE, force = FALSE)
{

  standardGeneric("xmlSource")
})



# Break up into methods for character and xmlInternalDocument


setMethod("xmlSource", c("character"),
function(url, ...,
          envir =globalenv(),
          xpath = character(),         
          ids = character(),
          omit = character(),
          ask = FALSE,
          example = NA,         
          fatal = TRUE, verbose = TRUE, echo = verbose, print = echo,
          xnodes = DefaultXMLSourceXPath,
          namespaces = DefaultXPathNamespaces,
          section = character(), eval = TRUE, init = TRUE, setNodeNames = FALSE, parse = TRUE, force = FALSE)
{
  doc = xmlTreeParse(url, ..., useInternalNodes = TRUE)
  xmlSource(doc, ..., envir = envir, xpath = xpath, ids = ids, omit = omit,
             ask = ask, example = example, fatal = fatal, verbose = verbose,
              print = print, xnodes = xnodes, namespaces = namespaces,
               section = section, eval = eval, init = init, setNodeNames = setNodeNames, parse = parse, force = force)
})


setMethod("xmlSource", c("XMLInternalDocument"),
function(url, ...,
          envir =globalenv(),
          xpath = character(),         
          ids = character(),
          omit = character(),
          ask = FALSE,
          example = NA,         
          fatal = TRUE, verbose = TRUE, echo = verbose, print = echo,
          xnodes = DefaultXMLSourceXPath,
          namespaces = DefaultXPathNamespaces,
          section = character(), eval = TRUE, init = TRUE, setNodeNames = FALSE, parse = TRUE, force = FALSE)
{
  doc = url
  if(inherits(verbose, "numeric"))
    verbose = verbose - 1

  if(!is.character(section))
    section = as.integer(section)

      #XXX use section when processing the examples
  if(length(example) && !all(is.na(example)))  {
    egs = getNodeSet(doc, "//r:example", namespaces)
    if(length(egs)) {
      ids = sapply(egs, xmlGetAttr, "id")
      if(length(example) == 1 && is.na(example)) {
        cat("Select an example\n")
        example = ids[w <- menu(ids)]
      } 

      if(inherits(example, "numeric")) {
        i = example
      } else {
        i = pmatch(example, ids)
        if(all(is.na(i)))
          stop("no example named ", example)
      }

       # find any r:init nodes which are not inside an example.
      init = getNodeSet(doc, "//r:init[not(ancestor::r:example)]",
                           c(r = "http://www.r-project.org"))
      if(length(init)) {
        xmlSource(init, envir = envir, omit = omit, verbose = verbose, namespaces = namespaces, eval = eval, force = force)
        cat("Done doc-level init", length(init), "\n")
      }

      ans = sapply(i, function(x) {
                        nodes = getNodeSet(egs[[x]], paste(xnodes, collapse = "|"), namespaces)
                        if(verbose) 
                          cat("Example", ids[x], "\n")

                             #XXX put the correct ids in her.
                        xmlSource(nodes, envir = envir, omit = omit, verbose = verbose, namespaces = namespaces, eval = eval, setNodeNames = setNodeNames, parse = parse, force = force)
                        
                      })
      return(ans)
   }
  }

  
#  if(length(section) && is.character(section))
#    section = paste("@id", ddQuote(section), sep = "=")
  
  if(length(xpath)) {
       # do an XPath query and then look inside the resulting nodes
       # for the xnodes of interest.

    if(length(section)) {
          # XXX assumes just one section. What about c(1, 2, 4)
       xpath =paste("//section[", section, "]", xpath, sep = "")
    }
    
    nodes = getNodeSet(doc, xpath, namespaces)
    v =
      unlist(lapply(nodes, function(n) {
                     unlist(lapply(xnodes,
                                   function(p)
                                      getNodeSet(n, p, namespaces)),
                             recursive = FALSE)
                  }), recursive = FALSE)
  } else {

    functions = limitXPathToSection(section, "//r:function[not(@eval = 'false') and not(ancestor::ignore)]")
    xnodes = limitXPathToSection(section, xnodes)
        # Do we need to ensure the order for the functions first?
    v = getNodeSet(doc, paste(c(functions, xnodes), collapse = "|"), namespaces)
#   v = getNodeSet(doc, functions, namespaces)
#   w = getNodeSet(doc, xnodes, namespaces)
#   v = c(v, w)    
  }

  if(is.null(v))
    stop("No matching nodes in the document found")
  
  class(v) <- "XMLNodeSet"



    # deal with a top-level node r:codeIds which is of the form
    #   abc
    #   def
    #   ghi
    # i.e. a single entry on each line which identifies the nodes that are to be read.
  if(missing(ids) && missing(xnodes) && length(ids <- getNodeSet(doc, "/*/r:codeIds|/*/invisible/r:codeIds",
                                               namespaces = c(r = "http://www.r-project.org")))) {

     if(length(ids) > 1) {
       warning("more than one r:codeIds node. Using the first one")
     }
#     txt = paste(sapply(ids, xmlValue))
     ids = strsplit(xmlValue(ids[[1]]), "\\\n")[[1]]
     ids = unique(ids)
     ids = ids[ids != ""]
   }

  xmlSource(v, ids = ids, omit = omit, ask = ask, fatal = fatal, verbose = verbose,  envir = envir,
            section = if(!is.character(section)) section else character(),
            eval = eval, setNodeNames = setNodeNames, parse = parse, force = force)
})


limitXPathToSection =
  #
  # limitToSection(1:3)
  # limitToSection(letters[1:3])
  # limitToSection(letters[1:3], "//r:plot")
function(section, xpath = c("//r:code", "//r:func", "//r:plot", "//r:expr"))
{
  if(length(section) == 0)
    return(paste(xpath, collapse = "|"))
  
  if(is.character(section))
    section = paste("@id=", sQuote(section), sep = "")
  
  paste(outer(section, xpath,
                         function(sect, xp)
                           paste("//section[", sect, "]", xp, sep = "")),
        collapse = "|")
}  




setMethod("xmlSource", "XMLNodeSet",
function(url, ..., envir =globalenv(),
          xpath = character(),         
          ids = character(),
          omit = character(),
          ask = FALSE,
          example = NA,         
          fatal = TRUE, verbose = TRUE, echo = verbose, print = echo,
          xnodes = c("r:function[not(@val='false')]", "r:init[not(@eval='false')]", "r:code[not(@eval='false')]", "//r:plot[not(@eval='false')]"),         
          namespaces =  DefaultXPathNamespaces, section = character(), eval = TRUE, init = TRUE, setNodeNames = FALSE, parse = TRUE, force = FALSE)
{
  if(ask) {
     doc = as(url[[1]], "XMLInternalDocument")     #XXXX  no doc here now.
     v = getNodeSet(doc, "//r:function|//r:init|//r:code|//r:plot", namespaces)
     funs = sapply(v, xmlName) == "function"
     if(any(funs)) {
#XXX
     }
   }

  ans = sapply(url, evalNode, envir = envir, verbose = verbose, ids = ids,
                omit = omit, echo = echo, print = print, ask = ask, eval = eval, parse = parse, force = force)

  if(setNodeNames)
     names(ans) = sapply(url, getRCodeNodeName)
  else
     names(ans) = sapply(url, xmlName, full = TRUE)
  invisible(ans)
})


evalNode = 
function(node, envir = globalenv(), ids = character(), verbose = FALSE, echo = verbose, omit = character(),
         namespaces = c(r = "http://www.r-project.org"), print = echo, ask = FALSE, eval = TRUE, parse = TRUE, force = FALSE)
{
            #XXX check all ancestors. Ideally exclude them in the XPath query
   if(!force && (xmlName(xmlParent(node)) == "ignore" ||
          length(getNodeSet(node, "./ancestor::section[@r:eval='false']|./ancestor::para[@r:eval='false']",
                              c(r = "http://www.r-project.org"))) > 0))
     return(FALSE)

   tmp = xmlGetAttr(node, "id", NA)
   if(is.na(tmp) && length(ids) > 0 && !("" %in% ids))
     return()
   
   if(!is.na(tmp)) {
     if(length(omit) > 0 && tmp %in% omit) {
       if(verbose)
         warning("skipping id ", tmp)
       return()
     } else if(length(ids) > 0 && !(tmp %in% ids)) {
       if(verbose)
         warning("ignoring id ", tmp)
       return()
     }
   }

   tmp = xmlGetAttr(node, "ignore", NA, converter = as.logical)
   if(!is.na(tmp) && tmp) {
     if(verbose)
        warning("ignoring node", as(node, "character"))
     return()
   }

     # go through the node and see if there are any r:code nodes
     # with a ref attribute
     # and go fetch the corresponding <r:code id='@ref'> node.
   txt = paste(getRCode(node, namespaces),  collapse = "\n")

   if(!parse)
     return(txt)
   
#   txt = xmlValue(node)
   if(verbose)
     cat("*************\nEvaluating node\n")
   cmd = parse(text = txt)
   if(echo)
     print(cmd)

   if(eval) {
     if(ask) {
       w = menu(c("evaluate", "skip", "terminate"))
       if(w == 2)
         return(NULL)
       else if(w == 3)
         stop("User terminated the xmlSource")
     }

     eval(cmd, envir)     
   } else
     cmd
}  


getRCode =
function(node, namespaces = c(r = "http://www.r-project.org"), recursive = TRUE)
{
 tmp = xmlSApply(node, function(x) {

  if(inherits(x, c("XMLInternalCommentNode", "XMLInternalPINode"))) {

  } else if(inherits(x, "XMLInternalElementNode") && xmlName(x, full = TRUE) %in% c("r:code", "r:frag")) {
     ref = xmlGetAttr(x, "ref", NA)
     if(!is.na(ref)) {
         v = getNodeSet(as(x, "XMLInternalDocument"),
                        paste(sapply(c("code", "frag"),
                                      function(x) paste("//r:", x , "[@id='", ref, "']", sep = "")), collapse = "|"),
                        namespaces)
         if(length(v) == 0)
           stop("No code block/fragment named ", ref)
         else if(length(v) > 1)
           stop("More than 1 code block/fragment named ", ref)         
         else
            if(recursive)
              getRCode(v[[1]], namespaces, recursive = TRUE)
            else
              xmlValue(v[[1]])
     } else {
         if(recursive)
              getRCode(x, namespaces, recursive = TRUE)
            else
              xmlValue(x)
     }
  } else if(inherits(x, "XMLInternalElementNode") && xmlName(x, full = TRUE) %in% c("r:error", "r:output")) {
  }  else
     xmlValue(x)
 })
 paste(tmp, collapse = "\n")
}



setClass("XMLCodeFile", contains = "character")
setClass("XMLCodeDoc", contains = "XMLInternalDocument")
setAs("XMLCodeFile", "XMLCodeDoc",
        function(from) {
           new("XMLCodeDoc", xmlParse(from))
        })

setAs("character", "XMLCodeFile",
        function(from) {
           xmlCodeFile(from)
        })

setAs("character", "XMLCodeDoc",
        function(from) {
           xmlCodeFile(from, TRUE)
        })

xmlCodeFile =
function(f, parse = FALSE)
{
  if(parse)
     new("XMLCodeDoc", xmlParse(f))
  else
     new("XMLCodeFile", f)
}


tmp.source =
function (file, local = FALSE, echo = verbose, print.eval = echo, 
    verbose = getOption("verbose"), prompt.echo = getOption("prompt"), 
    max.deparse.length = 150, chdir = FALSE, encoding = getOption("encoding"), 
    continue.echo = getOption("continue"), skip.echo = 0, 
    keep.source = getOption("keep.source"))
   {
      if(length(verbose) == 0)
        verbose = FALSE
      
      if(chdir) {
        cwd = getwd()
        on.exit(setwd(cwd))
        setwd(dirname(file))
      }

      xmlSource(file, verbose = verbose)
   }

if(compareVersion(as.character(getRversion()), "2.8.0") < 0) {
 cat("Fixing source definition\n")
 formals(tmp.source) =  formals(tmp.source)[ - length(formals(tmp.source)) ]
}

setMethod("source", "XMLCodeFile", tmp.source)


setMethod("[[", "XMLCodeFile",
          function(x, i, j, ...,  env = globalenv()) {
            doc = as(x, "XMLCodeDoc")
            n = getNodeSet(doc, paste("//*[@id=", sQuote(i), "]"))
            if(length(n) == 0) {
              # This needs code from ptoc to determine the name of an "element"
             doc = updateIds(doc, save = x)
            }

            eval(parse(text = xmlValue(n[[1]])), envir = env)
          })

updateIds =
function(doc)
{
   nodes = getNodeSet(doc,
                      "//r:function[not(@id) and not(@eval = 'false')]|//r:code[not(@id) and not(@eval = 'false')]",
                        c("r" = "http://www.r-project.org"))
   sapply(nodes, getCodeVar)
}

getCodeVar =
function(node)
{
     e = parse(text = getRCode(node))
     e = e[[length(e)]]
     # This should use the code in ptoc in RTools.
     id = if(class(e) %in% c("=", "<-"))
             as.character(e[[2]])
          else
             NA

     if(!is.na(id))
        addAttributes(node, id = id)
     id
}

#
# f = xmlCodeFile("~/Classes/stat242-08/Code/FormulaWeights/rpartScope.xml")
# source(f)
# f[["rpart.model"]]


setGeneric("xmlSourceFunctions",
            function(doc, ids = character(), parse = TRUE, ...) {
                  standardGeneric("xmlSourceFunctions")
            })


setMethod("xmlSourceFunctions", "character",
  #
  # evaluate the r:function nodes, or restricted to @id from ids.
  #
function(doc, ids = character(), parse = TRUE, ...)
{
  invisible(xmlSourceFunctions(xmlParse(doc), ids, parse = parse, ...))
})


sQuote =
function(x)
  sprintf("'%s'", as.character(x))

setMethod("xmlSourceFunctions", "XMLInternalDocument",
  #
  # evaluate the r:function nodes, or restricted to @id from ids.
  #
function(doc, ids = character(), parse = TRUE, setNodeNames = FALSE, ...)
{

  if(length(ids))
     nodes = getNodeSet(doc, paste("//r:function[", paste("@id", sQuote(ids), sep = "=", collapse = " or " ), "]"), c(r = "http://www.r-project.org"))
  else
     nodes = getNodeSet(doc, "//r:function", c(r = "http://www.r-project.org"))

  if(parse == FALSE)
     return(nodes)

  ans = xmlSource(nodes, ...)
  if(setNodeNames)
    names(ans) = sapply(nodes, getRCodeNodeName)
  
  invisible(ans)
})

getRCodeNodeName =
function(node)
{
   xmlGetAttr(node, "name", xmlGetAttr(node, "idx", getTaskId(node)))
}


################


setGeneric("xmlSourceSection",
           function(doc, ids = character(),
                      xnodes = c(".//r:function", ".//r:init[not(@eval='false')]", ".//r:code[not(@eval='false')]", ".//r:plot[not(@eval='false')]"), namespaces = DefaultXPathNamespaces, ...)
              standardGeneric("xmlSourceSection"))

setMethod("xmlSourceSection", "character",
           function(doc, ids = character(), xnodes = c(".//r:function", ".//r:init[not(@eval='false')]", ".//r:code[not(@eval='false')]", ".//r:plot[not(@eval='false')]"), namespaces = DefaultXPathNamespaces, ...)
              xmlSourceSection(xmlParse(doc), ids, xnodes, namespaces, ...))

setMethod("xmlSourceSection", "XMLInternalDocument",
           function(doc, ids = character(),
                      xnodes = c(".//r:function", ".//r:init[not(@eval='false')]", ".//r:code[not(@eval='false')]", ".//r:plot[not(@eval='false')]"),
                      namespaces = DefaultXPathNamespaces, ...) {
              nodes = getNodeSet(doc, "//section")
              
              aids = sapply(nodes, xmlGetAttr, "id", NA)
              m = pmatch(ids, aids)
              if(any(is.na(m))) {
                       # for those ids the caller gave us that didn't match
                       # compare these to the titles.
                 i = which(is.na(m))
                 tmp = ids[i]
                 j = pmatch(ids[i], sapply(nodes, function(x) {
                                                   tmp = getNodeSet(x, "./title")
                                                   if(length(tmp))
                                                     xmlValue(tmp[[1]])
                                                   else
                                                     ""
                                                  }))
                 m[i [!is.na(j)]] = j[!is.na(j)]
              }

              if(any(is.na(m)))
                stop("cannot match section id or title for ", paste(m[is.na(m)], collapse = ", "))

              lapply(nodes[m], evalSection, xnodes, namespaces, ...)
           })

evalSection =
function(node, xnodes, namespaces = DefaultXPathNamespaces, envir = globalenv(), ...)
{
    # Or use xnodes by stripping away any [] and .//
  if(xmlName(node, TRUE) %in% c("r:function", "r:plot", "r:code", "r:graphics"))
    return(evalNode(node, envir, ...))
  
  xpath = paste(xnodes, collapse = "|")
  nodes = getNodeSet(node, xpath, namespaces)
  sapply(nodes, evalNode, envir, ...)
}



#############################################

# r:code[@thread='name']|r:code[ancestor::*[@thread='name']]

setGeneric("xmlSourceThread",
           function(doc, id, envir = globalenv(), ...,
                    xnodes = c("r:function", "r:init", "r:code", "r:plot"))
              standardGeneric("xmlSourceThread"))

setMethod("xmlSourceThread", "character",
           function(doc, id, envir = globalenv(), ...,
                    xnodes = c("r:function", "r:init", "r:code", "r:plot"))                    
              xmlSourceThread(xmlParse(doc), id, envir, ..., xnodes = xnodes)
          )

setMethod("xmlSourceThread", "list",
           function(doc, id, envir = globalenv(), ...,
                    xnodes = c("r:function", "r:init", "r:code", "r:plot"))                                        
             sapply(doc, evalNode, envir = envir, ..., xnodes = xnodes))

if(FALSE)
setMethod("xmlSourceThread", "XMLNodeList",
           function(doc, id, envir = globalenv(), ...,
                    xnodes = c(".//r:function", ".//r:init[not(@eval='false')]", ".//r:code[not(@eval='false')]", ".//r:plot[not(@eval='false')]"))
             sapply(doc, evalNode, envir = envir, ...))


setMethod("xmlSourceThread", "XMLInternalDocument",
           function(doc, id, envir = globalenv(), ...,
                    xnodes = c("r:function", "r:init", "r:code", "r:plot"))  {

                 # all the nodes that are "under" this thread.
              xp = sprintf("//*[@thread='%s']", id)
              anc = sprintf("//%s[not(ancestor::*[@thread]) and not(ancestor::altApproach)]", xnodes)
              xp = paste(c(xp, anc), collapse = " | ")

              nodes = getNodeSet(doc, xp)

              sapply(nodes, evalSection, envir = envir, ..., xnodes = xnodes)
           })

setGeneric("xmlSourceTask",
            function(doc, id, envir = globalenv(), ...) {
              standardGeneric("xmlSourceTask")
            }
           )

#                nodes = c("r:code", "r:plot", "r:expr")
#                fmt = paste(nodes, "[@thread='%s']", sep = "")
#                xp = paste(sprintf(fmt, id), collapse = " | ")
#                paste(nodes[]
#                getNodeSet
#                "r:code[@thread='%s']|r:plot[@thread='%s']|r:expr['
#                paste(getNode
#           })


# See tangle.R and xmlTangle.
# Fix this up. Just thrown down one morning.
xmlToCode = tangle =
function(doc, file = stdout())
{
  e = xmlSourceFunctions(doc, eval = FALSE) # want to avoid parsing.
  if(!is(file, "connection"))
    con = file(file, "w")
  sapply(e, function(x) cat(x, "\n", file = file))
  file
}



getTaskId =
function(node) {
  els = getNodeSet(node, ".//ancestor::task")
  if(length(els))
    xmlGetAttr(els[[1]], "id")
  else
    ""
}
