# Parse all nodes in <rootNode> 
# (i.e., all nodes of a certain type),
# and write them to a list
parseGeneAttrs <- function(rootNode)
{
    i <- 0
    res <- xmlApply(rootNode,function(child)
        {
            i <<- i + 1
            attrs <- xmlAttrs(child)
            
            name <- unname(attrs["name"])
            if (name == "" | is.na(name))
              name <- paste(xmlName(child),i,sep="_")
            
            # find the "simulationLogic" node, which contains
            # the necessary logic information
            logic <- xmlFindNode(child,"simulationLogic")
            
            initVal <- NULL
            if (is.null(logic))
            {
                warning(paste("Gene ",name," does not specify a simulation logic, assuming AND!",sep=""))
                logic <- "AND"
            }
            else
            {
              # try to find initial values in the simulation parameters
              logic <- logic[[1]]
              params <- xmlFindNode(logic,"simParameters")
              if (!is.null(params))
              {
                for (param in xmlChildren(params[[1]]))
                {
                  pa <- xmlAttrs(param)
                  if (pa["name"] == "initVal")
                  {
                    initVal <- as.numeric(pa["value"])
                    if (initVal != 0 & initVal != 1)
                    # only 0 and 1 are allowed for the initialization
                    {
                       warning("Initial value for gene ",name," is neither 0 nor 1 and is ignored!",sep="")
                       initVal <- NULL
                    }
                    break
                  }
                }
              }
              
              # get the Boolean function of the node
              logic <- unname(xmlAttrs(logic)["type"])
            }
            
            list(id=unname(attrs["id"]),
                 name=adjustGeneNames(name),logic=logic,initVal=initVal,inputs=NULL)
        })     
}

# Parse all types of nodes (genes, signals, etc)
parseAllGenes <- function(rootNode)
{
  genes <- unlist(xmlApply(rootNode,parseGeneAttrs),recursive=FALSE)
  names(genes) <- lapply(genes,function(gene)gene$id)
  return(genes)
}

# Parse the links in the children of <node> and insert them
# into the input lists of the genes
parseAllLinks <- function(rootNode,geneList)
{
  for (link in xmlChildren(rootNode))
  {
      attrs <- xmlAttrs(link)
      
      if (is.na(attrs["sign"]))
      # a link with unspecified direction (inhibition or activation) was found
      {
        warning(paste("The link between \"",geneList[[attrs["src"]]]$name,
                      "\" and \"",geneList[[attrs["targ"]]]$name,
                      "\" is not specified as an enhancer or suppressor and is therefore ignored!",sep=""))
        
      }
      else
      {
        geneList[[attrs["targ"]]]$input[[length(geneList[[attrs["targ"]]]$input) + 1]] <- 
            list(input=unname(attrs["src"]),sign=attrs["sign"])
      }    
  }
  return(geneList)
}

# Load a BioTapestry file (*.btp) from <file>
# and convert it to a Boolean network
loadBioTapestry <- function(file, symbolic=FALSE)
{
    doc <- xmlRoot(xmlTreeParse(file))
    
    # detect the "genome" XML node which holds the nodes and links of the network
    genome <- xmlFindNode(doc,"genome",throwError=TRUE)[[1]]
    
    # find the "nodes" node that contains all genes/inputs/etc
    nodes <- xmlFindNode(genome,"nodes",throwError=TRUE)[[1]]
    
    # read in genes
    geneList <- parseAllGenes(nodes)
    
    # find the "links" node that contains the links/dependencies
    links <- xmlFindNode(genome,"links",throwError=TRUE)[[1]]
    
    # add links to gene list
    geneList <- parseAllLinks(links,geneList)
    
    # build up network structure
    genes <- unname(sapply(geneList,function(gene)gene$name))
    
    geneIds <- names(geneList)
    
    fixed <- rep(-1L,length(genes))
    
    i <- 0
    interactions <- lapply(geneList,function(gene)
        {
            i <<- i + 1
            input <- sapply(gene$input,function(inp)
            {
                which(geneIds == inp$input) 
            })
            
            if (length(input) == 0)
            {
                if (!is.null(gene$initVal))
                # input-only gene without fixed value (=> depending on itself)
                {
                    input <- 0
                    func <- gene$initVal
                    fixed[i] <<- func
                    expr <- func
                } 
                else
                # input-only gene with fixed value
                {
                    input <- i
                    func <- c(0,1)
                    expr <- gene$name
                }
            }    
            else
            # dependent gene
            {
                # determine signs/negations of genes
                inputSigns <- sapply(gene$input,function(inp)inp$sign)  
        
                if (!symbolic || gene$logic == "XOR")
                {
                  tt <- as.matrix(allcombn(2,length(input)) - 1)
                              
                  func <- as.integer(switch(gene$logic,
                      AND = {
                              # calculate truth table for AND
                              apply(tt,1,function(assignment)
                              {
                                  res <- 1
                                  for (i in seq_along(assignment))
                                  {
                                      if (inputSigns[i] == "+")
                                      res <- res & assignment[i]
                                      else
                                      res <- res & !assignment[i]
                                      if (!res)
                                      break
                                  }
                                  res
                              })
                          },
                      OR = {
                              # calculate truth table for OR
                              apply(tt,1,function(assignment)
                              {
                                  res <- 0
                                  for (i in seq_along(assignment))
                                  {
                                      if (inputSigns[i] == "+")
                                      res <- res | assignment[i]
                                      else
                                      res <- res | !assignment[i]
                                      if (res)
                                      break
                                  }
                                  res
                              })
                          },
                      XOR = {
                              # calculate truth table for XOR
                              apply(tt,1,function(assignment)
                              {
                                  res <- assignment[1]
                                  for (i in 2:length(assignment))
                                  {
                                      res <- xor(res,assignment[i])
                                  }
                                  res
                              })
                          },
                      stop(paste("Unknown Boolean operator \"",
                                gene$logic,"\"!",sep=""))
                  ))
                }                
                # get string representation of the input gene literals
                literals <- mapply(function(gene,sign)
                                    {
                                    if (sign == "+")
                                        gene
                                    else
                                        paste("!",gene,sep="")
                                    }, genes[input], inputSigns)
                
                # get string representation of function                    
                expr <- switch(gene$logic,
                    AND = {
                            paste(literals,collapse=" & ")
                        },
                    OR =  {
                            paste(literals,collapse=" | ")
                        },
                    XOR = {
                            getDNF(func,genes[input]) 
                          }      
                        )
            }
            
            if (symbolic)
              return(parseBooleanFunction(expr,varNames=genes))
            else
              return(list(input=input,func=func,expression=expr))
        })
        
    names(interactions) <- genes
    fixed <- as.integer(fixed)
    names(fixed) <- genes
        
    net <- list(genes=genes,interactions=interactions,fixed=fixed)
    
    if (symbolic)
    {
      net$timeDelays <- apply(sapply(net$interactions,maxTimeDelay,genes=net$genes),1,max)
      net$internalStructs <- .Call("constructNetworkTrees_R",net)
      class(net) <- "SymbolicBooleanNetwork"
    }
    else    
      class(net) <- "BooleanNetwork"
    return(net)
}
