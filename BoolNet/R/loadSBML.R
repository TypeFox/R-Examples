# Check whether N is an integer using
# a regular expression.
check.integer <- function(N)
{
    !length(grep("[^[:digit:]]", as.character(N)))
}

# Parse an SBML species list in <rootNode>.
# Returns a list containing a map of gene ids to gene names,
# a vector specifying which genes are fixed
# and a vector of initial levels.
parseSBMLSpecies <- function(rootNode)
{
  genes <- c()
                  
  fixed <- c()
  
  initialLevels <- c()

  # iterate over species
  for (species in xmlChildren(rootNode))
  {
    attrs <- xmlAttrs(species)
    
    # if available, use the "name" attribute as the gene name
    # (and replace unsuitable characters by underscores)
    # otherwise use the gene id.
    if (is.na(attrs["name"]) || attrs["name"] == "")
      # for better compatibility with CoLoMoTo, 
      # remove "s_" at the beginning of identifiers
      name <- sub("^s_","",attrs["id"])
    else
      name <- gsub("[^a-zA-Z0-9_]+","_",attrs["name"])
     
    # gene names must be unique 
    if (name %in% genes)
    {
      suffix <- 1
      while (paste(name,suffix,sep="_") %in% genes)
        suffix <- suffix + 1
        
      warning(paste("Duplicate gene \"",name,"\", renaming to ",
                    name,"_",suffix,"!",sep=""))
                    
      name <- paste(name,suffix,sep="_")
    }
    
    # Reject logical networks with more than two values for a gene
    if (!is.na(attrs["maxLevel"]) && as.integer(attrs["maxLevel"]) > 1)
      stop(paste("BoolNet supports only binary genes, but gene \"",
                  name,"\" has a maximum level of ",attrs["maxLevel"],"!", sep=""))
    
    # build a lookup table id -> gene name  
    genes[attrs["id"]] <- adjustGeneNames(name)
    
    if (!is.na(attrs["constant"]) && tolower(attrs["constant"]) == "true")
    # if the gene is constant, save its initial level in the "fixed" vector
    {
      #if (is.na(attrs["initialLevel"]))
        #stop(paste("Gene \"", name, "\" is constant, but no initial level is supplied!", sep=""))
      #  warning(paste("Gene \"", name, "\" is constant, but no initial level is supplied! Assuming an input!", sep=""))
      
      fixed[name] <- TRUE 
    }
    else
    # this gene is not constant
      fixed[name] <- FALSE
          
    initialLevels[name] <- as.integer(attrs["initialLevel"])
  }
  return(list(genes = genes, fixed = fixed, initialLevels = initialLevels))
}

# Parse a list of transitions in <rootNode>.
# Here, <genes> is the result of parseSBMLSpecies()
# containing a map of gene ids/assignments, a vector
# specifying fixed genes and a vector of initial levels.
# Returns a list of interactions in the format 
# of class BooleanNetwork
parseSBMLTransitions <- function(rootNode, genes, symbolic)
{ 
  # iterate over all transitions 
  transitions <- xmlApply(rootNode,function(transition)
  {
    # parse inputs
    inputList <- xmlFindNode(transition, "listOfInputs")
    if (is.null(inputList))
    {
      inputs <- c()
    }
    else
    {
      inputList <- inputList[[1]]
      inputs <- c()
      inputThresholds <- c()
      
      # iterate over inputs
      for (input in xmlChildren(inputList))
      {
        attrs <- xmlAttrs(input)
        id <- attrs["qualitativeSpecies"]

        # verify gene name in species list
        if (is.na(genes$genes[id]))
          stop(paste("Unknown input \"",id,"\"!",sep=""))
        else
          inputs <- union(inputs, id)
        
        # check whether attributes of the input are compatible 
        # with Boolean logic
        if (tolower(attrs["transitionEffect"]) != "none")
          stop(paste("Transition effect for input gene \"",id,"\" is \"",
               attrs["transitionEffect"],"\", expected \"none\"!", sep=""))

        if (!is.na(attrs["thresholdLevel"]) && !is.na(attrs["id"]))
        # if a threshold level has been specified, save it with the 
        # corresponding ID for later use in the MathML terms
        {
          inputThresholds[attrs["id"]] <- as.integer(attrs["thresholdLevel"])
          
          if (inputThresholds[attrs["id"]] > 1)
            stop("Threshold levels must be 0 or 1!")        
        }
          
      }
    }
    
    # parse outputs
    outputList <- xmlFindNode(transition, "listOfOutputs")
    outputs <- c()
    
    if (!is.null(outputList))
    {
      outputList <- outputList[[1]]
      for (output in xmlChildren(outputList))
      # iterate over outputs
      {
        attrs <- xmlAttrs(output)
        id <- attrs["qualitativeSpecies"]
        
        # verify gene list in species list
        if (is.na(genes$genes[id]))
          stop(paste("Unknown output ",id,"!",sep=""))
        else  
          outputs <- union(outputs, id)
        
        # check whether attributes of the output are compatible 
        # with Boolean logic
        if (tolower(attrs["transitionEffect"]) != "assignmentlevel")
          stop(paste("Transition effect for output gene \"",id,"\" is \"",
               attrs["transitionEffect"],"\", expected \"assignmentLevel\"!", sep=""))
               
        if (!is.na(attrs["outputLevel"]))
          stop("Output levels for transitions are not supported in Boolean models!")
      }
    }
    
    # parse function terms
    functionTermList <- xmlFindNode(transition, "listOfFunctionTerms", throwError=TRUE)[[1]]
    transitionFunction <- parseSBMLFunctionTerms(functionTermList, 
                                                 genes$genes[inputs],
                                                 inputThresholds)
    
    return(list(inputs=inputs, outputs=outputs, transitionFunction=transitionFunction))
  })
  
  # now convert the read data to the BoolNet interaction format
  interactions <- lapply(names(genes$genes), function(gene)
  {
    # identify the transitions linked to each gene
    linkedTransitions <- which(sapply(transitions, 
                                      function(transition)gene %in% transition$output))
    
    if (length(linkedTransitions) == 0)
    # no transitions are assigned to this gene
    {
      if (!genes$fixed[genes$gene[gene]] || is.na(genes$initialLevels[genes$gene[gene]]))
      {
        if (is.na(genes$initialLevels[genes$gene[gene]]))
        {
          # Assume an input if the gene has no transition function and no initial value
          warning(paste("There is no transition and no initial level for gene \"",
                     gene,"\"! Assuming an input!",sep=""), call.=FALSE)
                     
        if (symbolic)
          return(parseBooleanFunction(genes$gene[gene]))
        else
          return(list(input=which(names(genes$genes) == gene), 
                      func=c(0,1), 
                      expression=genes$gene[gene]))
        }
        else
        if (!genes$fixed[genes$gene[gene]])
          warning(paste("There is no transition for the non-constant gene \"",
                     gene,"\"! Setting it to a constant ", 
                     genes$initialLevels[genes$gene[gene]], "!" ,sep=""), call.=FALSE)
      }
      
      # build a constant interaction
      
      if (symbolic)
        return(parseBooleanFunction(as.character(genes$initialLevels[genes$gene[gene]])))
      else
        return(list(input=0, 
                    func=genes$initialLevels[genes$gene[gene]], 
                    expression=as.character(genes$initialLevels[genes$gene[gene]])))
    }
    else
    {
    
      if (length(linkedTransitions) > 1)
      # multiple transitions per gene are not allowed, as these may be conflicting
      {
        stop(paste("Gene \"",gene,"\" is affected by multiple transitions!",sep=""))
      }
      
      if (genes$fixed[genes$gene[gene]])
      # a constant gene should not be the output of a transition
        stop(paste("Gene \"",gene,"\" has been specified as constant, but there is a transition!",sep=""))

      if (symbolic)
      {
        # parse the Boolean expression, and construct a symbolic expression tree
        return(parseBooleanFunction(transitions[[linkedTransitions]]$transitionFunction, genes$genes))
      }
      else
      {
        # parse the Boolean expression, and generate an interaction with 
        # the corresponding truth table
        return(generateInteraction(transitions[[linkedTransitions]]$transitionFunction, 
                                   #genes$genes[transitions[[linkedTransitions]]$input], 
                                   genes$genes))
      }
    }
  })
  names(interactions) <- genes$genes
  return(interactions)
}

# Parse a list of function terms in <rootNode>, where <genes> specifies
# the assignment of gene identifiers to gene names,
# and <inputThresholds> specifies the assignment of threshold identifiers
# to values (see also parseMathML()).
# Returns a single character string representing the function term 
# as an R expression.
parseSBMLFunctionTerms <- function(rootNode, genes, inputThresholds)
{
  # iterate over function terms
  functionTerms <- xmlApply(rootNode, function(term)
  {
    attrs <- xmlAttrs(term)
    
    outputLevel <- as.integer(attrs["resultLevel"])
    
    if (outputLevel > 1)
      stop("The result level of a function term must be 0 or 1!")
    
    if (tolower(xmlName(term)) == "defaultterm")
    # this is the default term => no expression
    {
      return(list(term="",outputLevel=outputLevel))
    }
    else
    # parse the MathML expression in the function term    
    {
      math <- xmlFindNode(term, "math", throwError=TRUE)[[1]][[1]]
              
      return(list(term=parseMathML(math, genes, inputThresholds), 
                  outputLevel=outputLevel))
    }
  })
  
  # build lists of terms with result 0 (negative terms)
  # and terms with result 1 (positive term)
  posTerms <- c()
  negTerms <- c()
  defaultValue <- NA
  
  for (term in functionTerms)
  {
    if (term$term != "")
    {
      if (term$outputLevel == 0)
      {
        negTerms <- c(negTerms, term$term)
      }
      else
      {
        posTerms <- c(posTerms, term$term)
      }
    }
    else
      defaultValue <- term$outputLevel

  }
  
  if (is.na(defaultValue))
    stop("Missing default term in transition!")
  
  if (defaultValue == 0)
  # if the default is 0, the result is a DNF of the positive terms
  {
    if (length(posTerms) > 0)
      totalTerm <- paste(posTerms, collapse=" | ")
    else
      totalTerm <- "0"
    
    if (length(negTerms) > 0)
    # if the default value is 0, additional negative terms are ignored, 
    # as they should be part of the default value 
    # (otherwise this is a contradiction to the positive terms)
      warning("Potentially contradictory terms in a transition have been ignored!")
  }
  else
  # if the default value is 1, the result is a negated DNF of the negative terms
  {
    if (length(negTerms) > 0)
      totalTerm <- paste("!(",paste(negTerms, collapse=" | "),")",sep="")
    else
      totalTerm <- "1"
      
    if (length(posTerms) > 0)
    # if the default value is 1, additional positive terms are ignored, 
    # as they should be part of the default value 
    # (otherwise this is a contradiction to the negative terms)    
      warning("Potentially contradictory terms in a transition have been ignored!")
  }
  return(totalTerm)
}

# Recursively parse the MathML expression in <rootNode>.
# Here, valid identifiers are the gene names in <names(genes)>
# and the input thresholds in <names(inputThresholds)>,
# which are replaced by the corresponding values.
# Returns an R expression string representing the MathML expression.
parseMathML <- function(rootNode, genes, inputThresholds)
{
  name <- xmlName(rootNode)
  if (name == "apply")
  # a bracket
  {
    operator <- xmlName(xmlChildren(rootNode)[[1]])
    
    children <- sapply(xmlChildren(rootNode)[-1],parseMathML, genes, inputThresholds)

    if (operator == "and" || operator == "times")
    {
      # treat "and" and "times" equally, but warn    
      if (operator == "times")
        warning("Interpreting \"times\" operator as a logical \"and\"!")
        
      return(paste("(",paste(children,
                   collapse = " & "),")",sep=""))
    }
    else
    if (operator == "or" || operator == "plus")
    {
      # treat "or" and "plus" equally, but warn    
      if (operator == "plus")
        warning("Interpreting \"plus\" operator as a logical \"or\"!")
        
      return(paste("(",paste(children,
                   collapse = " | "),")",sep=""))
    }
    else
    if (operator == "xor")
    {
      # convert XOR to a DNF by pasting all odd entries in the truth table
      tt <- allcombn(2, length(children)) - 1
      tt <- apply(tt,1,function(x)sum(x) %/% 2  == 1)
      return(paste("(",getDNF(tt, children),")",sep=""))
    }
    else
    if (operator %in% c ("eq", "neq", "gt", "lt", "geq", "leq"))
    {
      # comparison operators have to be converted to Boolean logic
      
      if (length(children) != 2)
        stop(paste("Operator \"",operator,"\" requires two operands!",sep=""))

      # check which of the child expressions are constant
      isConst <- sapply(children,function(x)
                        {
                          check.integer(x)
                        })
                        
      if (all(isConst))
      # two constants are compared (this does not really make sense!)
      {
        children <- as.integer(children)
        return(as.integer(switch(operator,
          eq = {children[1] == children[2]},
          neq = {children[1] != children[2]},
          gt = {children[1] > children[2]},
          lt = {children[1] < children[2]},
          geq = {children[1] >= children[2]},
          leq = {children[1] <= children[2]}
        )))
      }
      else
      if (any(isConst))
      # one constant and one variable are compared
      {
        constChild <- as.integer(children[isConst])
        varChild  <- children[!isConst]
        return(switch(operator,
          eq = 
          {
            if (constChild == 1)
              varChild
            else
              paste("!",varChild,sep="")
          },
          neq = 
          {
            if (constChild == 0)
              varChild
            else
              paste("!",varChild,sep="")
          },
          gt = 
          {
            if (constChild == 1)
              "0"
            else
              varChild
          },
          lt = 
          {
            if (constChild == 0)
              "0"
            else
              paste("!",varChild,sep="")
          },
          geq = 
          {
            if (constChild == 0)
              "1"
            else
              varChild
          },
          leq = 
          {
            if (constChild == 1)
              "1"
            else
              paste("!",varChild,sep="")
          }
        ))
      }
      else
      # two variables are compared
      {
        return(switch(operator,
          eq = 
          {
            paste("((",children[1]," & ",children[2],") | ",
                  "(!",children[1]," & !",children[2],"))", sep="")
          },
          neq = 
          {
            paste("((",children[1]," & !",children[2],") | ",
                  "(!",children[1]," & ",children[2],"))", sep="")
          },
          gt = 
          {
            paste("(",children[1]," & !",children[2],")", sep="")
          },
          lt = 
          {
            paste("(!",children[1]," & ",children[2],")", sep="")
          },
          geq = 
          {
            paste("(",children[1], " | !",children[2],")", sep="")
          },
          leq = 
          {
            paste("(!",children[1], " | ",children[2],")", sep="")
          }
        ))
      }
    }
    else
    if (operator == "not")
    {
      if (length(children) > 2)
        stop("Multiple arguments supplied to unary operator \"neg\"!")
      return(paste("!",children,
                   sep=""))
    }
    else
    # an unsupported symbol has been specified
      stop(paste("Unsupported math symbol: ",
                 operator,"!",sep=""))
  }
  else
  if (name == "ci")
  # this is a gene identifier or a threshold level
  {
    id <- trim(xmlValue(rootNode))
    if (!(id %in% names(genes)))
    {
      if (!(id %in% names(inputThresholds)))
      # neither threshold identifier nor gene
        stop(paste("Unspecified input \"",id,"\" in transition function!",sep=""))
      else
      # this is a threshold identifier
        return(inputThresholds[id])
    }
    # this is a gene
    return(genes[id])  
  }
  else
  if (name == "cn")
  # this is a constant
  {
    # convert value and ensure it is Boolean
    attrs <- xmlAttrs(rootNode)
    
    if (!is.null(attrs) && !is.na(attrs["type"]) && tolower(attrs["type"]) != "integer")
      stop("\"cn\" nodes must be of type \"integer\"!")
      
    val <- trim(xmlValue(rootNode))
    
    if (!check.integer(val) || !(as.integer(val) %in% c(0,1)))
      stop("\"cn\" nodes must be 0 or 1!")
    
    return(val)  
  }
  else
  if (name == "true")
    return("1")
  else
  if (name == "false")
    return("0")
  else
  # an unsupported symbol has been specified
    stop(paste("Unsupported math symbol: ",
               name,"!",sep=""))
}

# Import the sbml-qual document <file>
loadSBML <- function(file, symbolic=FALSE)
{
  # load XML document  
  doc <- xmlRoot(xmlParse(file))
  
  # remove comments from the document
  suppressWarnings(comments <- getNodeSet(doc,"//comment()"))
  
  if (length(comments) > 0)
    removeNodes(comments)
  
  # do various checks to ensure this is an sbml-qual document
  if (xmlName(doc) != "sbml")
    stop("Not an SBML document!")
    
  if (as.integer(xmlAttrs(doc)["level"]) > 3 || 
      as.integer(xmlAttrs(doc)["version"]) > 1)
      warning("This import is designed for documents up to SBML level 3 version 1!")
  
  if (is.null(xmlNamespaces(doc)$qual))
    stop("This document does not import the sbml-qual package!")
  
  model <- xmlFindNode(doc, "model", throwError=TRUE)[[1]]
  
  # compartments are ignored
  
  # parse species
  species <- xmlFindNode(model, "listOfQualitativeSpecies", throwError=TRUE)[[1]]
  genes <- parseSBMLSpecies(species)
  
  # parse transitions  
  transitions <- xmlFindNode(model, "listOfTransitions", throwError=TRUE)[[1]]
  interactions <- parseSBMLTransitions(transitions, genes, symbolic)
  
  if (symbolic)
  {
    delays <- apply(sapply(interactions,maxTimeDelay,genes=genes$genes),1,max)

    fixed <- as.integer(rep(-1L,length(genes$genes)))
    names(fixed) <- genes$genes
    
    res <- list(genes = genes$genes, interactions=interactions, fixed=fixed)
    res$internalStructs <- .Call("constructNetworkTrees_R",res)
    res$timeDelays <- delays
    
    class(res) <- "SymbolicBooleanNetwork"
  }
  else
  {
    # create BooleanNetwork structure
    res <- list(genes = genes$genes,
                fixed = sapply(interactions,function(i)
                              {
                                if (i$input[1] == 0)
                                  i$func[1]
                                else
                                  -1
                              }),
                interactions = interactions)
    class(res) <- "BooleanNetwork"
  }
  return(res)              
}
