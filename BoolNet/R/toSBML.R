# Indent <string> using <count> tabs.
indent <- function(string,count)
{
  if (count == 0)
    return(string)
  return(paste(paste(rep("\t",count),collapse=""), string, sep=""))
}

# Export a Boolean network <network> to an sbml-qual file <fileName>.
# If <generateDNFs>, a new symbolic representation for the interactions
# is generated on the basis of the truth tables (in disjunctive normal form).
# Otherwise, the $expression elements of the interactions are parsed.
# If <saveFixed> is true, constant transition functions are exported for fixed genes
# instead of their true transition functions
toSBML <- function(network, file, generateDNFs=FALSE, saveFixed = TRUE)
{
  symbolic <- inherits(network,"SymbolicBooleanNetwork")
  stopifnot(inherits(network,"BooleanNetwork") || symbolic)
    
  if (symbolic)
  {
    if (any(network$timeDelays > 1))
      stop("SBML does not support networks with time delays!")
      
    parseTrees <- network$interactions
  }
  else
  {
    parseTrees <- NULL
    
    if (generateDNFs == FALSE)
    # Check whether all interactions have suitable string representations
    {
      tryCatch(
      {
         # parse the string representations
         parseTrees <- lapply(network$interactions,
                             function(int)parseBooleanFunction(int$expression))
      },
      error=function(e)
      {
        warning(paste("The transition functions of this network did not contain valid symbolic expressions!",
                "Generating DNF representations from the truth tables!"))

        # There was an error parsing a string representation => generate DNFs              
        generateDNFs <<- TRUE
      })
    }
    
    if (generateDNFs != FALSE)
    # build new representations of the functions in disjunctive normal form
    {
      network$interactions <- lapply(network$interactions, function(interaction)
      {
        table <- allcombn(2, length(interaction$input)) - 1
        interaction$expression <- getDNF(interaction$func, 
                                         network$genes[interaction$input],
                                         generateDNFs)
        return(interaction)
      })
      
      # parse the DNF representations
      parseTrees <- lapply(network$interactions,
                          function(int)parseBooleanFunction(int$expression))
    }
    
    names(parseTrees) <- network$genes
  }  
  # generate a network identifier from the file name
  id <- sub(".sbml", "", basename(file), fixed=TRUE)
  id <- gsub("[^a-zA-Z0-9_]+","_",id)
  
  # open a string connection
  output <- NULL
  f <- textConnection("output",  encoding="UTF-8", open="w", local=TRUE)
  
  # write document header
  cat(file=f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
  cat(file=f, "<sbml xmlns=\"http://www.sbml.org/sbml/level3/version1/core\" level=\"3\" version=\"1\" xmlns:qual=\"http://www.sbml.org/sbml/level3/version1/qual/version1\" qual:required=\"true\">\n")
  cat(file=f, "\t<model id=\"", id, "\">\n", sep="")
  
  # write default compartment
  cat(file=f, "\t\t<listOfCompartments>\n")
  cat(file=f, "\t\t\t<compartment id=\"default\" constant=\"true\"/>\n")
  cat(file=f, "\t\t</listOfCompartments>\n")
  
  # write genes
  cat(file=f, "\t\t<qual:listOfQualitativeSpecies>\n")
  for (gene in network$genes)
  {
    if ((saveFixed && network$fixed[gene] != -1) || 
        (!symbolic && network$interactions[[gene]]$input[1] == 0))
    {
      if (saveFixed && network$fixed[gene] != -1)
        level <- network$fixed[gene]
      else
        level <- network$interactions[[gene]]$func[1]
        
      cat(file=f, "\t\t\t<qual:qualitativeSpecies qual:compartment=\"default\"",
                  " qual:constant=\"true\" qual:id=\"", gene,
                  "\" qual:name=\"", gene,
                  "\" qual:initialLevel=\"", level,
                  "\" qual:maxLevel=\"", level, 
                  "\"/>\n", sep="")
    }
    else
      cat(file=f, "\t\t\t<qual:qualitativeSpecies qual:compartment=\"default\"",
                  " qual:constant=\"false\" qual:id=\"", gene, 
                  "\" qual:name=\"", gene, "\" qual:maxLevel=\"1\"/>\n", sep="")
  }
  cat(file=f, "\t\t</qual:listOfQualitativeSpecies>\n")
  
  # write transition functions
  cat(file=f, "\t\t<qual:listOfTransitions>\n")
  
  for (gene in network$genes)
  {
    if ((!saveFixed || network$fixed[[gene]] == -1)  && 
        (symbolic || network$interactions[[gene]]$input[1] != 0))
    {
      cat(file=f, "\t\t\t<qual:transition qual:id=\"tr_", gene,
          "\" qual:name=\"Interactions targeting ", gene, "\">\n", sep="")
      cat(file=f, "\t\t\t\t<qual:listOfInputs>\n")
      
      if (symbolic)
        inputs <- getInputs(network$interactions[[gene]])
      else
        inputs <- network$genes[network$interactions[[gene]]$input]
      
      for (input in inputs)
        cat(file=f, "\t\t\t\t\t<qual:input qual:qualitativeSpecies=\"", input,
            "\" qual:transitionEffect=\"none\"/>\n", sep="")
      cat(file=f, "\t\t\t\t</qual:listOfInputs>\n")
      cat(file=f, "\t\t\t\t<qual:listOfOutputs>\n")
      cat(file=f, "\t\t\t\t\t<qual:output qual:qualitativeSpecies=\"", gene,
          "\" qual:transitionEffect=\"assignmentLevel\"/>\n", sep="")
      cat(file=f, "\t\t\t\t</qual:listOfOutputs>\n")
      
      cat(file=f, "\t\t\t\t<qual:listOfFunctionTerms>\n")
      cat(file=f, "\t\t\t\t\t<qual:functionTerm qual:resultLevel=\"1\">\n")
      cat(file=f, "\t\t\t\t\t\t<math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n")

      parseTree <- parseTrees[[gene]]
      cat(file=f, MathMLFromParseTree(parseTree, indentLevel=7))

      cat(file=f, "\t\t\t\t\t\t</math>\n")
      cat(file=f, "\t\t\t\t\t</qual:functionTerm>\n")
      cat(file=f, "\t\t\t\t\t<qual:defaultTerm qual:resultLevel=\"0\"/>\n")
      cat(file=f, "\t\t\t\t</qual:listOfFunctionTerms>\n")
      cat(file=f, "\t\t\t</qual:transition>\n")
    }
  }
  # finish document
  cat(file=f, "\t\t</qual:listOfTransitions>\n")
  cat(file=f, "\t</model>\n")
  cat(file=f, "</sbml>\n")
  close(f)
  
  # open file and write the complete XML string
  f <- file(file, encoding="UTF-8", open="w")
  cat(file=f,output,sep="\n")
  close(f)  
}

# Build a MathML representation of a parse tree <tree>
# that represents a symbolic Boolean expression.
# Indentation of the XML nodes starts with <indentLevel>.
# Returns a string with MathML code.
MathMLFromParseTree <- function(tree,indentLevel=0)
{
  res <- switch(tree$type,
    operator = 
    {
      if (tree$operator %in% c("timeis","timegt","timelt"))
        stop(sprintf("Operator \"%s\" not supported in SBML!", tree$operator));
      if (tree$operator %in% c("maj","sumis","sumlt","sumgt"))
      # convert special predicates to general Boolean formulae
        tree <- expandCountPredicate(tree)
        
      if (tree$negated)
        paste(indent("<apply>\n", indentLevel),
              indent("<not/>\n", indentLevel+1),
              indent("<apply>\n" ,indentLevel+1),
              indent({if (tree$operator=="|" || tree$operator=="any")
              {"<or/>\n"} else{"<and/>\n"}}, indentLevel+2),
              paste(sapply(tree$operands,MathMLFromParseTree, indentLevel+2), 
                    collapse=""),
              indent("</apply>\n", indentLevel+1), 
              indent("</apply>\n", indentLevel), sep="")
      else
        paste(indent("<apply>\n", indentLevel),
              indent({if (tree$operator=="|" || tree$operator=="any")
              {"<or/>\n"} else{"<and/>\n"}}, indentLevel+1),
              paste(sapply(tree$operands, MathMLFromParseTree, indentLevel+1), 
                    collapse=""),
              indent("</apply>\n", indentLevel), sep="")
    },
    atom =  
    {
      if ((tree$name == "0" && !tree$negated) || (tree$name == "1" && tree$negated))
        indent("<cn type=\"integer\">0</cn>\n", indentLevel)
      else
      if ((tree$name == "1" && !tree$negated) || (tree$name == "0" && tree$negated))
        indent("<cn type=\"integer\">1</cn>\n", indentLevel)
      else      
        paste(indent("<apply>\n",indentLevel),
              indent("<eq/>\n", indentLevel+1),
              indent(paste("<ci>",tree$name,"</ci>\n",sep=""),indentLevel+1),
              indent(paste("<cn type=\"integer\">", {if (tree$negated) 0 else 1},
                     "</cn>\n", sep=""), indentLevel+1),
              indent("</apply>\n",indentLevel), sep="")
    })
  return(res)
}
