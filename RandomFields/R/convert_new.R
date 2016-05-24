
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Alexander Malinowski, Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



### Diese Datei wandelt RMmodel in eine Liste um

### Kleinigkeiten geaendert: Jan 2015, M. Schlather


reps <- function(n, sign=",") paste(rep(sign, n), collapse="")




## @FUNCTION-STARP***********************************************************************************
# @NAME		parseModel
# @PARAM		$model - list, formula
# @RETURN		list
# @REQUIRE	$model is a linear mixed model in list or formula syntax
# @ENSURE		$listModel is a linear mixed model in list syntax
# @SEE		RMmodel, RFsimulate, devel-doc
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
parseModel <- function(model, ..., x=NULL) {

#Print(model)
  
  ## check whether $model is already in list syntax
  if (is.list(model)) return(model)
	
  ## check whether $model is already in list syntax
  if (isModel(model)) return(buildCovList(model, x=x))
	
  ## check whether $model has correct formula syntax
  if (!isFormulaModel(model)) stop(syntaxError)

  
  ## extract tokens
  summands <- extractSummands(model)
  #  Print(summands)
  
  listModel <- list()	
  if (length(summands) == 1) {
    listModel <- buildFactorList(summands[[1]], ..., x=x)#, last=TRUE)
  } else {
    listModel <- list(ZF_SYMBOLS_PLUS)		
    ##last <- getLastCovIndex(summands)
    for (i in 1:length(summands)) {
       listModel <- c(listModel, list(buildFactorList(summands[[i]], ..., x=x)))
    }
  }

#  Print(listModel); kkk
  
  return(listModel)
}



# @FUNCTION-STARP***********************************************************************************
# @NAME		isListModel
# @PARAM		$model - any r-object
# @RETURN		TRUE, FALSE
# @REQUIRE	none
# @ENSURE		it is confirmed that either $model has a correct model syntax or not
# @SEE		devel-doc
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
isModel <- function(model)
{
	return(is(model, ZF_MODEL))
}




# @FUNCTION-STARP***********************************************************************************
# @NAME		isFormulaModel
# @PARAM		$model - any r-object
# @RETURN		TRUE, FALSE
# @REQUIRE	none
# @ENSURE		it is confirmed that either $model has a correct model syntax or not
# @SEE		devel-doc
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
isFormulaModel <- function(model)
{
        if (missing(model) || is.null(model)) return(FALSE)
        if (!is(model, "formula")) return(FALSE)
		
	# ensure the @ operator is just bivariate
	if (regexpr("[[:alnum:]_]+@[[:alnum:]_]*@",
                    tail(as.character(model), 1)) != -1)
		return(FALSE)

	return(TRUE)
}


# @FUNCTION-STARP***********************************************************************************
# @NAME		buildCovList
# @PARAM		$model - RMmodel
# @RETURN		list
# @REQUIRE	none
# @ENSURE		isListModel(output) == TRUE
# @SEE		RMmodel, devel-doc
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
#                        modified 2015 by Martin Schlather
## @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
buildCovList <- function(model, x=NULL) {
  if (is.atomic(model) || is.list(model) ||
      is.language(model) || is.environment(model))
    return(model) ## for recursive calling
  
  if (!is(model, ZF_MODEL))
    stop('model must be of class ZF_MODEL') 

  if (model@name==ZF_COORD) model@name <- ZF_MIXED[1]
  
  li <- c(list(model@name),
          lapply(model@par.model[model@par.model != ZF_DEFAULT_STRING],
                 FUN=buildCovList, x=x),
          lapply(model@submodels,
                 FUN=buildCovList, x=x)
          )
  if (li[[1]] == ZF_PLUS[1]) li[[1]] <- ZF_SYMBOLS_PLUS
  if (li[[1]] == ZF_MULT[1]) li[[1]] <- ZF_SYMBOLS_MULT
  if (li[[1]] == ZF_COVARIATE && length(x) > 0 &&
      all(names(li) != "x")) {
    li$x <- x
  }
  
  
  ##  par.general.is.default <-
  ##    unlist(lapply(model@par.general, FUN=function(x) x==ZF_DEFAULT_STRING))
  if (length(model@par.general)>0 &
      !all(model@par.general==ZF_DEFAULT_STRING)) {
    li <- c(DOLLAR[1],
            lapply(model@par.general[model@par.general != ZF_DEFAULT_STRING],
                   FUN=buildCovList, x=x),
            list(li))
    if (length(pos <- which(names(li)=="Aniso")) > 0)
      ## in c-level, parameter is called 'A'
      names(li)[pos] <- "A"
  }
  
  return(li)           
}


# @FUNCTION-STARP***********************************************************************************
# @NAME		syntaxError
# @PARAM		none
# @RETURN		string
# @REQUIRE	none
# @ENSURE		none
# @SEE		none
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
syntaxError <- "Malformed model expression -- maybe you have used a wrong or obsolete definition, or just used an incorrect option name. See ?RMmodel for the model definition. Check manual for further information (RMmodel, RFsimulate)"




# @FUNCTION-STARP***********************************************************************************
# @NAME		extractSummands
# @PARAM		$model - formula
# @RETURN		list[string]
# @REQUIRE	isFormulaModel(model) == TRUE
# @ENSURE		the output is a complete list of all summands
# @SEE		RMModel, devel-doc
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
extractSummands <- function(model) {
  tmpList <- list()

  ## ignore rest of the formula
  rightSide <- tail(as.character(model), 1)
  chars <- strsplit(rightSide, "")[[1]]
	
  ## toggles parenthesis, eg. whether we confront a toplevel plus or not
  parToggle <- 0
  
  token <- ""
  for (char in chars) {
    if (char == ZF_SYMBOLS_PLUS && parToggle == 0) {
      tmpList <- c(tmpList, token)
      token <- ""			
    } else {
      if (char != " ") token <- paste(token, char, sep="")		
      if (char == ZF_SYMBOLS_L_PAR) parToggle <- parToggle+ 1
      if (char == ZF_SYMBOLS_R_PAR) parToggle <- parToggle- 1
    }
  }

  tmpList <- c(tmpList, token)
  tokenList <- list()
  for (token in tmpList) {
    tokenList <- c(tokenList, removeParenthesis(token))
  }
        
  if (length(tokenList)==1) 
    return(tokenList)  ## sonst steht unten paste(NULL), was "" gibt
        
  iscov <- unlist(lapply(tokenList, FUN=isGenuineCovModel))

#  Print(iscov)
  
  tokenList <- c(tokenList[!iscov],
                 list(paste(unlist(tokenList[iscov]),
                            collapse=ZF_SYMBOLS_PLUS)))

  return(tokenList)
}


# @FUNCTION-STARP***********************************************************************************
# @NAME		removeParenthesis
# @PARAM		$string - string
# @RETURN		string
# @REQUIRE	none
# @ENSURE		The returned string is one not enclosed by parenthesis
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
removeParenthesis <- function(string)
{
	# split the string
	chars <- strsplit(string, "")[[1]]

	while (head(chars, 1) == ZF_SYMBOLS_L_PAR &&
               tail(chars, 1) == ZF_SYMBOLS_R_PAR)
	{
		chars <- chars[-1]
		chars <- chars[-length(chars)]
	}

	# rejoin the string
	string <- ""
	for (char in chars)
	{
		string <- paste(string, char, sep="")
	}
	#print(string)
	return (string)
}

# @FUNCTION-STARP***********************************************************************************
# @NAME		buildFactorList
# @PARAM		$summand - string
# @RETURN		list
# @REQUIRE	$summand is one of the form: string"@"string
#			$last - an indicator telling the function if the covariance function is the last of
#			its kind, for the rule exists to encapsulate all covariance functions in RMmixed
#			except the last
# @ENSURE		the output is a correct list
# @SEE		RMModel, devel-doc
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
## @DATE		29.08.2011
## changed by Martin Schlather 2015
# @FUNCTION-END*************************************************************************************
buildFactorList <- function(summand, ..., x=x) { #, last)

  factorList <- strsplit(summand, ZF_SYMBOLS_AT)[[1]]  
  factorsNr <- length(factorList)

  ## remove parenthesis
  factorA <- removeParenthesis(factorList[[1]])
  if (factorsNr == 2) {
    factorB <- removeParenthesis(factorList[[2]])
    if (isFormalCovModel(factorA))
      stop(paste(factorA, "must NOT be a covariance model"))
  }

  ##  Print(summand, factorsNr, isFormalCovModel(factorA), factorA)

  ## do we have a mixed model
  X <- catch(factorA, ...)    
  if (!(factorsNr == 1 && isFormalCovModel(factorA))) {# && last))
    if (factorsNr == 1) {
      #str(X)
      #     Print(X, is.factor(X))
      if (is.factor(X)) {        
        lev <- levels(X)
        m <- matrix(nrow=length(X), ncol=length(lev) - 1)
        for (i in 2:length(lev)) m[, i - 1] <- as.numeric(X == lev[i])
        tmpList <- list(ZF_COVARIATE)         
        tmpList[[COVARIATE_C_NAME]] <- m
        tmpList[[COVARIATE_X_NAME]] <- x
        tmpList[[COVARIATE_ADDNA_NAME]] <- TRUE
      } else if (is.vector(X)) {
        if (length(X) == 1) {
          tmpList <- list(ZF_SYMBOLS_CONST)
          tmpList[[CONST_A_NAME]] <-
            if (is.finite(X) && X == 1) NA else as.numeric(X)
        } else {
          tmpList <- list(ZF_COVARIATE)         
          tmpList[[COVARIATE_C_NAME]] <- as.numeric(X)
          tmpList[[COVARIATE_X_NAME]] <- x
          tmpList[[COVARIATE_ADDNA_NAME]] <- TRUE
        }
      } else {
        tmpList <- list(ZF_MIXED[1])
        tmpList[[MIXED_X_NAME]] <- X
        tmpList[[MIXED_BETA_NAME]] <- NA
      }
    } else {
      tmpList <- list(ZF_MIXED[1])
      tmpList[[MIXED_X_NAME]] <- X
      if (isGenuineCovModel(factorB)) {
        model <- catch(factorB, ...)
        if (model@name==ZF_COORD) {
          tmpList[["coord"]] <- model@par.model$coord
          tmpList[["dist"]] <- model@par.model$dist
          tmpList[["cov"]] <-
            (buildCovList(model@submodels[[1]], x=x))
        } else {
          tmpList[["cov"]] <-
            (buildCovList(model, x=x))
        }
      } else {
        tmpList[["b"]] <- extractFixed(factorB, ...)
      }
    }
   # Print(tmpList)
    return(tmpList)
  } else {
    tmpList <- buildCovList(X, x=x)
    return(tmpList)
  }
}

# @FUNCTION-STARP***********************************************************************************
# @NAME		isCovModel: isFormalCovModel, isGenuineCovModel
# @PARAM		$name - string
# @RETURN		TRUE, FALSE
# @REQUIRE	none
# @ENSURE		The $name argument is a function that returns an RMmodel object
# @SEE		RFModel, devel-doc
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		29.08.2011
# @FUNCTION-END*************************************************************************************
isFormalCovModel <- function(name)
{
  ## Martin: habe hier RFtrend hinzugefuegt!
  if (is(try(tmp <- eval(parse(text=name)), silent=TRUE),ZF_MODEL))
    return(TRUE)
          
  ## has signature, eg. funName"("funPar")"
  if (regexpr("^[[:alnum:]_]+\\([[:print:]]*\\)$", name) != 1)
    return(FALSE)
  
  fun <- strsplit(name, "\\(")[[1]][1]
  
  return(exists(fun) && is(get(fun), ZF_MODEL_FACTORY))
}

isGenuineCovModel <- function(name)
{
  ## Martin: habe hier RFtrend hinzugefuegt!
  if (substr(name, 1, length(ZF_TRENDFCT)) == ZF_TRENDFCT) return(FALSE)
  if (is(try(tmp <- eval(parse(text=name)), silent=TRUE), ZF_MODEL))
    return(TRUE)
          
  ## has signature, eg. funName"("funPar")"
  if (regexpr("^[[:alnum:]_]+\\([[:print:]]*\\)$", name) != 1)
    return(FALSE)
  
  fun <- strsplit(name, "\\(")[[1]][1]
  
  return(exists(fun) && is(get(fun), ZF_MODEL_FACTORY))
}

# @FUNCTION-STARP***********************************************************************************
# @NAME		extractFixed
# @PARAM		$factor - string
# @RETURN		vector
# @REQUIRE	$factor has format RMfixed(<b=>r-vector)
# @ENSURE		none
# @SEE		RMmodel
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
extractFixed <- function(factor, ...)
{
  
	# has signature, eg. funName"("funPar")"
	if (regexpr(paste("^",ZF_FIXED,"\\([[:print:]]*\\)$",sep=""),
                    factor) != 1) {
          stop(paste("Second factor is not a cov model AND does not start with",
                     ZF_FIXED, "\n", syntaxError))
        }

	# extract the argument of RMfixed

  first_par <- regexpr("{1}\\(", factor)
  last_par <- regexpr("\\)$", factor)
  if (first_par == -1 || last_par == -1)
    arg <- factor
  else
    arg <- substr(factor, first_par+1, last_par-1)
  
  #arg <- strsplit(factor,"\\)")[[1]]
  #arg <- strsplit(arg, "\\(")[[1]][2]

  beta <- catch(arg, ...)
  
  return(beta)
}


# @FUNCTION-STARP***********************************************************************************
# @NAME		catch
# @PARAM		$expr - expression
#			$handler - the error handler function
# @RETURN		r-object
# @REQUIRE	none
# @ENSURE		stops if the expression fails
# @SEE		
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
catch <- function(expr, handler=function(res){stop(res)}, ...)
{
  #   Print("catch", expr)
  tmpENV <- new.env(parent=.GlobalEnv)
  dots <- list(...)
  assign("dots", dots, envir=tmpENV)
  if (length(dots)>0) {
    text <- paste(names(dots), "<- dots[[", 1:length(dots), "]]", collapse=";")
    eval(envir=tmpENV, parse(text=text))
  }
  res <-try(eval(envir=tmpENV, parse(text=expr)), silent=TRUE)
  
  if (class(res) == "try-error")
    handler(res)
  
  return(res)
}


# @FUNCTION-STARP***********************************************************************************
# @NAME		prepareData
# @REQUIRE	$model is of type formula
#			$simObj is an RFsp object
# @ENSURE		
# @SEE		
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
selectDataAccordingFormula <- function(simObj, model) {
  varNames <- extractVarNames(model)
  #Print(varNames)

  ## are there any variablenames given
  if (is.null(varNames)) {
    ##warning("'model' is not given as a formula or model formula contains an empty left side --> all colums of the data matrix in 'data' are selected")
    coord <- RFoptions()$coord
    if (is.na(coord$varnames[1])) return (simObj)
    varnames <- coord$varnames
  }

  cleanNames <- (if (simObj@.RFparams$n == 1) colnames(simObj@data) else 
                 sapply(colnames(simObj@data), FUN=cleanse))
    
    mymatch <- match(cleanNames, varNames)#, dup=TRUE)
    if (!all(varNames %in% cleanNames))
      stop("response variable names could not be found in colnames of data object")

  #Print(simObj, mymatch)
  
    simObj <- simObj[!is.na(mymatch)]
    
    simObj@.RFparams$vdim <- length(varNames)
   
    return(simObj)
}
  
selectAccordingFormula <- function(data, model)
{
  varNames <- extractVarNames(model)
  #Print(varNames)

  if (is.null(varNames)) {
    coord <- RFoptions()$coord
    if (is.na(coord$varnames[1])) return(NULL)
    varNames <- coord$varnames
  }
  
  ## are there any variablenames given
  if (is.matrix(data) || is.data.frame(data)) {
    cleanNames <- colnames(data)
    if (!all(varNames %in% cleanNames))
      stop("response variable names could not be found in colnames of data object")
    mymatch <- !is.na(match(cleanNames, varNames))#, dup=TRUE)
    return(mymatch)
  }
}



# @FUNCTION-STARP***********************************************************************************
# @NAME		extractVarNames
# @REQUIRE	$model is a formula or a RMmodel
# @ENSURE		the return value is a vector with the names of the response variables of the formula 
#			or NULL if no left side given
# @SEE		
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011; 2014 Martin Schlather modified
# @FUNCTION-END*************************************************************************************
extractVarNames <- function(model)
{
	if (missing(model) ||
            length(model) <= 2 || !is(model, "formula")) return (NULL)
			
	tmp <- as.character(model)[2]
	tmp <- strsplit(tmp, "c\\(")[[1]]
	tmp <- paste(tmp, sep="", collapse="")
	tmp <- strsplit(tmp, "\\)")[[1]]
	tmp <- paste(tmp, sep="", collapse="")
	varNames <- strsplit(tmp, ", ")[[1]]

	## ignore numeric formated varNames
	i<- 1
	while (i <= length(varNames))
	{
		if (regexpr("^[[:digit:]]", varNames[i]) > 0)
			varNames<- varNames[-i]
		else
			i<- i+ 1
	}
	if (length(varNames) == 0)
		return (NULL)
	
	return (varNames)
}

# @FUNCTION-STARP***********************************************************************************
# @NAME		cleanse
# @REQUIRE	$x a character string
# @ENSURE		if $x ends with '.n[:digit:]+' this part is cut off
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @FUNCTION-END*************************************************************************************
cleanse <- function(x)
{
	return (strsplit(x, "\\.n[[:digit:]]+$")[[1]][1])
}

# @FUNCTION-STARP***********************************************************************************
# @NAME		cutoffArray
# @REQUIRE	$data is a numeric array
#			$grid is a sequence matrix and each $data point has a corresponding $grid point
#			in particular the expanded $grid has the same size as $data
#			$len is a numeric vector sorted in ascending order
#			all three parameters start with 0
# @ENSURE		$data is returned and the original array with each point lying on an coordinate axis
#			having a maximum distance less or equal the last entry of $len
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		14.09.2011
# @FUNCTION-END*************************************************************************************
cutoffArray <- function(data, grid, len)

  ## Was macht diese Funktion? (Frage von Martin)
  ## arrays abschneiden, sodass es in keiner orthogonalen raumrichtung punkte gibt
  ## die einen abstand groesser len haben
{
	len <- tail(len, 1)
	
	## cut off in each direction
	for (dir in c(1:dim(grid)[2]))
	{
		## retrieve the number of points within bin length
		nrPts <- ceiling(len/grid[2,dir])+ 1		
				
          text <- paste("data[", reps(dir-1), "1:", nrPts,
                        reps(dim(grid)[2]-dir), ",drop=FALSE]")                
                              
		data <- eval(parse(text=text))
	}
	
	return (data)
}




