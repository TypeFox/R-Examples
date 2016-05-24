################################################################################################
##    Copyright (C) 2015,  Constantinos Tsirogiannis and Brody Sandel.  
##
##    Email: analekta@gmail.com and brody.sandel@bios.au.dk
##
##    This file is part of PhyloMeasures.
##
##    PhyloMeasures is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    PhyloMeasures is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with PhyloMeasures.  If not, see <http://www.gnu.org/licenses/>
################################################################################################

require(ape)

#################################################
#################################################
########## Matrix Query Functions ###############
#################################################
#################################################

pd.query = function(tree, matrix, is.standardised = FALSE)
{
  matrix.data = unlist(t(matrix))

  if(any(is.na(matrix.data)) || any(matrix.data < 0))
    stop("The input matrix contains negative or undefined values.\n\n")


  if(any(as.logical(matrix.data) != matrix.data))
    cat("\n Warning: the input matrix contains values which are not 0 or 1.\n")

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species = colnames(matrix)
  matrix.rows = nrow(matrix)
  matrix.columns = ncol(matrix)

  output = vector(mode="numeric",matrix.rows)

  error.message = as.character(paste0(rep(" ",length = 150), collapse=""))
  message.size = 0

  res <-.C("pd_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
           as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species), 
           as.integer(matrix.rows), as.integer(matrix.columns), as.integer(matrix.data), 
           as.logical(is.standardised), as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[13]] > 0)
    stop(substr(res[[12]],1,res[[13]]))

  return (res[[11]])

} # pd.query = function(...)


mpd.query = function(tree, matrix, is.standardised = FALSE)
{
  matrix.data = unlist(t(matrix))

  if(any(is.na(matrix.data)) || any(matrix.data < 0))
    stop("The input matrix contains negative or undefined values.\n\n")

  if(any(as.logical(matrix.data) != matrix.data))
    cat("\n Warning: the input matrix contains values which are not 0 or 1.\n")

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species = colnames(matrix)
  matrix.rows = nrow(matrix)
  matrix.columns = ncol(matrix)

  output = vector(mode="numeric",matrix.rows)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  res <-.C("mpd_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
           as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species), 
           as.integer(matrix.rows), as.integer(matrix.columns), as.integer(matrix.data), 
           as.logical(is.standardised), as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[13]] > 0)
    stop(substr(res[[12]],1,res[[13]]))


  return (res[[11]])

} # mpd.query = function(...)


mntd.query = function(tree, matrix, is.standardised = FALSE)
{
  matrix.data = unlist(t(matrix))

  if(any(is.na(matrix.data)) || any(matrix.data < 0))
    stop("The input matrix contains negative or undefined values.\n\n")

  if(any(as.logical(matrix.data) != matrix.data))
    cat("\n Warning: the input matrix contains values which are not 0 or 1.\n")

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species = colnames(matrix)
  matrix.rows = nrow(matrix)
  matrix.columns = ncol(matrix)

  output = vector(mode="numeric",matrix.rows)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  res <-.C("mntd_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
           as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species), 
           as.integer(matrix.rows), as.integer(matrix.columns), as.integer(matrix.data), 
           as.logical(is.standardised), as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[13]] > 0)
    stop(substr(res[[12]],1,res[[13]]))

  return (res[[11]])

} # mntd.query = function(...)

cac.query = function(tree, matrix, chi, is.standardised = FALSE)
{
  matrix.data = unlist(t(matrix))

  if(any(is.na(matrix.data)) || any(matrix.data < 0))
    stop("The input matrix contains negative or undefined values.\n\n")


  if(any(as.logical(matrix.data) != matrix.data))
    cat("\n Warning: the input matrix contains values which are not 0 or 1.\n")

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species = colnames(matrix)
  matrix.rows = nrow(matrix)
  matrix.columns = ncol(matrix)

  output = vector(mode="numeric",matrix.rows)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  res <-.C("cac_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
           as.numeric(tree$edge.length), as.character(tree$tip.label), as.numeric(chi), as.character(matrix.species), 
           as.integer(matrix.rows), as.integer(matrix.columns), as.integer(matrix.data), 
           as.logical(is.standardised), as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[14]] > 0)
    stop(substr(res[[13]],1,res[[14]]))

  return (res[[12]])

} # cac.query = function(...)


cbl.query = function(tree, matrix.a, matrix.b = NULL, query.matrix=NULL, is.standardised = FALSE)
{
  matrix.rows.a = nrow(matrix.a)
  matrix.data.a = unlist(t(matrix.a))

  matrix.rows.b = 0
  matrix.data.b = NULL

  query.matrix.data = NULL

  if(any(is.na(matrix.data.a)) || any(matrix.data.a < 0))
    stop("One of the input matrices contains negative or undefined values.\n\n")


  if(any(as.logical(matrix.data.a) != matrix.data.a))
    cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")

  if(is.null(matrix.b) == FALSE)
  {
    matrix.rows.b = nrow(matrix.b)
    matrix.data.b = unlist(t(matrix.b))

    if(any(is.na(matrix.data.b)) || any(matrix.data.b < 0))
      stop("One of the input matrices contains negative or undefined values.\n\n")

    if(any(as.logical(matrix.data.b) != matrix.data.b))
      cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")
  }

  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.data = unlist(t(query.matrix))

    if(any(is.na(query.matrix.data)) || any(query.matrix.data < 1))
      stop("The matrix with the query pairs contains undefined, or out of range values.\n\n")

    col.a = query.matrix[,1]
    col.b = query.matrix[,2]

    if(any(col.a > matrix.rows.a))
      stop("The matrix with the query pairs contains out of range values.\n\n")

    temp.size.b = matrix.rows.a

    if(is.null(matrix.b) == FALSE)
      temp.size.b = matrix.rows.b

    if(any(col.b > temp.size.b))
      stop("The matrix with the query pairs contains out of range values. \n\n")

    
  } # if(is.null(query.matrix) == FALSE)


  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species.a = colnames(matrix.a)
  matrix.columns.a = ncol(matrix.a)

  matrix.species.b = NULL
  matrix.columns.b = 0

  query.matrix.rows = 0

  if(is.null(matrix.b) == FALSE)
  {
    matrix.species.b = colnames(matrix.b)
    matrix.columns.b = ncol(matrix.b)
  } 

  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.rows = nrow(query.matrix)
  }

  output.size = 0

  if(!is.null(query.matrix))
  {
    output.size = query.matrix.rows

  } else {

    if(is.null(matrix.b))   
    {
      output.size = matrix.rows.a * matrix.rows.a;
    } else {
      output.size = matrix.rows.a * matrix.rows.b;
    }

  }

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("cbl_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows), as.integer(query.matrix.data),
             as.logical(is.standardised), as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

    if(res[[19]] > 0)
      stop(substr(res[[18]],1,res[[19]]))

    return (res[[17]])  
  }

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("cbl_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows), as.integer(NULL), 
             as.logical(is.standardised), as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

    if(res[[19]] > 0)
      stop(substr(res[[18]],1,res[[19]]))

    results = t(matrix(res[[17]], nrow=matrix.rows.b, ncol=matrix.rows.a))   
    return (results)  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("cbl_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(query.matrix.data), as.logical(is.standardised), 
             as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

    if(res[[19]] > 0)
      stop(substr(res[[18]],1,res[[19]]))

    return (res[[17]])  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("cbl_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(NULL), as.logical(is.standardised), 
             as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

    if(res[[19]] > 0)
      stop(substr(res[[18]],1,res[[19]]))

    results = t(matrix(res[[17]], nrow=matrix.rows.a, ncol=matrix.rows.a))   
    return (results)  
  }

} # cbl.query = function(...)


cd.query = function(tree, matrix.a, matrix.b = NULL, query.matrix=NULL, is.standardised = FALSE)
{
  matrix.rows.a = nrow(matrix.a)
  matrix.data.a = unlist(t(matrix.a))

  matrix.rows.b = 0
  matrix.data.b = NULL

  query.matrix.data = NULL

  if(any(is.na(matrix.data.a)) || any(matrix.data.a < 0))
    stop("One of the input matrices contains negative or undefined values. \n\n")


  if(any(as.logical(matrix.data.a) != matrix.data.a))
    cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")

  if(is.null(matrix.b) == FALSE)
  {
    matrix.rows.b = nrow(matrix.b)
    matrix.data.b = unlist(t(matrix.b))

    if(any(is.na(matrix.data.b)) || any(matrix.data.b < 0))
      stop("One of the input matrices contains negative or undefined values. \n\n")


    if(any(as.logical(matrix.data.b) != matrix.data.b))
      cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")
  }


  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.data = unlist(t(query.matrix))

    if(any(is.na(query.matrix.data)) || any(query.matrix.data < 1))
      stop("The matrix with the query pairs contains undefined, or out of range values. \n\n")


    col.a = query.matrix[,1]
    col.b = query.matrix[,2]

    if(any(col.a > matrix.rows.a))
      stop("The matrix with the query pairs contains out of range values. \n\n")

    temp.size.b = matrix.rows.a

    if(is.null(matrix.b) == FALSE)
      temp.size.b = matrix.rows.b

    if(any(col.b > temp.size.b))
      stop("The matrix with the query pairs contains out of range values. \n\n")

    
  } # if(is.null(query.matrix) == FALSE)


  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species.a = colnames(matrix.a)
  matrix.rows.a = nrow(matrix.a)
  matrix.columns.a = ncol(matrix.a)

  matrix.species.b = NULL
  matrix.columns.b = 0
  query.matrix.rows = 0

  if(is.null(matrix.b) == FALSE)
  {
    matrix.species.b = colnames(matrix.b)
    matrix.columns.b = ncol(matrix.b)
  } 

  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.rows = nrow(query.matrix)
  }

  output.size = 0

  if(!is.null(query.matrix))
  {
    output.size = query.matrix.rows

  } else {

    if(is.null(matrix.b))   
    {
      output.size = matrix.rows.a * matrix.rows.a;
    } else {
      output.size = matrix.rows.a * matrix.rows.b;
    }

  }

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("cd_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(query.matrix.data), as.logical(is.standardised),
             as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

    if(res[[19]] > 0)
      stop(substr(res[[18]],1,res[[19]]))

    return (res[[17]])  
  }

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("cd_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(NULL), as.logical(is.standardised), as.numeric(output), 
             as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

    if(res[[19]] > 0)
      stop(substr(res[[18]],1,res[[19]]))

    results = t(matrix(res[[17]], nrow=matrix.rows.b, ncol=matrix.rows.a))   
    return (results)  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("cd_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(NULL), as.integer(query.matrix.rows),
             as.integer(query.matrix.data), as.logical(is.standardised), 
             as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

    if(res[[19]] > 0)
      stop(substr(res[[18]],1,res[[19]]))

    return (res[[17]])  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("cd_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(NULL), as.logical(is.standardised), 
             as.numeric(output), as.character(error.message), as.integer(message.size),
           PACKAGE = "PhyloMeasures")

    if(res[[19]] > 0)
      stop(substr(res[[18]],1,res[[19]]))

    results = t(matrix(res[[17]], nrow=matrix.rows.a, ncol=matrix.rows.a))   
    return (results)  
  }


} # cd.query = function(...)


cdnt.query = function(tree, matrix.a, matrix.b = NULL, query.matrix=NULL)
{
  matrix.rows.a = nrow(matrix.a)
  matrix.data.a = unlist(t(matrix.a))

  matrix.rows.b = 0
  matrix.data.b = NULL

  query.matrix.data = NULL

  if(any(is.na(matrix.data.a)) || any(matrix.data.a < 0))
    stop("One of the input matrices contains negative or undefined values. \n\n")

  if(any(as.logical(matrix.data.a) != matrix.data.a))
    cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")

  if(is.null(matrix.b) == FALSE)
  {
    matrix.rows.b = nrow(matrix.b)
    matrix.data.b = unlist(t(matrix.b))

    if(any(is.na(matrix.data.b)) || any(matrix.data.b < 0))
      stop("One of the input matrices contains negative or undefined values. \n\n")


    if(any(as.logical(matrix.data.b) != matrix.data.b))
      cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")
  }


  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.data = unlist(t(query.matrix))

    if(any(is.na(query.matrix.data)) || any(query.matrix.data < 1))
      stop("The matrix with the query pairs contains undefined, or out of range values. \n\n")


    col.a = query.matrix[,1]
    col.b = query.matrix[,2]

    if(any(col.a > matrix.rows.a))
      stop("The matrix with the query pairs contains out of range values. \n\n")


    temp.size.b = matrix.rows.a

    if(is.null(matrix.b) == FALSE)
      temp.size.b = matrix.rows.b

    if(any(col.b > temp.size.b))
      stop("The matrix with the query pairs contains out of range values. \n\n")

    
  } # if(is.null(query.matrix) == FALSE)

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species.a = colnames(matrix.a)
  matrix.rows.a = nrow(matrix.a)
  matrix.columns.a = ncol(matrix.a)

  matrix.species.b = NULL
  matrix.columns.b = 0
  query.matrix.rows = 0

  if(is.null(matrix.b) == FALSE)
  {
    matrix.species.b = colnames(matrix.b)
    matrix.columns.b = ncol(matrix.b)
  } 

  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.rows = nrow(query.matrix)
  }

  output.size = 0

  if(!is.null(query.matrix))
  {
    output.size = query.matrix.rows

  } else {

    if(is.null(matrix.b))   
    {
      output.size = matrix.rows.a * matrix.rows.a;
    } else {
      output.size = matrix.rows.a * matrix.rows.b;
    }

  }

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("cdnt_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(query.matrix.data), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    return (res[[16]])  
  }

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("cdnt_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(NULL), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    results = t(matrix(res[[16]], nrow=matrix.rows.b, ncol=matrix.rows.a))   
    return (results)  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("cdnt_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(query.matrix.data), 
             as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    return (res[[16]])  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("cdnt_query",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(NULL), as.numeric(output), 
             as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    results = t(matrix(res[[16]], nrow=matrix.rows.a, ncol=matrix.rows.a))   
    return (results)  
  }

} # cdnt.query = function(...)


cdnt.averaged.query = function(tree, matrix.a, matrix.b = NULL, query.matrix=NULL)
{
  matrix.rows.a = nrow(matrix.a)
  matrix.data.a = unlist(t(matrix.a))

  matrix.rows.b = 0
  matrix.data.b = NULL

  query.matrix.data = NULL

  if( any(is.na(matrix.data.a)) || any(matrix.data.a < 0) )
    stop("One of the input matrices contains negative or undefined values. \n\n")

  if(any(as.logical(matrix.data.a) != matrix.data.a))
    cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")

  if(is.null(matrix.b) == FALSE)
  {
    matrix.rows.b = nrow(matrix.b)
    matrix.data.b = unlist(t(matrix.b))

    if(any(is.na(matrix.data.b)) || any(matrix.data.b < 0))
      stop("One of the input matrices contains negative or undefined values. \n\n")

    if(any(as.logical(matrix.data.b) != matrix.data.b))
      cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")
  }


  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.data = unlist(t(query.matrix))

    if(any(is.na(query.matrix.data)) || any(query.matrix.data < 1))
      stop("The matrix with the query pairs contains undefined, or out of range values. \n\n")


    col.a = query.matrix[,1]
    col.b = query.matrix[,2]

    if(any(col.a > matrix.rows.a))
      stop("The matrix with the query pairs contains out of range values. \n\n")


    temp.size.b = matrix.rows.a

    if(is.null(matrix.b) == FALSE)
      temp.size.b = matrix.rows.b

    if(any(col.b > temp.size.b))
      stop("The matrix with the query pairs contains out of range values. \n\n")

    
  } # if(is.null(query.matrix) == FALSE)

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species.a = colnames(matrix.a)
  matrix.rows.a = nrow(matrix.a)
  matrix.columns.a = ncol(matrix.a)

  matrix.species.b = NULL
  matrix.columns.b = 0
  query.matrix.rows = 0

  if(is.null(matrix.b) == FALSE)
  {
    matrix.species.b = colnames(matrix.b)
    matrix.columns.b = ncol(matrix.b)
  } 

  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.rows = nrow(query.matrix)
  }

  output.size = 0

  if(!is.null(query.matrix))
  {
    output.size = query.matrix.rows

  } else {

    if(is.null(matrix.b))   
    {
      output.size = matrix.rows.a * matrix.rows.a;
    } else {
      output.size = matrix.rows.a * matrix.rows.b;
    }

  }

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("cdnt_averaged_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(query.matrix.data), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    return (res[[16]])  
  }

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("cdnt_averaged_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(NULL), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    results = t(matrix(res[[16]], nrow=matrix.rows.b, ncol=matrix.rows.a))
    return (results)  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("cdnt_averaged_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(query.matrix.data), 
             as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    return (res[[16]])  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("cdnt_averaged_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(NULL), as.numeric(output), 
             as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    results = t(matrix(res[[16]], nrow=matrix.rows.a, ncol=matrix.rows.a))   
    return (results)  
  }

} # cdnt.averaged.query = function(...)


cdnt.directed.query = function(tree, matrix.a, matrix.b = NULL, query.matrix=NULL)
{
  matrix.rows.a = nrow(matrix.a)
  matrix.data.a = unlist(t(matrix.a))

  matrix.rows.b = 0
  matrix.data.b = NULL

  query.matrix.data = NULL

  if(any(is.na(matrix.data.a)) || any(matrix.data.a < 0))
    stop("One of the input matrices contains negative or undefined values. \n\n")


  if(any(as.logical(matrix.data.a) != matrix.data.a))
    cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")

  if(is.null(matrix.b) == FALSE)
  {
    matrix.rows.b = nrow(matrix.b)
    matrix.data.b = unlist(t(matrix.b))

    if(any(is.na(matrix.data.b)) || any(matrix.data.b < 0))
      stop("One of the input matrices contains negative or undefined values. \n\n")


    if(any(as.logical(matrix.data.b) != matrix.data.b))
      cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")
  }


  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.data = unlist(t(query.matrix))

    if(any(is.na(query.matrix.data)) || any(query.matrix.data < 1))
      stop("The matrix with the query pairs contains undefined, or out of range values. \n\n")

    col.a = query.matrix[,1]
    col.b = query.matrix[,2]

    if(any(col.a > matrix.rows.a))
      stop("The matrix with the query pairs contains out of range values. \n\n")


    temp.size.b = matrix.rows.a

    if(is.null(matrix.b) == FALSE)
      temp.size.b = matrix.rows.b

    if(any(col.b > temp.size.b))
      stop("The matrix with the query pairs contains out of range values. \n\n")

    
  } # if(is.null(query.matrix) == FALSE)

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species.a = colnames(matrix.a)
  matrix.rows.a = nrow(matrix.a)
  matrix.columns.a = ncol(matrix.a)

  matrix.species.b = NULL
  matrix.columns.b = 0
  query.matrix.rows = 0

  if(is.null(matrix.b) == FALSE)
  {
    matrix.species.b = colnames(matrix.b)
    matrix.columns.b = ncol(matrix.b)
  } 

  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.rows = nrow(query.matrix)
  }

  output.size = 0

  if(!is.null(query.matrix))
  {
    output.size = 2*query.matrix.rows

  } else {

    if(is.null(matrix.b))   
    {
      output.size = matrix.rows.a * matrix.rows.a;
    } else {
      output.size = 2*matrix.rows.a * matrix.rows.b;
    }

  } # else of if(!is.null(query.matrix))

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0
  
  res = NULL

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == FALSE)
    res <-.C("cdnt_directed_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(query.matrix.data), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == TRUE)
    res <-.C("cdnt_directed_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), 
             as.integer(matrix.data.a), as.character(matrix.species.b), 
             as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(NULL), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == FALSE)
    res <-.C("cdnt_directed_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), 
             as.integer(matrix.data.a), as.character(NULL), 
             as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(NULL), as.integer(query.matrix.rows),
             as.integer(query.matrix.data), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == TRUE)
    res <-.C("cdnt_directed_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), 
             as.integer(matrix.data.a), as.character(NULL), 
             as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(NULL), as.integer(query.matrix.rows),
             as.integer(NULL), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

  if(res[[18]] > 0)
    stop(substr(res[[17]],1,res[[18]]))

    
  result.a.to.b = NULL
  result.b.to.a = NULL

  if(is.null(query.matrix) == TRUE)
  {
    if(is.null(matrix.b) == TRUE)
    {
      result.a.to.b = t(matrix(res[[16]][1:(matrix.rows.a*matrix.rows.a)], nrow=matrix.rows.a, ncol=matrix.rows.a))
      result = list("result.A.to.B" = result.a.to.b)
      return (result)
    } 
    else 
    {

      result.a.to.b = t(matrix(res[[16]][1:(matrix.rows.a*matrix.rows.b)], nrow=matrix.rows.b, ncol=matrix.rows.a))
      result.b.to.a = t(matrix(tail(res[[16]],(matrix.rows.a*matrix.rows.b)),
                                       nrow=matrix.rows.a, ncol=matrix.rows.b))
    }
    
  } 
  else 
  {
   result.a.to.b= head(res[[16]],query.matrix.rows)
   result.b.to.a= tail(res[[16]],query.matrix.rows)
  } # else of if(is.null(query.matrix) == TRUE) 

  result = list("result.A.to.B" = result.a.to.b, "result.B.to.A" = result.b.to.a)
  
  return (result)

} # cdnt.directed.query = function(...)


phylosor.query = function(tree, matrix.a, matrix.b = NULL, query.matrix=NULL)
{
  matrix.rows.a = nrow(matrix.a)
  matrix.data.a = unlist(t(matrix.a))

  matrix.rows.b = 0
  matrix.data.b = NULL

  query.matrix.data = NULL

  if( any(is.na(matrix.data.a)) || any(matrix.data.a < 0) )
    stop("One of the input matrices contains negative or undefined values. \n\n")

  if(any(as.logical(matrix.data.a) != matrix.data.a))
    cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")

  if(is.null(matrix.b) == FALSE)
  {
    matrix.rows.b = nrow(matrix.b)
    matrix.data.b = unlist(t(matrix.b))

    if(any(is.na(matrix.data.b)) || any(matrix.data.b < 0))
      stop("One of the input matrices contains negative or undefined values. \n\n")

    if(any(as.logical(matrix.data.b) != matrix.data.b))
      cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")
  }


  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.data = unlist(t(query.matrix))

    if(any(is.na(query.matrix.data)) || any(query.matrix.data < 1))
      stop("The matrix with the query pairs contains undefined, or out of range values. \n\n")


    col.a = query.matrix[,1]
    col.b = query.matrix[,2]

    if(any(col.a > matrix.rows.a))
      stop("The matrix with the query pairs contains out of range values. \n\n")


    temp.size.b = matrix.rows.a

    if(is.null(matrix.b) == FALSE)
      temp.size.b = matrix.rows.b

    if(any(col.b > temp.size.b))
      stop("The matrix with the query pairs contains out of range values. \n\n")

    
  } # if(is.null(query.matrix) == FALSE)

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species.a = colnames(matrix.a)
  matrix.rows.a = nrow(matrix.a)
  matrix.columns.a = ncol(matrix.a)

  matrix.species.b = NULL
  matrix.columns.b = 0
  query.matrix.rows = 0

  if(is.null(matrix.b) == FALSE)
  {
    matrix.species.b = colnames(matrix.b)
    matrix.columns.b = ncol(matrix.b)
  } 

  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.rows = nrow(query.matrix)
  }

  output.size = 0

  if(!is.null(query.matrix))
  {
    output.size = query.matrix.rows

  } else {

    if(is.null(matrix.b))   
    {
      output.size = matrix.rows.a * matrix.rows.a;
    } else {
      output.size = matrix.rows.a * matrix.rows.b;
    }

  }

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("phylosor_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(query.matrix.data), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    return (res[[16]])  
  }

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("phylosor_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(NULL), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    results = t(matrix(res[[16]], nrow=matrix.rows.b, ncol=matrix.rows.a))
    return (results)  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("phylosor_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(query.matrix.data), 
             as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    return (res[[16]])  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("phylosor_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(NULL), as.numeric(output), 
             as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    results = t(matrix(res[[16]], nrow=matrix.rows.a, ncol=matrix.rows.a))   
    return (results)  
  }

} # phylosor.query = function(...)


unifrac.query = function(tree, matrix.a, matrix.b = NULL, query.matrix=NULL)
{
  matrix.rows.a = nrow(matrix.a)
  matrix.data.a = unlist(t(matrix.a))

  matrix.rows.b = 0
  matrix.data.b = NULL

  query.matrix.data = NULL

  if( any(is.na(matrix.data.a)) || any(matrix.data.a < 0) )
    stop("One of the input matrices contains negative or undefined values. \n\n")

  if(any(as.logical(matrix.data.a) != matrix.data.a))
    cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")

  if(is.null(matrix.b) == FALSE)
  {
    matrix.rows.b = nrow(matrix.b)
    matrix.data.b = unlist(t(matrix.b))

    if(any(is.na(matrix.data.b)) || any(matrix.data.b < 0))
      stop("One of the input matrices contains negative or undefined values. \n\n")

    if(any(as.logical(matrix.data.b) != matrix.data.b))
      cat("\n Warning: one of the input matrices contains values which are not 0 or 1.\n")
  }


  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.data = unlist(t(query.matrix))

    if(any(is.na(query.matrix.data)) || any(query.matrix.data < 1))
      stop("The matrix with the query pairs contains undefined, or out of range values. \n\n")


    col.a = query.matrix[,1]
    col.b = query.matrix[,2]

    if(any(col.a > matrix.rows.a))
      stop("The matrix with the query pairs contains out of range values. \n\n")


    temp.size.b = matrix.rows.a

    if(is.null(matrix.b) == FALSE)
      temp.size.b = matrix.rows.b

    if(any(col.b > temp.size.b))
      stop("The matrix with the query pairs contains out of range values. \n\n")

    
  } # if(is.null(query.matrix) == FALSE)

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  number.of.leaves = length(tree$tip.label)

  matrix.species.a = colnames(matrix.a)
  matrix.rows.a = nrow(matrix.a)
  matrix.columns.a = ncol(matrix.a)

  matrix.species.b = NULL
  matrix.columns.b = 0
  query.matrix.rows = 0

  if(is.null(matrix.b) == FALSE)
  {
    matrix.species.b = colnames(matrix.b)
    matrix.columns.b = ncol(matrix.b)
  } 

  if(is.null(query.matrix) == FALSE)
  {
    query.matrix.rows = nrow(query.matrix)
  }

  output.size = 0

  if(!is.null(query.matrix))
  {
    output.size = query.matrix.rows

  } else {

    if(is.null(matrix.b))   
    {
      output.size = matrix.rows.a * matrix.rows.a;
    } else {
      output.size = matrix.rows.a * matrix.rows.b;
    }

  }

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("unifrac_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(query.matrix.data), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    return (res[[16]])  
  }

  if(is.null(matrix.b) == FALSE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("unifrac_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(matrix.species.b), as.integer(matrix.rows.b), as.integer(matrix.columns.b), 
             as.integer(matrix.data.b), as.integer(query.matrix.rows),
             as.integer(NULL), as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    results = t(matrix(res[[16]], nrow=matrix.rows.b, ncol=matrix.rows.a))
    return (results)  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == FALSE)
  {
    res <-.C("unifrac_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(query.matrix.data), 
             as.numeric(output), as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    return (res[[16]])  
  }

  if(is.null(matrix.b) == TRUE &&  is.null(query.matrix) == TRUE)
  {
    res <-.C("unifrac_query",as.integer(length(tree$edge.length)), 
             as.integer(number.of.leaves), as.integer(tree$edge), 
             as.numeric(tree$edge.length), as.character(tree$tip.label), as.character(matrix.species.a), 
             as.integer(matrix.rows.a), as.integer(matrix.columns.a), as.integer(matrix.data.a), 
             as.character(NULL), as.integer(matrix.rows.b), as.integer(matrix.columns.b), as.integer(NULL), 
             as.integer(query.matrix.rows), as.integer(NULL), as.numeric(output), 
             as.character(error.message), as.integer(message.size),
             PACKAGE = "PhyloMeasures")

    if(res[[18]] > 0)
      stop(substr(res[[17]],1,res[[18]]))

    results = t(matrix(res[[16]], nrow=matrix.rows.a, ncol=matrix.rows.a))   
    return (results)  
  }

} # unifrac.query = function(...)


#################################################
#################################################
############# Moments Functions #################
#################################################
#################################################

pd.moments = function(tree, sample.sizes, comp.expectation = TRUE, comp.deviation = TRUE)
{
  number.of.leaves = length(tree$tip.label)

  if(any(is.na(sample.sizes))||any(sample.sizes < 0)||any(sample.sizes > number.of.leaves))
    stop("The list with the sample sizes contains values which are undefined, or out of range. \n\n")


  if(any(as.integer(sample.sizes) != sample.sizes))
    stop("The list with the sample sizes contains non-integer values.\n\n")
 

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  vector.size = length(sample.sizes)

  output.size = 0 

  if(comp.expectation == TRUE)
    output.size = output.size + vector.size

  if(comp.deviation == TRUE)
    output.size = output.size + vector.size

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  res <-.C("pd_moments",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), 
           as.integer(tree$edge), as.numeric(tree$edge.length), as.character(tree$tip.label), 
           as.integer(vector.size), as.integer(sample.sizes), as.logical(comp.expectation), 
           as.logical(comp.deviation), as.numeric(output),as.character(error.message),as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[12]] > 0)
    stop(substr(res[[11]],1,res[[12]]))

  result = NULL

  if( comp.expectation == TRUE && comp.deviation == TRUE)
    result = matrix(res[[10]], nrow=output.size/2, ncol=2)
  else
    result = res[[10]]

  return (result)

} # pd.moments = function(...)


mpd.moments = function(tree, sample.sizes, comp.expectation = TRUE, comp.deviation = TRUE)
{
  number.of.leaves = length(tree$tip.label)

  if(any(is.na(sample.sizes))||any(sample.sizes < 0)||any(sample.sizes > number.of.leaves))
    stop("The list with the sample sizes contains values which are undefined, or out of range. \n\n")


  if(any(as.integer(sample.sizes) != sample.sizes))
    stop("The list with the sample sizes contains non-integer values.\n\n")


  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  vector.size = length(sample.sizes)

  output.size = 0 

  if(comp.expectation == TRUE)
    output.size = output.size + vector.size

  if(comp.deviation == TRUE)
    output.size = output.size + vector.size

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  res <-.C("mpd_moments",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), 
           as.integer(tree$edge), as.numeric(tree$edge.length), as.character(tree$tip.label), 
           as.integer(vector.size), as.integer(sample.sizes), as.logical(comp.expectation), 
           as.logical(comp.deviation), as.numeric(output),as.character(error.message),as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[12]] > 0)
    stop(substr(res[[11]],1,res[[12]]))

  result = NULL

  if( comp.expectation == TRUE && comp.deviation == TRUE)
    result = matrix(res[[10]], nrow=output.size/2, ncol=2)
  else
    result = res[[10]]

  return (result)

} # mpd.moments = function(...)


mntd.moments = function(tree, sample.sizes, comp.expectation = TRUE, comp.deviation = TRUE)
{
  number.of.leaves = length(tree$tip.label)

  if(any(is.na(sample.sizes))||any(sample.sizes < 0)||any(sample.sizes > number.of.leaves))
    stop("The list with the sample sizes contains values which are undefined, or out of range. \n\n")


  if(any(as.integer(sample.sizes) != sample.sizes))
    stop("The list with the sample sizes contains non-integer values.\n\n")


  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  vector.size = length(sample.sizes)

  output.size = 0 
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  if(comp.expectation == TRUE)
    output.size = output.size + vector.size

  if(comp.deviation == TRUE)
    output.size = output.size + vector.size

  output = vector(mode="numeric",output.size)

  res <-.C("mntd_moments",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), 
           as.integer(tree$edge), as.numeric(tree$edge.length), as.character(tree$tip.label), 
           as.integer(vector.size), as.integer(sample.sizes), as.logical(comp.expectation), 
           as.logical(comp.deviation), as.numeric(output),as.character(error.message),as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[12]] > 0)
    stop(substr(res[[11]],1,res[[12]]))

  result = NULL

  if( comp.expectation == TRUE && comp.deviation == TRUE)
    result = matrix(res[[10]], nrow=output.size/2, ncol=2)
  else
    result = res[[10]]

  return (result)

} # mntd.moments = function(...)


cac.moments <- function(tree, chi, sample.sizes, k=2)
{
  number.of.leaves = length(tree$tip.label)

  if(any(is.na(sample.sizes))||any(sample.sizes < 0)||any(sample.sizes > number.of.leaves))
    stop("The list with the sample sizes contains values which are undefined, or out of range. \n\n")


  if(any(as.integer(sample.sizes) != sample.sizes))
    stop("The list with the sample sizes contains non-integer values.\n\n")


  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  vector.size = length(sample.sizes)

  output = vector(mode="numeric",k*vector.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  res <-.C("cac_moments",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), 
           as.integer(tree$edge), as.numeric(tree$edge.length), as.character(tree$tip.label),
           as.numeric(chi), as.integer(vector.size), as.integer(sample.sizes), as.integer(k), 
           as.numeric(output),as.character(error.message),as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[12]] > 0)
    stop(substr(res[[11]],1,res[[12]]))

  result = matrix(res[[10]], nrow=vector.size, ncol=k)

  return (result)

} # cac.moments = function(...)

cbl.moments = function(tree, sample.sizes, comp.expectation = TRUE, comp.deviation = TRUE)
{
  size.data = unlist(sample.sizes)
  number.of.leaves = length(tree$tip.label)

  if(any(is.na(size.data)) || any(size.data < 0) || any(size.data > number.of.leaves) )
    stop("The list with the sample sizes contains values which are undefined, or out of range. \n\n")


  if(any(as.integer(size.data) != size.data))
    stop("The list with the sample sizes contains non-integer values.\n\n")


  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  vector.size = nrow(sample.sizes)

  output.size = 0 

  if(comp.expectation == TRUE)
    output.size = output.size + vector.size

  if(comp.deviation == TRUE)
    output.size = output.size + vector.size

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  res <-.C("cbl_moments",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), 
           as.integer(tree$edge), as.numeric(tree$edge.length), as.character(tree$tip.label), 
           as.integer(2*vector.size), as.integer(size.data), as.logical(comp.expectation), 
           as.logical(comp.deviation), as.numeric(output),as.character(error.message),as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[12]] > 0)
    stop(substr(res[[11]],1,res[[12]]))

  result = NULL

  if( comp.expectation == TRUE && comp.deviation == TRUE)
    result = matrix(res[[10]], nrow=output.size/2, ncol=2)
  else
    result = res[[10]]

  return (result)

} # cbl.moments = function(...)


cd.moments = function(tree, sample.sizes, comp.expectation = TRUE, comp.deviation = TRUE)
{
  size.data = unlist(sample.sizes)
  number.of.leaves = length(tree$tip.label)

  if(any(is.na(size.data)) || any(size.data < 0) || any(size.data > number.of.leaves) )
    stop("The list with the sample sizes contains values which are undefined, or out of range. \n\n")


  if(any(as.integer(size.data) != size.data))
    stop("The list with the sample sizes contains non-integer values.\n\n")

  #dyn.load(paste("PhyloMeasures", .Platform$dynlib.ext, sep = ""))

  number.of.edges = length(tree$edge)/2
  vector.size = nrow(sample.sizes)

  output.size = 0 

  if(comp.expectation == TRUE)
    output.size = output.size + vector.size

  if(comp.deviation == TRUE)
    output.size = output.size + vector.size

  output = vector(mode="numeric",output.size)
  error.message = paste0(rep(" ",length = 150), collapse="")
  message.size = 0

  res <-.C("cd_moments",as.integer(length(tree$edge.length)), as.integer(number.of.leaves), 
           as.integer(tree$edge), as.numeric(tree$edge.length), as.character(tree$tip.label), 
           as.integer(2*vector.size), as.integer(size.data), as.logical(comp.expectation), 
           as.logical(comp.deviation), as.numeric(output),as.character(error.message),as.integer(message.size),
           PACKAGE = "PhyloMeasures")

  if(res[[12]] > 0)
    stop(substr(res[[11]],1,res[[12]]))

  result = NULL

  if( comp.expectation == TRUE && comp.deviation == TRUE)
    result = matrix(res[[10]], nrow=output.size/2, ncol=2)
  else
    result = res[[10]]

  return (result)

} # cd.moments = function(...)
