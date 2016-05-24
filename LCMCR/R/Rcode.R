 # 
 # Copyright (C) 2012-2016 Daniel Manrique-Vallier
 # 
 # This program is free software; you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation; either version 2 of the License, or (at
 # your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful, but
 # WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 # General Public License for more details.
 # 
 # You should have received a copy of the GNU General Public License
 # along with this program; if not, write to the Free Software
 # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 #

##########################################
# Functions and objects definitions for 
#  manipulationg MCMC environments.
# (implementation with reference classes)

#this is a base class. Has to be used by inheritance. NEED TO SET THE POINTER!!!
#To do: VALIDATION!

setRefClass(
  Class="MCMCenviron",
  fields=list(
    pointer='externalptr', 
    blobsize = 'numeric',
	local_seed ='numeric'
  ), 
  methods=list(
    initialize = function(){
      #.self$initialized = FALSE
		usingMethods('Init_Model', 
		'Update',
		'Get_Iteration',
		'Get_Trace_List', 
		'Get_Param_List',
		'Get_Trace',
		'Get_Param',
		'Get_Trace_Size',
		'Set_Trace',
		'Reset_Traces', 
		'Activate_Tracing', 
		'Deactivate_Tracing', 
		'Change_SubSamp', 
		'Change_Trace_Length', 
		'Set_Seed') 
    },
    Init_Model = function(output = TRUE, seed = c('auto','r.seed')){
	"Initializes the sampler. Set output = FALSE to suppress console output. Seed can be an integer, 'auto' or 'r.seed'."
	  #seed the internal RNG
	  if (!is.numeric(seed)){
			seed <- match.arg(seed)
	  }
	  if (seed == 'r.seed'){
		runif(1)
		.self$local_seed <- Reduce(bitwXor, .Random.seed)
		.self$Set_Seed(.self$local_seed)
	  } else if (is.numeric(seed)){
		.self$local_seed <- as.integer(seed)
		.self$Set_Seed(.self$local_seed)
	  } #default behavior: do nothing. C++ class will seed automatically from clock.
	  #Control messages
	  if (!interactive()) output = FALSE
	  if (output){
		tmp <- .Call('R_Activate_Chain_Messages', .self$pointer)
	  } else {
		tmp <- .Call('R_Deactivate_Chain_Messages', .self$pointer)
	  }
      #SEXP R_Init_Model(SEXP p);
      tmp <- .Call('R_Init_Model', .self$pointer)
    },
    Update = function(iter, output = TRUE){
	"Runs num_iter iterations of the sampler. Set output = FALSE to suppress console output."
      #SEXP R_Update_Model(SEXP p, SEXP int_iter)
	  if (!interactive()) output = FALSE
	  if (output){
		tmp <- .Call('R_Activate_Updating_Output', .self$pointer)
	  } else {
		tmp <- .Call('R_Deactivate_Updating_Output', .self$pointer)
	  }
      tmp <- .Call('R_Update_Model', .self$pointer, as.integer(iter))
    },
    Get_Iteration = function(){
	"Retrieves the current number of iterations the sampler."
      return(.Call('R_Get_Iteration', .self$pointer))
    },
    Get_Trace_List = function(){
	"Retrieves the names of the parameters being currently traced."
      tmp <- .Call('R_Get_Trace_List', .self$pointer)
      return(tmp)
    },
    Get_Param_List = function(){
	"Retrieves the names of the parameters of the model."
      tmp <- .Call('R_Get_Param_List', .self$pointer)
      return(tmp)
    },
    Get_Trace = function(param){
	"Retrieves the tracing buffer for parameter param. Returns an array whose first dimension is the sample index."
      tmp <- .Call('R_Get_Trace', .self$pointer, as.character(param))
      return(aperm(tmp, length(dim(tmp)):1))
    },
    Get_Param = function(param){
	"Retrieves the current value of the parameter param."
      tmp <- .Call('R_Get_Param', .self$pointer, as.character(param))
      return(aperm(tmp, length(dim(tmp)):1))
    },
   Get_Trace_Size = function(){
   "Retrieves the size (in iterations) of the trace buffer."
      tmp <- .Call('R_Get_Trace_Size', .self$pointer)
      return(tmp)
    },
    Set_Trace = function(params){
	"Adds parameters to the tracer. To list the available parameters for tracing use the Get_Param_List() method."
      for (tra in params){
        tmp <- .Call('R_Set_Trace', .self$pointer, as.character(tra))
      }
    },
    Reset_Traces = function(){
	"Deletes the content of the tracing buffer."
      tmp <- .Call('R_Reset_Traces', .self$pointer)
    },
    Activate_Tracing = function(){
	"Activates the tracing buffer."
      tmp <- .Call('R_Activate_Tracing', .self$pointer)
    },
    Deactivate_Tracing = function(){
	"Deactivates the tracing buffer."
      tmp <- .Call('R_Deactivate_Tracing', .self$pointer)
    },
    Change_SubSamp = function(new_subsamp = 1){
	"Changes the sub-sampling period (thinning) of the tracing buffer."
      tmp <- .Call('R_Change_SubSamp', .self$pointer, as.integer(new_subsamp))
    },
    Change_Trace_Length = function(new_length){
	"Changes the size (in number of samples) of the tracing buffer."
      tmp <- .Call('R_Change_Trace_Size', .self$pointer, as.integer(new_length))
    },
	Set_Seed = function(seed){
	"Seeds the internal random number generator. Doesn't affect R's internal RNG."
		#SEXP R_Set_Seed(SEXP p, SEXP seed);
		tmp <- .Call('R_Set_Seed', .self$pointer, as.integer(seed))
		.self$local_seed <- seed
	},
	Get_Status = function(){
	#SEXP R_Get_Status(SEXP p)
	#Returns a vector c(<iteration>, <initialized>, <buffer size>, <buffer used>, <tracer activated>, <thinning>)  		
		tmp <- .Call('R_Get_Status', .self$pointer)
		return(list(
			iteration = tmp[1],
			initialized = (tmp[2] == 1),
			buffer_size = tmp[3],
			buffer_used = tmp[4],
			tracing = (tmp[5] == 1),
			thinning = tmp[6]
			)
		)	
	},
	show = function(){
		tmp <- .self$Get_Status()
		cat('\tInitialized =', tmp$initialized,'\n')
		cat('\tCurrent iteration =', tmp$iteration,'\n')
		cat('Tracer:\n')
		cat('\tActivated =', tmp$tracing, '\n')
		cat('\tCapacity =',tmp$buffer_size,'samples.\n')
		cat('\tUsed  = ', tmp$buffer_used, ' (', round(tmp$buffer_used/tmp$buffer_size*100, digits = 2), '%)\n', sep ='');
		cat('\tThining =', tmp$thinning,'\n')
		cat('\tCurrently Tracing:', paste(.self$Get_Trace_List(), collapse=', '), '\n');
	})
  )


###############################################################################
#Functions for manipulating matrices and dataframes with categorical values
# (c) Daniel Manrique-Vallier 2012
###############################################################################

fn_df_nlevels <- function(d){
    sapply(names(d), FUN = function(x)nlevels(d[,x]))
}

fn_df_discretize <- function(d, cols = 1:NCOL(d)){
  for (i in cols){
    d[,i] <- factor(d[,i])
  }
  return(d)
}
fn_df_levels <- function(d, cols = 1:NCOL(d)){

  res <- list()
  for (i in cols){
    res[[names(d)[i]]] <- levels(d[1,i])
  }
  return(res)
}

fn_numeric_to_factor <- function(
  x, list_levls = rep(list(1:2), NCOL(x)), cols = 1:NCOL(x), missing = -1, offset = 0){
  x <- as.data.frame(x)
  for (c in cols){
    lv <- list_levls[[c]]
    nlv <- length(lv)
    x[,c] <- factor(x[,c], levels = 1:nlv + offset)
    levels(x[,c]) <- lv
  }
  return(x)
}

fn_apply_levels_from <- function(dest, from, cols = 1:NCOL(from)){
  dest <- as.data.frame(dest)
  lev <- sapply(cols, function(x)levels(from[1,x]), simplify=FALSE)
  for (i in cols){
    dest[,i] <- factor(dest[,i], levels = lev[[i]])
  }
  names(dest[,cols]) <- names(from[,cols])
  return(dest)
}

fn_dataframe2num_matrix <- function(d, offset = -1, missing = -1, C_style= FALSE){
  offset <- as.integer(offset); missing <- as.integer(missing)
  r <- data.matrix(d) + offset
  r[is.na(r)] <- missing
  if(C_style){
    r <- t(r)
  }
  class(r) <- c('matrix', 'DMVmatrix')
  attr(r, 'levels') <- fn_df_levels(d)
  attr(r, 'na.code') <- missing
  attr(r, 'offset') <- offset
  attr(r, 'C_style') <- C_style
  return(r)
}

fn_df_2_C_Data <- function(df, prefix = 'array', filename = 'console'){
  #Use this function for generating test data. Works with arrays if factors.
  breakdown <- fn_dataframe2num_matrix(df, offset = -1, missing=-1, C_style=T)
  n <- NROW(df)
  J <- NCOL(df)
  l1 <- paste('int ',prefix,'_data_raw[]= {', paste(breakdown, collapse=', '), '};\n', sep='')
  l2 <- paste('int ', prefix,'_levels[]= {', paste(fn_df_nlevels(df), collapse=', '), '};\n', sep ='')
  l3 <- paste('int ',prefix, '_n_glob = ', n,';\n', 'int ', prefix,'_J_glob=', J,';\n', sep='')
  if(filename != 'console'){
    fl = file(description = filename, open = 'w')
    writeLines(strwrap(paste(l1, l2, l3, sep ='\n')), con = fl)
    close(fl)
  } else {
    writeLines(strwrap(paste(l1, l2, l3, sep ='\n')))
  }
}
######################################################################
# Functions for creating and manipulating NP_LCM objects
# Depends on functions in: 
# - ../GeneralCode/ArrayUtils.R
# - ../GeneralMCMC/MCMCenv_refClass.R
# (c) Daniel Manrique-Vallier 2013
######################################################################

#prototype:
#SEXP R_Create_LCM_CR(SEXP x_flat, SEXP J, SEXP n, SEXP covI, SEXP cov_nlevels, SEXP K, SEXP Nmis_max,
#						SEXP a_alpha, SEXP b_alpha, 
#						SEXP len_buffer, SEXP subsamp)

lcm_CR_Basic_generator <- setRefClass(
  Class = "lcm_CR_Basic",
  fields = list(
	J = 'numeric',
	K = 'numeric',
	n = 'numeric',
	Captures = 'data.frame'
    ),
  methods = list(
    initialize =   function(data_captures, K, a_alpha, b_alpha, in_list_symbol = '1',
			len_buffer, subsamp){
      callSuper()
	  zeros <- apply(as.matrix(data_captures != in_list_symbol), MARGIN = 1, FUN = prod)==1
	  if (sum(zeros) > 0){
		data_captures <- data_captures[!zeros,]
		warning('Inconsistent rows. ', sum(zeros), ' "no capture" rows eliminated.')
	  }
      .self$J = NCOL(data_captures)
      .self$K = K
      .self$n = NROW(data_captures)
      .self$Captures = data_captures
	  #SEXP R_Create_LCM_CR_Basic(SEXP x_flat, SEXP J, SEXP n, SEXP K, SEXP Nmis_max, 
		#				SEXP a_alpha, SEXP b_alpha, 
		#				SEXP len_buffer, SEXP subsamp);
      tmp <- .Call('R_Create_LCM_CR_Basic', 
				   as.integer(fn_factor2CR(data_captures, in_list_symbol)), 
				   as.integer(.self$J), 
                   as.integer(.self$n), 
				   as.integer(K), 
                   as.double(a_alpha), 
				   as.double(b_alpha), 
                   as.integer(len_buffer), as.integer(subsamp)
      )
      .self$pointer = tmp
    }
  ),
  contains = "MCMCenviron"
)
.onUnload <- function (libpath) {
  gc(verbose = FALSE)
  library.dynam.unload("LCMCR", libpath)
}
lcmCR <- function(captures, tabular = FALSE, in_list_label = '1', not_in_list_label = '0', K = 5, a_alpha = 0.25, b_alpha=0.25, buffer_size=10000, thinning = 10, seed = 'auto', verbose = TRUE){
	if(is.matrix(captures)){
		captures <- fn_CRmatrix2dataframe(captures, in_list_label=in_list_label, not_in_list_label = not_in_list_label, tabulated = tabular)
	}
	if (tabular){
		captures <- fn_CRtabular2indiv(captures)
	}
	o <- lcm_CR_Basic_generator(data_captures = captures, K, a_alpha, b_alpha, in_list_symbol = in_list_label, len_buffer = buffer_size, subsamp = thinning)
	o$Init_Model(output=verbose, seed = seed)
	return(o)
}

lcmCR_PostSampl <- function(object, burnin=10000, samples = 1000, thinning = 10, clear_buffer = FALSE, output = TRUE){
  object$Update(burnin, output)
  object$Change_SubSamp(thinning)
  object$Change_Trace_Length(samples)
  if (!('n0' %in% object$Get_Trace_List())){
    object$Set_Trace('n0')
  }
  if (clear_buffer){
    object$Reset_Traces()
  }
  object$Activate_Tracing()
  object$Update(samples * thinning, output)
  N <- object$Get_Trace('n0') + object$n
  return(as.numeric(N))
}

fn_factor2CR <- function(df.factors, sym.capture = '1', missing = -1, offset = -1, C_style = TRUE){
  #check that all variables have exactly 2 levels
  if (sum(fn_df_nlevels(df.factors) != 2) > 0 ) stop("Variable doesn't have 2 levels")
  if ( sum(sapply(fn_df_levels(df.factors), FUN = function(x) ! (sym.capture %in% x), simplify = T))  > 0)
    stop(paste( "'",sym.capture,"'", " not a level in all variables.", sep = ''  ))
  offset <- as.integer(offset); missing <- as.integer(missing)
  r <- ifelse(df.factors == sym.capture,as.integer(1),as.integer(0)) + offset + 1
  r[is.na(r)] <- missing
  if(C_style){
    r <- t(r)
  }
  class(r) <- c('matrix', 'DMVmatrix')
  attr(r, 'levels') <- fn_df_levels(df.factors)
  attr(r, 'na.code') <- missing
  attr(r, 'offset') <- offset
  attr(r, 'C_style') <- C_style
  return(r)  
}
fn_CRfact2OnesZeros <- function(data_vector, in_list_label = 'Yes'){
  if(!is.factor(data_vector)) stop("Data vector is not a factor.")
  if(!(in_list_label %in% levels(data_vector))) stop(in_list_label," is not a level in data_vector")
  return(ifelse(data_vector == in_list_label, 1, 0))
}
fn_CRdataframe2matrix <- function(data, in_list_label = 'Yes'){
  if(!is.data.frame(data)) stop('Object "data" is not a dataframe')
  lastindx <- NCOL(data)
  res <- matrix(NA, ncol = NCOL(data), nrow = NROW(data))
  if(names(data)[lastindx] != 'Freq'){
    cols <- 1:lastindx
  } else {
    cols <- 1:(lastindx - 1)
    res[,lastindx] <- data[,lastindx] 
  }
  for (i in cols){
    res[, i] <- fn_CRfact2OnesZeros(data[,i], in_list_label = in_list_label)
  }
  return(res)
}
fn_CRtabular2indiv <- function(tabular){
  indiv <- tabular[unlist(sapply(1:NROW(tabular), FUN = function(x)rep(x, tabular$Freq[x]))),1:(NCOL(tabular) - 1)]
  return(indiv)
}
fn_CRmatrix2dataframe <- function(data, in_list_label = 'Yes', not_in_list_label = 'No', tabulated = FALSE){
  #
  if(!is.matrix(data)) stop('Object "data" is not a matrix')
  res <- as.data.frame(data)
  lastindx <- NCOL(data)
  if(!tabulated){
    cols <- 1:lastindx
  } else {
    cols <- 1:(lastindx - 1)
    if (!is.numeric(res[,lastindx])) stop('Tabulation column is not numeric')
    res[,lastindx] <- data[,lastindx]
    names(res)[lastindx] <- 'Freq'
  }
  if (!setequal(unique(data[,cols]), c(0,1)) ) stop("Elements of matrix are not exclusively zeros and ones.")
  for (i in cols){
    res[, i] <- factor(data[, i], levels = c(0, 1), labels = c(not_in_list_label, in_list_label) )
  }
  return(res)
}
