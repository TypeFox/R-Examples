## 
# @file  runifInterface.R
# @brief R file for runif interface
#
# @author Petr Savicky
#
#
# Copyright (C) 2009, Petr Savicky, Academy of Sciences of the Czech Republic.
# All rights reserved.
#
# The new BSD License is applied to this software.
# Copyright (c) 2009 Petr Savicky. 
# All rights reserved.
#
#      Redistribution and use in source and binary forms, with or without
#      modification, are permitted provided that the following conditions are
#      met:
#      
#          - Redistributions of source code must retain the above copyright
#          notice, this list of conditions and the following disclaimer.
#          - Redistributions in binary form must reproduce the above
#          copyright notice, this list of conditions and the following
#          disclaimer in the documentation and/or other materials provided
#          with the distribution.
#          - Neither the name of the Academy of Sciences of the Czech Republic
#          nor the names of its contributors may be used to endorse or promote 
#          products derived from this software without specific prior written
#          permission.
#     
#      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#      "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#      LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#      A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#      OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#      SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#      LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#      DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#      THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
#
#############################################################################
### runif interface
###
###			R functions
### 


set.generator <- function(name=c("congruRand", "WELL", "MersenneTwister", "default"), parameters=NULL, seed=NULL, ...,
	only.dsc=FALSE)
{
	name <- match.arg(name)
	dots <- list(...)
	if (name == "congruRand")
	{
		if (is.null(parameters))
			parameters <- c(mod=dots$mod, mult=dots$mult, incr=dots$incr)
		if (length(parameters) == 0)
			parameters <- c(mod="2147483647", mult="16807", incr="0")
		if (!identical(names(parameters), c("mod", "mult", "incr")))
		{
			param.names <- paste(names(parameters),collapse=" ")
			stop("parameter list \"", param.names, "\" is not correct for congruRand")
		}
		if (is.null(seed))
			seed <- floor(as.double(parameters["mod"]) * runif(1))
		if (is.numeric(parameters))
			parameters <- formatC(parameters, format="f", digits=0)
		if (is.numeric(seed))
			seed <- formatC(seed, format="f", digits=0)
		state <- c(seed=seed)
		description <- list(name=name, parameters=parameters, state=state)
	} else if (name == "WELL")
	{
		if (is.null(parameters))
		{
			if (is.null(dots$order)) dots$order <- ""
			if (is.null(dots$version)) dots$version <- ""
			if (dots$order != "" & nchar(dots$version) != 1) {
				stop("unsupported parameters order=", dots$order, ", version=", dots$version," for WELL")
			}
			version.name <- paste(dots$order, dots$version, sep="")
			order <- substr(version.name, 1, nchar(version.name) - 1)
			version <- substr(version.name, nchar(version.name), nchar(version.name))
			parameters <- c(order=order, version=version)
		}
		if (!identical(names(parameters), c("order", "version")))
		{
			param.names <- paste(names(parameters),collapse=" ")
			cat("parameters required for WELL: order, version\n")
			cat("parameters provided: ", param.names, "\n")
			stop("parameter list is not correct for WELL")
		}
		if (! paste(parameters, collapse="") %in% c("512a", "521a", "521b", "607a", "607b", "800a", "800b", "1024a", "1024b",
			"19937a", "19937b", "19937c", "21701a", "23209a", "23209b", "44497a", "44497b"))
			stop("unsupported parameters order=", parameters["order"], ", version=", parameters["version"]," for WELL")
		if (is.null(seed))
			seed <- floor(2^31 * runif(1))
		size <- ceiling(as.numeric(parameters["order"])/32)
		state <- doinitMT2002(seed, size, size)[[3]]
#		state <- .C("initMT2002",
#					as.integer(seed),
#					as.integer(size),
#					integer(size),
#					PACKAGE="rngWELL")[[3]]
		description <- list(name=name, parameters=parameters, state=state)
	} else if (name == "MersenneTwister")
	{
		if (is.null(parameters))
			parameters <- c(initialization=dots$initialization, resolution=dots$resolution)
		if (!identical(names(parameters), c("initialization", "resolution")))
		{
			param.names <- paste(names(parameters),collapse=" ")
			stop("parameter list \"", param.names, "\" is not correct for MersenneTwister")
		}
		type <- match(parameters["initialization"], c("init2002", "array2002"), nomatch=0)
		if (type == 0)
			stop("initialization ", parameters["initialization"], " is not in c(\"init2002\", \"array2002\")")
		if ( ! parameters["resolution"] %in% c("32", "53"))
			stop("resolution \"", parameters["resolution"], "\" is not in c(\"32\", \"53\")")
		if (is.null(seed))
			seed <- floor(2^31 * runif(1))
		state <- .C("initMersenneTwister",
					as.integer(type),
					length(seed),
					as.integer(seed),
					state=integer(625),
					PACKAGE="randtoolbox")$state
		description <- list(name=name, parameters=parameters, state=state)
	} else if (name == "default")
	{
		RNGkind("default")
		if (!is.null(seed))
			set.seed(seed)
		return(invisible(NULL))
	} else
		stop("unsupported generator: ", name)
	if (only.dsc)
		return(description)
	put.description(description)
	invisible(NULL)
}

put.description <- function(description)
{
	name <- description$name
	parameters <- description$parameters
	state <- description$state
	if (name == "congruRand")
	{
		aux <- .C("put_state_congru",
			parameters,
			state,
			err = integer(1),
			PACKAGE="randtoolbox")
		if (aux$err != 0)
			stop("check congruRand error: ", aux$err)
		if (RNGkind()[1] != "user-supplied")
		{
			.C("set_noop", PACKAGE="randtoolbox")
			RNGkind("user-supplied")
			aux <- .C("put_state_congru",
				parameters,
				state,
				err = integer(1),
				PACKAGE="randtoolbox")
			if (aux$err != 0)
				stop("check congruRand error: ", aux$err)
		}
	} else if (name == "WELL")
	{
		.C("set_noop", PACKAGE="randtoolbox")
		RNGkind("user-supplied")
		doputRngWELL(parameters["order"], parameters["version"], state)
#		.C("putRngWELL",
#			as.integer(parameters["order"]),
#			match(parameters["version"], c("a", "b", "c"), nomatch=0),
#			as.integer(state),
#			PACKAGE="rngWELL")
	} else if (name == "MersenneTwister")
	{
		.C("set_noop", PACKAGE="randtoolbox")
		RNGkind("user-supplied")
		.C("putMersenneTwister",
			match(parameters["initialization"], c("init2002", "array2002"), nomatch=0),
			as.integer(parameters["resolution"]),
			as.integer(state),
			NAOK=TRUE,
			PACKAGE="randtoolbox")
	} else 
		stop("unsupported generator: ", name)
	invisible(NULL)
}

get.description <- function()
{
	if (RNGkind(NULL)[1] != "user-supplied")
		stop("For R base generators, use .Random.seed, not get.description()")
	generator <- .C("current_generator",
		integer(1),
		PACKAGE="randtoolbox")[[1]]
	if (generator == 1)
	{
		name <- "congruRand"
		outspace <- "18446744073709551616" # 2^64
		aux <- .C("get_state_congru",
			parameters=rep(outspace, times=3),
			seed=outspace,
			PACKAGE="randtoolbox")
		parameters <- aux$parameters
		names(parameters) <- c("mod", "mult", "incr")
		seed <- aux$seed
		state <- c(seed=aux$seed)
		if(parameters[1] == "4294967296" && parameters[2] == "1664525" && parameters[3] == "1013904223")
			literature <- "Knuth - Lewis"
		else if(parameters[1] == "281474976710656" && parameters[2] == "31167285" && parameters[3] == "1")
			literature <- "Lavaux - Jenssens"
		else if(parameters[1] == "18446744073709551616" && parameters[2] == "636412233846793005" && parameters[3] == "1")
			literature <- "Haynes"
		else if(parameters[1] == "4294967296" && parameters[2] == "69069" && parameters[3] == "0") 
			literature <- "Marsaglia"
		else if(parameters[1] == "4294967295" && parameters[2] == "16807" && parameters[3] == "0") 
			literature <- "Park - Miller"
		else 
			literature <- "Unknown"
	} else if (generator == 2)
	{
		name <- "WELL"
		tmp <- dogetRngWELL(1, 1, 2000)
#		tmp <- .C("getRngWELL",
#			order = integer(1),
#			version = integer(1),
#			state = integer(2000),
#			PACKAGE="rngWELL")
		order <- as.character(tmp$order)
		print(tmp)
		version <- letters[tmp$version]
		parameters <- c(order=order, version=version)
		size <- ceiling(tmp$order/32)
		state <- tmp$state[1:size]
		literature <- "Panneton - L'Ecuyer - Matsumoto"
	} else if (generator == 3)
	{
		name <- "MersenneTwister"
		tmp <- .C("getMersenneTwister",
			initialization = integer(1),
			resolution = integer(1),
			state = integer(625),
			PACKAGE="randtoolbox")
		initialization <- c("init2002", "array2002")[tmp$initialization]
		resolution <- as.character(tmp$resolution)
		parameters <- c(initialization=initialization, resolution=resolution)
		state <- tmp$state
		literature <- "M. Matsumoto, T. Nishimura, 1998"
	} else
		stop("internal error of randtoolbox")
	list(name=name, parameters=parameters, state=state, authors=literature)
}

