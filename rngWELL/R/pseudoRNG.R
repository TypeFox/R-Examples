## 
# @file  pseudoRNG.R
# @brief R file for all pseudo RNGs
#
# @author Christophe Dutang
# @author Petr Savicky 
#
#
# Copyright (C) 2009, Christophe Dutang, 
# Petr Savicky, Academy of Sciences of the Czech Republic. 
# All rights reserved.
#
# The new BSD License is applied to this software.
# Copyright (c) 2009 Christophe Dutang, Petr Savicky. 
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
### pseudo random generation
###
###			R functions
### 


### set the seed ###

setSeed4WELL <- function(seed)
	invisible( .Call("doSetSeed4WELL", seed) )


### pseudo random generation ###
 
WELL2test <- function(n, dim = 1, order = 512, temper = FALSE, version = "a")
{
    if(n <0 || is.array(n))
        stop("invalid argument 'n'")
    if(dim < 0 || length(dim) >1)
            stop("invalid argument 'dim'")
    if(!is.numeric(order))
            stop("invalid argument 'order'")
    if( !(order %in% c(512, 521, 607, 800, 1024, 19937, 21701, 23209, 44497) ) )
            stop("'order' must be in {512, 521, 607, 800, 1024, 19937, 21071, 23209, 44497}.")
    if( !(version %in% c("a", "b") ) )
            stop("'version' must be either 'a' or 'b'.")

    if(!is.logical(temper))
        stop("invalid argument 'temper'")
    if(temper && order %in% c(512, 521, 607, 1024))
        stop("tempering impossible")
    
    zeversion <- 0
    if(version == "a")
        zeversion <- 1
    if(version == "b")
        zeversion <- 2
    if(zeversion == 0)
        stop("wrong version for WELL RNG")
    if(version == "b" && order %in% c(512,  21701) ) 
        stop("this WELL RNG does not have a 'b' version")
    
    if(length(n) > 1)
        res <- .Call("doWELL", length(n), dim, order, temper, zeversion)
    else
        res <- .Call("doWELL", n, dim, order, temper, zeversion)	
    
	cat("Warning: you should use the function WELL from the randtoolbox package.\n")
    
    if(dim == 1)
        as.vector(res)
    else
        as.matrix(res)
}

doinitMT2002 <- function(seed, n, state)
{
	if(n <= state)
		.C("initMT2002", as.integer(seed), as.integer(n), integer(state))
	else
		NA
}

doputRngWELL <- function(order, version, state)
	.C("putRngWELL", as.integer(order), match(version, c("a", "b", "c"), nomatch=0), as.integer(state))

dogetRngWELL <- function(order, version, state)
	.C("getRngWELL", order=integer(order), version=integer(version), state=integer(state))



