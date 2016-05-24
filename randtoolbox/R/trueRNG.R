## 
# @file  trueRNG.R
# @brief R file using the random package
#
# @author Christophe Dutang
#
#
# Copyright (C) 2009, Christophe Dutang, 
# All rights reserved.
#
# The new BSD License is applied to this software.
# Copyright (c) 2009 Christophe Dutang. 
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
#          - Neither the name of its contributors may be used to endorse or promote 
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
### true random generation
###
###			R functions
### 


### true random generation ###

#	trueRNG <- function(n, dim = 1)
#	{
#		if(n <0 || is.array(n))
#			stop("invalid argument 'n'")
#		if(dim < 0 || length(dim) >1)
#			stop("invalid argument 'dim'")
#
#		#constraints on calls to www.random.org
#		maxsequence <- 10000
#		nbcall <- ceiling(n / maxsequence)
#		
#		res <- matrix(0, nrow=n, ncol=dim)
#		if(nbcall > 1)
#		{
#			for(i in 1:(nbcall-1))
#			{
#				subset <- ( maxsequence * (i-1) + 1):(maxsequence * i)
#	#			print(range(subset))
#				options(show.error.messages=FALSE, warn=-1)
#				testtry <- try(res[subset, ] <- randomNumbers(length(subset), min=0, max=2^16-1, col=dim, base=10, check=FALSE)/2^16)	
#				options(show.error.messages=TRUE, warn=0)
#				
#				if(class(testtry) == "try-error")
#				{
#					stop(paste("www.random.org service no longer available at ",i,"th call.\n", sep=""))
#				}
#			}
#		}else
#		{ 
#			i <- 1
#		}
#		
#		#final subset to be filled up
#		subset <- ( maxsequence * (i-1) + 1):(n)
#	#			print(range(subset))
#		options(show.error.messages=FALSE, warn=-1)
#		testtry <- try(res[subset, ] <- randomNumbers(length(subset), min=0, max=2^16-1, col=dim, base=10, check=FALSE)/2^16)	
#		options(show.error.messages=TRUE, warn=0)
#
#		if(class(testtry) == "try-error")
#		{
#			stop(paste("www.random.org service no longer available at ",i,"th call.\n", sep=""))
#		}
#
#		
#		#return result
#		if(dim == 1)
#			as.vector(res)
#		else
#			as.matrix(res)
#	}
#

