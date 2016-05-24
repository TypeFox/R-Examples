## 
# @file  randtoolboxEnv.R
# @brief R file for environment of randtoolbox
#
# @author Yohann Chalabi
# @author Diethelm Wuertz 
#
#
# Copyright (C) 2008, Yohann Chalabi, Diethelm Wuertz, ETH Zurich. 
# All rights reserved.
#
# The new BSD License is applied to this software.
# Copyright (c) 2008 Yohann Chalabi, Diethelm Wuertz. 
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
#          - Neither the name of the ETH Zurich nor the names of its 
#          contributors may be used to endorse or promote 
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
#                                                                                                                                         #
#############################################################################
### environment functions
###
###			R functions
### 

.randtoolboxEnv <- new.env(hash = TRUE)

.setrandtoolboxEnv <-
    function(...)
{
    x <- list(...)
    nm <- names(x)
     if (is.null(nm) || "" %in% nm)
        stop("all arguments must be named")
    sapply(nm, function(nm) assign(nm, x[[nm]],
                                 envir = .randtoolboxEnv))
    invisible()
}

.getrandtoolboxEnv <-
    function(x = NULL, unset = "")
{
    if (is.null(x))
        x <- ls(all.names = TRUE, envir = .randtoolboxEnv)
###     unlist(mget(x, envir = .randtoolboxEnv, mode = "any",
###                 ifnotfound = as.list(unset)), recursive = FALSE)
    get(x, envir = .randtoolboxEnv, mode = "any")
}

