## Copyright (c) 2014, Pioneer Hi-Bred International, Inc.

## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## Redistributions of source code must retain the above copyright
## notice, this list of conditions and the following disclaimer.

##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.

##     Neither the name of Pioneer Hi-Bred International, Inc. nor the
##     names of its contributors may be used to endorse or promote
##     products derived from this software without specific prior written
##     permission.

## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


globalVariables(c('ps'))

evalProb <- function(node.prob, x = 0, y = 0, z = 0, chunk.size = min(length(node.prob),75)){

    ## remove redundancies
    tab.prob <- table(node.prob)
    node.prob <- paste(tab.prob, names(tab.prob), sep = "*")
    
    ## in parts
    prob.split <- split(node.prob, ceiling(seq_along(node.prob)/chunk.size))
    
    value <- foreach(ps = 1:length(prob.split), .combine = "+") %do% {
        raw <- paste(prob.split[[ps]], collapse = "+")
        eval(parse(text = raw))
    }

    return(value)
    
}

