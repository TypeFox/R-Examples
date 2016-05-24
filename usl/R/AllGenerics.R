# Copyright (c) 2013 Stefan Moeding
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.


setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setGeneric("print", function(x, ...) standardGeneric("print"))
setGeneric("predict", function(object, ...) standardGeneric("predict"))
setGeneric("summary", function(object, ...) standardGeneric("summary"))
setGeneric("confint", function(object, parm, level) standardGeneric("confint"))

#
# USL specific methods
#
setGeneric("scalability",
           function(object, sigma, kappa) standardGeneric("scalability"))
setGeneric("peak.scalability",
           function(object, sigma, kappa) standardGeneric("peak.scalability"))
setGeneric("efficiency",
           function(object) standardGeneric("efficiency"))
setGeneric("overhead",
           function(object, newdata) standardGeneric("overhead"))
