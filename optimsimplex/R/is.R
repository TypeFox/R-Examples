# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#

is.optimsimplex <- function(x=NULL){
  
  inherits(x,'optimsimplex')
  
}

is.simplex <- function(x=NULL){
  
  inherits(x,'simplex')
  
}

is.vertex <- function(x=NULL){
  
  inherits(x,'vertex')
  
}