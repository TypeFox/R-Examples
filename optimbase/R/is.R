# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#

is.optimbase <- function(x=NULL){
  
  inherits(x,'optimbase')
  
}

is.optimbase.functionargs <- function(x=NULL){
  
  inherits(x,'optimbase.functionargs')
  
}

is.optimbase.outputargs <- function(x=NULL){
  
  inherits(x,'optimbase.outputargs')
  
}