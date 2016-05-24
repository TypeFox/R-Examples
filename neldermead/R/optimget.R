# Copyright (C) 2008-2009 - INRIA - Michael Baudin
# Copyright (C) 2009-2010 - DIGITEO - Michael Baudin
# Copyright (C) 2010-2015 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#
# This source code is a R port of the neldermead component
# originally written by Michael Baudin for Scilab :
# "Nelder-Mead User's Manual", 2010, Consortium Scilab - Digiteo,
# Michael Baudin, http://wiki.scilab.org/The_Nelder-Mead_Component

optimget <- function(options=NULL, key=NULL, value=NULL){

  if (is.null(options)){
    stop('optimget: No option list provided.',
         call.=FALSE)
  }

  if (is.null(key)){
    stop('optimget: No option key provided.',
         call.=FALSE)
  }

  if (length(key)>1){
    stop('optimget: More than one key provided.',
         call.=FALSE)
  }

  # Search the field by index
  fields <- names(options)

  # Search for the given key in the list of available fields.
  # Use a regexp which ignores the case
  r <- grep(toupper(key),toupper(fields))
  opsize <- size(r,2)

  if (opsize==0){
    cat(sprintf('optimget: No match found between key (%s) and options elements (%s).\n',
                key, paste(names(options),collapse=', ')))
    return(value)
  } else if (opsize!=1){
    matching <- paste(r[1:opsize], collapse=' ')
    stop(sprintf('optimget: Ambiguous property name matches several fields %s\n',
                matching),
         call.=FALSE)
  }
  
  # Get the matching field
  name <- fields[r[1]]
  val <- options[[name]]
  
  # When the value is given and field is empty, return the value
  if (!is.null(value) & is.null(val)){
    val <- value
  }

  return(val)

}
