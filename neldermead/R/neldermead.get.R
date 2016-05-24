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

neldermead.get <- function(this=NULL,key=NULL){
  
  if (!any(key==c('method','coords0','simplex0method','simplex0length',
                  'simplex0deltausual','simplex0deltazero','rho','chi',
                  'gamma','sigma','tolsimplexizeabsolute','tolsimplexizerelative',
                  'tolsimplexizemethod','toldeltafv','tolssizedeltafvmethod',
                  'restartmax','restarteps','restartstep','kelleystagnationflag',
                  'kelleynormalizationflag','kelleystagnationalpha0',
                  'restartflag','restartdetection','restartsimplexmethod',
                  'checkcostfunction','boxnbpoints','boxineqscaling',
                  'scalingsimplex0','guinalphamin',
                  'boxtermination','boxtolf','boxnbmatch','boxreflect',
                  'mymethod','myterminate','myterminateflag','tolvarianceflag',
                  'tolabsolutevariance','tolrelativevariance','greedy',
                  'historysimplex','simplexopt','simplex0','restartnb'))){
    # Delegate to the optimization object
    value <- optimbase.get(this=this$optbase,key=key)
  } else {
    if (key=='historysimplex'){
      storehistory <- optimbase.get(this=this$optbase,key='storehistory')
      if (!storehistory){
        stop('neldermead.get: History disabled ; turn on -storehistory option.',
          call.=FALSE)
      }
    }
    
    value <- this[[key]]
  }

  return(value)

}

