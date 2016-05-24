# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.

#' @title Place Legend
#' @description This function decides where the legend should be put (top left or bottom left)
#' @param wv_1 A \code{double} that indicates the first value of \code{wv.empir}
#' @param low_n A \code{doble} that indicates the last value of \code{ci.low}
#' @param high_n A \code{dobule} that indicates the last value of \code{ci.high}
#' @return A numeric vector containing 4 elements. The first two elements indicate legend.justification, the last two elements indicate legend.position (see \code{?theme}).
#' @keywords internal
placeLegend = function(wv_1, low_n, high_n){
  if(log10(wv_1) > ( log10(low_n) + log10(high_n) )/2 ){
    # legend should be placed in bottom left
    legend.justification = c(0,0)
    legend.position = c(0,0)
    #x = wv_1[1]/xlim_length
    #y = high_n/ylim_length
  }
  else{
    # legend should be placed in top left
    legend.justification = c(0,1)
    legend.position = c(0,1)
    #x = wv_1[1]/xlim_length
    #y = low_n/ylim_length
  }
  return( c(legend.justification, legend.position) )
  
}

#' @title Frequent Graph Setting for Paper
#' @description This function sets some parameters such as plot.margin.
#' @return A ggplot2 panel containing the frequent graph setting for paper.
#' @keywords internal
paperSetting = function(){
  p = theme(axis.title.y=element_text(vjust=3.5), axis.title.x=element_text(vjust=-2.5), 
            plot.title = element_text(vjust=2.5), plot.margin=unit(c(0.7,0.1,0.7,0.7),"cm"))
  return(p)
}

#' @title Emulate ggplot2 default color palette
#' @description Autogenerate a colors according to the ggplot selection mechanism. 
#' @param n An \code{integer} indicating how many colors user wants.
#' @return A \code{vector} containing \code{n} colors
#' @author John Colby
#' @keywords internal
ggColor <- function(n) {
  hues = seq(15, 375, length=n+1)
  rev(hcl(h=hues, l=70, c=100)[1:n])
}

#' @title Get the model in a \code{gmwm} object
#' @description Extracts and formats the model string.
#' @param object A \code{gmwm} object
#' @return A \code{string} containing the model
#' @keywords internal
getModel.gmwm = function(object){
  if( !is(object, 'gmwm') ){
    stop('It must be a gmwm object')
  }
  model.desc = object$model.hat$desc
  count.map = count_models(model.desc)
  all.model = names(count.map)
  
  model = ''
  for (i in 1:length(all.model)){
    if(length(model)==1){
      if(count.map[all.model[i]] == 1){
        model = bquote(.(all.model[i]))
      }else if(count.map[all.model[i]]>1){
        model = bquote(.(count.map[all.model[i]])%*%.(all.model[i]))
      }
    }else{
      if(count.map[all.model[i]] == 1){
        model = bquote(.(model) * "+" *.(all.model[i]))
      }else if(count.map[all.model[i]]>1){
        model = bquote(.(model) * "+" *.(count.map[all.model[i]])%*%.(all.model[i]) )
      }
    }
  }
  return(as.expression(model) )
}


#' @title Order the Model
#' @description Orders the model and changes it to the correct format
#' @param models A vector of \code{string} that specifies the models
#' @details If the \code{models} are c("AR1", "WN", "AR1", "WN", "AR1+WN+AR1+WN"), it will be converted to 
#' c("AR1-1", "WN-1", "AR1-2", "WN-2", "AR1+WN+AR1+WN").
#' 
#' This function is used in \code{gen.lts()}
#' @keywords internal
#' @examples 
#' models = c("AR1", "WN", "AR1", "WN", "AR1+WN+AR1+WN")
#' new.models = orderModel(models)
#' new.models
#' 
#' models = c('AR1', 'QN', 'WN', 'AR1+QN+WN')
#' new.models = orderModel(models)
#' new.models
orderModel = function(models){
  count = table(models)
  if( any(count>1)){
    multi.models = names( count[count>1] )
    
    for(model in multi.models){
      num = count[model]
      models[models == model] = paste( rep(model, num), rep('-', num), 1:num, sep = '' )
    }
    
    return(models)
  }else{
    return(models)
  }
}
