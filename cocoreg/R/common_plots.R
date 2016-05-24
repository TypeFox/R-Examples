#Copyright (c) 2016 Jussi Korpela (Finnish Institute of Occupational Healt, jussi.korpela@ttl.fi) and Andreas Henelius (Finnish Institute of Occupational Healt, andreas.henelius@iki.fi)
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

## -----------------------------------------------------------------------------
## Plotting tools for the data formats used by cocoreg  
# dataset = a data.frame with several variables
# data collection = a list of datasets i.e. a list of data.frames
## -----------------------------------------------------------------------------

#' Plotting data collections using ggplot
#'  
#'  @param dc_lst A list of data collections i.e. a list of lists of data.frames (see examples)
#'  @param ylim (optional) y-axis limits as [1,2] numeric, passed on to dfplot() as 'ylim'
#'  @param titleArr (optional)  Plot column titles as [1, length(dc_lst)] string array
#'  @param legendMode (optional) Where to put legend, allowed values c('none','first','all')
#'  @param dfplot (optional) Function used to plot a data.frame (one panel in final plot)
#'  
#'  @return Produces a plot to the active graphics device
#'  
#' @examples
#' \dontrun{
#' dc <- create_syn_data_toy()
#' ccr <- cocoreg(dc$data)
#' ggplot_dclst(list(d1 = dc$data, d2 = ccr$data, dn = dc$data))
#' }
#'  
#'  @importFrom gridExtra grid.arrange
#'  @import ggplot2
#'  
#'  @export
ggplot_dclst <- function(dc_lst,
                         ylim=NULL,
                         titleArr = names(dc_lst),
                         legendMode = 'none',
                         dfplot = ggplot_df){
  
#   if (is.null(dfplot)){
#     dfplot = function(df, ylim=NULL, titlestr=NULL){
#       p <- ggplot_df(df, ylim=ylim, titlestr=titlestr)
#       p} 
#   }
  
  if (is.null(titleArr)){
    titleArr = paste0('dc', 1:length(dc_lst))
  }
  
  p_lst <- vector("list", length(dc_lst))
  for (i in 1:length(dc_lst)){
    #p_lst[i] <- list( ggplot_dflst(dc_lst[[i]], plotfun=dfplot,
    #                              plot=F, ylim=ylim, titlestr="test"))
    M <- length(dc_lst[i])
    tmpPLst <- c( lapply(dc_lst[[i]][1], dfplot, ylim=ylim, titlestr=titleArr[i]),
                  lapply(dc_lst[[i]][2:length(dc_lst[[i]])], dfplot, ylim=ylim) )
    
    # Switch legends on/off
    tmpPLst <- switch(legendMode,
                      all = tmpPLst,
                      first = if (i > 1){ lapply(tmpPLst, function(p){p <- p + theme(legend.position="none")}) } else {tmpPLst} ,
                      none = lapply(tmpPLst, function(p){p <- p + theme(legend.position="none")}) 
    )
    
    #     if (i > 1){
    #       tmpPLst <- lapply(tmpPLst, function(p){p <- p + theme(legend.position="none")})
    #     }
    #browser()
    p_lst[i] <- list(tmpPLst)
  }
  
  p_lst <- nplst_reorder_grid(p_lst, length(dc_lst))
  do.call(gridExtra::grid.arrange, c(p_lst, list(ncol=length(dc_lst))) )  
}


#' Compare data collections variable by variable
#' 
#' @param dclst A (named) list of data collections i.e. a list of lists of data.frames (see examples)
#'  
#' @return Returns a ggplot object (which is by default printed if not assigned to variable)
#' 
#' @import ggplot2
#'  
#' @export
ggcompare_dclst <- function(dclst){
  
  if (is.null(names(dclst))){
    names(dclst) <- sprintf('DC %d', 1:length(dclst))
  }
  
  pd <- data.frame()
  for (i in 1:length(dclst)){
    dtmp <- dflst2dfmelt(dclst[[i]])
    dtmp$dc <- names(dclst)[[i]]
    pd <- rbind(pd, dtmp)
  }
  
  pd$vards <- as.character(pd$variable)
  pd$vards <- gsub('\\_.*','', pd$vards)
  pd$varno <- as.character(pd$variable)
  pd$varno <- gsub('x.*\\_','', pd$varno)
  # Lines above do roughly the same as: 
  #pd <- tidyr::separate(pd, variable, c('vards','varno'))
  
  pd$varname <- sprintf('variable %s', pd$varno)
  
  p <- ggplot2::ggplot(data = pd, aes_string(group = 'dc', color = 'dc'))
  p <- p + ggplot2::geom_line(aes_string(x = 'obs', y = 'value'))
  p <- p + ggplot2::facet_grid(dataset ~ varname)
  p
}


#' Plot a list of data.frames using ggplot2
#' 
#'  @param dflst A list of datasets as a list of data.frames
#'  @param ncol (optional) Number of columns in final plot
#'  @param plot (optional) Plot or not: if TRUE produces a plot else returns a list of ggplot objects
#'  @param plotfun (optional) Function used to plot a data.frame (one panel in final plot)
#'  @param ... (optional) Additional parameters passed on to plotfun
#'  
#'  @return Produces a plot to the active graphics device or returns a list of ggplot objects
#'  
#' @examples
#' \dontrun{
#' dc <- create_syn_data_toy()
#' ggplot_dflst(dc$data)
#' }
#'  @importFrom gridExtra grid.arrange
#'
#'  @export
ggplot_dflst <- function(dflst, ncol = 1, plot = T, plotfun = ggplot_df, ...){
  plst <- lapply(dflst, plotfun, ...)
  if (plot){
    do.call(grid.arrange, c(plst, list(ncol=ncol)))
  } else {
    plst
  }
}


#' Plotting data.frame using ggplot
#' 
#' @param df A data.frame to plot
#' @param titlestr (optional) Title of plot as string
#' @param ylabstr (optional) Y-axis label as string
#' @param ylim (optional) y-axis limits as [1,2] numeric, passed on to dfplot() as 'ylim'
#' @param color (optional) Input for manual color scale
#' @param linetype (optional) Input for manual linetype scale
#' @param logy (optional) Should y-axis be logarithmic? A boolean value.
#'  
#' @return Returns a ggplot2 object
#'  
#' @examples
#' \dontrun{
#' dc <- create_syn_data_toy()
#' ggplot_df(dc$data[[1]])
#' }
#' 
#' @import ggplot2
#' 
#' @export
ggplot_df <- function(df, titlestr=NULL, ylabstr=NULL, ylim=NULL, color=NULL,
                      linetype=NULL, logy=F) {
  n_views <- length(df)
  
  if (is.null(ylim)){
    ylim <- c(min(df), max(df))
  }

  df <- df_ggplot_melt(df) #melt into ggplot compatible format
  
  p <- ggplot(data=df)
  p <- p + geom_line(aes_string(x = 't', y = 'value', group = 'variable',
                                color = 'variable', linetype = 'variable'))
  if (logy){
    p <- p + scale_y_log10(limits=ylim)
  } else {
    p <- p + scale_y_continuous(limits=ylim)  
  }
  
  if (!is.null(color)){
    p <- p + scale_color_manual(values=color)   
  } 
  if (!is.null(linetype)){
    p <- p + scale_linetype_manual(values=linetype)   
  } else {
    p <- p + scale_linetype_manual(values=rep(1,length(unique(df$variable))) )  
  }
  p <- p + theme_bw()
  if (!is.null(titlestr)){
    p <- p + ggtitle(titlestr)
  }
  if (!is.null(ylabstr)){
    p <- p + ylab(ylabstr)
  }
  p
}


######################################################################################
## Tools needed in the plotting functions
######################################################################################

#' Reorders a nested list of ggplots
#' 
#' @description
#' Reorders a nested list of ggplots to ncol columns prior to calling grid.arrange()
#' Note: p_list is a list of lists of ggplots.
#' p_list = list(p_list1, p_list2,...)
#' 
#' @param p_list A list of lists of ggplots
#' @param ncol Target number of columns, integer value
#'  
#' @return  A reordered and flattened version of input list as a list of ggplot2 objects
#'  
#' @export
nplst_reorder_grid <- function(p_list, ncol){
  p_list_out <- vector("list", length(p_list)*length(p_list[[1]]))
  for(i in 1:ncol){
    p_list_out[seq(i, length(p_list_out), ncol)] <- p_list[[i]]  
  }
  p_list_out
}


#' Melt data.frame into ggplottable format
#' 
#' @description 
#' Melts a data.frame into format that is suitable for use with ggplot2.
#' Creates the time variable 't' used by plotting functions.
#' 
#' @param df A data.frame
#' 
#' @return A ggplot2 compatible data.frame with time variable
#' 
#' @examples
#' \dontrun{
#' dc <- create_syn_data_toy()
#' df <- dc$data[[1]]
#' str(df)
#' str(df_ggplot_melt(df))
#' }
#' 
#' @importFrom reshape melt melt.data.frame
#' @export
df_ggplot_melt <- function(df){
  if (!is.null( attr(df,'orig_dimnames') )){
    names(df) <- attr(df,'orig_dimnames')[[2]]
    timevec <- attr(df,'orig_dimnames')[[1]]
    if (!is.numeric(timevec)){
      timevec <- as.numeric(timevec)
    }
  } else {
    timevec <- as.numeric(dimnames(df)[[1]])
    if (sum(is.na(timevec))>1){
      timevec <- 1:nrow(df)
    }
  }
  def_names <- sprintf("X%d",1:ncol(df))
  names(df) <- paste0(def_names,": ",names(df))
  
  df$t <- timevec
  df <- reshape::melt(df, id.vars = "t")
} 

