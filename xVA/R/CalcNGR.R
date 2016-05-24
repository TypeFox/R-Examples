#' Calculates the Net/Gross ratio used under the CEM regulatory framework
#' @title Calculates the Net/Gross ratio (NGR)
#' @param MtM_Vector A vector containing the trades to be netted
#' @return The Net-Gross ratio (NGR)
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#'
CalcNGR = function(MtM_Vector)
{
  if (all(MtM_Vector== 0)) 
  { NGR = 1  
  } else if(all(MtM_Vector< 0))
  { NGR = 0
  } else
  {   NGR = max(sum(MtM_Vector),0)/sum(max(MtM_Vector,0))
  }
  return(NGR)
}

