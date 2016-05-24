#' Calculates the addon amount after splitting the trades into 'basis swap', 'volatility' and 'normal' transactions.
#' The corresponding penalty factors are applied to the supervisory factors for each trade group.
#' @title Calculates the Addon amount after handling basis and Volatility trades
#' @param trades The full list of the Trade Objects
#' @return The aggregate amount of the addon summed up for all the asset classes
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>
#' @references Basel Committee: The standardised approach for measuring counterparty credit risk exposures
#' http://www.bis.org/publ/bcbs279.htm
#'
HandleBasisVol <- function(trades)  {

  Addon_aggregate = 0
  vol_trades        = list()
  basis_swap_trades = list()
  normal_trades     = list()

  for (i in 1:length(trades))
  {
    if(is(trades[[i]],"Vol"))
      vol_trades[[length(vol_trades)+1]]=trades[[i]]
    else if(is(trades[[i]],"Swap"))
    {
      if(trades[[i]]$isBasisSwap())
      {
        basis_swap_trades[length(basis_swap_trades)+1] = trades[[i]]
      } else
      {
        normal_trades[length(normal_trades)+1] = trades[[i]]
      }
    }
    else
      normal_trades[length(normal_trades)+1] = trades[[i]]
  }
  if(length(basis_swap_trades)!=0)
  {
    swap_pairs = unique(lapply(basis_swap_trades, function(x) paste(x$pay_leg_ref,x$rec_leg_ref)))
    for (i in 1:length(swap_pairs))
    {
      split_pair = strsplit(swap_pairs[[i]]," ")
      group_trades <- basis_swap_trades[sapply(basis_swap_trades, function(x) x$pay_leg_ref==split_pair[[1]][1]&&x$rec_leg_ref==split_pair[[1]][2])]
      Addon_aggregate <- Addon_aggregate + CalcAddon(trades=group_trades,factor_mult=1.5)
    }
  }
  if(length(vol_trades)!=0)
  {
    vOl_refs = unique(lapply(vol_trades, function(x) x$reference))
    for (i in 1:length(vOl_refs))
    {
      #picking up the trades belonging to this specific hedging set
      group_trades <- vol_trades[sapply(vol_trades, function(x) x$reference==vOl_refs[i])]
      Addon_aggregate <- Addon_aggregate + CalcAddon(trades=group_trades,factor_mult=5)
    }
  }
  if(length(normal_trades)!=0)
    Addon_aggregate <- Addon_aggregate + CalcAddon(trades=normal_trades)

  return(Addon_aggregate)
}
