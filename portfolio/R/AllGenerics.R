################################################################################
##
## $Id: AllGenerics.R 387 2007-01-10 04:14:02Z enos $
##
## All generic functions for the portfolio class.
##
################################################################################

if(!isGeneric("create"))
  setGeneric("create", function(object, ...) standardGeneric("create"))

if(!isGeneric("scaleWeights"))
  setGeneric("scaleWeights", function(object, ...) standardGeneric("scaleWeights"))

if(!isGeneric("balance"))
  setGeneric("balance", function(object, in.var, ...) standardGeneric("balance"))

if(!isGeneric("exposure"))
  setGeneric("exposure", function(object, exp.var, ...) standardGeneric("exposure"))

if(!isGeneric("performance"))
  setGeneric("performance", function(object, ...) standardGeneric("performance"))

if(!isGeneric("totalReturn"))
  setGeneric("totalReturn", function(object, ...) standardGeneric("totalReturn"))

if(!isGeneric("portfolioDiff"))
  setGeneric("portfolioDiff", function(object, x, ...) standardGeneric("portfolioDiff"))

if(!isGeneric("contribution"))
  setGeneric("contribution", function(object, contrib.var, ...) standardGeneric("contribution"))

if(!isGeneric("securityInfo"))
  setGeneric("securityInfo", function(object, id, ...) standardGeneric("securityInfo"))

## Class portfolio only.

if(!isGeneric("calcWeights"))
  setGeneric("calcWeights", function(object, ...) standardGeneric("calcWeights"))

if(!isGeneric("calcShares"))
  setGeneric("calcShares", function(object, ...) standardGeneric("calcShares"))

if(!isGeneric("mvLong"))
  setGeneric("mvLong", function(object, ...) standardGeneric("mvLong"))

if(!isGeneric("mvShort"))
  setGeneric("mvShort", function(object, ...) standardGeneric("mvShort"))

if(!isGeneric("sizeLong"))
  setGeneric("sizeLong", function(object, ...) standardGeneric("sizeLong"))

if(!isGeneric("sizeShort"))
  setGeneric("sizeShort", function(object, ...) standardGeneric("sizeShort"))

if(!isGeneric("updatePrices"))
  setGeneric("updatePrices", function(object, id, price, ...) standardGeneric("updatePrices"))

if(!isGeneric("matching"))
  setGeneric("matching", function(object, ...) standardGeneric("matching"))

if(!isGeneric("getYahooData"))
  setGeneric("getYahooData", function(object, symbol.var, ...) standardGeneric("getYahooData"))

if(!isGeneric("expandData"))
  setGeneric("expandData", function(object, ...) standardGeneric("expandData"))

if(!isGeneric("expose"))
  setGeneric("expose", function(object, trades, ...) standardGeneric("expose"))


## Class tradelist

## Main methods

if(!isGeneric("calcCandidates"))
  setGeneric("calcCandidates", function(object, orig, target, ...) standardGeneric("calcCandidates"))

if(!isGeneric("calcRanks"))
  setGeneric("calcRanks", function(object, ...) standardGeneric("calcRanks"))

if(!isGeneric("calcChunks"))
  setGeneric("calcChunks", function(object, ...) standardGeneric("calcChunks"))

if(!isGeneric("calcSwaps"))
  setGeneric("calcSwaps", function(object, ...) standardGeneric("calcSwaps"))

if(!isGeneric("calcSwapsActual"))
  setGeneric("calcSwapsActual", function(object, ...) standardGeneric("calcSwapsActual"))

if(!isGeneric("calcChunksActual"))
  setGeneric("calcChunksActual", function(object, ...) standardGeneric("calcChunksActual"))

if(!isGeneric("calcActual"))
  setGeneric("calcActual", function(object, ...) standardGeneric("calcActual"))

if(!isGeneric("calcFinal"))
  setGeneric("calcFinal", function(object, ...) standardGeneric("calcFinal"))

## Utility methods

if(!isGeneric("candidatesCols"))
  setGeneric("candidatesCols", function(object, ...) standardGeneric("candidatesCols"))

if(!isGeneric("ranksCols"))
  setGeneric("ranksCols", function(object, ...) standardGeneric("ranksCols"))

if(!isGeneric("actualCols"))
  setGeneric("actualCols", function(object, ...) standardGeneric("actualCols"))

if(!isGeneric("finalCols"))
  setGeneric("finalCols", function(object, ...) standardGeneric("finalCols"))

if(!isGeneric("chunksCols"))
  setGeneric("chunksCols", function(object, ...) standardGeneric("chunksCols"))

if(!isGeneric("restrictedCols"))
  setGeneric("restrictedCols", function(object, ...) standardGeneric("restrictedCols"))

if(!isGeneric("trimSide"))
  setGeneric("trimSide", function(object, side, value, ...) standardGeneric("trimSide"))

if(!isGeneric("dummyChunks"))
  setGeneric("dummyChunks", function(object, side, num, quality, ...) standardGeneric("dummyChunks"))

if(!isGeneric("securityInfo"))
  setGeneric("securityInfo", function(object, id, ...) standardGeneric("securityInfo"))

if(!isGeneric("mapMarket"))
  setGeneric("mapMarket", function(object, ...) standardGeneric("mapMarket"))
