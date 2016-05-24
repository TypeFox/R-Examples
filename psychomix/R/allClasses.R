## define own class for raschmix and stepRaschmix objects
setClass("FLXMCrasch", contains = "FLXMC")

setClass("raschmix",
         representation(scores = "character",
                        restricted = "logical",
                        deriv = "character",
                        extremeScoreProbs = "numeric",
                        rawScoresData = "table",
                        flx.call = "call",
                        nobs = "numeric",
                        identified.items = "factor"),
         contains = "flexmix")

setClass("stepRaschmix",
         representation(flx.call = "call"),
         contains = "stepFlexmix")

setClass("FLXRoptimRasch",
         representation(scores = "ANY"),
         contains="FLXRoptim")


## define own class for btmix and stepBTmix objects
setClass("FLXMCbt", contains = "FLXMC")

setClass("btmix",
         representation(flx.call = "call",
                        nobs = "numeric",
			labels = "character",
			mscale = "numeric",
                        undecided = "logical",
                        ref = "ANY",
                        type = "character"),
         contains = "flexmix")

setClass("stepBTmix",
         representation(flx.call = "call",
	                labels = "character",
			mscale = "numeric"),
         contains = "stepFlexmix")


