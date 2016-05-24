## PST
## Sequence objects as produced by TraMineR are S3 classes
setOldClass("stslist")

setClass("cprobd",
	representation(
		context="character"
	),
	contains="matrix"
)

setClass("cprobd.list",
	representation(
		alphabet="character",
		labels="character",
		cpal="character"
	),
	contains="list"
)

setClass("PSTr",
	representation(
		alphabet="character",
		labels="character",
		cpal="character",
		index="matrix",
		counts="matrix",
    		n="matrix",
		prob="matrix",
		path="character",
		order="integer",
		leaf="matrix",
		pruned="matrix"),
	contains="list"
)

setClass("PSTf",
	representation(
		data="stslist",
		cdata="stslist",
		alphabet="character",
		labels="character",
		cpal="character",
		segmented="logical",
		group="factor",
		call="call",
		logLik="numeric"),
	contains="list"
)

setClass("PST.summary",
	representation(
		alphabet="character",
		labels="character",
		cpal="character",
		ns="integer",
		depth="integer",
		nodes="integer",
		leaves="integer",
		freepar="integer")
)




