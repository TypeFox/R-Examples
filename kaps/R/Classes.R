## S4 classes for kaps package
## Soo-Heang Eo, 2013-08-28
setOldClass("Formula")
setOldClass("kapsOptions")
setOldClass("Surv")

##########################################
### A class for adaptive partitioning
setClass(Class = "kapsOptions",
	representation = representation(
		V 			= 	"numeric",
		pre.pt 		=	"list",
		scope		=	"list",
		lower.limit	=	"numeric",
		upper.limit	=	"numeric",
		N.perm		=	"numeric",
		N.boot		=	"numeric",
		alpha 		=	"numeric",
		rho			=	"numeric",
		fold		=	"logical",
		ncl			=	"integer",
		splits 		=	"character",
		sel 		=  	"character",
		shortcut	=	"logical",
		correct		=	"character",
		p.adjust.methods	= 	"character"
	)
)

setClass(Class = "kaps", 
	representation = representation(
		call				= "language",
		formula				= "Formula",
		data				= "data.frame",
		groupID				= "vector",  
		index				= "integer",
		X					= "numeric",
		Z 		 			= "numeric",
		pair				= "numeric",
		split.var		 	= "character",
		split.pt		 	= "numeric",
		mindat				= "numeric",
		test.stat			= "matrix",
		over.stat.sample	= "matrix",
		pair.stat.sample	= "matrix",
		groups				= "vector",
		results				= "list",
		Options				= "kapsOptions"
	)
)

#####################################
### A class for data by recursive binary splits
#setClass(Class = "dataset",
#	representation = representation(
#		Y = "Surv",
#		X = "data.frame",
#		Z = "data.frame",
#		resid = "vector",
#		resid.sign = "vector"
#	)
#)

# END by Soo-Heang Eo