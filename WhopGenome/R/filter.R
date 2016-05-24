###
###
###		Data filtering setup and management
###
###

#TODO : removefilter (move following filters one up)
#

##
##
.vcf.column.name.to.idx <- function( cname )
{
	cname = toupper( cname )
	return( as.integer( switch(		EXPR=cname,
							CHROM = 0,
							POS = 1,
							ID = 2,
							REF = 3,
							ALT = 4,
							QUAL = 5,
							FILTER = 6,
							INFO = 7,
							FORMAT = 8,
							-1
						)#switch
					)#as.integer
		)#return
}


##
##
.vcf.comparison.name.to.idx <- function( comparison )
{
	comparison = toupper( comparison )
	return( as.integer( switch( EXPR = comparison,
								DOES_EXIST = 0,
									EXISTS = 0,
									HASKEY = 0,
									"HAS KEY" = 0,
								INT_CMP=1,
									"INT=" =1,
									"INT=="=1,
								INT_CMP_OO=2,
									"INT()"=2,
								INT_CMP_OC=3,
									"INT(]"=3,
								INT_CMP_CO=4,
									"INT[)"=4,
								INT_CMP_CC=5,
									"INT[]"=5,
								FLT_CMP=6,
									"FLT=" =6,
									"FLT=="=6,
								FLT_CMP_OO=7,
									"FLT()"=7,
								FLT_CMP_OC=8,
									"FLT(]"=8,
								FLT_CMP_CO=9,
									"FLT[)"=9,
								FLT_CMP_CC=10,
									"FLT[]"=10,
								-1
							)#switch
						)#as.integer
			)#return
}


##
##
.vcf.action.name.to.idx <- function( action )
{
	action = toupper( action )
	return( as.integer( switch( EXPR = action,
								NOP=0,
								SKIP=1,
									DROP=1,			# DROP as alternative name for SKIP
								KEEP=2,
								SKIP_NOT=3,
									DROP_NOT=3,
									SKIP_IF_NOT=3,
									DROP_IF_NOT=3,
								KEEP_NOT=4,
								NOP_DISABLED	=	0x80,
								SKIP_DISABLED	=	0x81,	#DISABLED : skip if false, dont test further
									DROP_DISABLED	=	0x81,
								KEEP_DISABLED	=	0x82,	#DISABLED : keep if true, dont test further
								SKIP_NOT_DISABLED=	0x83,	#DISABLED : skip if true, dont test further
									DROP_NOT_DISABLED=	0x83,
								KEEP_NOT_DISABLED=	0x84,	#DISABLED : keep if false, dont test further
								-1
							)#switch
						)#as.integer
			)#return
}




##
##
##	Add a filter
##
##
vcf_addfilter <- function( vcf, columnnam, fieldnam, cmptype, cmpvalue1, cmpvalue2=0, action )
{
	
	#0. turn column name into column index
	#
	columnidx = .vcf.column.name.to.idx( columnnam )
	if( columnidx < 0 )
	{
		stop("(!!) Column name not valid!");
	}
	
	
	#1. add fieldnam to ruleset of vcf
	#	error if ruleset full; name might already be in there
	fieldnam = toupper(fieldnam)
	
	#2. set cmptype
	#
	cmptype = .vcf.comparison.name.to.idx( cmptype )
	if( cmptype < 0 )	stop("Illegal comparison type!");

	#3. set cmp parameters
	#

		#
		#
	if( (cmptype >= .vcf.comparison.name.to.idx("INT_CMP")) && (cmptype <= .vcf.comparison.name.to.idx("INT_CMP_CC")) )
	{
		if( !is.numeric(cmpvalue1) || ( (cmptype > .vcf.comparison.name.to.idx("INT_CMP")) && !is.numeric( cmpvalue2 )) )
			stop("cmpvalue1 and cmpvalue2 must be numeric for use in integer comparisons!");
		cmp1 = as.integer( cmpvalue1 )
		cmp2 = as.integer( cmpvalue2 )
	}
	# for all FLT_CMP types
	else if( (cmptype >= .vcf.comparison.name.to.idx("FLT_CMP")) && (cmptype <= .vcf.comparison.name.to.idx("FLT_CMP_CC")) )
	{
		# check whether ref1 (and for range-comparisons, also ref2) are invalid
		#
		if( !is.numeric(cmpvalue1) || ( (cmptype > .vcf.comparison.name.to.idx("FLT_CMP")) && !is.numeric( cmpvalue2 )) )
			stop("cmpvalue1 and cmpvalue2 must be numeric for use in floating-point comparisons!");
		cmp1 = as.double( cmpvalue1 )
		cmp2 = as.double( cmpvalue2 )
	}
	else	# for DOES_EXIST
	{
		cmp1 <- 0
		cmp2 <- 0
	}
	
	#4. set action on cmp-result
	#
	cmpaction = .vcf.action.name.to.idx( action )
	if( cmpaction == -1 )
		stop("Unknown action selected!");
	
	#	add rule parameters
	#
	return(
		.Call("VCF_addFilter",vcf,columnidx, fieldnam,cmptype,cmpaction,cmp1,cmp2,PACKAGE="WhopGenome")
		)

	#
}#...addFilter

##
##	Remove all filters
##
##
##
vcf_clearfilters <- function( vcffh )	.Call("VCF_clearFilters" , vcffh , PACKAGE="WhopGenome")

##
##	Print an understandable description of the current filering rules
##
##
##
vcf_describefilters <- function( vcffh )	.Call("VCF_describeFilterConfig" , vcffh , PACKAGE="WhopGenome")


##
##	
##
##
##
vcf_rule.setaction <- function( vcffh, ruleidx, action )
{
	ruleidx <- as.integer( ruleidx )
	action <- .vcf.action.name.to.idx( action )
	if( action < 0 ) stop("unknown action in vcf_setruleaction!");
	.Call("VCF_setRuleAction" , vcffh , ruleidx , action , PACKAGE="WhopGenome")
}


##
##	Enable a rule
##
##
##
vcf_rule.enable <- function( vcffh, ruleidx )
{
	ruleidx <- as.integer( ruleidx )
	.Call("VCF_setRuleEnabled" , vcffh , ruleidx ,  PACKAGE="WhopGenome")
}


##
##	Disable a rule
##
##
##
vcf_rule.disable <- function( vcffh, ruleidx )
{
	ruleidx <- as.integer( ruleidx )
	.Call("VCF_setRuleDisabled" , vcffh , ruleidx ,  PACKAGE="WhopGenome")
}


##
##	Set which column should be checked
##
##
##
vcf_rule.setcolumn <- function( vcffh, ruleidx, column )
{
	column <- .vcf.column.name.to.idx( column )
	ruleidx <- as.integer( ruleidx )
	.Call("VCF_setRuleColumn" , vcffh , ruleidx , column , PACKAGE="WhopGenome")
}


##
##	Set the reference values for the comparison operation of the given rule
##
##
##
vcf_rule.setrefvalues <- function( vcffh, ruleidx, ref1, ref2 )
{
	ruleidx <- as.integer( ruleidx )
	stopifnot( length(ref1) == 1 )
	stopifnot( length(ref2) == 1 )
	stopifnot( mode(ref1) == "numeric" )
	stopifnot( mode(ref2) == "numeric" )
	.Call("VCF_setRuleRefValues" , vcffh , ruleidx , ref1 , ref2 , PACKAGE="WhopGenome")
}

#	TODO:
#VCF_getRuleRefValues
#vcf_rule.setleftref 
#vcf_rule.setrightref


##
##	Set which field (or key) in a column should be checked
##
##
##
vcf_rule.setfield <- function( vcffh, ruleidx, field )
{
	ruleidx <- as.integer( ruleidx )
	stopifnot( mode(field) == "character" )
	stopifnot( length(field) == 1 )
	.Call("VCF_setRuleField" , vcffh , ruleidx , field , PACKAGE="WhopGenome")
}


##
##	Set the kind of comparison to perform
##
##	TODO : make validity checks so that the comparison value types do not suddenly become invalid!
##
vcf_rule.setcomparison <- function( vcffh, ruleidx, cmpop )
{
	ruleidx <- as.integer( ruleidx )
	cmpop <- .vcf.comparison.name.to.idx( cmpop )
	stopifnot( cmpop >= 0 )
	.Call("VCF_setRuleCmpOp" , vcffh , ruleidx , cmpop , PACKAGE="WhopGenome")
}

#
#TODO:
#
#
#	removefilter: vcffh, filterentrynum
#
#	replacefilter: = addfilter + filternum-parameter
#
#	getfilters := get a structure describing all filters
#
#
#
#
