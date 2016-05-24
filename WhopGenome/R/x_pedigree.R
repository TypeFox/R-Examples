##
##
##
##	manage populations in VCF files
##
##
##
##
##
##

#
#	NOTE:	Gender	1=male	2=female
#

#"/home/ami/Downloads/1000Genomes/G1K_samples_20111108.ped.txt"
#

#
#
whop.ped.load	<-	function( filename )
{
	pops = read.delim(filename,header=T,sep="\t")
	names(pop1)[ names(pop1) == "Gender" ] <- "Sex"
	pops[ pops[,"Sex"] == "M" , "Sex" ] <- 1
	pops[ pops[,"Sex"] == "F" , "Sex" ] <- 2
	return( pops )
}



#
#
whop.ped.save	<-	function( p , filename )
{
	write.table( x=p , file=filename )
}


#
#	returns the Individual.IDs from the given table
#
whop.ped.entriesOf	<-	function( p , invids )
{
	p[ p[,"Individual.ID"] %in% invids , ]
}


#
#	returns the Individual.IDs from the given table
#
whop.ped.names	<-	function( p )
{
	as.character( unlist( p[,"Individual.ID"] ) )
}


#
#	returns the fathers from the given table
#
whop.ped.fathers	<-	function( p  )
{
	p[,"Paternal.ID"]
}


#
#	returns the mothers from the given table
#
whop.ped.mothers	<-	function( p  )
{
	p[,"Maternal.ID"]
}


#
#	Returns all entries with a population ID that it also found in <popids>
#
whop.ped.fromPop	<-	function( p , popids )
{
	p[ which( p[,"Population"] %in% popids ) , ]
}


#
#	Returns all males in the population-matrix
#
whop.ped.males	<-	function( p )
{
	p[ p[,"Sex"] == 1 , ]
}


#
#	Returns all females in the population-matrix
#
whop.ped.females	<-	function( p )
{
	p[ p[,"Sex"] == 2 , ]
}


#
#
#
whop.ped.parentsOf	<-	function( p , invids )
{
	children = whop.ped.entriesOf( p , invids )
	parentids = as.character(unlist( children[, c("Paternal.ID","Maternal.ID") ] ))
	whop.ped.entriesOf( p , parentids )
}


#
#
whop.ped.familyOf	<-	function( p , lis )
{
	#	get siblings' names
	#
	s1 = whop.ped.siblingsOf( p, lis )
	sibs = whop.ped.names( s1 )
	
	#	get parents' names
	#
	p1 = whop.ped.parentsOf( p, lis )
	parents = whop.ped.names( p1 )
	
	#	get entries with both parents' and siblings' names
	#
	res = whop.ped.entriesOf( p, c(sibs,parents) )
	return( res )
}


#
#	Returns the entries of sons and daughters of the individuals in <lis>
#
whop.ped.siblingsOf	<-	function( p , lis )
{
	res <- c()
	z <- 0
	for( z in lis )
	{
		sib1 = as.character( p[ p[,"Paternal.ID"] %in% z , "Individual.ID" ] )
		sib2 = as.character( p[ p[,"Maternal.ID"] %in% z , "Individual.ID" ] )
		sib1 = sib1[ !is.na(sib1) ]
		sib2 = sib2[ !is.na(sib2) ]
		res <- c( res, sib1, sib2 )
	}
	whop.ped.entriesOf( p , res )
}




#
#	Returns a vector of Individual.IDs of sons and daughters 
#
whop.ped.daughtersOf	<-	function( p , lis )
{
	whop.ped.females( whop.ped.siblingsOf( p, lis ) )
}


#
#	Returns a vector of Individual.IDs of sons and daughters 
#
whop.ped.sonsOf	<-	function( p , lis )
{
	whop.ped.males( whop.ped.siblingsOf( p, lis ) )
}





