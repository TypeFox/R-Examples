###
###
###
###
###
###


#-----
#Notes:
#
#	MAPCOUNTS : list of numbers of mapped keys per table
#
#-----


###
###		org.<species>.eg database 
###
###
###



#
#	Available org.SPEC.eg.db databases in BioconductoR at end of June 2012 :
#		(excerpt of http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData )
#
#	Package:			Maintainer:			Title:
#	org.Ag.eg.db		Biocore Data Team	Genome wide annotation for Anopheles
#	org.At.tair.db		Biocore Data Team	Genome wide annotation for Arabidopsis
#	org.Bt.eg.db		Biocore Data Team	Genome wide annotation for Bovine
#	org.Ce.eg.db		Biocore Data Team	Genome wide annotation for Worm
#	org.Cf.eg.db		Biocore Data Team	Genome wide annotation for Canine
#	org.Dm.eg.db		Biocore Data Team	Genome wide annotation for Fly
#	org.Dr.eg.db		Biocore Data Team	Genome wide annotation for Zebrafish
#	org.EcK12.eg.db		Biocore Data Team	Genome wide annotation for E coli strain K12
#	org.EcSakai.eg.db	Biocore Data Team	Genome wide annotation for E coli strain Sakai
#	org.Gg.eg.db		Biocore Data Team	Genome wide annotation for Chicken
#	org.Hs.eg.db		Biocore Data Team	Genome wide annotation for Human
#	org.Hs.ipi.db		Hong Li				A data package containing annotation data for org.Hs.ipi.db
#	org.Mm.eg.db		Biocore Data Team	Genome wide annotation for Mouse
#	org.Mmu.eg.db		Biocore Data Team	Genome wide annotation for Rhesus
#	org.Pf.plasmo.db	Biocore Data Team	Genome wide annotation for Malaria
#	org.Pt.eg.db		Biocore Data Team	Genome wide annotation for Chimp
#	org.Rn.eg.db		Biocore Data Team	Genome wide annotation for Rat
#	org.Sc.sgd.db		Biocore Data Team	Genome wide annotation for Yeast
#	org.Sco.eg.db		Roxane Legaie		Genome wide annotation for Streptomyces coelicolor
#	org.Ss.eg.db		Biocore Data Team	Genome wide annotation for Pig
#	org.Tgondii.eg.db	Olivier Lucas		Genome wide annotation for Toxoplasma gondii
#	org.Xl.eg.db		Biocore Data Team	Genome wide annotation for Xenopus
#


#
#	Get <organism> for the org.<organism>.eg.eb-style database names, given one of the
#		commonly used and plausible names for the wanted organism (see list above)
#
whop.eg.abbrevForOrganism <- function( organismname )
{
	organismname <- gsub( "\\W|_","", organismname )
	org <- NULL
	try( org <- .WHOPGENenv$eg[["organisms_list"]][[ tolower(organismname) ]] , silent=T)
	if( is.null( org ) )
	{
		cat( paste("Unrecognised organism:",organismname,"\n") )
	}
	return( org )
}

#
#	Is org.<organism>.eg.db loaded ?
#
whop.eg.orgdb_loaded <- function( organismname )
{
	paste(sep="","package:org.",whop.eg.abbrevForOrganism( organismname ),".eg.db") %in% search()
}

#
#
#
whop.eg.installdb <- function( organismname )
{
	orglibnam <- NULL ; rm(orglibnam)
	biocLite <- NULL ; rm(biocLite)
	abbrev <- whop.eg.abbrevForOrganism( organismname )
	if( ! is.null( abbrev ) )
	{		
		source("http://bioconductor.org/biocLite.R")
		orglibnam <- paste(sep="","org.",abbrev,".eg.db" )
		biocLite( orglibnam )

		if( library( orglibnam, character.only=TRUE ) )	return(TRUE)
			
		print( .WHOPGENenv$pp("Could not load library '",orglibnam,"'!") )
	}

	#
	return(FALSE)
}


##########
#
#	Load an org.<organism>.eg.db and, if needed, install it first from BioconductoR
#
##########
whop.eg.load_orgdb <- function( organismname, install.if.missing=F )
{
	# find abbrev used in the bioconductor package names
	orgdblibnam <- NULL
	abbrev <- whop.eg.abbrevForOrganism( organismname )
	
	# if the abbreviation was found,...
	if( ! is.null( abbrev ) )
	{		
		# try to load the database
		if( ! whop.eg.orgdb_loaded( abbrev ) )
		{
			orgdblibnam <- paste(sep="","org.",abbrev,".eg.db" )
			print( .WHOPGENenv$pp("Loading ", orgdblibnam,"...") )
			if( ! library( orgdblibnam, character.only=TRUE ) )
			{
				
				if( ! install.if.missing )
				{
					print( .WHOPGENenv$pp("Could not load library '",orgdblibnam,"'!") )
					return( FALSE )
				}
				
				return( whop.eg.installdb( organismname ) )

			}#..if( could not load org.eg.db package )

		}#...if( org.eg.db for organism not loaded )
		
		return( TRUE )

	}#...if ( could translate organismname to abbreviation )
	
	#
	return( FALSE );
}


##########
#
#	Select an organism
#
##########
whop.eg.selectOrganism <- function( organismname, dontload=FALSE, install.if.missing=F )
{
	# find abbrev used in the bioconductor package names
	abbrev <- whop.eg.abbrevForOrganism( organismname )
	
	# if the abbreviation was found,...
	if( ! is.null( abbrev ) )
	{
		# set the organism abbreviation currently used
		.WHOPGENenv$eg[["organism_name"]] <- abbrev
		res = TRUE
		
		# also try to load the database
		if( (!dontload) && (! whop.eg.orgdb_loaded( abbrev )) )
		{
			res = whop.eg.load_orgdb( abbrev , install.if.missing )
		}
		return( res )
	}
	return(FALSE)
}


##########
#
#	Given a common gene alias (like LCT, MTLA1) into the numeric identifier used in eg
#
##########

	#	
	#
whop.eg.eg_lookupAll <- function( id, subdbname, db=.WHOPGENenv$eg[["databasename"]] ) mget( id , get( paste(sep="","org.",.WHOPGENenv$eg[["organism_name"]],".eg",subdbname) ) )

	#	map a single Entrez gene id to some information in a given table
	#
whop.eg.eg_lookupSingle <- function( id, subdbname, db=.WHOPGENenv$eg[["databasename"]] ) mget( id , get( paste(sep="","org.",.WHOPGENenv$eg[["organism_name"]],".eg",subdbname) ) )[[1]]

	#	map each of a vector of Entrez gene ids to some information in a given table
	#
whop.eg.eg_lookup <- function( ids, subdbname, db=.WHOPGENenv$eg[["databasename"]] )
{
	res<-c();
	x <- 0
	for( x in ids )
	{
		res <- c( res, mget( x , get( paste(sep="","org.",.WHOPGENenv$eg[["organism_name"]],".eg",subdbname) ) ) )
	}
	return (res);
}
	#	map each of a vector of Entrez gene ids to some information in a given table
	#
whop.eg.eg_RevLookup <- function( ids, subdbname, db=.WHOPGENenv$eg[["databasename"]] )
{
	res<-c();
	y <- 0
	for( y in ids )
	{
		res <- c( res, mget( y , 
							AnnotationDbi::revmap( 
								get( paste(sep="","org.",.WHOPGENenv$eg[["organism_name"]],".eg",subdbname) ) 
								) 
							) 
				)
	}
	return (res);
}

#################################################################
#
#	Other Identifiers -> EG
#
#################################################################

	#	From
	#

whop.eg.fromAlias <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "ALIAS2EG", db )
whop.eg.fromAccnum <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "ACCNUM2EG", db )
whop.eg.fromEnsembl <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "ENSEMBL2EG", db )
whop.eg.fromEnsemblProt <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "ENSEMBLPROT2EG", db )
whop.eg.fromEnsemblTrans <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "ENSEMBLTRANS2EG", db )
whop.eg.fromEnzyme <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "ENZYME2EG", db )
whop.eg.toEnzyme <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "ENZYME", db )
whop.eg.fromGO <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "GO2EG", db )
whop.eg.fromGO2AllEgs <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "GO2ALLEGS", db )
whop.eg.fromOmim <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "OMIM2EG", db )
whop.eg.fromPath <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "PATH2EG", db )
whop.eg.fromPmid <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "PMID2EG", db )
whop.eg.fromRefseq <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "REFSEQ2EG", db )
whop.eg.fromUnigene <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "UNIGENE2EG", db )
whop.eg.fromUniprot <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "UNIPROT2EG", db )

	#	To
	#

whop.eg.toAlias <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "SYMBOL2EG", db )
whop.eg.toAccnum <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "ACCNUM2EG", db )
whop.eg.toEnsembl <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "ENSEMBL2EG", db )
whop.eg.toEnsemblProt <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "ENSEMBLPROT2EG", db )
whop.eg.toEnsemblTrans <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "ENSEMBLTRANS2EG", db )
#whop.eg.toEnzyme <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "ENZYME2EG", db )
whop.eg.toGO <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "GO2EG", db )
#whop.eg.fromGO2AllEgs <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "GO2ALLEGS", db )
whop.eg.toOmim <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "OMIM2EG", db )
whop.eg.toPath <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "PATH2EG", db )
whop.eg.toPmid <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "PMID2EG", db )
whop.eg.toRefseq <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "REFSEQ2EG", db )
whop.eg.toUnigene <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "UNIGENE2EG", db )
whop.eg.toUniprot <- function( id , db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_RevLookup( id, "UNIPROT2EG", db )

#PFAM?
#	> org.Hs.egPFAM$`1`
#	IPI00022895 IPI00644018 IPI00902880 
#         NA          NA          NA 
#	> names(org.Hs.egPFAM$`1`)
#	[1] "IPI00022895" "IPI00644018" "IPI00902880"
#
#PROSITE?
#
#	same pattern as with PFAM
#

#################################################################
#
#	EG -> Other Identifiers
#
#################################################################

#org.Hs.egACCNUM           org.Hs.egENSEMBLTRANS     org.Hs.egMAP              org.Hs.egPFAM             org.Hs.egUCSCKG
#org.Hs.egACCNUM2EG        org.Hs.eg_dbfile          org.Hs.egENSEMBLTRANS2EG  org.Hs.egMAP2EG           org.Hs.egPMID             org.Hs.egUNIGENE
#org.Hs.egALIAS2EG         org.Hs.eg_dbInfo          org.Hs.egENZYME           org.Hs.egMAPCOUNTS        org.Hs.egPMID2EG          org.Hs.egUNIGENE2EG
#org.Hs.egCHR              org.Hs.eg_dbschema        org.Hs.egENZYME2EG        org.Hs.egOMIM             org.Hs.egPROSITE          org.Hs.egUNIPROT
#org.Hs.egCHRLENGTHS       org.Hs.egENSEMBL          org.Hs.egGENENAME         org.Hs.egOMIM2EG          org.Hs.egREFSEQ           org.Hs.eg.db::
#org.Hs.egCHRLOC           org.Hs.egENSEMBL2EG       org.Hs.egGO               org.Hs.egORGANISM         org.Hs.egREFSEQ2EG        
#org.Hs.egCHRLOCEND        org.Hs.egENSEMBLPROT      org.Hs.egGO2ALLEGS        org.Hs.egPATH             org.Hs.egSYMBOL

whop.eg.toAccnum <- function( id , db=.WHOPGENenv$eg[["databasename"]]  ) whop.eg.eg_lookup( as.character(id) ,"ACCNUM",db)

#################################################################
#
#	Other handy lookup functions
#
#################################################################

	##########
	#	gene descriptive name for egID
	#
whop.eg.genename <- function( id, db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "GENENAME", db )

	##########
	#	enzyme nomenclature identifiers for egID
	#
whop.eg.enzyme <- function( id, db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "ENZYME", db )

	##########
	#	chromosomal region for egID
	#
whop.eg.region <- function( id, db=.WHOPGENenv$eg[["databasename"]] ) c( whop.eg.eg_lookup( id, "CHRLOC", db ), whop.eg.eg_lookup( id, "CHRLOCEND", db ) )

	##########
	#	chromosome on which gene entrez <id> lies
	#
whop.eg.chromosome <- function( id, db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "CHR", db )

	##########
	#	GO identifiers relating to gene with entrez <id>
	#
whop.eg.goIds <- function( id, db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "GO", db )

#MAP2EG	(cytogenetic location) e.g. 19q13.4

#################################################################
#
#	Miscelleaneous functions
#
#################################################################


whop.eg.Organism <- function() get( paste(sep="","org.",.WHOPGENenv$eg[["organism_name"]],".egORGANISM") )

#################################################################
#
#	KEGG related (pathways)
#
#################################################################


whop.eg.keggpathways <- function( id, db=.WHOPGENenv$eg[["databasename"]] ) whop.eg.eg_lookup( id, "PATH", db )



	
	#	print URL pointing to the KEGG pathway webpage for a given KEGG pathway ID
	#
	#NOTE: tested to work correctly with a vector of pathwayids
whop.kegg.pathway_url <- function( pathwayids ) print( .WHOPGENenv$pp( "http://www.genome.jp/dbget-bin/www_bget?pathway+",.WHOPGENenv$eg[["organism_prefixed"]][ tolower(.WHOPGENenv$eg[["organism_name"]]) ],pathwayids) )
