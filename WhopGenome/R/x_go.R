
##TODO : 
##	Using and citing a mapping
##
##	If you have used a mapping in a publication or presentation, please ensure
##	that you cite both the GO project and the source of the mapping (detailed below).
##	See the GO citation guide < http://www.geneontology.org/GO.cite.shtml > for citing the GO project.
##
##

#
#
#


#
#
#
whop.go.load <- function( filename=NA )
{
	#
	if( is.na( filename ) )
	{
		filename <- .WHOPGENenv$go[["default.filename"]]
	}
	
	#
	.WHOPGENenv$go[["gotable"]] <- read.delim( filename ,comment.char=c("!"),col.names=c("prim.term","sec.terms","name","type","isObs"))
	
	#TODO : return whether any rows have been read
}

#
#
#
whop.go.match <- function( tofind )
{
	if( is.null(.WHOPGENenv$go[["gotable"]]) )
		stop("whop.go.match : this call relies on a GO table loaded from file. Please use whop.go.load!")
	matching = grep( tofind , unlist( .WHOPGENenv$go[["gotable"]]["name"] ))
	.WHOPGENenv$go[["gotable"]][matching,]
}

#
#	find any GO terms where the ID is similar to <idmatch>
#
whop.go.goid_like <- function( idmatch )
{
	if( is.null(.WHOPGENenv$go[["gotable"]]) )
		stop("whop.go.match : this call relies on a GO table loaded from file. Please use whop.go.load!")
	#
	.WHOPGENenv$go[["gotable"]][ grep( idmatch , unlist(.WHOPGENenv$go[["gotable"]]["prim.term"]) ), ]
}

#############################################################################################
#############################################################################################
#############################################################################################


#
#	Establish a connection to the AmiGO database server
#
whop.go.connect <- function( althost=NA, altport=NA, altuser=NA, altpass=NA, altdb=NA, altdbdrivername=NA, dbdrvpkgnam=NA )
{
	if( ! is.na( dbdrvpkgnam ) )
		library( dbdrvpkgnam ,character.only=TRUE )
		
	#TODO : sever any existing connection to an AmiGO db server
	#

	#	set default values for GO database
	#
	.WHOPGENenv$go$defaultAmiGO()
	
	#	check for overrides to default values for EBI AmiGO servers
	#
	if( !is.na( althost ) )	.WHOPGENenv$go[["db_host"]] <- althost
	if( !is.na( altuser ) )	.WHOPGENenv$go[["db_user"]] <- altuser
	if( !is.na( altpass ) )	.WHOPGENenv$go[["db_pass"]] <-	altpass
	if( !is.na( altport ) )	.WHOPGENenv$go[["db_port"]] <-	altport
	if( !is.na( altdb ) )	.WHOPGENenv$go[["db_godb"]] <- altdb
	if( !is.na( altdbdrivername ) )	.WHOPGENenv$go[["db_drivername"]] <- altdbdrivername
	
	#
	print( "Currently configured AmiGO DB server:" )
	print( paste(sep="", "sql://",.WHOPGENenv$go[["db_user"]],":",.WHOPGENenv$go[["db_pass"]],"@",.WHOPGENenv$go[["db_host"]],":",.WHOPGENenv$go[["db_port"]],"/",.WHOPGENenv$go[["db_godb"]]) )
	
	#	instantiate DBMS driver
	#
	.WHOPGENenv$go[["db_driver"]] <- DBI::dbDriver(.WHOPGENenv$go[["db_drivername"]])
	if( is.null( .WHOPGENenv$go[["db_driver"]] ) )
		stop( paste(sep="","Could not instantiate DBMS driver ",.WHOPGENenv$go[["db_drivername"]],"!") );
	
	#	connect to database
	#
	.WHOPGENenv$go[["db_connection"]] <- DBI::dbConnect( 
		drv=.WHOPGENenv$go[["db_driver"]],
		host=.WHOPGENenv$go[["db_host"]],
		user=.WHOPGENenv$go[["db_user"]],
		pass=.WHOPGENenv$go[["db_pass"]],
		port=.WHOPGENenv$go[["db_port"]],
		dbname=.WHOPGENenv$go[["db_godb"]] 
	 )

}#...whop.go.connect()

#
#	find GO terms whose name is similar to <tomatch>
#
whop.go.terms_match <- function( tomatch )
{
	#TODO check single string in tomatch
	#dbGetQuery(conn=.WHOPGENenv$go[["db_connection"]],statement=paste(sep="",
	.WHOPGENenv$go$goquery(
		"select * from term where name like '%",tomatch,"%'"
		)
	#	) )
}

#
#	for the GO-term name <tomatch> find other GO terms synonymous with
#
whop.go.term_synonyms <- function( tomatch )
{
	#TODO check single string in tomatch
	.WHOPGENenv$go$goquery( "SELECT * FROM term INNER JOIN term_synonym ON (term.id=term_synonym.term_id) WHERE term_synonym LIKE '%",tomatch,"%';")
}

#
#	find ancestors of GO terms LIKE <tomatch>
#
whop.go.term_ancestors_similar <- function( tomatch )
{
	#TODO check single string in tomatch
	.WHOPGENenv$go$goquery(
			"SELECT ancestor.*,graph_path.distance,graph_path.term1_id AS ancestor_id",
			" FROM term child, graph_path, term ancestor",
			" WHERE child.id=graph_path.term2_id AND ancestor.id=graph_path.term1_id",
			" AND child.name LIKE '%",tomatch,"%';"
		)
}

#
#	find ancestors of GO term <tomatch>
#
whop.go.term_ancestors <- function( tomatch )
{
	#TODO check single string in tomatch
	.WHOPGENenv$go$goquery(
			"SELECT ancestor.*,graph_path.distance,graph_path.term1_id AS ancestor_id",
			" FROM term child, graph_path, term ancestor",
			" WHERE child.id=graph_path.term2_id AND ancestor.id=graph_path.term1_id",
			" AND child.name='",tomatch,"';"
		)
}

#
#	find children of GO term named <tomatch>
#
whop.go.term_children <- function( tomatch )
{
	#TODO check single string in tomatch
	.WHOPGENenv$go$goquery(
			"SELECT DISTINCT descendant.acc, descendant.name, descendant.term_type",
			" FROM term",
			" INNER JOIN graph_path ON (term.id=graph_path.term1_id)",
			" INNER JOIN term AS descendant ON (descendant.id=graph_path.term2_id)",
			" WHERE term.name='",tomatch,"' AND distance <> 0 ;"
		)
}

#
#	find out if GO term LIKE <tomatch> is obsolete
#
whop.go.is_obsolete_byname <- function( tomatch )
{
	#TODO check single string in tomatch
	.WHOPGENenv$go$goquery(
			"SELECT * FROM term WHERE is_obsolete=1 AND name LIKE '",tomatch,"';"
		)
}

#
#	find out if GO term LIKE <idmatch> is obsolete
#
whop.go.is_obsolete_byid <- function( idmatch )
{
	#TODO check single string in tomatch
	.WHOPGENenv$go$goquery(
			"SELECT * FROM term WHERE is_obsolete=1 AND acc LIKE '",idmatch,"';"
		)
}

#
#	returns all genes annotated to certain GO term
#
whop.go.all_genes_for_term <- function( tomatch )
{
	#TODO check single string in tomatch
	.WHOPGENenv$go$goquery(
			"SELECT association.is_not, term.name, term.acc, term.term_type, gene_product.symbol AS gp_symbol,",
			" gene_product.symbol AS gp_full_name, dbxref.xref_dbname AS gp_dbname, dbxref.xref_key AS gp_acc,",
			" species.genus, species.species, species.common_name, species.ncbi_taxa_id, association.assocdate,",
			" db.name AS assigned_by, db.fullname",
			" FROM term INNER JOIN association ON (term.id=association.term_id)",
			" INNER JOIN gene_product ON (association.gene_product_id=gene_product.id)",
			" INNER JOIN species ON (gene_product.species_id=species.id)",
			" INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id)",
			" INNER JOIN db ON (association.source_db_id=db.id)",
			" WHERE term.name = '",tomatch,"';"
		)
 }




