#  ZZZ_whopgenome.R
#
#  Copyright (C) 2012-2013 Ulrich Wittelsbuerger, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: ulrich.wittelsbuerger@uni-duesseldorf.de
#
#  This file is part of WhopGenome.
#
#  WhopGenome contains slightly modified source code of Zlib and
#	Tabix. Modified parts of their sources are marked and commented.
#
#  WhopGenome is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  WhopGenome is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with WhopGenome.  If not, see <http://www.gnu.org/licenses/>.


.WHOPGENenv <- new.env()

.onLoad <- function(lib, pkg) {

    # -------------------------------------------------------------- #

	#
	#
	#		Gene Ontology - via AmiGO or file database
	#
	#
    .WHOPGENenv$go <- list(
		default.filename = "http://www.geneontology.org/doc/GO.terms_alt_ids",
		gotable = NULL,
		limit = 100,
		
		db_drivername = "MySQL",
		db_host <- "mysql.ebi.ac.uk",
		db_user <- "go_select",
		db_pass <-	"amigo",
		db_port <-	4085,
		db_godb <- "go_latest",
		db_connection = NULL,
		db_driver = NULL,
		
		##
		##
		defaultAmiGO = function(){
			.WHOPGENenv$go[["db_drivername"]] = "MySQL"
			.WHOPGENenv$go[["db_host"]] <- "mysql.ebi.ac.uk"
			.WHOPGENenv$go[["db_user"]] <- "go_select"
			.WHOPGENenv$go[["db_pass"]] <-	"amigo"
			.WHOPGENenv$go[["db_port"]] <-	4085
			.WHOPGENenv$go[["db_godb"]] <- "go_latest"
			db_connection = NULL
			db_driver = NULL
		},
		
		##
		##
		getConnection = function(){
			if( is.null( .WHOPGENenv$go[["db_connection"]] ) )
			{
				.WHOPGENenv$go$connect()
			}
			return( .WHOPGENenv$go[["db_connection"]] )
		},

		##
		##	connect to AmiGO SQL server
		##
		connect = function()
		{
			.WHOPGENenv$go[["db_connection"]] <- DBI::dbConnect( DBI::dbDriver( .WHOPGENenv$go[["db_drivername"]] ), user=.WHOPGENenv$go[["db_user"]], pass=.WHOPGENenv$go[["db_pass"]], host=.WHOPGENenv$go[["db_host"]], port=.WHOPGENenv$go[["db_port"]] )
		},

		##
		##	disconnect from AmiGO SQL server
		##
		disconnect = function()
		{
			DBI::dbDisconnect( conn=.WHOPGENenv$go[["db_connection"]] )
			.WHOPGENenv$go[["db_connection"]] = NULL
		},

		##
		##	disconnect from AmiGO SQL server
		##
		goquery = function(...)
		{
			DBI::dbGetQuery( conn=.WHOPGENenv$go$getConnection(),
				statement=.WHOPGENenv$pp(...)
				)
		}
		
		
    )#GO

    # -------------------------------------------------------------- #

	###
	###
	###		Bioconductor ORG.eg.Db
	###
	###
    .WHOPGENenv$eg <- list(
		databasename = NULL,
		organism_name <- "Hs",		#.WHOPGENenv$eg[["organism_name"]]
		organisms_list <- list( 
				homosapiens="Hs", human = "Hs", hs="Hs", hsap="Hs", hsapiens = "Hs",
				mouse="Mm", mm="Mm", mmus="Mm", mmusculus="Mm"
				# TODO : extend
				#	Ag		Anopheles gambiae		Anopheles
				#	At		Arabidopsis thaliana	Arabidopsis
				#	Bt		Bos tauris				Bovine
				#	Ce		Caenorhabditis elegans	Worm
				#	Cf		Canis f??				Canine
				#	Dm		Drosophila melanogaster	Fly
				#	Dr		Danio rerio				Zebrafish
				#	EcK12	Escherichia coli K12
				#	EcSakai	Escherichia coli Sakai
				#	Gg		Gallus gallus			Chicken
				#	org.Hs.ipi.db	not actually related? IPI protein database information relating to humans
				#	Mmu		??						Rhesus
				#	Pf		Plasmodium falciparum	Malaria
				#	Pt		Pan troglodytes			Chimp
				#	Rn		Rattus norvegicus		Rat
				#	Sc		Saccharomyces cerevisiae Yeast
				#	Sco		Streptomyces coelicolor
				#	Ss		Sus scrofa				Pig
				#	Tgondii	Toxoplasma gondii
				#	Xl		Xenopus laevus			Frog
			),
			
		organism_prefixed = list(
					"ag"="aga",	"at"="ath",
					"bt"="bta",
					"ce"="cel", "cf"="cfa",
					"dm"="dme", "dr"="dre",
					"ec"="eco","eck12"="eco","ecsakai"="eco",
					"gg"="gga",
					"hs"="hsa",
					"mm"="mmu","mmu"="mcc" ,#Mm= mus musculus = mouse BUT mmu=macaca mulatta=rhesus monkey
					"pf"="pfa","pt"="ptr",
					"rn"="rno",
					"sc"="sce","sco"="sco","ss"="ssc",
					"tgondii"="tgo",
					"xl"="xla"
			)
    )#ORG.eg.Db
    
    # -------------------------------------------------------------- #

	.WHOPGENenv$pp = function( ... ) return( paste( ... , sep="") )
	.WHOPGENenv$pps = function( ... ) return( paste( ... , sep=" ") )

    # -------------------------------------------------------------- #


	##
	##
	##			UCSC
	##
	##
    .WHOPGENenv$ucsc <- list(
		limit = 100,
		db_host = "genome-mysql.cse.ucsc.edu",
		db_username = "genome",
		db_password = "",

		#
		#
		connection = NULL,
		
		#
		#
		getConnection = function(){
			if( is.null( .WHOPGENenv$ucsc[["connection"]] ) )
			{
				.WHOPGENenv$ucsc$connect()
			}
			return( .WHOPGENenv$ucsc[["connection"]] )
		},

		##
		##	connect to UCSC SQL server
		##
		connect = function()
		{
			.WHOPGENenv$ucsc[["connection"]] <- DBI::dbConnect( DBI::dbDriver("MySQL"), user=.WHOPGENenv$ucsc[["db_username"]], pass=.WHOPGENenv$ucsc[["db_password"]], host=.WHOPGENenv$ucsc[["db_host"]] )
		},

		##
		##	disconnect from UCSC SQL server
		##
		disconnect = function()
		{
			DBI::dbDisconnect(conn=.WHOPGENenv$ucsc[["connection"]])
			.WHOPGENenv$ucsc[["connection"]] = NULL
		}
		
		##
		##
    )#UCSC

}#.onLoad



