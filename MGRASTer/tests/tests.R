##############################################################################
#
#  Although this package is essentially a wrapper for a distinct API,
#  the design is intended to ensure stability with respect to the R platform,
#  independent of changes to that API.
#
#  The MG-RAST API is developed in a dynamic environment.  Our goal here is
#  for package functionality to handle well in the face of significant API changes,
#  patches, hotfixes, etc.
#
#  (Of course, this only applies to interim periods between package versions,
#  which do regularly incorporate latest API updates.)
#
#  In particular, these package tests test the package, not the API.  They are 
#  written with the thought that the package will still "work" in a formal sense,
#  even when the API does not.  This refers to buggy behavior, s as well as 
#  occasional server downtime.
#
#  The effort made to accommodate API updates dynamically is par of this approach, 
#  comprised of the build.MGRAST() and load.MGRAST() functions.  They can only 
#  provide so much resilience, of course.
#
#  One of the challenges of this project has been that, although the API
#  is formally specified, best results in practice have required experimentation
#  and been modeled on examples.
#
#  Accordingly, these URLs from API documentation form the backbone of testing. 
#  Some are also copied to the doc examples of call.MGRAST() and parse.MGRAST().
#  Updates here may require updates there.  Last updated April 2014.
#
#  Per correspondence with CRAN, the policy in these tests is not to actually
#  attempt communicaton with the MG-RAST API server.  Such tests are present 
#  below in comments, however.
#
#
#  http://api.metagenomics.anl.gov/annotation/sequence/mgm4447943.3?evalue=10&type=organism&source=SwissProt
#  http://api.metagenomics.anl.gov/annotation/similarity/mgm4447943.3?identity=80&type=function&source=KO
#  http://api.metagenomics.anl.gov/compute/alphadiversity/mgm4447943.3?level=order
#  http://api.metagenomics.anl.gov/download/mgm4447943.3?file=350.1
#  http://api.metagenomics.anl.gov/download/mgm4447943.3?stage=650
#  http://api.metagenomics.anl.gov/library?limit=20&order=name
#  http://api.metagenomics.anl.gov/library/mgl52924?verbosity=full
#  http://api.metagenomics.anl.gov/m5nr/ontology?source=Subsystems&min_level=level3
#  http://api.metagenomics.anl.gov/m5nr/taxonomy?filter=Bacteroidetes&filter_level=phylum&min_level=genus
#  http://api.metagenomics.anl.gov/m5nr/sources
#  http://api.metagenomics.anl.gov/m5nr/accession/YP_003268079.1
#  http://api.metagenomics.anl.gov/m5nr/alias/IPR001478
#  http://api.metagenomics.anl.gov/m5nr/md5/000821a2e2f63df1a3873e4b280002a8?source=InterPro
#  http://api.metagenomics.anl.gov/m5nr/function/sulfatase?source=GenBank
#  http://api.metagenomics.anl.gov/m5nr/organism/akkermansia?source=KEGG
#  http://api.metagenomics.anl.gov/m5nr/sequence/MAGENHQWQGSIL?source=TrEMBL
#  http://api.metagenomics.anl.gov/matrix/organism?id=mgm4447943.3&id=mgm4447192.3&id=mgm4447102.3&id=mgm4447103.3&group_level=family&source=RefSeq&result_type=abundance&evalue=15
#  http://api.metagenomics.anl.gov/matrix/function?id=mgm4447943.3&id=mgm4447192.3&id=mgm4447102.3&id=mgm4447103.3&group_level=level3&source=Subsystems&result_type=abundance&identity=80&filter_level=phylum&filter=Firmicutes
#  http://api.metagenomics.anl.gov/matrix/feature?id=mgm4447943.3&id=mgm4447192.3&id=mgm4447102.3&id=mgm4447103.3&source=KEGG&result_type=evalue&length=25
#  http://api.metagenomics.anl.gov/metadata/template
#  http://api.metagenomics.anl.gov/metadata/cv
#  http://api.metagenomics.anl.gov/metadata/export/mgp128
#  http://api.metagenomics.anl.gov/metagenome?limit=20&order=name
#  http://api.metagenomics.anl.gov/metagenome/mgm4447943.3?verbosity=metadata
#  http://api.metagenomics.anl.gov/project?limit=20&order=name
#  http://api.metagenomics.anl.gov/project/mgp128?verbosity=full
#  http://api.metagenomics.anl.gov/sample?limit=20&order=name
#  http://api.metagenomics.anl.gov/sample/mgs25823?verbosity=full
#  http://api.metagenomics.anl.gov/validation/template/
#  http://api.metagenomics.anl.gov/validation/data/?template=
#
##############################################################################


library(MGRASTer)

#-----------------------------------------------------------------------------
#  test package environment & local API representation.
#  also, test access and slicing of documentation.
#  we would like a simpler way to list controlled vocabularies
#  than the last two examples, but that must wait.
#-----------------------------------------------------------------------------

ls (.MGRAST)
get ("server", envir=.MGRAST)
get ("this.package", envir=.MGRAST)
get ("API.version", envir=.MGRAST)
get ("API", envir=.MGRAST)
API <- get ("API", envir=.MGRAST)

#-----------------------------------------------------------------------------
#  test docs
#-----------------------------------------------------------------------------

doc.MGRAST()
doc.MGRAST(2)
doc.MGRAST(3)

doc.MGRAST (head="matrix")
doc.MGRAST (2, head="matrix")
doc.MGRAST (3, head="matrix")

doc.MGRAST (stratum="parameters")
doc.MGRAST (2, stratum="parameters")
doc.MGRAST (3, stratum="parameters")

doc.MGRAST (stratum="options")
doc.MGRAST (2, stratum="options")
doc.MGRAST (3, stratum="options")

doc.MGRAST (stratum="required")
doc.MGRAST (2, stratum="required")
doc.MGRAST (3, stratum="required")

doc.MGRAST (stratum="attributes")
doc.MGRAST (2, stratum="attributes")
doc.MGRAST (3, stratum="attributes")

doc.MGRAST (2, head = c('mat','orga','param','opt','group_level'))
doc.MGRAST (2, head = c('mat','func','param','opt','group_level'))

#-----------------------------------------------------------------------------
#  test URL construction from arguments, without issuing any calls.
#
#  the first set are calls from API documentation, adapted by hand, April 2014.
#  some are copied in the doc examples of call.MGRAST() and parse.MGRAST().
#  if these are updated, so should those be.  (see above for the original URLs.)
#
#  the second set are test cases put together by hand.
#  calls intended to produce errors are "try"d.
#-----------------------------------------------------------------------------

call.MGRAST ('an', 'se', id=4447943.3, ev=10, ty='or', so='Sw', destfile="no.file", issue=FALSE)
call.MGRAST ('an', 'si', id=4447943.3, iden=80, ty='fu', so='KO', destfile="no.file", issue=FALSE)
call.MGRAST ('co', 'al', id=4447943.3, le='or', issue=FALSE)
call.MGRAST ('do', 'ins', id=4447943.3, fi=350.1, destfile="no.file", issue=FALSE)
call.MGRAST ('do', 'se', id=4447943.3, st=650, issue=FALSE)
call.MGRAST ('li', 'qu', li=20, or='na', issue=FALSE)
call.MGRAST ('li', 'ins', id=52924, ve='fu', issue=FALSE)
call.MGRAST ('m5', 'on', so='Subsystems', mi='level3', issue=FALSE)
call.MGRAST ('m5', 'ta', filter='Bacteroidetes', filter_l='phylum', mi='genus', issue=FALSE)
call.MGRAST ('m5', 'so', issue=FALSE)
call.MGRAST ('m5', 'ac', id='YP_003268079.1', issue=FALSE)
call.MGRAST ('m5', 'al', text='IPR001478', issue=FALSE)
call.MGRAST ('m5', 'md', id='000821a2e2f63df1a3873e4b280002a8', so='InterPro', issue=FALSE)
call.MGRAST ('m5', 'fu', text='sulfatase', so='GenBank', issue=FALSE)
call.MGRAST ('m5', 'or', text='akkermansia', so='KEGG', issue=FALSE)
call.MGRAST ('m5', 'se', text='MAGENHQWQGSIL', so='TrEMBL', issue=FALSE)
call.MGRAST ('ma', 'or', id=c(4447943.3, 4447192.3, 4447102.3, 4447103.3), gro='family', so='Ref', resu='ab', ev=15, issue=FALSE)
call.MGRAST ('ma', 'fu', id=c(4447943.3, 4447192.3, 4447102.3, 4447103.3), gro='level3', so='Sub', resu='ab', iden=80, filter_l='phylum', filter='Firmicutes', issue=FALSE)
call.MGRAST ('ma', 'fe', id=c(4447943.3, 4447192.3, 4447102.3, 4447103.3), so='KE', resu='ev', len=25, issue=FALSE)
call.MGRAST ('metadata', 'te', issue=FALSE)
call.MGRAST ('metadata', 'cv', issue=FALSE)
call.MGRAST ('metadata', 'ex', id=128, issue=FALSE)
call.MGRAST ('metageno','qu', li=20, ord='name', issue=FALSE)
call.MGRAST ('metageno','ins', id=4447943.3, ve='me', issue=FALSE)
call.MGRAST ('proj', 'qu', li=20, or='na', issue=FALSE)
call.MGRAST ('proj', 'ins', id=128, ve='fu', issue=FALSE)
call.MGRAST ('sa', 'qu', li=20, or='na', issue=FALSE)
call.MGRAST ('sa', 'ins', id=25823, ve='fu', issue=FALSE)
call.MGRAST ('va', 'te', issue=FALSE)
call.MGRAST ('va', 'da', issue=FALSE)

try (call.MGRAST ("ann", issue=FALSE))
try (call.MGRAST ("ann", "seque", issue=FALSE))
call.MGRAST ("ann", "seque", destfile="no.file", issue=FALSE)
call.MGRAST ("compu", "alpha", issue=FALSE)
call.MGRAST ("downloa", "setl", issue=FALSE)
call.MGRAST ("downloa", "setl", id=5, issue=FALSE)
try (call.MGRAST ("downloa", "sel", id=5, issue=FALSE))
call.MGRAST ("download", "setlist", st=5, id=10, i=11, issue=FALSE)
call.MGRAST ("download", "setlist", stage=5, id=10, id=11, issue=FALSE)
call.MGRAST ("mat", "org", issue=FALSE)
try (call.MGRAST ("mat", "func", fil="a", issue=FALSE))
call.MGRAST ("mat", "func", filter="a", issue=FALSE)
try (call.MGRAST ("mat", "func", filter="a", filter_="b", issue=FALSE))
try (call.MGRAST ("mat", "func", filter="a", filter_l="b", filter_s="c", issue=FALSE))

#-----------------------------------------------------------------------------
#  test calls retrieving info pages
#-----------------------------------------------------------------------------
# 
# for (xx in names (API)) {
# 	if (!inherits (try (call.MGRAST (xx, 'info', verify=FALSE)), "try-error")) {
# 		message ("\'info\' call successful for: ", xx)
# 	} else
# 		message ("\'info\' call failed for: ", xx)
# 	}

#-----------------------------------------------------------------------------
#  test parsing of URLs from API documentation, without issuing calls.
#-----------------------------------------------------------------------------

flat.API <- unlist (get ("API", .MGRAST))
ee <- flat.API [grep("example", names(flat.API), fixed=T)]
ee <- unname(ee) [substr(ee, 1, 4) == "http"]

for (xx in ee) {
	message ("Parsing test URL: ", xx)
	str (try (parse.MGRAST (xx)))
	}

#-----------------------------------------------------------------------------
#  test parsing of URLS from API documentation [and the calls themselves]
#  for simplicity, we send all results to file,
#  rather than bothering to identify cases where that's definitely required.
#-----------------------------------------------------------------------------

ff <- list()
for (xx in ee) {
	li <- parse.MGRAST (xx)
	ff [[xx]] <- tempfile (li$resource)
	li [c ('quiet', 'destfile', 'issue')] <- list (TRUE, ff [[xx]], FALSE)
	res <- try (do.call (call.MGRAST, li))	
	if (xx == res) {
		message ('Matched URL: ', xx)
	} else
		message('Unmatched URL:\n\t', xx, '\n\t', res)

#	li$issue <- TRUE
#	try (do.call (call.MGRAST, li))

	}
print (unlist (unname (ff)))
sapply (ff, unlink)



