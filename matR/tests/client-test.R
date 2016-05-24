library(matR)

# #-----------------------------------------------------------------------------------------
# #  NOT OK FOR CRAN
# #-----------------------------------------------------------------------------------------
# 	
# #-----------------------------------------------------------------------------------------
# #  biomRequest() and biom()
# #-----------------------------------------------------------------------------------------
# tt <- tempfile()
# ff <- demoSets() [1]
# gg <- demoSets () [4]
# 
# xx <- readSet (ff)							#  7 metagenomes; IDs only --- vector
# yy <- readSet (gg)							#  32 metagenomes; with metadata --- data.frame
# 
# biomRequest (file=ff)
# biomRequest (file=gg)
# 
# biomRequest(xx)
# biomRequest(xx, request="organism")
# biomRequest(xx, block=3)
# biomRequest(xx, block=3, request="organism")
# 
# biomRequest(xx, quiet=TRUE)
# biomRequest(xx, quiet=TRUE, request="organism")
# biomRequest(xx, quiet=TRUE, block=3)
# biomRequest(xx, quiet=TRUE, block=3, request="organism")
# 
# biomRequest(xx, wait=FALSE) -> zz ; biom(zz, wait=TRUE)
# biomRequest(xx, wait=FALSE, request="organism") -> zz ; biom(zz, wait=TRUE)
# biomRequest(xx, wait=FALSE, block=3) -> zz ; biom(zz, wait=TRUE)
# biomRequest(xx, wait=FALSE, block=3, request="organism") -> zz ; biom(zz, wait=TRUE)
# 
# biomRequest(xx, wait=FALSE, quiet=TRUE) -> zz ; biom(zz, wait=TRUE)
# biomRequest(xx, wait=FALSE, quiet=TRUE, request="organism") -> zz ; biom(zz, wait=TRUE)
# biomRequest(xx, wait=FALSE, quiet=TRUE, block=3) -> zz ; biom(zz, wait=TRUE)
# biomRequest(xx, wait=FALSE, quiet=TRUE, block=3, request="organism") -> zz ; biom(zz, wait=TRUE)
# 
# biomRequest(xx, wait=FALSE) -> zz ; biom(zz, wait=FALSE)
# biomRequest(xx, wait=FALSE, request="organism") -> zz ; biom(zz, wait=FALSE)
# biomRequest(xx, wait=FALSE, block=3) -> zz ; biom(zz, wait=FALSE)
# biomRequest(xx, wait=FALSE, block=3, request="organism") -> zz ; biom(zz, wait=FALSE)
# 
# biomRequest(xx, wait=FALSE, quiet=TRUE) -> zz ; biom(zz, wait=FALSE)
# biomRequest(xx, wait=FALSE, quiet=TRUE, request="organism") -> zz ; biom(zz, wait=FALSE)
# biomRequest(xx, wait=FALSE, quiet=TRUE, block=3) -> zz ; biom(zz, wait=FALSE)
# biomRequest(xx, wait=FALSE, quiet=TRUE, block=3, request="organism") -> zz ; biom(zz, wait=FALSE)
# 
# biomRequest(yy)
# biomRequest(yy, request="organism", source="Greengenes")
# biomRequest(yy, block=10)
# biomRequest(yy, block=10, request="organism", source="Greengenes")
# 
# biomRequest(yy, quiet=TRUE)
# biomRequest(yy, quiet=TRUE, request="organism", source="Greengenes")
# biomRequest(yy, quiet=TRUE, block=10)
# biomRequest(yy, quiet=TRUE, block=10, request="organism", source="Greengenes")
# 
# biomRequest(yy, wait=FALSE) -> zz ; biom(zz, wait=TRUE)
# biomRequest(yy, wait=FALSE, request="organism", source="Greengenes") -> zz ; biom(zz, wait=TRUE)
# biomRequest(yy, wait=FALSE, block=10) -> zz ; biom(zz, wait=TRUE)
# biomRequest(yy, wait=FALSE, block=10, request="organism", source="Greengenes") -> zz ; biom(zz, wait=TRUE)
# 
# biomRequest(yy, wait=FALSE, quiet=TRUE) -> zz ; biom(zz, wait=TRUE)
# biomRequest(yy, wait=FALSE, quiet=TRUE, request="organism", source="Greengenes") -> zz ; biom(zz, wait=TRUE)
# biomRequest(yy, wait=FALSE, quiet=TRUE, block=10) -> zz ; biom(zz, wait=TRUE)
# biomRequest(yy, wait=FALSE, quiet=TRUE, block=10, request="organism", source="Greengenes") -> zz ; biom(zz, wait=TRUE)
# 
# biomRequest(yy, wait=FALSE) -> zz ; biom(zz, wait=FALSE)
# biomRequest(yy, wait=FALSE, request="organism", source="Greengenes") -> zz ; biom(zz, wait=FALSE)
# biomRequest(yy, wait=FALSE, block=10) -> zz ; biom(zz, wait=FALSE)
# biomRequest(yy, wait=FALSE, block=10, request="organism", source="Greengenes") -> zz ; biom(zz, wait=FALSE)
# 
# biomRequest(yy, wait=FALSE, quiet=TRUE) -> zz ; biom(zz, wait=FALSE)
# biomRequest(yy, wait=FALSE, quiet=TRUE, request="organism", source="Greengenes") -> zz ; biom(zz, wait=FALSE)
# biomRequest(yy, wait=FALSE, quiet=TRUE, block=10) -> zz ; biom(zz, wait=FALSE)
# biomRequest(yy, wait=FALSE, quiet=TRUE, block=10, request="organism", source="Greengenes") -> zz ; biom(zz, wait=FALSE)
# 
# 
# 
# 
# 
# 
# biomRequest (xx, quiet=TRUE)
# biomRequest (xx, blocking=3)
# biomRequest (xx, outfile=tt)
# biomRequest (xx, blocking=3, outfile=tt)
# unlink(tt)
# 
# biomRequest (xx, request="function")
# biomRequest (xx, request="function", group_level="level1")
# biomRequest (xx, request="function", group_level="level2")
# biomRequest (xx, request="function", group_level="level3")
# biomRequest (xx, request="function", group_level="function")
# biomRequest (xx, request="function", group_level="level1", evalue=1)
# biomRequest (xx, request="function", group_level="level1", evalue=1, length=20)
# biomRequest (xx, request="function", group_level="level1", evalue=1, length=20, identity=85)
# biomRequest (xx, request="function", group_level="level1", evalue=1, length=20, identity=85, filter_source="NOG")
# 
# biomRequest (xx, request="organism")
# biomRequest (xx, request="organism", group_level="domain")
# biomRequest (xx, request="organism", group_level="phylum")
# biomRequest (xx, request="organism", group_level="species")
# biomRequest (xx, request="organism", group_level="strain")
# biomRequest (xx, request="organism", group_level="domain", evalue=1)
# biomRequest (xx, request="organism", group_level="domain", evalue=1, length=20)
# biomRequest (xx, request="organism", group_level="domain", evalue=1, length=20, filter_source="Greengenes")
# 
# biomRequest (xx, request="organism", hit_type="all")
# biomRequest (xx, request="organism", hit_type="single")
# 
# biomRequest (xx, request="organism", result_type="abundance")
# biomRequest (xx, request="organism", result_type="identity")
# 
# ticket <- biomRequest (xx, wait=FALSE)
# # ... here you can go for a coffee break, or do some other calculations; then later:
# xx <- biom (ticket, wait=TRUE)
# 
# 
# #-----------------------------------------------------------------------------------------
# #  metadata(detail=NULL)			...which means just lookup ids
# #-----------------------------------------------------------------------------------------
# metadata ("mgp21")
# metadata ("mgp21 mgp24")
# metadata ("mgp21 mgp24 mgp30")				# example
# 
# metadata("mgm4440066.3")
# metadata("mgm4440066.3 mgm4440062.3")
# metadata("mgm4440066.3 mgm4440062.3 mgm4440055.3")
# metadata("mgm4440066.3 mgm4441681.3")
# metadata("mgm4440066.3 mgm4440062.3 mgm4441681.3")
# metadata("mgm4440066.3 mgm4440062.3 mgm4440055.3 mgm4441681.3")
# metadata("mgm4440066.3 mgm4441681.3 mgm4441682.3")
# metadata("mgm4440066.3 mgm4440062.3 mgm4441681.3 mgm4441682.3")
# metadata("mgm4440066.3 mgm4440062.3 mgm4440055.3 mgm4441681.3 mgm4441682.3")
# metadata("mgm4440066.3 mgm4441681.3 mgm4441682.3 mgm4440463.3")
# metadata("mgm4440066.3 mgm4440062.3 mgm4441681.3 mgm4441682.3 mgm4440463.3")
# metadata("mgm4440066.3 mgm4440062.3 mgm4440055.3 mgm4441681.3 mgm4441682.3 mgm4440463.3")
# metadata("mgm4440066.3 mgm4441681.3 mgm4440463.3 mgm4440464.3")
# metadata("mgm4440066.3 mgm4440062.3 mgm4441681.3 mgm4440463.3 mgm4440464.3")
# metadata("mgm4440066.3 mgm4440062.3 mgm4440055.3 mgm4441681.3 mgm4440463.3 mgm4440464.3")   # example
# 
# #-----------------------------------------------------------------------------------------
# #  metadata(detail=TRUE)			...which just means verbosity="minimal"
# #-----------------------------------------------------------------------------------------
# metadata ("mgp21", detail=TRUE)
# metadata ("mgp21 mgp24", detail=TRUE)
# metadata ("mgp21 mgp24 mgp30", detail=TRUE)			# example
# metadata ("mgm4440463.3", detail=TRUE)
# metadata ("mgm4440463.3 mgm4440464.3", detail=TRUE)
# metadata ("mgm4440463.3 mgm4440464.3 mgm4441679.3", detail=TRUE)
# metadata ("mgm4440066.3 mgm4440062.3 mgm4440055.3 mgm4441681.3 mgm4440463.3 mgm4440464.3", detail=TRUE) # example
# 
# #-----------------------------------------------------------------------------------------
# #  metadata(detail=		c("minimal","verbose","full")				for projects
# #						c("minimal","metadata","stats","full")		for metagenomes
# #  ...relayed directly as "verbosity" to call.MGRAST()
# #-----------------------------------------------------------------------------------------
# metadata ("mgp21", detail="verbose")
# metadata ("mgp21 mgp24", detail="verbose")
# metadata ("mgp21 mgp24 mgp30", detail="verbose")	# example; show names() of
# metadata ("mgp21", detail="full")
# metadata ("mgp21 mgp24", detail="full")
# metadata ("mgp21 mgp24 mgp30", detail="full")
# 
# metadata ("mgm4440463.3", detail="metadata")
# metadata ("mgm4440463.3 mgm4440464.3", detail="metadata")
# metadata ("mgm4440463.3 mgm4440464.3 mgm4441679.3", detail="metadata")
# metadata ("mgm4440066.3 mgm4440062.3 mgm4440055.3 mgm4441681.3 mgm4440463.3 mgm4440464.3", detail="metadata")    # example; show names() of
# metadata ("mgm4440463.3", detail="stats")										# bad idea, don't show
# metadata ("mgm4440463.3 mgm4440464.3", detail="stats")							# bad idea, etc
# metadata ("mgm4440463.3 mgm4440464.3 mgm4441679.3", detail="stats")				# bad idea, etc
# metadata ("mgm4440463.3", detail="full")										# bad idea..
# metadata ("mgm4440463.3 mgm4440464.3", detail="full")							# bad idea..
# metadata ("mgm4440463.3 mgm4440464.3 mgm4441679.3", detail="full")				# bad idea..
# 
# #-----------------------------------------------------------------------------------------
# #  metadata(file=)
# #-----------------------------------------------------------------------------------------
# tt <- tempfile()
# writeLines ("mgm4440463.3", file=tt) 
# metadata (map=TRUE, file=tt)
# writeLines (c("mgm4440463.3", "mgm4440464.3"), file=tt) 
# metadata (map=TRUE, file=tt)
# writeLines (c("mgm4440463.3", "mgm4440464.3", "mgm4441679.3"), file=tt) 
# metadata (map=TRUE, file=tt)
# writeLines ("mgm4440463.3", file=tt) 
# metadata (map=FALSE, file=tt)
# writeLines (c("mgm4440463.3", "mgm4440464.3"), file=tt) 
# metadata (map=FALSE, file=tt)
# writeLines (c("mgm4440463.3", "mgm4440464.3", "mgm4441679.3"), file=tt) 
# metadata (map=FALSE, file=tt)
# unlink(tt)
# 
# tt <- tempfile()
# writeLines ("mgp21", file=tt); metadata (map=TRUE, file=tt)
# writeLines (c ("mgp21", "mgp24"), file=tt) ; metadata (map=TRUE, file=tt)
# writeLines (c ("mgp21", "mgp24", "mgp30"), file=tt) ; metadata (map=TRUE, file=tt)
# writeLines ("mgp21", file=tt); metadata (map=FALSE, file=tt)
# writeLines (c ("mgp21", "mgp24"), file=tt); metadata (map=FALSE, file=tt)
# writeLines (c ("mgp21", "mgp24", "mgp30"), file=tt) ; metadata (map=FALSE, file=tt)
# unlink(tt)
# 
# #-----------------------------------------------------------------------------------------
# #  dir.MGRAST()
# #-----------------------------------------------------------------------------------------
# dir.MGRAST()
# dir.MGRAST (1, 50)
# dir.MGRAST (1, 100)
# dir.MGRAST (1, 200)
# dir.MGRAST (1, 500)
# dir.MGRAST (1, 1000)
# dir.MGRAST (100, 150)
# dir.MGRAST (100, 200)
# dir.MGRAST (100, 300)
# dir.MGRAST (100, 600)
# dir.MGRAST (100, 1100)
# dir.MGRAST (500, 550)
# dir.MGRAST (500, 600)
# dir.MGRAST (500, 700)
# dir.MGRAST (500, 1000)
# dir.MGRAST (100, len=50)
# dir.MGRAST (100, len=100)
# dir.MGRAST (100, len=200)
# dir.MGRAST (100, len=500)
# dir.MGRAST (100, len=1000)
# 
# dir.MGRAST (1, 50, order="id")
# dir.MGRAST (1, 100, order="id")
# dir.MGRAST (1, 200, order="id")
# dir.MGRAST (1, 500, order="id")
# dir.MGRAST (1, 1000, order="id")
# dir.MGRAST (100, 150, order="id")
# dir.MGRAST (100, 200, order="id")
# dir.MGRAST (100, 300, order="id")
# dir.MGRAST (100, 600, order="id")
# dir.MGRAST (100, 1100, order="id")
# dir.MGRAST (500, 550, order="id")
# dir.MGRAST (500, 600, order="id")
# dir.MGRAST (500, 700, order="id")
# dir.MGRAST (500, 1000, order="id")
# dir.MGRAST (100, len=50, order="id")
# dir.MGRAST (100, len=100, order="id")
# dir.MGRAST (100, len=200, order="id")
# dir.MGRAST (100, len=500, order="id")
# dir.MGRAST (100, len=1000, order="id")
# 
# dir.MGRAST (1, 50, verbosity="verbose")
# dir.MGRAST (1, 100, verbosity="verbose")
# dir.MGRAST (1, 200, verbosity="verbose")
# dir.MGRAST (1, 500, verbosity="verbose")
# dir.MGRAST (1, 1000, verbosity="verbose")
# dir.MGRAST (100, 150, verbosity="verbose")			# lot of columns named NA here.  debug that sometime
# dir.MGRAST (100, 200, verbosity="verbose")
# dir.MGRAST (100, 300, verbosity="verbose")
# dir.MGRAST (100, 600, verbosity="verbose")
# dir.MGRAST (100, 1100, verbosity="verbose")
# dir.MGRAST (500, 550, verbosity="verbose")
# dir.MGRAST (500, 600, verbosity="verbose")
# dir.MGRAST (500, 700, verbosity="verbose")
# dir.MGRAST (500, 1000, verbosity="verbose")
# dir.MGRAST (100, len=50, verbosity="verbose")
# dir.MGRAST (100, len=100, verbosity="verbose")
# dir.MGRAST (100, len=200, verbosity="verbose")
# dir.MGRAST (100, len=500, verbosity="verbose")
# dir.MGRAST (100, len=1000, verbosity="verbose")
# 
# dir.MGRAST (1, 50, verbosity="full")
# dir.MGRAST (1, 100, verbosity="full")
# dir.MGRAST (1, 200, verbosity="full")
# # dir.MGRAST (1, 500, verbosity="full")			# hangs
# # dir.MGRAST (1, 1000, verbosity="full")			# didn't try
# dir.MGRAST (100, 150, verbosity="full")
# dir.MGRAST (100, 200, verbosity="full")
# dir.MGRAST (100, 300, verbosity="full")
# # dir.MGRAST (100, 600, verbosity="full")			# didn't try
# # dir.MGRAST (100, 1100, verbosity="full")		# didn't try
# dir.MGRAST (500, 550, verbosity="full")
# dir.MGRAST (500, 600, verbosity="full")
# dir.MGRAST (500, 700, verbosity="full")
# # dir.MGRAST (500, 1000, verbosity="full")		# didn't try
# dir.MGRAST (100, len=50, verbosity="full")
# dir.MGRAST (100, len=100, verbosity="full")
# dir.MGRAST (100, len=200, verbosity="full")
# # dir.MGRAST (100, len=500, verbosity="full")		# didn't try
# # dir.MGRAST (100, len=1000, verbosity="full")	# didn't try
# 
# dir.MGRAST (offset=0, limit=50)
# dir.MGRAST (offset=0, limit=100)
# dir.MGRAST (offset=0, limit=200)
# dir.MGRAST (offset=0, limit=500)
# dir.MGRAST (offset=0, limit=1000)
# dir.MGRAST (offset=99, limit=50)
# dir.MGRAST (offset=99, limit=100)
# dir.MGRAST (offset=99, limit=200)
# dir.MGRAST (offset=99, limit=500)
# dir.MGRAST (offset=99, limit=1000)
# dir.MGRAST (offset=499, limit=50)
# dir.MGRAST (offset=499, limit=100)
# dir.MGRAST (offset=499, limit=200)
# dir.MGRAST (offset=499, limit=500)
# 
# #-----------------------------------------------------------------------------------------
# #  search.MGRAST()
# #-----------------------------------------------------------------------------------------
# 
# 
# 
# 
# 
# #-----------------------------------------------------------------------------------------
# #  API resources
# #-----------------------------------------------------------------------------------------
# 
# names(.MGRAST$API [[c('matrix','organism','parameters','options')]])
# names(.MGRAST$API [[c('matrix','function','parameters','options')]])
# names(.MGRAST$API [[c('matrix','feature','parameters','options')]])
# # 
# #  $ organism.parameters.options:List of 15
# #   ..$ asynchronous : chr [1:2] "boolean" "if true return process id to query status resource for results, default is false"
# #   ..$ source       :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 14
# #   ..$ result_type  :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 4
# #   ..$ filter       : chr [1:2] "string" "filter the return results to only include abundances based on genes with this function"
# #   ..$ group_level  :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 8
# #   ..$ taxid        : chr [1:2] "boolean" "if true, return annotation ID as NCBI tax id. Only for group_levels with a tax_id"
# #   ..$ grep         : chr [1:2] "string" "filter the return results to only include annotations that contain this text"
# #   ..$ length       : chr [1:2] "int" "value for minimum alignment length cutoff: default is 15"
# #   ..$ evalue       : chr [1:2] "int" "negative exponent value for maximum e-value cutoff: default is 5"
# #   ..$ identity     : chr [1:2] "int" "percent value for minimum % identity cutoff: default is 60"
# #   ..$ filter_source:List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 4
# #   ..$ hide_metadata: chr [1:2] "boolean" "if true do not return metagenome metadata in 'columns' object, default is false"
# #   ..$ id           : chr [1:2] "string" "one or more metagenome or project unique identifier"
# #   ..$ filter_level :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 4
# #   ..$ hit_type     :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 3
# # 
# #  $ function.parameters.options:List of 13
# #   ..$ asynchronous : chr [1:2] "boolean" "if true return process id to query status resource for results, default is false"
# #   ..$ source       :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 4
# #   ..$ result_type  :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 4
# #   ..$ filter       : chr [1:2] "string" "filter the return results to only include abundances based on genes with this organism"
# #   ..$ group_level  :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 4
# #   ..$ grep         : chr [1:2] "string" "filter the return results to only include annotations that contain this text"
# #   ..$ length       : chr [1:2] "int" "value for minimum alignment length cutoff: default is 15"
# #   ..$ evalue       : chr [1:2] "int" "negative exponent value for maximum e-value cutoff: default is 5"
# #   ..$ identity     : chr [1:2] "int" "percent value for minimum % identity cutoff: default is 60"
# #   ..$ filter_source:List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 14
# #   ..$ hide_metadata: chr [1:2] "boolean" "if true do not return metagenome metadata in 'columns' object, default is false"
# #   ..$ id           : chr [1:2] "string" "one or more metagenome or project unique identifier"
# #   ..$ filter_level :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 8
# # 
# #  $ feature.parameters.options :List of 11
# #   ..$ asynchronous   : chr [1:2] "boolean" "if true return process id to query status resource for results, default is false"
# #   ..$ source         :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 12
# #   ..$ result_type    :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 4
# #   ..$ filter         : chr [1:2] "string" "filter the return results to only include abundances based on genes with this organism"
# #   ..$ length         : chr [1:2] "int" "value for minimum alignment length cutoff: default is 15"
# #   ..$ evalue         : chr [1:2] "int" "negative exponent value for maximum e-value cutoff: default is 5"
# #   ..$ identity       : chr [1:2] "int" "percent value for minimum % identity cutoff: default is 60"
# #   ..$ hide_annotation: chr [1:2] "boolean" "if true do not return feature metadata in 'rows' object, default is false"
# #   ..$ id             : chr [1:2] "string" "one or more metagenome or project unique identifier"
# #   ..$ hide_metadata  : chr [1:2] "boolean" "if true do not return metagenome metadata in 'columns' object, default is false"
# #   ..$ filter_level   :List of 2
# #   .. ..$ : chr "cv"
# #   .. ..$ :List of 8
