
library(BIOM.utils)

str (jtxt, list.len=5)
str (smat, list.len=5)
str (dmat, list.len=5)
str (li1, list.len=5)
str (li2, list.len=5)
str (li3, list.len=5)
str (li4, list.len=5)

#-----------------------------------------------------------------------------
#  test matrix manipulations
#-----------------------------------------------------------------------------

onecol <- dmat [,1,drop=FALSE]
dense2sparse (onecol)
sparse2dense (dense2sparse (onecol))
matrix2list (onecol)

storage.mode(onecol) <- "character"
dense2sparse (onecol)
sparse2dense (dense2sparse (onecol))
matrix2list (onecol)

onerow <- dmat [1,,drop=FALSE]
dense2sparse (onerow)
sparse2dense (dense2sparse (onerow))
matrix2list (onerow)

storage.mode(onerow) <- "character"
dense2sparse (onerow)
sparse2dense (dense2sparse (onerow))
matrix2list (onerow)

dm <- dmat
dense2sparse (dm)
sparse2dense (dense2sparse (dm))
matrix2list (dm)

storage.mode(dm) <- "character"
dense2sparse (dm)
sparse2dense (dense2sparse (dm))
matrix2list (dm)

sm <- smat ; sm [,1:2] <- sm [,1:2] + 1
sparse2dense (sm)
dense2sparse (sparse2dense (sm))
matrix2list (smat)

storage.mode(sm) <- "character"
sparse2dense (sm)
dense2sparse (sparse2dense (sm))
matrix2list (smat)

#-----------------------------------------------------------------------------
#  constructing from matrix
#-----------------------------------------------------------------------------

#  from named matrix
applyBiomMethods (biom (dmat))

#  from nameless matrix
applyBiomMethods (biom (unname (dmat)))

#  from named matrix; type specified
applyBiomMethods (biom (dmat, type="Taxon"))

#  from nameless matrix; sparse not assumed
applyBiomMethods (biom (smat, type="Taxon"))

#  from nameless matrix; sparse specified via dims
#  note, nothing prevents declaring a sparse matrix bigger than its minimnum size
applyBiomMethods (biom (smat, type="Taxon", sparse=c(266,4)))
applyBiomMethods (biom (smat, type="Taxon", sparse=c(270,5)))
applyBiomMethods (biom (smat, type="Taxon", sparse=c(270,6)))

#  but it may not be smaller, so this is an error
try (biom (smat, type="Taxon", sparse=c(266,3)))

#  from nameless matrix; sparse specified via dimnames
dd <- list (as.character (200:465), c('A','B','C','D'))
applyBiomMethods (biom (smat, type="Taxon", sparse=dd))

#  from named matrix; ensure matrix_element_type is "int"
mm <- dmat
storage.mode(mm) <- "integer"
applyBiomMethods (biom(mm, type="Taxon"))

#  from named matrix; ensure matrix_element_type is "float"
storage.mode(mm) <- "double"
applyBiomMethods (biom(mm, type="Taxon"))

#  from named matrix; ensure matrix_element_type is "Unicode"
storage.mode(mm) <- "character"
applyBiomMethods (biom(mm, type="Taxon"))

#  from nameless matrix; sparse specified; ensure matrix_element_type is "int"
mm <- smat
storage.mode(mm) <- "integer"
applyBiomMethods (biom(mm, type="Taxon", sparse=c(266,4)))

#  from nameless matrix; sparse specified; ensure matrix_element_type is "float"
storage.mode(mm) <- "double"
applyBiomMethods (biom(mm, type="Taxon", sparse=c(266,4)))

#  from named matrix; ensure matrix_element_type is "Unicode"
storage.mode(mm) <- "character"
applyBiomMethods (biom(mm, type="Taxon", sparse=c(266,4)))

#-----------------------------------------------------------------------------
#  constructing from character (JSON text)
#-----------------------------------------------------------------------------

applyBiomMethods (biom (jtxt))
applyBiomMethods (biom (file=exampleBiomFile()))

#-----------------------------------------------------------------------------
#  constructing from list
#-----------------------------------------------------------------------------

#  an empty object
applyBiomMethods (biom (list()))

#  from list of:  named matrix (dense)
applyBiomMethods (biom (list(
	data=dmat)))

#  from list of:  nameless matrix (dense)
applyBiomMethods (biom (list(
	data=unname(dmat))))

#  from list of:  named matrix (dense), biom type
applyBiomMethods (biom (list(
	data=dmat, 
	type="Taxon")))

#  from example list:  unnamed matrix (dense), biom type
li <- li1
li$data <- unname(li$data)
applyBiomMethods (biom (li))

#  from list of:  named matrix (dense) and some annotations (but no metadata)
applyBiomMethods (biom (list(
	data=dmat, 
	type="Taxon", 
	matrix_type="dense",
	id="my first biom matrix",
	generated_by="science",
	comment="a very profound result")))

#  from list of:  unnamed matrix, sparse indicated
applyBiomMethods (biom (list(
	data=smat,
	matrix_type="sparse",
	shape=c(266,4))))

#  from list of:  named matrix, matrix_element_type
applyBiomMethods (biom (list (
	data=dmat,
	matrix_element_type = "int")))
applyBiomMethods (biom (list (
	data=dmat,
	matrix_element_type = "float")))
applyBiomMethods (biom (list (
	data=dmat,
	matrix_element_type = "unicode")))

#  from list of:  named matrix (dense), biom type, metadata (which will supercede dimnames)
applyBiomMethods (biom (li2))

#  from list of:  rowlist (dense)
applyBiomMethods (biom (li3 ["data"]))

#  from list of:  rowlist (dense), summary(biom type, metadata
applyBiomMethods (biom (li3))

#  from list of:  rowlist (dense), biom type, metadata for columns only
li <- li3
li$rows <- NULL
applyBiomMethods (biom(li))

#  from complete component list; data given as (sparse) entry list
applyBiomMethods (biom (li4))

#  from complete component list; data given as sparse matrix
li <- li4
li$data <- t (simplify2array (li$data))
applyBiomMethods (biom (li4))

#  from complete component list; data given as matrix; metadata supercedes dimnames
li <- li4
li$data <- sparse2dense (t (simplify2array (li$data) + c(1,1,0)))
dimnames (li$data) <- list (1:nrow(li$data), 1:ncol(li$data))
li$matrix_type <- "dense"
applyBiomMethods (biom (li))

#  same as above, but omit to reset matrix_type ... so, an error
li <- li4
li$data <- sparse2dense (t (simplify2array (li$data) + c(1,1,0)))
try (applyBiomMethods (biom (li)))

#  from complete component list; data given as dense row list
li <- li4
li$data <- matrix2list (sparse2dense (t (simplify2array (li$data) + c(1,1,0))))
li$matrix_type <- "dense"
applyBiomMethods (biom (li))

#-----------------------------------------------------------------------------
#  robustness of conversion
#
#  1 - return to original type via biom
#  2 - convert to another type via biom
#  3 - convert to biom via another type via biom
#-----------------------------------------------------------------------------

tt <- tempfile()
as.character (biom (jtxt))								# similar JSON text, prob not same
as.character (biom (file=exampleBiomFile()), file=tt)
as.matrix (biom (jtxt))								# sparse matrix
as.matrix (biom (jtxt), expand=TRUE)				# dense matrix
as.list (biom (jtxt))								# list of components
biom (as.matrix (biom (jtxt)))						# note, sparseness is lost
biom (as.matrix (biom (jtxt)), sparse=c(266,4))		# sparse retained
biom (as.matrix (biom (jtxt), expand=TRUE))			# sparseness correctly lost
biom (as.list (biom (jtxt)))						# drops metadata #####
unlink (tt)

as.matrix (biom (smat, sparse=c(266,4)))				# identity
as.matrix (biom (smat, sparse=c(266,4)), expand=TRUE)	# expansion
as.character (biom (smat, sparse=TRUE))				# JSON text
as.list (biom (smat, sparse=TRUE))					# list of components
biom (as.character (biom (smat, sparse=c(266,4))))		# same result from 1st and 2nd biom()

as.matrix (biom (dmat))				# identity
as.character (biom (dmat))			# JSON text
as.list (biom (dmat))				# list of components
biom (as.character (biom (dmat)))	# same result from 1st and 2nd biom()
biom (as.list (biom (dmat)))		# loses dimnames (via row.ids and column.ids) #####

as.list (biom (li1))			# list of components
as.matrix (biom (li1))			# just the data
as.character (biom (li))		# JSON text from short list
biom (as.matrix (biom (li1)))		# same result from 1st and 2nd biom()
biom (as.character (biom (li)))		# same result from 1st and 2nd biom()

as.list (biom (li4))					# slight modification of li4
as.matrix (biom (li4))					# strip sparse matrix
as.matrix (biom (li4), expand=TRUE)		# expand sparse matrix
as.character (biom (li4))						# JSON text
biom (as.matrix (biom (li4)), expand=TRUE, sparse=TRUE)		# metadata and annotations stripped
biom (as.character (biom (li4)))		# same result from 1st and 2nd biom()


#-----------------------------------------------------------------------------
# test passing in and out of a file
#-----------------------------------------------------------------------------

biom(as.character(biom(li4)))
tt <- tempfile()
as.character (biom (li4), file=tt)
biom (file=tt)

#-----------------------------------------------------------------------------
#  test specifying matrix_element_type implicitly with storage.mode in matrix construction
#-----------------------------------------------------------------------------

storage.mode(dmat) <- "integer"			# matrix_element_type "int"
biom(dmat)
as.character(biom(dmat))
as.character (biom (dmat), file=tt)
biom (file=tt)

storage.mode(dmat) <- "double"			# matrix_element_type "float"
biom(dmat)
as.character(biom(dmat))
as.character (biom (dmat), file=tt)
biom (file=tt)

# storage.mode(dmat) <- "character"		# matrix_element_type "unicode"
biom(dmat)
as.character(biom(dmat))
as.character (biom (dmat), file=tt)
biom (file=tt)

#-----------------------------------------------------------------------------
#  test specifying matrix_element_type in list construction
#-----------------------------------------------------------------------------

biom(li4, matrix_element_type="int")
as.character (biom (li4, matrix_element_type="int"))
as.character (biom (li4, matrix_element_type="int"), file=tt)
biom (file=tt)

biom(li4, matrix_element_type="float")
as.character(biom(li4, matrix_element_type="float"))
as.character (biom(li4, matrix_element_type="float"), file=tt)
biom (file=tt)

biom(li4, matrix_element_type="unicode")
as.character (biom (li4, matrix_element_type="unicode"))
as.character (biom (li4, matrix_element_type="unicode"), file=tt)
biom (file=tt)

unlink (tt)
