# Concept -------------------------------------------------------------------

# The package 'rstreams' uses S4 classes and methods to implement random
# streams. It allows to have multiple streams that can be handled 
# independently. However, there is an important difference to other R
# objects. Unlike other R objects Rstream objects store necessary data
# outside the R object. Thus they behave like environments which
# cannot be copied. If you assign the same Rstream object (or
# environment) to several symbols and change one, the others will
# change, too (see 'R Language Definition', Version 1.9.0,
# Sect. 2.1.10 [Environments])
#
# All data are changed during sampling from the Rstream and all data
# that can be modified by replacement methods must be stored outside
# the R object, in particular:
#  (*) the state of the generator of the Rstream;
#  (*) the name of the Rstream;
#  (*) the antithetic flag; and
#  (*) the increased precision flag.
# This can be done by means of an external pointer (as in
# "rstream.lecuyer"), or by variables in environments each associated
# to a particular instance of a Stream object (as in "rstream.runif"),
# or by both techniques.
#
# When adding a new source of random streams one should keep in mind
# that all 'get' methods only return default methods and that
# replacement methods are discabled. Thus if data should be subject to
# changes by the user the extended class must provide the
# corresponding methods.
#
# To make Rstream objects able to being copied they must be PACKED.
# Hence an pack/unpack method should be provided.
# When packing an Rstream object the name, antithetic and increased
# precision flag must be put into the pack list. Use the slots "name",
# "anti", and "incp" (see "rstream.runif" for an example). These will
# be used by the print command when the Rstream object is packed.
# Make sure that one cannot change a PACKED instance of an Rstream object.

