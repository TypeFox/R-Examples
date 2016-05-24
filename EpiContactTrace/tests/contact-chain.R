## Copyright 2013-2014 Stefan Widgren and Maria Noremark,
## National Veterinary Institute, Sweden
##
## Licensed under the EUPL, Version 1.1 or - as soon they
## will be approved by the European Commission - subsequent
## versions of the EUPL (the "Licence");
## You may not use this work except in compliance with the
## Licence.
## You may obtain a copy of the Licence at:
##
## http://ec.europa.eu/idabc/eupl
##
## Unless required by applicable law or agreed to in
## writing, software distributed under the Licence is
## distributed on an "AS IS" basis,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
## express or implied.
## See the Licence for the specific language governing
## permissions and limitations under the Licence.

library(EpiContactTrace)

##
## Check in- and outgoing contact chain methods
##

##
## Case 1
##

movements <- structure(list(source = 1:7, destination = c(4L, 5L, 5L,
6L, 8L, 8L, 8L), t = structure(c(14849, 14846, 14847, 14850, 14848,
14851, 14852), class = "Date")), .Names = c("source", "destination",
"t"), class = "data.frame", row.names = c(NA, -7L))
ct <- Trace(movements,
            root = 8L,
            inBegin = as.Date('2010-08-22'),
            inEnd = as.Date('2010-10-01'),
            outBegin = as.Date('2010-08-01'),
            outEnd = as.Date('2010-08-31'))
ct
stopifnot(identical(IngoingContactChain(ct)$ingoingContactChain, 7L))
stopifnot(identical(OutgoingContactChain(ct)$outgoingContactChain, 0L))

##
## Case 2
##

movements <- structure(list(source = c(1L, 2L, 3L, 3L), destination =
c(3L, 3L, 4L, 4L), t = structure(c(14834, 14838, 14836, 14841), class
= "Date"), individual = c(NA_character_, NA_character_, NA_character_,
NA_character_), n = c(NA_integer_, NA_integer_, NA_integer_,
NA_integer_)), .Names = c("source", "destination", "t", "individual",
"n"), row.names = c(NA, -4L), class = "data.frame")
ct <- Trace(movements,
            root = 4L,
            inBegin = as.Date('2010-07-22'),
            inEnd = as.Date('2010-08-21'),
            outBegin = as.Date('2010-08-01'),
            outEnd = as.Date('2010-08-31'))
ct
stopifnot(identical(IngoingContactChain(ct)$ingoingContactChain, 3L))
stopifnot(identical(OutgoingContactChain(ct)$outgoingContactChain, 0L))

##
## Case 3
##

movements <- structure(list(source = 1:2, destination = c(2L, 1L), t =
structure(c(14834, 14834), class = "Date"), individual =
c(NA_character_, NA_character_ ), n = c(NA_integer_, NA_integer_)),
.Names = c("source", "destination", "t", "individual", "n"), row.names
= c(NA, -2L), class = "data.frame")
ct <- Trace(movements,
            root = 1L,
            inBegin = as.Date('2010-08-02'),
            inEnd = as.Date('2010-09-01'),
            outBegin = as.Date('2010-09-01'),
            outEnd = as.Date('2010-10-01'))
ct
stopifnot(identical(IngoingContactChain(ct)$ingoingContactChain, 1L))
stopifnot(identical(OutgoingContactChain(ct)$outgoingContactChain, 0L))

##
## Case 4
##

movements <- structure(list(source = c(1L, 2L, 2L, 1L, 3L, 7L, 1L),
destination = c(2L, 5L, 6L, 3L, 7L, 8L, 4L), t = structure(c(14834,
14838, 14836, 14857, 14860, 14862, 14884), class = "Date"), individual
= c(NA_character_, NA_character_, NA_character_, NA_character_,
NA_character_, NA_character_, NA_character_), n = c(NA_integer_,
NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_,
NA_integer_)), .Names = c("source", "destination", "t", "individual",
"n"), row.names = c(NA, -7L ), class = "data.frame")
ct <- Trace(movements,
            root = 1L,
            inBegin = as.Date('2010-08-02'),
            inEnd = as.Date('2010-09-01'),
            outBegin = as.Date('2010-08-01'),
            outEnd = as.Date('2010-11-09'))
ct
stopifnot(identical(IngoingContactChain(ct)$ingoingContactChain, 0L))
stopifnot(identical(OutgoingContactChain(ct)$outgoingContactChain, 7L))

##
## Case 5
##

movements <- structure(list(source = 1:2, destination = c(2L, 1L), t =
structure(c(14834, 14834), class = "Date"), individual =
c(NA_character_, NA_character_ ), n = c(NA_integer_, NA_integer_)),
.Names = c("source", "destination", "t", "individual", "n"), row.names
= c(NA, -2L), class = "data.frame")
ct <- Trace(movements,
            root = 1L,
            inBegin = as.Date('2010-07-02'),
            inEnd = as.Date('2010-08-01'),
            outBegin = as.Date('2010-08-01'),
            outEnd = as.Date('2010-08-31'))
ct
stopifnot(identical(IngoingContactChain(ct)$ingoingContactChain, 0L))
stopifnot(identical(OutgoingContactChain(ct)$outgoingContactChain, 1L))

##
## Case 6
##

movements <- structure(list(source = c(1L, 2L, 1L, 2L, 1L, 3L, 1L),
destination = c(2L, 3L, 2L, 3L, 2L, 4L, 2L), t = structure(c(1L, 2L,
3L, 4L, 7L, 6L, 5L), .Label = c("2010-10-01", "2010-10-05",
"2010-10-10", "2010-10-15", "2010-10-20", "2010-10-25", "2010-10-30"),
class = "factor")), .Names = c("source", "destination", "t"), class =
"data.frame", row.names = c(NA, -7L))
ct <- Trace(movements,
            root = 1L,
            inBegin = as.Date('2010-10-10'),
            inEnd = as.Date('2010-10-20'),
            outBegin = as.Date('2010-10-10'),
            outEnd = as.Date('2010-10-20'))
ct
stopifnot(identical(IngoingContactChain(ct)$ingoingContactChain, 0L))
stopifnot(identical(OutgoingContactChain(ct)$outgoingContactChain, 2L))

##
## Case 7
##

movements <- structure(list(source = c(1L, 2L, 1L, 2L, 1L, 3L, 1L),
destination = c(2L, 3L, 2L, 3L, 2L, 4L, 2L), t = structure(c(1L, 2L,
3L, 4L, 7L, 6L, 5L), .Label = c("2010-10-01", "2010-10-05",
"2010-10-10", "2010-10-15", "2010-10-20", "2010-10-25", "2010-10-30"),
class = "factor")), .Names = c("source", "destination", "t"), class =
"data.frame", row.names = c(NA, -7L))

ns <- NetworkSummary(movements, root=1, tEnd="2010-10-20", days=10)

df <- structure(list(root = structure(1L, .Label = "1", class =
"factor"), inBegin = structure(14892, class = "Date"), inEnd =
structure(14902, class = "Date"), inDays = 10L, outBegin =
structure(14892, class = "Date"), outEnd = structure(14902, class =
"Date"), outDays = 10L, inDegree = 0L, outDegree = 1L,
ingoingContactChain = 0L, outgoingContactChain = 2L), .Names =
c("root", "inBegin", "inEnd", "inDays", "outBegin", "outEnd",
"outDays", "inDegree", "outDegree", "ingoingContactChain",
"outgoingContactChain"), row.names = c(NA, -1L), class = "data.frame")
ns
stopifnot(identical(ns, df))
