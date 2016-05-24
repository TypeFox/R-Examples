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
## Check in- and outgoing degree methods
##

##
## Case 1
##
movements <- structure(list(source = c(1L, 2L, 3L, 3L),
                            destination = c(3L, 3L, 4L, 4L),
                            t = structure(c(14834, 14838, 14836, 14841), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -4L), class = "data.frame")

ct <- Trace(movements,
            root = 4L,
            inBegin = as.Date('2010-08-02'),
            inEnd = as.Date('2010-09-01'),
            outBegin = as.Date('2010-08-01'),
            outEnd = as.Date('2010-08-31'))

stopifnot(identical(InDegree(ct)$inDegree, 1L))
stopifnot(identical(OutDegree(ct)$outDegree, 0L))

##
## Case 2
##
movements <- structure(list(source = c(1L, 2L, 3L, 3L),
                            destination = c(3L, 3L, 4L, 4L),
                            t = structure(c(14834, 14838, 14836, 14841), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -4L), class = "data.frame")

ns <- NetworkSummary(movements, root = 4, tEnd = '2010-09-01', days = 30)
stopifnot(identical(ns$inDegree, 1L))

##
## Case 3
##
movements <- structure(list(source = c(1L, 2L, 3L, 3L),
                            destination = c(3L, 3L, 4L, 4L),
                            t = structure(c(14834, 14838, 14836, 14841), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -4L), class = "data.frame")

ns <- NetworkSummary(movements, root = 4, tEnd = '2010-08-31', days = 30)
stopifnot(identical(ns$outDegree, 0L))

##
## Case 4
##
movements <- structure(list(source = c(1L, 2L, 3L, 3L),
                            destination = c(3L, 3L, 4L, 4L),
                            t = structure(c(14834, 14838, 14836, 14841), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -4L), class = "data.frame")

ct <- Trace(movements,
            root = 4L,
            inBegin = as.Date('2010-08-27'),
            inEnd = as.Date('2010-09-01'),
            outBegin = as.Date('2010-08-01'),
            outEnd = as.Date('2010-08-31'))

stopifnot(identical(InDegree(ct)$inDegree, 0L))
stopifnot(identical(OutDegree(ct)$outDegree, 0L))

##
## Case 5
##
movements <- structure(list(source = c(1L, 2L, 3L, 3L),
                            destination = c(3L, 3L, 4L, 4L),
                            t = structure(c(14834, 14838, 14836, 14841), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -4L), class = "data.frame")

ns <- NetworkSummary(movements, root = 4, tEnd = '2010-09-01', days = 5)
stopifnot(identical(ns$inDegree, 0L))

##
## Case 6
##
movements <- structure(list(source = c(1L, 2L, 3L, 3L),
                            destination = c(3L, 3L, 4L, 4L),
                            t = structure(c(14834, 14838, 14836, 14841), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -4L), class = "data.frame")

ns <- NetworkSummary(movements, root = 4, tEnd = '2010-08-31', days = 30)
stopifnot(identical(ns$outDegree, 0L))

##
## Case 7
##
movements <- structure(list(source = c("1", "1", "1"),
                            destination = c("2", "3", "4"),
                            t = structure(c(14834, 14838, 14836), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -3L), class = "data.frame")

ct <- Trace(movements,
            root = 1L,
            inBegin = as.Date('2010-08-02'),
            inEnd = as.Date('2010-09-01'),
            outBegin = as.Date('2010-08-01'),
            outEnd = as.Date('2010-08-31'))

stopifnot(identical(InDegree(ct)$inDegree, 0L))
stopifnot(identical(OutDegree(ct)$outDegree, 3L))

##
## Case 8
##
movements <- structure(list(source = c("1", "1", "1"),
                            destination = c("2", "3", "4"),
                            t = structure(c(14834, 14838, 14836), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -3L), class = "data.frame")

ns <- NetworkSummary(movements, root = 1, tEnd = '2010-09-01', days = 30)
stopifnot(identical(ns$inDegree,0L))

##
## Case 9
##
movements <- structure(list(source = c("1", "1", "1"),
                            destination = c("2", "3", "4"),
                            t = structure(c(14834, 14838, 14836), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -3L), class = "data.frame")

ns <- NetworkSummary(movements, root = 1, tEnd = '2010-08-31', days = 30)
stopifnot(identical(ns$outDegree, 3L))

##
## Case 10
##
movements <- structure(list(source = c("1", "1", "1"),
                            destination = c("2", "3", "4"),
                            t = structure(c(14834, 14838, 14836), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -3L), class = "data.frame")

ct <- Trace(movements,
            root = 1L,
            inBegin = as.Date('2010-08-02'),
            inEnd = as.Date('2010-09-01'),
            outBegin = as.Date('2010-08-01'),
            outEnd = as.Date('2010-08-16'))

stopifnot(identical(InDegree(ct)$inDegree, 0L))
stopifnot(identical(OutDegree(ct)$outDegree, 2L))

##
## Case 11
##
movements <- structure(list(source = c("1", "1", "1"),
                            destination = c("2", "3", "4"),
                            t = structure(c(14834, 14838, 14836), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -3L), class = "data.frame")

ns <- NetworkSummary(movements, root = 1, tEnd = '2010-09-01', days = 30)
stopifnot(identical(ns$inDegree, 0L))

##
## Case 12
##
movements <- structure(list(source = c("1", "1", "1"),
                            destination = c("2", "3", "4"),
                            t = structure(c(14834, 14838, 14836), class = "Date"),
                            individual = c(NA_character_, NA_character_, NA_character_),
                            n = c(NA_integer_, NA_integer_, NA_integer_)),
                       .Names = c("source", "destination", "t", "individual", "n"),
                       row.names = c(NA, -3L), class = "data.frame")

ns <- NetworkSummary(movements, root = 1, tEnd = '2010-08-16', days = 15)
stopifnot(identical(ns$outDegree, 2L))
