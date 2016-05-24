## rmatio, a R interface to the C library matio, MAT File I/O Library.
## Copyright (C) 2013-2014  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## rmatio is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(rmatio)

##
## Check write and read of cell arrays in MAT5 format:
## 1) without compression
## 2) with compression
##
## In rmatio, cell arrays are mapped to an unnamed list

##
## Note:
## If the list contains elements of differents lengths i.e.
## a14.in <- list(c("a", "bb"), c("c", "dd"))
## then the expected result is not identical, since each element
## is saved in a cell and the expected result of a14.in is therefore
## a14.exp <- list(list("a", "bb"), list("c", "dd"))
##

##
## cell: case-1
##
a1.exp <- list()
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a1.exp), filename=filename, compression=FALSE, version='MAT5')
a1.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a1.obs)
stopifnot(identical(a1.obs, a1.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a1.exp), filename=filename, compression=TRUE, version='MAT5')
    a1.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a1.zlib.obs)
    stopifnot(identical(a1.zlib.obs, a1.exp))
}

##
## cell: case-2
##
a2.exp <- list(complex(0), logical(0), character(0), numeric(0), integer(0))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a2.exp), filename=filename, compression=FALSE, version='MAT5')
a2.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a2.obs)
stopifnot(identical(a2.obs, a2.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a2.exp), filename=filename, compression=TRUE, version='MAT5')
    a2.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a2.zlib.obs)
    stopifnot(identical(a2.zlib.obs, a2.exp))
}

##
## cell: case-3
##
a3.exp <- list(list(array(c(1, 3, 2, 4), c(2, 2)),
                    array(c(5, 8, 6, 9, 7, 10), c(2,3)),
                    array(c(11, 15, 12, 16, 13, 17, 14, 18), c(2, 4))),
               list(array(c(19, 21, 20, 22), c(2, 2)),
                    array(c(23, 25, 27, 24, 26, 28), c(3L, 2L)),
                    array(c(29, 31, 33, 35, 30, 32, 34, 36), c(4, 2))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a3.exp), filename=filename, compression=FALSE, version='MAT5')
a3.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a3.obs)
stopifnot(identical(a3.obs, a3.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a3.exp), filename=filename, compression=TRUE, version='MAT5')
    a3.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a3.zlib.obs)
    stopifnot(identical(a3.zlib.obs, a3.exp))
}

##
## cell: case-4
##
a4.exp <- list(list(array(c(1L, 3L, 2L, 4L), c(2, 2)),
                   array(c(5L, 8L, 6L, 9L, 7L, 10L), c(2,3)),
                   array(c(11L, 15L, 12L, 16L, 13L, 17L, 14L, 18L), c(2, 4))),
              list(array(c(19L, 21L, 20L, 22L), c(2, 2)),
                   array(c(23L, 25L, 27L, 24L, 26L, 28L), c(3L, 2L)),
                   array(c(29L, 31L, 33L, 35L, 30L, 32L, 34L, 36L), c(4, 2))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a4.exp), filename=filename, compression=FALSE, version='MAT5')
a4.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a4.obs)
stopifnot(identical(a4.obs, a4.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a4.exp), filename=filename, compression=TRUE, version='MAT5')
    a4.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a4.zlib.obs)
    stopifnot(identical(a4.zlib.obs, a4.exp))
}

##
## cell: case-5
##
a5.exp <- list(list(triu(Matrix(1:20, nrow=4, ncol=5, sparse=TRUE)),
                    tril(Matrix(1:20, nrow=5, ncol=4, sparse=TRUE, byrow=TRUE))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a5.exp), filename=filename, compression=FALSE, version='MAT5')
a5.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a5.obs)
stopifnot(identical(a5.obs, a5.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a5.exp), filename=filename, compression=TRUE, version='MAT5')
    a5.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a5.zlib.obs)
    stopifnot(identical(a5.zlib.obs, a5.exp))
}

##
## cell: case-6
##
a6.exp <- list(array(c(1+21i, 0+0i, 0+0i, 0+0i, 5+25i,
                       6+26i, 0+0i, 0+0i, 9+29i, 10+30i, 11+31i, 0+0i,
                       13+33i, 14+34i, 15+35i, 16+36i, 17+37i, 18+38i,
                       19+39i, 20+40i), c(4,5)),
               array(c(1-21i, 5-25i, 9-29i, 13-33i, 17-37i,
                       0+0i, 6-26i, 10-30i, 14-34i, 18-38i, 0+0i, 0+0i,
                       11-31i, 15-35i, 19-39i, 0+0i, 0+0i, 0+0i,
                       16-36i, 20-40i), c(5,4)))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a6.exp), filename=filename, compression=FALSE, version='MAT5')
a6.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a6.obs)
stopifnot(identical(a6.obs, a6.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a6.exp), filename=filename, compression=TRUE, version='MAT5')
    a6.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a6.zlib.obs)
    stopifnot(identical(a6.zlib.obs, a6.exp))
}

##
## cell: case-7
##
a7.exp <- list(list("abcdefghijklmnopqrstuvwxyz",
                    "1234567890!@#$%^&*()-_=+`~"), #
               list("ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                    "[{]}\\|;:'\",<.>/?          "))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a7.exp), filename=filename, compression=FALSE, version='MAT5')
a7.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a7.obs)
stopifnot(identical(a7.obs, a7.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a7.exp), filename=filename, compression=TRUE, version='MAT5')
    a7.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a7.zlib.obs)
    stopifnot(identical(a7.zlib.obs, a7.exp))
}

##
## cell: case-8
##
a8.exp <- list(structure(list(), .Names = character(0)),
               list(),
               structure(list(field1 = numeric(0),
                              field2 = character(0)),
                         .Names = c("field1", "field2")))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a8.exp), filename=filename, compression=FALSE, version='MAT5')
a8.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a8.obs)
stopifnot(identical(a8.obs, a8.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a8.exp), filename=filename, compression=TRUE, version='MAT5')
    a8.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a8.zlib.obs)
    stopifnot(identical(a8.zlib.obs, a8.exp))
}

##
## cell: case-9
##
a9.exp <- list(list(structure(list(
    field1 = list(1, 14),
    field2 = list(
        structure(c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13),
                  .Dim = 3:4),
        structure(c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 ),
                  .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1, 14),
                       field2 = list(structure(c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13),
                           .Dim = 3:4),
                           structure(c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26),
                                     .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1, 14),
                       field2 = list( structure(c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13),
                           .Dim = 3:4),
                           structure(c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 ),
                                     .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1, 14),
                       field2 = list(structure(c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13),
                           .Dim = 3:4),
                           structure(c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1L, 14L),
                       field2 = list( structure(2:13,
                           .Dim = 3:4),
                           structure(15:26,
                                     .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1, 14),
                       field2 = list( structure(c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13),
                           .Dim = 3:4),
                           structure(c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 ),
                                     .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1L, 14L),
                       field2 = list(structure(2:13,
                           .Dim = 3:4),
                           structure(15:26,
                                     .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1L, 14L),
                       field2 = list( structure(2:13, .Dim = 3:4),
                           structure(15:26,
                                     .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1L, 14L),
                       field2 = list( structure(2:13, .Dim = 3:4),
                           structure(15:26,
                                     .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1L, 14L),
                       field2 = list( structure(2:13, .Dim = 3:4),
                           structure(15:26, .Dim = 3:4))),
                             .Names = c("field1", "field2"))),
              list(structure(list(
                  field1 = list(1+51i, 14+64i),
                  field2 = list(structure(c(2+52i, 3+53i, 4+54i, 5+55i, 6+56i,
                      7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                      13+63i),
                      .Dim = 3:4),
                      structure(c(15+65i, 16+66i,
                                  17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                  23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1+51i, 14+64i),
                       field2 = list( structure(c(2+52i, 3+53i, 4+54i, 5+55i,
                           6+56i, 7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                           13+63i),
                           .Dim = 3:4),
                           structure(c(15+65i, 16+66i,
                                       17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                       23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1+51i, 14+64i),
                       field2 = list( structure(c(2+52i, 3+53i, 4+54i, 5+55i,
                           6+56i, 7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                           13+63i), .Dim = 3:4), structure(c(15+65i, 16+66i,
                                        17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                        23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1+51i, 14+64i),
                       field2 = list( structure(c(2+52i, 3+53i, 4+54i, 5+55i,
                           6+56i, 7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                           13+63i), .Dim = 3:4), structure(c(15+65i, 16+66i,
                                        17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                        23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1+51i, 14+64i),
                       field2 = list( structure(c(2+52i, 3+53i, 4+54i, 5+55i,
                           6+56i, 7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                           13+63i), .Dim = 3:4), structure(c(15+65i, 16+66i,
                                        17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                        23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1+51i, 14+64i),
                       field2 = list( structure(c(2+52i, 3+53i, 4+54i, 5+55i,
                           6+56i, 7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                           13+63i), .Dim = 3:4), structure(c(15+65i, 16+66i,
                                        17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                        23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1+51i, 14+64i),
                       field2 = list( structure(c(2+52i, 3+53i, 4+54i, 5+55i,
                           6+56i, 7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                           13+63i), .Dim = 3:4), structure(c(15+65i, 16+66i,
                                        17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                        23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1+51i, 14+64i),
                       field2 = list( structure(c(2+52i, 3+53i, 4+54i, 5+55i,
                           6+56i, 7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                           13+63i), .Dim = 3:4), structure(c(15+65i, 16+66i,
                                        17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                        23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1+51i, 14+64i),
                       field2 = list( structure(c(2+52i, 3+53i, 4+54i, 5+55i,
                           6+56i, 7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                           13+63i), .Dim = 3:4), structure(c(15+65i, 16+66i,
                                        17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                        23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2")),
                   structure(list(
                       field1 = list(1+51i, 14+64i),
                       field2 = list( structure(c(2+52i, 3+53i, 4+54i, 5+55i,
                           6+56i, 7+57i, 8+58i, 9+59i, 10+60i, 11+61i, 12+62i,
                           13+63i), .Dim = 3:4), structure(c(15+65i, 16+66i,
                                        17+67i, 18+68i, 19+69i, 20+70i, 21+71i, 22+72i,
                                        23+73i, 24+74i, 25+75i, 26+76i), .Dim = 3:4))),
                             .Names = c("field1", "field2"))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a9.exp), filename=filename, compression=FALSE, version='MAT5')
a9.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a9.obs)
stopifnot(identical(a9.obs, a9.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a9.exp), filename=filename, compression=TRUE, version='MAT5')
    a9.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a9.zlib.obs)
    stopifnot(identical(a9.zlib.obs, a9.exp))
}

##
## cell: case-10
##
a10.exp <- list(list(field1=list(triu(Matrix(1:20, nrow=4,
                         ncol=5, sparse=TRUE))),
                     field2=list(tril(Matrix(1:20, nrow=5, ncol=4,
                         sparse=TRUE, byrow=TRUE)))),
                list(field1=list(array(c(1+21i, 0+0i,
                         0+0i, 0+0i, 5+25i, 6+26i, 0+0i, 0+0i, 9+29i,
                         10+30i, 11+31i, 0+0i, 13+33i, 14+34i, 15+35i,
                         16+36i, 17+37i, 18+38i, 19+39i, 20+40i),
                         c(4,5))),
                     field2=list(array(c(1-21i, 5-25i,
                         9-29i, 13-33i, 17-37i, 0+0i, 6-26i, 10-30i,
                         14-34i, 18-38i, 0+0i, 0+0i, 11-31i, 15-35i,
                         19-39i, 0+0i, 0+0i, 0+0i, 16-36i, 20-40i),
                         c(5,4)))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a10.exp), filename=filename, compression=FALSE, version='MAT5')
a10.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a10.obs)
stopifnot(identical(a10.obs, a10.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a10.exp), filename=filename, compression=TRUE, version='MAT5')
    a10.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a10.zlib.obs)
    stopifnot(identical(a10.zlib.obs, a10.exp))
}

##
## cell: case-11
##
a11.exp <- list(list(field1 = "abcdefghijklmnopqrstuvwxyz",
                     field2 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"),
                list(field1 = "1234567890!@#$%^&*()-_=+`~", #
                     field2 = "[{]}\\|;:'\",<.>/?          "))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a11.exp), filename=filename, compression=FALSE, version='MAT5')
a11.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a11.obs)
stopifnot(identical(a11.obs, a11.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a11.exp), filename=filename, compression=TRUE, version='MAT5')
    a11.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a11.zlib.obs)
    stopifnot(identical(a11.zlib.obs, a11.exp))
}

##
## cell: case-12
##
a12.exp <- list(structure(c(FALSE, TRUE, FALSE, TRUE, FALSE,
                            TRUE, FALSE, TRUE, FALSE, TRUE,
                            FALSE, TRUE, FALSE, TRUE, FALSE,
                            TRUE, FALSE, TRUE, FALSE, TRUE),
                          .Dim = 4:5),
                structure(c(TRUE, FALSE, TRUE,
                            FALSE, TRUE, FALSE, TRUE, FALSE,
                            TRUE, FALSE, TRUE, FALSE, TRUE,
                            FALSE, TRUE, FALSE, TRUE, FALSE,
                            TRUE, FALSE),
                          .Dim = 4:5),
                structure(c(TRUE, TRUE, TRUE,
                            TRUE, TRUE, FALSE, TRUE, TRUE,
                            TRUE, TRUE, FALSE, FALSE, TRUE,
                            TRUE, TRUE, FALSE, FALSE, FALSE,
                            TRUE, TRUE, FALSE, FALSE, FALSE,
                            FALSE, TRUE ),
                          .Dim = c(5L, 5L)),
                structure(c(TRUE, FALSE,
                            FALSE, FALSE, FALSE, TRUE, TRUE,
                            FALSE, FALSE, FALSE, TRUE, TRUE,
                            TRUE, FALSE, FALSE, TRUE, TRUE,
                            TRUE, TRUE, FALSE, TRUE, TRUE,
                            TRUE, TRUE, TRUE),
                          .Dim = c(5L, 5L)))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a12.exp), filename=filename, compression=FALSE, version='MAT5')
a12.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a12.obs)
stopifnot(identical(a12.obs, a12.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a12.exp), filename=filename, compression=TRUE, version='MAT5')
    a12.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a12.zlib.obs)
    stopifnot(identical(a12.zlib.obs, a12.exp))
}

##
## cell: case-13
##
a13.exp <- list(structure(list(),
                          .Names = character(0)),
                list(field1=list(), field2=list()),
                structure(list(field1 = numeric(0),
                               field2 = character(0)),
                          .Names = c("field1", "field2")))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a13.exp), filename=filename, compression=FALSE, version='MAT5')
a13.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a13.obs)
stopifnot(identical(a13.obs, a13.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a13.exp), filename=filename, compression=TRUE, version='MAT5')
    a13.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a13.zlib.obs)
    stopifnot(identical(a13.zlib.obs, a13.exp))
}

##
## cell: case-14
##
a14.in <- list(c("a", "bb"), c("c", "dd"))
a14.exp <- list(list("a", "bb"), list("c", "dd"))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a14.in), filename=filename, compression=FALSE, version='MAT5')
a14.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a14.obs)
stopifnot(identical(a14.obs, a14.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a14.in), filename=filename, compression=TRUE, version='MAT5')
    a14.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a14.zlib.obs)
    stopifnot(identical(a14.zlib.obs, a14.exp))
}

##
## cell: case-15
##
a15.in <- list(c("a", "bb"), list(c("d", "eee")))
a15.exp <- list(list("a", "bb"), list(list(list("d", "eee"))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a15.in), filename=filename, compression=FALSE, version='MAT5')
a15.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a15.obs)
stopifnot(identical(a15.obs, a15.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a15.in), filename=filename, compression=TRUE, version='MAT5')
    a15.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a15.zlib.obs)
    stopifnot(identical(a15.zlib.obs, a15.exp))
}

##
## cell: case-16
##
a16.in <- list(c("a", "bb"), Matrix(c(0, 0, 0, 0, 0, 0, 1, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0, 1, 0,
                                      0, 0, 0, 0, 0, 0, 0, 0, 1),
                                    nrow=3, ncol=9, byrow=TRUE, sparse=TRUE))
a16.exp <- list(list("a", "bb"), list(a16.in[[2]]))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a16.in), filename=filename, compression=FALSE, version='MAT5')
a16.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a16.obs)
stopifnot(identical(a16.obs, a16.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a16.in), filename=filename, compression=TRUE, version='MAT5')
    a16.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a16.zlib.obs)
    stopifnot(identical(a16.zlib.obs, a16.exp))
}

##
## cell: case-17
##
a17.in <- list(list(c("a", "bb")), list())
a17.exp <- list(list(list(list("a", "bb"))), list())
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a17.in), filename=filename, compression=FALSE, version='MAT5')
a17.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a17.obs)
stopifnot(identical(a17.obs, a17.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a17.in), filename=filename, compression=TRUE, version='MAT5')
    a17.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a17.zlib.obs)
    stopifnot(identical(a17.zlib.obs, a17.exp))
}
