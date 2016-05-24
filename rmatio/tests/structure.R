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
## Check write and read of structure arrays in MAT5 format:
## 1) without compression
## 2) with compression
##
## In rmatio, structure arrays are mapped to a named list

##
## Note:
## If the list contains elements of differents lengths i.e.
## a14.in <- list(y=c("a", "bb"), z=c(1, 2))
## then the expected result is not identical, since each element
## is saved in a cell and the expected result of a14.in is therefore
## a14.exp <- list(y = list("a", "bb"), z = list(c(1, 2)))
##

##
## structure: case-1 (Empty structure array)
##
a1.exp <- structure(list(), .Names = character(0))
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
## structure: case-2 (Empty structure array with fields)
##
a2.exp <- list(field1=list(), field2=list())
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
## structure: case-3 (Structure array with empty fields)
##
a3.exp <- list(field1=numeric(0), field2=character(0), field3=complex(0), filed4=integer(0), field5=logical(0))
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
## structure: case-4
##
a4.exp <- list(field1=list(1, 14),
               field2=list(array(as.numeric(2:13), c(3,4)),
                   array(as.numeric(15:26), c(3,4))))
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
## structure: case-5
##
a5.exp <- list(field1=list(1L, 14L),
               field2=list(array(2:13, c(3,4)),
                   array(15:26, c(3,4))))
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
## structure: case-6
##
a6.exp <- list(field1=list(1+51i, 14+64i),
               field2=list(array(c(2+52i, 3+53i, 4+54i, 5+55i, 6+56i, 7+57i,
                   8+58i, 9+59i, 10+60i, 11+61i, 12+62i, 13+63i),
                   c(3,4)),
                   array(c(15+65i, 16+66i, 17+67i, 18+68i, 19+69i, 20+70i,
                           21+71i, 22+72i, 23+73i, 24+74i, 25+75i, 26+76i),
                         c(3,4))))
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
## structure: case-7
##
a7.exp <- list(field1=list(triu(Matrix(1:20, nrow=4, ncol=5, sparse=TRUE))),
               field2=list(tril(Matrix(1:20, nrow=5, ncol=4, sparse=TRUE, byrow=TRUE))))
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
## structure: case-8
##
a8.exp <- list(field1=list(array(c(1+21i, 0+0i, 0+0i, 0+0i, 5+25i,
                   6+26i, 0+0i, 0+0i, 9+29i, 10+30i, 11+31i, 0+0i,
                   13+33i, 14+34i, 15+35i, 16+36i, 17+37i, 18+38i,
                   19+39i, 20+40i), c(4,5))),
               field2=list(array(c(1-21i, 5-25i, 9-29i, 13-33i, 17-37i,
                   0+0i, 6-26i, 10-30i, 14-34i, 18-38i, 0+0i, 0+0i,
                   11-31i, 15-35i, 19-39i, 0+0i, 0+0i, 0+0i,
                   16-36i, 20-40i), c(5,4))))
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
## structure: case-9
##
a9.exp <- list(field1 = c("abcdefghijklmnopqrstuvwxyz",
                   "1234567890!@#$%^&*()-_=+`~"), #
               field2 = c("ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                   "[{]}\\|;:'\",<.>/?          "))
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
## structure: case-10 (Structure array with empty fields)
##
a10.exp <- list(field1=numeric(0))
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
    stopifnot(identical(a10.zlib.obs, a10.exp))
}

##
## structure: case-11
##
a11.exp <- list(field1=list(1))
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
## structure: case-12
##
a12.exp <- structure(list(field1 = list(structure(c(FALSE, TRUE,
                           FALSE, TRUE, FALSE, TRUE, FALSE, TRUE,
                           FALSE, TRUE, FALSE, TRUE, FALSE, TRUE,
                           FALSE, TRUE, FALSE, TRUE, FALSE,
                           TRUE), .Dim = 4:5),
                           structure(c(TRUE, TRUE, TRUE, TRUE, TRUE,
                                       FALSE, TRUE, TRUE, TRUE, TRUE,
                                       FALSE, FALSE, TRUE, TRUE, TRUE,
                                       FALSE, FALSE, FALSE, TRUE, TRUE,
                                       FALSE, FALSE, FALSE, FALSE, TRUE),
                                     .Dim = c(5L, 5L))),
                       field2 = list(structure(c(TRUE, FALSE, TRUE, FALSE, TRUE,
                           FALSE, TRUE, FALSE, TRUE, FALSE,
                           TRUE, FALSE, TRUE, FALSE, TRUE,
                           FALSE, TRUE, FALSE, TRUE, FALSE),
                           .Dim = 4:5),
                           structure(c(TRUE, FALSE,
                                       FALSE, FALSE, FALSE, TRUE, TRUE,
                                       FALSE, FALSE, FALSE, TRUE, TRUE, TRUE,
                                       FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,
                                       FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
                                     .Dim = c(5L, 5L)))),
                  .Names = c("field1", "field2"))
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
## structure: case-13
##
a13.exp <- structure(list(X = structure(list(x = list(structure(c(1, 4,
                                              2, 5, 3, 6.2),
                                              .Dim = 2:3)),
                           y = list(list("Hello", "world!")),
                           z = list(structure(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
                               .Dim = c(14L, 14L)))),
                           .Names = c("x", "y", "z"))),
                  .Names = "X")
filename <- tempfile(fileext = ".mat")
write.mat(a13.exp, filename=filename, compression=FALSE, version='MAT5')
a13.obs <- read.mat(filename)
unlink(filename)
str(a13.obs)
stopifnot(identical(a13.obs, a13.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(a13.exp, filename=filename, compression=TRUE, version='MAT5')
    a13.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a13.zlib.obs)
    stopifnot(identical(a13.zlib.obs, a13.exp))
}

##
## structure: case-14
##
a14.in <- list(y=c("a", "bb"), z=c(1, 2))
a14.exp <- list(y = list("a", "bb"), z = list(c(1, 2)))
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
## structure: case-15
##
a15.in <- list(y=c("a", "bb"), z=c("c", "dd"))
a15.exp <- list(y = list("a", "bb"), z = list("c", "dd"))
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
## structure: case-16
##
a16.in <- list(y=c("a", "b"), z=c(1, 2))
a16.exp <- list(a = list(y = list("a", "b"), z = list(c(1, 2))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a16.in), filename=filename, compression=FALSE, version='MAT5')
a16.obs <- read.mat(filename)
unlink(filename)
str(a16.obs)
stopifnot(identical(a16.obs, a16.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a16.in), filename=filename, compression=TRUE, version='MAT5')
    a16.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a16.zlib.obs)
    stopifnot(identical(a16.zlib.obs, a16.exp))
}

##
## structure: case-17
##
a17.in <- list(y=c("a", "b"), z=c(1, 2, 3))
a17.exp <- list(a = list(y = list("a", "b"), z = list(c(1, 2, 3))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a17.in), filename=filename, compression=FALSE, version='MAT5')
a17.obs <- read.mat(filename)
unlink(filename)
str(a17.obs)
stopifnot(identical(a17.obs, a17.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a17.in), filename=filename, compression=TRUE, version='MAT5')
    a17.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a17.zlib.obs)
    stopifnot(identical(a17.zlib.obs, a17.exp))
}

##
## structure: case-18
##
a18.in <- list(y=c("a", "bb"), z=c(1, 2, 3))
a18.exp <- list(a = list(y = list("a", "bb"), z = list(c(1, 2, 3))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a18.in), filename=filename, compression=FALSE, version='MAT5')
a18.obs <- read.mat(filename)
unlink(filename)
str(a18.obs)
stopifnot(identical(a18.obs, a18.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a18.in), filename=filename, compression=TRUE, version='MAT5')
    a18.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a18.zlib.obs)
    stopifnot(identical(a18.zlib.obs, a18.exp))
}

##
## structure: case-19
##
a19.in <- list(y=c("a", "bb"), z=c(1, 2))
a19.exp <- list(y = list('a', 'bb'), z = c(1, 2))
filename <- tempfile(fileext = ".mat")
write.mat(a19.in, filename=filename, compression=FALSE, version='MAT5')
a19.obs <- read.mat(filename)
unlink(filename)
str(a19.obs)
stopifnot(identical(a19.obs, a19.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(a19.in, filename=filename, compression=TRUE, version='MAT5')
    a19.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a19.zlib.obs)
    stopifnot(identical(a19.zlib.obs, a19.exp))
}

##
## structure: case-20
##
a20.in <- list(y=c("a", "bb"), z=c("c", "dd"))
a20.exp <- list(y = list("a", "bb"), z = list("c", "dd"))
filename <- tempfile(fileext = ".mat")
write.mat(a20.in, filename=filename, compression=FALSE, version='MAT5')
a20.obs <- read.mat(filename)
unlink(filename)
str(a20.obs)
stopifnot(identical(a20.obs, a20.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(a20.in, filename=filename, compression=TRUE, version='MAT5')
    a20.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a20.zlib.obs)
    stopifnot(identical(a20.zlib.obs, a20.exp))
}

##
## structure: case-21
##
a21.exp <- list(y=c("a", "b"), z=c(1, 2))
filename <- tempfile(fileext = ".mat")
write.mat(a21.exp, filename=filename, compression=FALSE, version='MAT5')
a21.obs <- read.mat(filename)
unlink(filename)
str(a21.obs)
stopifnot(identical(a21.obs, a21.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(a21.exp, filename=filename, compression=TRUE, version='MAT5')
    a21.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a21.zlib.obs)
    stopifnot(identical(a21.zlib.obs, a21.exp))
}

##
## structure: case-22
##
a22.exp <- list(y=c("a", "b"), z=c(1, 2, 3))
filename <- tempfile(fileext = ".mat")
write.mat(a22.exp, filename=filename, compression=FALSE, version='MAT5')
a22.obs <- read.mat(filename)
unlink(filename)
str(a22.obs)
stopifnot(identical(a22.obs, a22.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(a22.exp, filename=filename, compression=TRUE, version='MAT5')
    a22.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a22.zlib.obs)
    stopifnot(identical(a22.zlib.obs, a22.exp))
}

##
## structure: case-23
##
a23.in <- list(y=c("a", "bb"), z=c(1, 2, 3))
a23.exp <- list(y = list("a", "bb"), z = c(1, 2, 3))
filename <- tempfile(fileext = ".mat")
write.mat(a23.in, filename=filename, compression=FALSE, version='MAT5')
a23.obs <- read.mat(filename)
unlink(filename)
str(a23.obs)
stopifnot(identical(a23.obs, a23.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(a23.in, filename=filename, compression=TRUE, version='MAT5')
    a23.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a23.zlib.obs)
    stopifnot(identical(a23.zlib.obs, a23.exp))
}

##
## structure: case-24
##
a24.in <- list(y=c("a", "bb"), z=list(c('d', 'eee')))
a24.exp <- list(a = list(y = list("a", "bb"), z = list(list(list("d", "eee")))))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a24.in), filename=filename, compression=FALSE, version='MAT5')
a24.obs <- read.mat(filename)
unlink(filename)
str(a24.obs)
stopifnot(identical(a24.obs, a24.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a24.in), filename=filename, compression=TRUE, version='MAT5')
    a24.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a24.zlib.obs)
    stopifnot(identical(a24.zlib.obs, a24.exp))
}

##
## structure: case-25
##
a25.in <- list(y=c("a", "bb"), z=list(c('d', 'eee')))
a25.exp <- list(y = list("a", "bb"), z = list(list("d", "eee")))
filename <- tempfile(fileext = ".mat")
write.mat(a25.in, filename=filename, compression=FALSE, version='MAT5')
a25.obs <- read.mat(filename)
unlink(filename)
str(a25.obs)
stopifnot(identical(a25.obs, a25.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(a25.in, filename=filename, compression=TRUE, version='MAT5')
    a25.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a25.zlib.obs)
    stopifnot(identical(a25.zlib.obs, a25.exp))
}

##
## structure: case-26
##
a26.in <- list(y=c("a", "bb"), z=Matrix(c(0, 0, 0, 0, 0, 0, 1, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 1, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 1),
                                   nrow=3, ncol=9, byrow=TRUE, sparse=TRUE))
a26.exp <- list(y=list("a", "bb"), z=a26.in$z)
filename <- tempfile(fileext = ".mat")
write.mat(a26.in, filename=filename, compression=FALSE, version='MAT5')
a26.obs <- read.mat(filename)
unlink(filename)
str(a26.obs)
stopifnot(identical(a26.obs, a26.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(a26.in, filename=filename, compression=TRUE, version='MAT5')
    a26.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a26.zlib.obs)
    stopifnot(identical(a26.zlib.obs, a26.exp))
}

##
## structure: case-27
##
a27.in <- list(y=c("a", "bb"), z=Matrix(c(0, 0, 0, 0, 0, 0, 1, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 1, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 1),
                                   nrow=3, ncol=9, byrow=TRUE, sparse=TRUE))
a27.exp <- list(y=list("a", "bb"),z=list(a27.in$z))
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a27.in), filename=filename, compression=FALSE, version='MAT5')
a27.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a27.obs)
stopifnot(identical(a27.obs, a27.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a27.in), filename=filename, compression=TRUE, version='MAT5')
    a27.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a27.zlib.obs)
    stopifnot(identical(a27.zlib.obs, a27.exp))
}

##
## structure: case-28
##
a28.in <- list(y=list(c("a", "bb")), z=list())
a28.exp <- list(y=list(list(list("a", "bb"))), z=list())
filename <- tempfile(fileext = ".mat")
write.mat(list(a=a28.in), filename=filename, compression=FALSE, version='MAT5')
a28.obs <- read.mat(filename)[['a']]
unlink(filename)
str(a28.obs)
stopifnot(identical(a28.obs, a28.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(list(a=a28.in), filename=filename, compression=TRUE, version='MAT5')
    a28.zlib.obs <- read.mat(filename)[['a']]
    unlink(filename)
    str(a28.zlib.obs)
    stopifnot(identical(a28.zlib.obs, a28.exp))
}

##
## structure: case-29
##
a29.in <- list(y=list(c("a", "bb")), z=list())
a29.exp <- list(y = list(list("a", "bb")), z = list())
filename <- tempfile(fileext = ".mat")
write.mat(a29.in, filename=filename, compression=FALSE, version='MAT5')
a29.obs <- read.mat(filename)
unlink(filename)
str(a29.obs)
stopifnot(identical(a29.obs, a29.exp))

## Run the same test with compression
if(rmatio:::have.zlib()) {
    filename <- tempfile(fileext = ".mat")
    write.mat(a29.in, filename=filename, compression=TRUE, version='MAT5')
    a29.zlib.obs <- read.mat(filename)
    unlink(filename)
    str(a29.zlib.obs)
    stopifnot(identical(a29.zlib.obs, a29.exp))
}
