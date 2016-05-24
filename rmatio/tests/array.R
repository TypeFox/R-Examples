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
## Check write and read of array in MAT5 format:
## 1) without compression
## 2) with compression
##

##
## array: case-1
##
a1.exp <- array(seq_len(32^3), c(32,32,32))
storage.mode(a1.exp) <- 'integer'
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
## array: case-2
##
a2.exp <- array(seq_len(32^3), c(32,32,32))
storage.mode(a2.exp) <- 'double'
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
## array: case-3
##
a3.exp <- array(c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE,
                 TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE,
                 FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
                 FALSE, FALSE, TRUE), c(5L, 5L))
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
## array: case-4
##
a4.exp <- array(seq_len(32^3), c(32,32,32));
storage.mode(a4.exp) <- 'double'
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
## array: case-5
##
a5.exp <- array(seq_len(32^3), c(32,32,32));
storage.mode(a5.exp) <- 'integer'
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
## array: case-6
##
a6.exp <- array(c(seq_len(32767), 32767), c(32,32,32));
storage.mode(a6.exp) <- 'integer'
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
## array: case-7
##
a7.exp <- array(complex(real=seq(1, 2*32^3, 2), imaginary=seq(2, 2*32^3, 2)), c(32,32,32))
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
