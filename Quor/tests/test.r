## File Quor/tests/test.r
##
## Quor package for R (http://www.R-project.org)
## Copyright (C) 2014 Adriano Polpo, Carlos A. de B. Pereira, Cassio P. de Campos.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## This code is used for testing purposes. The Quor library does not
## depend on it for any of its functionalities

############################################
source('loadlibs.r')
loadlibs(forcelocal=TRUE)

## Evaluateing confidence statement for more 
## than 2 groups with normal distributions.
set.seed(42)
n <- c(100,100,100,100)
x.m <- c(1,2,3,4)
data <- NULL
data$x1 <- sort(rnorm(n[1],x.m[1],1))
data$x2 <- sort(rnorm(n[2],x.m[2],1))
data$x3 <- sort(rnorm(n[2],x.m[3],1))
data$x4 <- sort(rnorm(n[2],x.m[4],1))
res=conf.statement(data)

############################################
## Evaluateing confidence statement for
## 2 groups from schizophrenia data set.
if(FALSE) {
    load('../../data/schizophrenia.rda') ## the license for this data set does not
                                         ## clearly allow us to give it together with 
                                         ## Quor package, so we prefer not to do it.
                                         ## Please refer to
                                         ## http://www.stanleygenomics.org and
                                         ## http://dx.doi.org/10.1186/1471-2164-7-70
    source('loadlibs.r')
    loadlibs(forcelocal=TRUE)
    data <- NULL
    data[[1]] <- t(as.matrix(a[a[,1]==1,2:20993]))
    data[[2]] <- t(as.matrix(a[a[,1]==2,2:20993]))
    res=conf.statement(data)
    exp(res$confidence[,1:100])
}
