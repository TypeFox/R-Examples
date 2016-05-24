######################################################################
## Copyright (C) 2006--2008, Roger D. Peng <rpeng@jhsph.edu>
##     
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA
#####################################################################

dumpEnv <- function(env, dbName) {
        keys <- ls(env, all.names = TRUE)
        dumpObjects(list = keys, dbName = dbName, envir = env)
}

dumpImage <- function(dbName = "Rworkspace", type = NULL) {
        dumpObjects(list = ls(envir = globalenv(), all.names = TRUE),
                    dbName = dbName, type = type, envir = globalenv())
}

dumpObjects <- function(..., list = character(0), dbName, type = NULL,
                        envir = parent.frame()) {
        names <- as.character(substitute(list(...)))[-1]
        list <- c(list, names)
        if(!dbCreate(dbName, type))
                stop("could not create database file")
        db <- dbInit(dbName, type)

        for(i in seq(along = list)) 
                dbInsert(db, list[i], get(list[i], envir))
        db
}

dumpDF <- function(data, dbName = NULL, type = NULL) {
        if(is.null(dbName))
                dbName <- as.character(substitute(data))
        dumpList(as.list(data), dbName = dbName, type = type)
}

dumpList <- function(data, dbName = NULL, type = NULL) {
        if(!is.list(data))
                stop("'data' must be a list")
        vnames <- names(data)
        
        if(is.null(vnames) || isTRUE("" %in% vnames))
                stop("list must have non-empty names")
        if(is.null(dbName))
                dbName <- as.character(substitute(data))
        
        if(!dbCreate(dbName, type))
                stop("could not create database file")
        db <- dbInit(dbName, type)

        for(i in seq(along = vnames))
                dbInsert(db, vnames[i], data[[vnames[i]]])
        db
}

