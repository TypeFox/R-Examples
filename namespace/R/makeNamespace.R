#  Copy of function defined inside of base::loadNamespace() in  
#  The R package source at src/library/base/R/namespace.R, 
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# [NB: Other comments and code omitted]

makeNamespace <- function(name, version = NULL, lib = NULL) {
    impenv <- new.env(parent = .BaseNamespaceEnv, hash = TRUE)
    attr(impenv, "name") <- paste("imports", name, sep=":")
    env <- new.env(parent = impenv, hash = TRUE)
    name <- as.character(as.name(name))
    version <- as.character(version)
    info <- new.env(hash = TRUE, parent = baseenv())
    assign(".__NAMESPACE__.", info, envir = env)
    assign("spec", c(name = name,version = version), envir = info)
    setNamespaceInfo(env, "exports", new.env(hash = TRUE, parent = baseenv()))
    dimpenv <- new.env(parent = baseenv(), hash = TRUE)
    attr(dimpenv, "name") <- paste("lazydata", name, sep=":")
    setNamespaceInfo(env, "lazydata", dimpenv)
    setNamespaceInfo(env, "imports", list("base" = TRUE))
    ## this should be an absolute path
    setNamespaceInfo(env, "path",
                     normalizePath(file.path(lib, name), "/", TRUE))
    setNamespaceInfo(env, "dynlibs", NULL)
    setNamespaceInfo(env, "S3methods", matrix(NA_character_, 0L, 3L))
    assign(".__S3MethodsTable__.",
           new.env(hash = TRUE, parent = baseenv()),
           envir = env)
    registerNamespace(name, env)
    env
}
