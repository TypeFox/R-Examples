#
# Copyright (c) 2011, 2014, Oracle and/or its affiliates. All rights reserved.
#
#    NAME
#      zzz.R - ROCI load and unload operations
#
#    DESCRIPTION
#      Load and unload code for ROCI package.
#
#    NOTES
#
#    MODIFIED   (MM/DD/YY)
#    rpingte     02/27/14 - use utils::globalVariables
#    demukhin    09/11/12 - add Extproc driver
#    demukhin    01/20/12 - cleanup
#    paboyoun    01/04/12 - minor code cleanup
#    demukhin    12/02/11 - add support for more methods
#    demukhin    10/19/11 - Creation
#

utils::globalVariables(c(".oci.GlobalEnv"))

.onLoad <- function(libname, pkgname)
{
  # create ROracle environment
  assign(".oci.GlobalEnv", new.env(parent = emptyenv()), envir = topenv())

  # create driver singleton
  hdl <- .Call("rociDrvAlloc", PACKAGE = "ROracle")
  drv <- new("OraDriver", handle = hdl)
  assign("ora.driver", drv, envir = .oci.GlobalEnv)

  # create extproc driver singleton
  hdl <- .Call("rociDrvAlloc", PACKAGE = "ROracle")
  drv <- new("ExtDriver", handle = hdl)
  assign("ext.driver", drv, envir = .oci.GlobalEnv)
}

# end of file zzz.R
