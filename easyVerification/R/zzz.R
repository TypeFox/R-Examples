# zzz.R Check Package Version for Potential Updates
#
#     Copyright (C) 2016 MeteoSwiss
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

.onAttach <- function(...) {
  pkgname <- "easyVerification"
  lib <- system.file(package = pkgname)
  ver <- utils::packageDescription(pkgname)$Version
  desturl <- "https://raw.githubusercontent.com/MeteoSwiss/easyVerification/master/DESCRIPTION"
  con <- tryCatch(RCurl::getURL(desturl, ssl.verifypeer = FALSE), error = function(er) {
    er <- NULL
    return(er)
  })
  if (!is.null(con)) {
    b <- readLines(textConnection(con))
    latest.ver <- package_version(gsub("Version: ", "", b[grep("Version", b)]))
    if (ver < latest.ver) {
      ver.mess1 <- paste("WARNING: Your current version of", pkgname,"is not up-to-date")
      ver.mess <- paste("Get the latest version", latest.ver, 'using install_github("MeteoSwiss/easyVerification")')      
      packageStartupMessage(ver.mess1)
      packageStartupMessage(ver.mess)
    }
  }   
} 
# End
