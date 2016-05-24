# Copyright (C) Tal Galili
#
# This file is part of dendextendRcpp.
#
# dendextend is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# dendextend is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


# get("old_get_branches_heights", envir = asNamespace("dendextend"), inherits = T)
# dendextend:::old_get_branches_heights
# dendextendRcpp:::old_get_branches_heights
# sys.frames(asNamespace("dendextend"))
# 
# environment(old_get_branches_heights)=as.environment("package:dendextendRcpp")
# dendextendRcpp:::old_get_branches_heights

assign_dendextendRcpp_to_dendextend <- function() {
   # assigns the FASTER dendextendRcpp functions to override
   # the dendextend functions....
   
#   if(TRUE) { #  (suppressWarnings(require(dendextend))) {   # no need to check this since dendextendRcpp Depends on dendextend!
      # This wouldn't work since it will only assign
      # the faster function in the current env      
      #       get_branches_heights <- dendextendRcpp::get_branches_heights
      #       heights_per_k.dendrogram <- dendextendRcpp::heights_per_k.dendrogram
      # for getting the functions "into" dendextend, we need to run this:
      
      # create a backup of these functions in order to later
      # compare them using benchmark (their kept invisible - but can be accessed)
      
      # I just ended up making a manual copy of these objects...
#       assign("old_get_branches_heights", get_branches_heights, envir=as.environment("package:dendextendRcpp"))
#       assign("old_heights_per_k.dendrogram", heights_per_k.dendrogram, envir=as.environment("package:dendextendRcpp"))
#       assign("old_cut_lower_fun", cut_lower_fun, envir=as.environment("package:dendextendRcpp"))

      # Using this:
      # assign("old_cut_lower_fun", dendextend::cut_lower_fun, envir=as.environment("package:dendextendRcpp"))
      # makes it impossible to lose the <environment: namespace:dendextend>
      # While it is actually NOT available through namespace:dendextend!!!
      # I guess it is a bug (maybe).

# 
# environment(old_cut_lower_fun)=.GlobalEnv
# 
#    as.environment("package:dendextendRcpp")
# dendextendRcpp::old_cut_lower_fun
# old_cut_lower_fun

      # require(utils) # doesn't help really...
      # but this does: (!)
      # http://stackoverflow.com/questions/13595145/overriding-a-package-function-inherited-by-another-package
      # 	  get("assignInNamespace", envir=asNamespace("utils"))
      # Using only "::" instead of ":::" will crash many tests...

      
#       # instead of assigining to the NAMESPACE, I'm moving to have the functions
#       # reside in the "options" world...
#       assignInNamespace(
#          x= "get_branches_heights",
#          value = dendextendRcpp:::get_branches_heights,
#          ns = "dendextend"
#       )   
#       assignInNamespace(
#          x= "heights_per_k.dendrogram",
#          value = dendextendRcpp:::heights_per_k.dendrogram,
#          ns = "dendextend"
#       )
#       assignInNamespace(
#          x= "cut_lower_fun",
#          value = dendextendRcpp:::cut_lower_fun,
#          ns = "dendextend"
#       )

   # options(dendextend_get_branches_heights = dendextendRcpp::get_branches_heights)
   # options(dendextend_heights_per_k.dendrogram = dendextendRcpp::heights_per_k.dendrogram)
   # options(dendextend_cut_lower_fun = dendextendRcpp::cut_lower_fun)

   dendextend::dendextend_options("get_branches_heights" , dendextendRcpp::dendextendRcpp_get_branches_heights)
   dendextend::dendextend_options("heights_per_k.dendrogram" , dendextendRcpp::dendextendRcpp_heights_per_k.dendrogram)
   dendextend::dendextend_options("cut_lower_fun" , dendextendRcpp::dendextendRcpp_cut_lower_fun)
   dendextend::dendextend_options("labels.dendrogram" , dendextendRcpp::dendextendRcpp_labels.dendrogram)
   
   # dendextend_options()
   
   # http://stackoverflow.com/questions/21921304/r-package-development-overriding-a-function-from-one-package-with-a-function-fr
   # http://stackoverflow.com/questions/15910113/what-where-are-the-attributes-of-a-function-object
   # http://stackoverflow.com/questions/8501906/define-a-function-in-a-specific-namespace
   
   
   
   
# dendextend:::cut_lower_fun
# Good! this works...

      ## p.s:
      # doing the following is a BAD IDEA!
      # This will not allow us to use labels.dendrogram when our Rcpp version fails...
      # assignInNamespace(
      #    x= "labels.dendrogram",
      #    value = dendextendRcpp:::labels.dendrogram,
      #    ns = "stats"
      #    )
      
      
#       
#    } else {
# #       warning("
# #               The 'dendextend' package runs 
# #               MUCH faster when you also have the dendextendRcpp package installed.
# #               Please consider running:
# #               install.packages('dendextendRcpp')
# #               and then re-load dendextend.
# #               ")
#    }
   
   }






remove_dendextendRcpp_options <- function() { 
   # assigns the functions which could later be replaced by the FASTER dendextendRcpp functions 
#    options(dendextend_get_branches_heights = NULL)
#    options(dendextend_heights_per_k.dendrogram = NULL)
#    options(dendextend_cut_lower_fun = NULL)
   dendextend::assign_dendextend_options()
}
   
   
   

.onLoad <- function(libname, pkgname){
   # Thanks for Romain: http://stackoverflow.com/questions/4369334/first-lib-idiom-in-r-packages
   
   # adding and removing menus from the Rgui when loading and detaching the library
   # setHook(packageEvent("installr", "attach"), {function(pkgname, libpath) {add.installr.GUI()}  } )
   # setHook(packageEvent("installr", "detach"), {function(pkgname, libpath) {remove.installr.GUI()}  } )

   setHook(packageEvent("dendextendRcpp", "detach"), {function(pkgname, libpath) {remove_dendextendRcpp_options()}  } )
   
}

# menus are added and removed as needed: !!


.onAttach <- function(lib, pkg,...){
   ####
   # packageStartupMessage("Consider running: labels.dendrogram <- dendextendRcpp_labels.dendrogram if you wish to gain extra speed on your own dendrogram-related functions")
   
   # look at the top of this file.
   assign_dendextendRcpp_to_dendextend()
   
   
}



# ' @import dendextend
# ' @import Rcpp




# devtools::use_travis()









# detach( 'package:dendextendRcpp', unload=TRUE )
# require( 'dendextendRcpp' )






# Steps:
# http://r.789695.n4.nabble.com/vignettes-problems-with-PDF-compaction-td4664909.html
# 1) install gs - http://www.ghostscript.com/download/gsdnld.html
# 2) find where it is, and update R_GSCMD:
# Sys.setenv(R_GSCMD="C:\\Program Files\\gs\\gs9.10\\bin\\gswin64c.exe")
# Sys.getenv("R_GSCMD")
# 3) Check that it works: 
# system2(Sys.getenv("R_GSCMD"), args="--version") 
# 4) use:
#### tools::compactPDF("inst\\doc\\dendextend-tutorial.pdf", gs_quality="ebook") 
# tools::compactPDF("inst\\doc\\dendextend-tutorial.pdf") 
###   compacted 'dendextend-tutorial.pdf' from 961Kb to 737Kb

# For checking:
# 1) get qpdf
#     http://sourceforge.net/projects/qpdf/files/
# 2) put it somewhere
# 3) set R_QPDF
#  Sys.setenv(R_QPDF="C:\\Rtools\\qpdf-5.1.1\\bin\\qpdf.exe")
#  Sys.which(Sys.getenv("R_QPDF", "qpdf"))

# Also, make sure to add:
# options(repos=c("http://cran.rstudio.com", "http://www.stats.ox.ac.uk/pub/RWin" ))
# to D:\R\R-devel\etc\Rprofile.site

# when a function is renamed, its document in man must be removed - otherwise it may cause problems with the built check (it will try to run the code in the example, and will fail.)
# When all is done, run:
# require(devtools)
# check(document=FALSE)
# browseURL(tempdir())
### http://www.rstudio.com/ide/docs/packages/build_options
# check(build_args="--no-build-vignettes --no-manual", args = "--no-examples --no-build-vignettes --no-manual",  cran = FALSE, cleanup = FALSE)
# check(build_args="--no-build-vignettes ", args = "--no-build-vignettes",  cran = FALSE, cleanup = FALSE)
# 
# check(args="--as-cran",document=FALSE) # I need to not check the documents since this seem to force a NAMESPACE change...
# require(devtools)
# check("D:/Dropbox/aaaa good R code/AA - My packages/dendextendRcpp", args="--as-cran",document=FALSE)
#                 Thanks to: http://stackoverflow.com/questions/10017702/r-cmd-check-options-for-more-rigorous-testing-2-15-0


# file.copy("NEWS", "NEWS.md")
# shell('git log --graph --stat --date=short --pretty=format:"%ad(%an) %s |%h" > ChangeLog', intern = TRUE)




#
#
# require(devtools)
# check(args="--as-cran",document=FALSE) # I need to not check the documents since this seem to force a NAMESPACE change...
# devtools::build_win(version="R-devel")
# release(check = FALSE)
# release()
