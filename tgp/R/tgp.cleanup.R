#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************


## tgp.cleanup
##
## gets called when the C-side is aborted by the R-side and enables
## the R-side to clean up the memory still allocaed to the C-side,
## as well as whatever files were left open on the C-side

"tgp.cleanup" <-
  function(message="INTERRUPT", verb, rmfiles=TRUE)
{
  .C("tgp_cleanup", PACKAGE = "tgp")

  ## remove the trace (and other) files?
  if(rmfiles) {
  
    if(file.exists(paste("./", "best_parts_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed best_parts_1.out\n", sep=""))
      unlink("best_parts_1.out")
    } 
    
    if(file.exists(paste("./", "tree_m0_posts.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed tree_m0_posts.out\n", sep=""))
      unlink("tree_m0_posts.out")
    }
    
    if(file.exists(paste("./", "trace_parts_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_parts_1.out\n", sep=""))
      unlink("trace_parts_1.out")
    } 
    
    if(file.exists(paste("./", "trace_post_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_post_1.out\n", sep=""))
      unlink("trace_post_1.out")
    }
    
    if(file.exists(paste("./", "trace_wlambda_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_wlambda_1.out\n", sep=""))
      unlink("trace_wlambda_1.out")
    } 
    
    if(file.exists(paste("./", "trace_hier_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_hier_1.out\n", sep=""))
      unlink("trace_hier_1.out")
    } 
    
    if(file.exists(paste("./", "trace_linarea_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_linarea_1.out\n", sep=""))
      unlink("trace_linarea_1.out")
    } 
    
    if(file.exists(paste("./", "trace_XX_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_XX_1.out\n", sep=""))
      unlink("trace_XX_1.out")
    } 

    if(file.exists(paste("./", "trace_Zp_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_Zp_1.out\n", sep=""))
      unlink("trace_Zp_1.out")
    }

    if(file.exists(paste("./", "trace_Zpkm_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_Zpkm_1.out\n", sep=""))
      unlink("trace_Zpkm_1.out")
    }
    
    if(file.exists(paste("./", "trace_Zpks2_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_Zpks2_1.out\n", sep=""))
      unlink("trace_Zpks2_1.out")
    }
    
    if(file.exists(paste("./", "trace_ZZ_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_ZZ_1.out\n", sep=""))
      unlink("trace_ZZ_1.out")
    }

    if(file.exists(paste("./", "trace_ZZkm_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_ZZkm_1.out\n", sep=""))
      unlink("trace_ZZkm_1.out")
    }
    
    if(file.exists(paste("./", "trace_ZZks2_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_ZZks2_1.out\n", sep=""))
      unlink("trace_ZZks2_1.out")
    }
    
    if(file.exists(paste("./", "trace_improv_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_improv_1.out\n", sep=""))
      unlink("trace_improv_1.out")
    }

    if(file.exists(paste("./", "trace_Ds2x_1.out", sep=""))) {
      if(verb >= 1) cat(paste(message, ": removed trace_Ds2x_1.out\n", sep=""))
      unlink("trace_Ds2x_1.out")
    }

    ## get all of the names of the tree files
    tree.files <- list.files(pattern="tree_m0_[0-9]+.out")
    
    ## for each tree file
    if(length(tree.files > 0)) {
      for(i in 1:length(tree.files)) {
        if(verb >= 1) cat(paste(message, ": removed ", tree.files[i], "\n", sep=""))
        if(rmfiles) unlink(tree.files[i])
      }
    }
  }
    
  if(verb >= 1 && message == "INTERRUPT") cat("\n")
}
