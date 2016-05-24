### ggm.test.edges  (2008-11-15)
###
###   Compute p-values, q-values and posterior probabilities for network edges
###
### Copyright 2003-07 Juliane Schaefer, Rainer Opgen-Rhein and Korbinian Strimmer
###
###
### This file is part of the `GeneNet' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


##############################

ggm.list.edges = function(r.mat)
{
   pcor = sm2vec(r.mat)
   indexes = sm.index(r.mat)
   colnames(indexes) = c("node1", "node2")

   result = cbind(pcor, indexes)

   # sort according to magnitude of correlation 
   sort.idx = order(-abs(result[,1]))
   result = result[sort.idx,]
   
   # return as data frame
   return(as.data.frame(result))
}


# assign p-values, q-values and posterior probabilities to each edge
network.test.edges = function(r.mat, fdr=TRUE, direct=FALSE, plot=TRUE, ...)
{
   pcor = sm2vec(r.mat)
   indexes = sm.index(r.mat)
   colnames(indexes) = c("node1", "node2")
   w = cbind(pcor, indexes)
   
   if(fdr==TRUE) 
   {
     # fit null distribution to partial correlations and compute pvalues etc
     
     cat("Estimate (local) false discovery rates (partial correlations):\n")
     fdr.out = fdrtool(w[,1], statistic="correlation", plot=plot, ...)
     pval = fdr.out$pval
     qval = fdr.out$qval
     prob = 1-fdr.out$lfdr
   }
   else
   {
      pval = rep(NA, length(w[,1]))
      qval = pval
      prob = pval
   }
   
   result = cbind(w, pval, qval, prob)

   ###########

   if(direct==TRUE)
   {  
     spvar=attr(r.mat, "spv")
     if(is.null(spvar)) # if not yet available compute standardized partial variances
     {
       r.mat.cor = pcor2cor(r.mat)
       spvar=1/diag(solve(r.mat.cor))
     }
        
     p = length(spvar)
     r.spvar=(t(spvar%*%t(rep(1,p)))/(spvar%*%t(rep(1,p))))

     log.spvar = log(sm2vec(r.spvar))
    
     if(fdr==TRUE) 
     {
       if(plot==TRUE)
       {
         dev.new()
       }

       cat("Estimate (local) false discovery rates (log ratio of spvars):\n")
       fdr.out = fdrtool(log.spvar, statistic="normal", plot=plot, ...)
       pval.dir = fdr.out$pval
       qval.dir = fdr.out$qval
       prob.dir = 1 - fdr.out$lfdr
     }
     else
     {
       pval.dir = rep(NA, length(w[,1]))
       qval.dir = pval.dir
       prob.dir = pval.dir
     }

     result = cbind(result, log.spvar, pval.dir, qval.dir, prob.dir)
   }

   ###########

   sort.idx = order(-abs(result[, 1]))
   result = as.data.frame(result[sort.idx, ])
  
   # return as data frame
   return(result)
}

