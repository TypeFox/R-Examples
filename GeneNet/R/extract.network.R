### extract.network  (2007-05-24)
###
###     Extract significant network elements
###
### Copyright 2007 Rainer Opgen-Rhein and Korbinian Strimmer
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




extract.network = function(network.all, method.ggm=c("prob", "qval","number"), 
      cutoff.ggm=0.8, method.dir=c("prob","qval","number", "all"), 
      cutoff.dir=0.8,verbose=TRUE)
{
  method.ggm = match.arg(method.ggm)
  method.dir = match.arg(method.dir)

  ## edges

  if(method.ggm=="prob")
  {
    if (cutoff.ggm<0|cutoff.ggm>1) stop("edges: cutoff for prob must be between 0 and 1")
    edges=network.all[network.all$prob > cutoff.ggm,]
  }
  else if(method.ggm=="qval")
  {
    if (cutoff.ggm<0|cutoff.ggm>1) stop("edges: cutoff for qval must be between 0 and 1")
    edges=network.all[network.all$qval < cutoff.ggm,]
  }
  else if(method.ggm=="number")
  {
    if (cutoff.ggm%%1!=0) stop("edges: cutoff for \"number\" must be integer")
    edges=network.all[1:cutoff.ggm,]
  }  

  if( ncol(network.all) == 10 ) # list includes directions  (6 colums without directions)
  {
    directions=rep("undirected",nrow(edges))
    if(method.dir=="prob")
    {
      if (cutoff.dir<0|cutoff.dir>1)
        stop("directions: cutoff for prob must be between 0 and 1")
      directions[edges$prob.dir>cutoff.dir&edges$log.spvar>0]="1to2"
      directions[edges$prob.dir>cutoff.dir&edges$log.spvar<0]="2to1"
      sig.dir.all=sum(network.all$prob.dir>cutoff.dir)
    }
    else if(method.dir=="qval")
    {
      print(1)	
      if (cutoff.dir<0|cutoff.dir>1)
        stop("directions: cutoff for qval must be between 0 and 1")

      directions[edges$qval.dir<cutoff.dir&edges$log.spvar>0]="1to2"
      directions[edges$qval.dir<cutoff.dir&edges$log.spvar<0]="2to1"
      sig.dir.all=sum(network.all$qval<cutoff.dir)
    }
    else if(method.dir=="number")
    {
      if (cutoff.dir%%1!=0)
        stop("directions: cutoff for \"number\" must be integer")
 
      sort.idx=order(-abs(edges$log.spvar))
      directions[sort.idx<=cutoff.dir&edges$log.spvar>0]="1to2"
      directions[sort.idx<=cutoff.dir&edges$log.spvar<0]="2to1"
    }
    else if(method.dir=="all")
    {
      directions[edges$log.spvar>0]="1to2"
      directions[edges$log.spvar<0]="2to1"
    }  
    n.dir=sum(directions!="undirected")

    network=cbind(edges, directions) 
  
    # some output if wanted

    if(verbose==TRUE)
    {
      cat("\nSignificant edges: ", nrow(edges),"\n")
      cat("    Corresponding to ", round(nrow(edges)/nrow(network.all),4)*100,"%  of possible edges \n")
      if(method.dir=="prob"|method.dir=="qval")
      {
        cat("\nSignificant directions: ", sig.dir.all,"\n")
        cat("    Corresponding to ", round(sig.dir.all/nrow(network.all),4)*100,"%  of possible directions \n")
      }
      cat("Significant directions in the network: ", n.dir,"\n")
      cat("    Corresponding to ", round(n.dir/nrow(edges),4)*100,"%  of possible directions in the network \n")
    }  
  }
  else
  {
    if(verbose==TRUE)
    {
      cat("\nSignificant edges: ", nrow(edges),"\n")
      cat("    Corresponding to ", round(nrow(edges)/nrow(network.all),4)*100,"%  of possible edges \n")
    }
    network=edges
  }

  return(network)

}
