# Multisensi R package ; file graph.bar.r (last modified: 2016-04-18) 
# Copyright INRA 2011-2015 
# Authors: C. Bidot, M. Lamboni, H. Monod
# MaIAGE, INRA, Univ. Paris-Saclay, 78350 Jouy-en-Josas, France
#
# More about multisensi in http://cran.r-project.org/web/packages/multisensi/
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
#===========================================================================
graph.bar <-function(x, col=1, nb.plot=15,xmax=NULL, beside=TRUE, xlab=NULL, ...)
#===========================================================================
{
    ##x :          GSI objects
    ##col :          A number describing the column for S and St
    ##nb.plot:       A number decribing the max number of factor bars to be ploted
    ##xmax:          A number that precise the indice limit

    out <- data.frame(x$tSI[,col],x$mSI[,col])
    rownames(out) <-  rownames(x$tSI)

    if(is.null(xlab)){xlab=colnames(x$tSI)[col]}

    out <- out[rev(order(out[,1])),]

    outgraph <- out[rev(1:min(nb.plot,nrow(out))),]
    #outgraph <- outgraph/100
    if (is.null(xmax)==TRUE){
        xmax <- 1.1*max(outgraph)
    } else {
        if (xmax>1) {xmax=1}
    }
    outnames <- rownames(outgraph)
    outnames <- gsub("([ ])","",outnames)
    if(beside==FALSE){
        outgraph.2 <- rbind(outgraph[,2],outgraph[,1]-outgraph[,2])
        invisible(
                  barplot(outgraph.2,
                          horiz=TRUE,
                          names.arg=outnames,
                          beside=beside,
                          las=2,
#                          cex.axis=2,
#                          cex.names=1,#1.5,
                          xlim=c(0,xmax),
                          xlab=xlab,
                          ...))
    } else {
        outgraph.2 <- rbind(outgraph[,1],outgraph[,2])
        invisible(
                  barplot(outgraph.2,
                          horiz=TRUE,
                          names.arg=outnames,
                          beside=beside,
                          las=2,
#                          cex.axis=2,
#                          cex.names=1,#1.5,
                          xlim=c(0,xmax),
                          xlab=xlab,
                          ...))
    }
}

