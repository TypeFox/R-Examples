# Multisensi R package ; file plot.gsi.r (last modified: 2016-04-18) 
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
plot.gsi <- function(x, nb.plot=10, nb.comp=3, graph=1:3, xmax=NULL, beside=TRUE, cor.plot=FALSE, xtick=TRUE, type="l",...)
#===========================================================================
{

    ## x:         GSI object

    nb.comp <- min(ncol(x$H), nb.comp)
    nb.plot <- min(nb.plot, nrow(x$tSI))

    ##if(is.null(graph)==FALSE && graph>3) {graph <- NULL}

    if(1 %in% graph){
        if(dev.cur() == 1 | dev.interactive()) dev.new();
        old.par <- par(no.readonly = TRUE)
        ## Correlation graph and main and total sensitivity indices bars on PCs
        graph.pc(x,
                 nb.plot=nb.plot,
                 nb.comp=nb.comp,
                 xmax=xmax,
                 cor.plot=cor.plot,
                 xtick=xtick,
                 beside=beside,
                 type=type,
                 ...)
        par(old.par)
    }

    if(2 %in% graph & colnames(x$mSI)[ncol(x$mSI)]=="GSI"){#x$call.info$analysis=="anova"
        # there is no GSI for sensitivity analysis
        if(dev.cur() == 1 | dev.interactive()) dev.new();
        old.par <- par(no.readonly = TRUE)
        ## Generalized main and total sensitivity indices bars
        graph.bar(x,
                  ncol(x$mSI),
                  nb.plot,
                  xmax=xmax,
                  beside=beside,
#                  xlab="GSI",passe dans graph.bar
                  ...)
        par(old.par)
    }

    if(3 %in% graph & !is.null(x$Rsquare)){
        # there is no Rsquare for sensitivity analysis
        if(dev.cur() == 1 | dev.interactive()) dev.new();
        old.par <- par(no.readonly = TRUE)
        ## Dynamic coefficient of determination
        plot(x$Rsquare,
             ylim=c(0,1),
             ylab="Rsquare",
             type=type,
             ...)
        par(old.par)
    }
}
