##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################


# =========================================================== =
#  Check that passing the wrong values to the pmaxT function  =
#  returns the value FALSE and is not craching                =
# =========================================================== =

test.wrong_args <- function() {

    size_of_rows <- 10

    # Load test data and class label
    data(golub)
    classlabel <- golub.cl

    # Get a random set of 100 rows
    smallgd <- golub[1:size_of_rows,]


    # -------
    # Check X
    # -------

    # ncol(X) != length(classlabel)
    checkTrue(suppressWarnings(identical(pmaxT(smallgd[,1:3], classlabel), FALSE)))
    #checkTrue(identical(pmaxT(smallgd[,1:3], classlabel), FALSE))

    # X ! numeric
    te <- rbind(rep(c('y', 'n'), 19), rep(c('y', 'n'), 19))
    checkTrue(suppressWarnings(identical(pmaxT(te, classlabel), FALSE)))
    #checkTrue(identical(pmaxT(as.matrix(te), classlabel), FALSE))


    # ----------------
    # Check classlabel
    # ----------------

    # column_length != length(classlabel)
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel=c(0,1)), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel=c(0,1)), FALSE))

    # classlabel=c('y', 'n')
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel=c('y','n')), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel=c('y','n')), FALSE))

    # test='pairt' && (length(classlabel) != even)
    checkTrue(suppressWarnings(identical(pmaxT(smallgd[,1:15], classlabel=rep(c(0,1,1),5), test="pairt"), FALSE)))
    #checkTrue(identical(pmaxT(smallgd[,1:15], classlabel=rep(c(0,1,1),5), test="pairt"), FALSE))

    # test='pairt' && (max(classlabel) > 1) => we support only two groups (one pair)
    checkTrue(suppressWarnings(identical(pmaxT(smallgd[,1:15], classlabel=rep(c(0,1,2),5), test="pairt"), FALSE)))
    #checkTrue(identical(pmaxT(smallgd[,1:15], classlabel=rep(c(0,1,2),5), test="pairt"), FALSE))



    # ----------
    # Check test
    # ----------

    # test='wrong'
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, test="wrong"), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, test="wrong"), FALSE))

    # 'test' (test=10)
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, test=10), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, test=10), FALSE))

    # (test=c('pairt', 'blockf'))
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, test=c("pairt", "blockf")), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, test=c("pairt", "blockf")), FALSE))


    # ----------
    # Check side
    # ----------

    # 'side' (side='wrong')
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, side="wrong"), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, side="wrong"), FALSE))

    # 'side' (side=10)
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, side=10), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, side=10), FALSE))

    # side=c('upper', 'abs'))
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, side=c("upper", "abs")), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, side=c("upper", "abs")), FALSE))


    # -------------------------
    # Check fixed.seed.sampling
    # -------------------------

    # 'fixed.seed.sampling' (fixed.seed.sampling='wrong')
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, fixed.seed.sampling="wrong"), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, fixed.seed.sampling="wrong"), FALSE))

    # 'fixed.seed.sampling' (fixed.seed.sampling=20)
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, fixed.seed.sampling=20), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, fixed.seed.sampling=20), FALSE))


    # -------
    # Check B
    # -------

    # 'B' ( B > 2^30 )
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, B=2^30+1), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, B=2^30+1), FALSE))

    # 'B' ( B < 0 )
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, B=-1), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, B=-1), FALSE))

    # 'B' ( B=c(10, 30) )
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, B=c(10, 30)), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, B=c(10, 30)), FALSE))


    # --------
    # Check na
    # --------
    # 'na' ( na=c(10, 30) )
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, na=c(10, 30)), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, na=c(10, 30)), FALSE))

    # 'na' ( na='wrong' )
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, na="wrong"), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, na="wrong"), FALSE))


    # ------------
    # Check nopara
    # ------------

    # 'nonpara' (nonpara='wrong')
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, nonpara="wrong"), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, nonpara="wrong"), FALSE))

    # 'nonpara' (nonpara=10)
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, nonpara=10), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, nonpara=10), FALSE))

    # 'nonpara' (nonpara=c('y', 'n'))
    checkTrue(suppressWarnings(identical(pmaxT(smallgd, classlabel, nonpara=c("y", "n")), FALSE)))
    #checkTrue(identical(pmaxT(smallgd, classlabel, nonpara=c("y", "n")), FALSE))

}

