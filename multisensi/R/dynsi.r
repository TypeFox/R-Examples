# Multisensi R package ; file dynsi.r (last modified: 2015-12-07) 
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
dynsi <- function(formula, model, factors, cumul=FALSE, simulonly=FALSE, nb.outp=NULL, Name.File=NULL, ...)
#===========================================================================
{
    ##INPUTS
    ## formula      : ANOVA formula like "A+B+c+A:B"   OR  The max interaction
    ##                order like 2 for example.
    ## model        : Data.frame output matrix OR The name of the R-function
    ##                which decribes the model. This function
    ##                must take only a vector corresponding to the input factors values
    ## factors      : Data.frame design if model is data.frame OR a list of
    ##                factors levels: factor<- list(A=c(0,1),B=c(0,1,4))
    ## cumul        : logical value. If TRUE the PCA will be done on the cumalative outputs
    ## simulonly    : logical value.  If TRUE the program simulates the
    ##                model outputs and stops
    ## Name.File    : Name of file containing the R-function model.
    ##                E.g  "exc.ssc"
    ## nb.outp      : number of the first output to be considered if it is not null
    ## ...          : possible fixed parameters of the model function

    ## OUTPUTS
    ## Objet de classe dynsi contenant
    ## X            : data.frame design of experiment (input sample)
    ## Y            : data.frame of model ouput output matrix (response)
    ## SI           : data.frame of first order, two ... Sensitivity Indices (SI) on model outputs
    ## tSI          : data.frame of total SI on model outputs
    ## mSI          : data.frame of principal SI on model outputs
    ## iSI          : data.frame of interaction SI on model outputs
    ## Att          : 

  if(is.null(dim(factors))){
    # factors is a list then we need to build a design
    multisensi.design=planfact.as
    d.args=factors
  }else{
    multisensi.design=factors
    d.args=list()
  }

  result <- multisensi(design=multisensi.design, model=model, reduction=NULL, dimension=nb.outp, center=FALSE, scale=FALSE, analysis=analysis.anoasg, cumul=cumul, simulonly=simulonly, Name.File=Name.File, design.args=d.args, analysis.args=list(formula=formula,keep.ouputs=FALSE), ...)
  cat("Warning : dynsi function will not be maintained in future version of multisensi package.\n")
  cat("You may use multisensi function instead, like this :\n")
  print(result$call.info$call)

  return(result)


}
