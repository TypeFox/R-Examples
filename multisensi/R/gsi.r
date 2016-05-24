# Multisensi R package ; file gsi.r (last modified: 2015-12-07) 
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
gsi <- function(formula, model, factors, inertia=0.95, normalized=TRUE, cumul=FALSE, simulonly=FALSE, Name.File=NULL, ...)
#===========================================================================
{
  ##INPUTS
  ## formula      : ANOVA formula like "A+B+c+A:B"   OR  The max interaction
  ##                order like 2 for example.
  ## model        : Data.frame output matrix OR The name of the R-function
  ##                which decribes the model. This function
  ##                must take only a vector corresponding to the input factors values
  ## factors      : Data.frame design if model is data.frame OR a list of factors
  ##                levels: factor<- list(A=c(0,1),B=c(0,1,4))
  ## inertia      : Inertia proportion account by Principal components <1 (0.95 default )
  ##                OR number of PCs to be used (E.g 3)
  ## normalized   : logical value. TRUE (default) computes a normalized Principal
  ##                Component analysis.
  ## cumul        : logical value. If TRUE the PCA will be done on the cumulative outputs
  ## simulonly    : logical value.  If TRUE the program simulates only the model outputs
  ##                and stops
  ## Name.File    : Name of file containing the R-function model.
  ##                E.g  "exc.ssc"
  ## ...          : possible fixed parameters of the model function
  
  ## OUTPUTS
  ##GSI objet de classe gsi contient
  ##
  ## X            : data.frame design of experiment (input sample)
  ## Y            : data.frame output matrix (response)
  ## H            :
  ## L            :
  ## lambda       :
  ## inertia      : vector of inertia per PCs and Global criterion
  ## cor          : data.frame of correlation between PCs and outputs
  ## SI           : data.frame of first, two ... order Sensitivity Indices (SI) on PCs and
  ##                        first, two...  order Generalized SI (GSI)
  ## mSI          : data.frame of principal SI on PCs and principal GSI
  ## tSI          : data.frame of total SI on PCs and total GSI
  ## iSI          : data.frame of interaction SI on PCs and interaction GSI
  ## pred         :
  ## residuals    :
  ## Rsquare      : vector of dynamic coefficient of determination
  ## Att          : matrice 0-1 de correspondance facteurs*termes-du-modele
  ## normalized   : logical value used for normalized
  ## cumul        : logical value used for cumul
 

  if(is.null(dim(factors))){
    # factors is a list then we need to build a design
    multisensi.design=planfact.as
    d.args=factors
  }else{
    multisensi.design=factors
    d.args=list()
  }

  result <- multisensi(design=multisensi.design, model=model, reduction=basis.ACP, dimension=inertia, center=TRUE, scale=normalized, analysis=analysis.anoasg, cumul=cumul, simulonly=simulonly, Name.File=Name.File, design.args=d.args, basis.args=list(), analysis.args=list(formula=formula,keep.ouputs=FALSE), ...)
  cat("Warning : gsi function will not be maintained in future version of multisensi package.\n")
  cat("You may use multisensi function instead, like this :\n")
  print(result$call.info$call)

  return(result)
}
