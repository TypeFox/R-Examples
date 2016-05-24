# Multisensi R package ; file biomasse.r (last modified: 2010-05-20) 
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
biomasse <- function(input , climdata, annee=3)
#===========================================================================
  {
    if(class(input)=="numeric"){input <- as.data.frame(as.list(input))}
    Eb <-    input[1,1]
    Eimax <- input[1,2]
    K <-     input[1,3]
    Lmax <-  input[1,4]
    A <-     input[1,5]
    B <-     input[1,6]
    TI <-    input[1,7]
    if (is.null(annee)){
           annee <- input[1,8]
    }


    #Calcul de PAR et de ST a partir des fichiers climatiques
    PAR<-0.5*0.01*climdata$RG[climdata$ANNEE==annee]
    Tmoy<-(climdata$Tmin[climdata$ANNEE==annee]+
           climdata$Tmax[climdata$ANNEE==annee])/2
    Tmoy[Tmoy<0]<-0
    ST<-Tmoy
    for (i in (2:length(Tmoy)))
      {
        ST[i]<-ST[i-1]+Tmoy[i]
      }

    #Calcul de LAI
    Tr<-(1/B)*log(1+exp(A*TI))
    LAI<-Lmax*((1/(1+exp(-A*(ST-TI))))-exp(B*(ST-(Tr))))
    LAI[LAI<0]<-0

    #Calcul de la biomasse (g/m2)
    U<-Eb*Eimax*(1-exp(-K*LAI))*PAR
    BIOMASSE<-sum(U)

   # BIOMASSE
    U <- cumsum(U)
    U
  }
#===========================================================================
#biomasse.nominal <- c(1.85,0.94,0.7,7.5,0.0065,0.00205,900)
#===========================================================================
#biomasse.input <- list(Eb=c(0.9,2.8),Eimax=c(0.9,0.99),K=c(0.6,0.8),Lmax=c(3,12),A=c(0.0035,0.01),B=c(0.0011,0.0025),TI=c(700,1100),Clim=c(1:14))
#attributes(biomasse.input)$qualitatif <- c(rep(F,7),T)
#===========================================================================

