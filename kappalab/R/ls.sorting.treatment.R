##############################################################################
#
# Copyright © 2005 Michel Grabisch, Ivan Kojadinovic, and Patrick Meyer
#
# Ivan.Kojadinovic@polytech.univ-nantes.fr
#
# This software is a package for the statistical system GNU R:
# http://www.r-project.org 
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
##############################################################################

## Least squares sorting capacity identification: analysis of the classification

##############################################################################

ls.sorting.treatment <- function(P, cl.proto, a, A, cl.orig.A = NULL){
     ## P : the prototypes
     ## cl.proto : the classes of the prototypes
     ## a : the Mobius.capacity object
     ## A : the alternatives to be classified
     ## cl.orig.A : the classes of the alternatives of A, in case
     ##                            they are known

     ## check P
     if (!(is.matrix(P) && dim(P)[2] == a@n)) 
        stop("wrong criteria matrix")

     ## check cl.proto
     if (!(is.vector(cl.proto) && length(cl.proto) == dim(P)[1])) 
        stop("wrong vector of classes")
     
     ## check A
     if (!(is.matrix(A) && dim(A)[2] == a@n)) 
        stop("wrong criteria matrix")

     ## check cl.proto
     if (!(is.vector(cl.orig.A) && length(cl.orig.A) == dim(A)[1] || is.null(cl.orig.A)))
        stop("wrong vector of classes")
     
     ## the min and the max of the classes
     min.y <- min(cl.proto)
     max.y <- max(cl.proto)

     ## calculation of the Choquet integral of the prototypes, with
     ## respect to a

     Choquet.proto <- numeric(length(cl.proto))
     for (i in 1:length(cl.proto))
       Choquet.proto[i] <- Choquet.integral(a,P[i,])

     ## the following is not used later, but can be interesting information
     ## it is the determination of the min and the max of each class
     ## according to the prototypes

     ## the Choquet.proto are assigned to their class
     counter<-1
     classes.choquet<-list()
     for (i in min.y:max.y){
       classes.choquet[[counter]]<-Choquet.proto[cl.proto==i]
       counter<-counter+1
     }

     ## minmax represents the min and the max of each class
     minmax <- matrix(numeric(2*(counter-1)),2,counter-1)
     for (i in 1:dim(minmax)[2]){
     	minmax[1,i] <- min(classes.choquet[[i]])
	minmax[2,i] <- max(classes.choquet[[i]])
     }

     ## end of the useless but possibly interesting part

     ## affect.proto is a matrix of 3 lines which containes first the
     ## indexes of the protos
     ## then the original classes of the protos
     ## and finally the Choquet integral of the protos
     ## All is sorted according to increasing choquet integral values 

     affect.proto <- matrix(numeric(3*length(cl.proto)),3,length(cl.proto))
     temp.rank<-rank(Choquet.proto)

     for (i in 1:(length(cl.proto))){
	affect.proto[3,temp.rank[i]]<-Choquet.proto[i]
	affect.proto[2,temp.rank[i]]<-cl.proto[i]
	affect.proto[1,temp.rank[i]]<-i
      }

     ## Calculation of the choquet integral of the elements of A

     Choquet.A <- numeric(dim(A)[1])
     for (i in 1:(dim(A)[1]))
       Choquet.A[i] <- Choquet.integral(a,A[i,])

     ## affectation of the interval of classes to each
     ## alternative of A
     ## (saved in class.A)
     
     class.A <- matrix(numeric(),2,dim(A)[1])

     for (i in 1:(dim(A)[1])){
	## on teste d'abord si Choquet.A[i] est plus petit que la plus petite eval des proto
	## si oui, left et right valent la plus petite classe
	if (Choquet.A[i] < affect.proto[3,1]) {
		left = min.y
		right = min.y
	} else {	
		## de meme on teste si Choquet.A[i] est plus grand que la plus grande eval des proto
        	## dans ce cas, right vaut la plus grande classe
		
		if (Choquet.A[i] > affect.proto[3,length(cl.proto)]){
			right = max.y
			left = max.y
		} else {
			## maintenant on teste les cas normaux
			## pour commencer il faut determiner la position
			## de Choquet.A[i] dans affect.proto[3,].
			## d'abord le max des indices pour lesquels affect.proto[3,] <= Choquet.A[i]
			max.temp <- max( c(1:length(cl.proto))[affect.proto[3,] <= Choquet.A[i]] )
			left <- max(affect.proto[2,1:max.temp])
			min.temp <- min(c(1:length(cl.proto))[affect.proto[3,] >= Choquet.A[i]])
			right <- min(affect.proto[2,min.temp:length(cl.proto)])
			}
		}
	class.A[1,i] <- min(left,right)
	class.A[2,i] <- max(left,right)
	}	

     ## if no cl.orig.A has been given, then the following code is
     ## not interpreted
     
     
     ## counting of the errors (i.e. comparison of class.A and cl.orig.A)
     ## correct.A : the matrix which counts how the elements of A are
     ## classed compared to cl.orig.A
     ## first line: correct of degree 0
     ## second line: correct of degree 1 (an ambiguity)
     ## ...
     ## last line: classification error
     
     if (!is.null(cl.orig.A)){
       correct.A <- matrix(numeric(dim(A)[1]*(max.y-min.y+2)), max.y-min.y+2, dim(A)[1])
       
       for (i in 1:length(cl.orig.A)){
         if ((cl.orig.A[i]<=class.A[2,i]) & (cl.orig.A[i]>=class.A[1,i])){                                        # alors on sait deja que la classe originale se trouve dans l'intervalle des classes du modele
           if (class.A[1,i] == class.A[2,i]){
                  # alors on sait que le modele a affecte cette alternative de A a une seule classe
                  # qui correspond a la classe originale
             correct.A[1,i] = 1
           } else {
		  # dans ce cas, le modele a affecte cette alternative de A a plusieurs classes
	     correct.A[(class.A[2,i]-class.A[1,i]+1), i] = 1
           }
         } else {
		# on a une erreur de classif
           correct.A[max.y-min.y+2,i] = 1
         }
      }

      ## the matrix of correct affectations (by affectation type)
       eval.correct <- numeric(dim(correct.A)[1])

       for (i in 1:dim(correct.A)[1]){
         eval.correct[i] <- sum(correct.A[i,])/(dim(A)[1])
       }
     }
     else {
       ## no class has been given to check correctness of class.A
       ## eval.correct and correct.A become NULL
       eval.correct <- NULL
       correct.A <- NULL
     }

     return(list(correct.A = correct.A,
                 class.A = class.A,
                 eval.correct = eval.correct,
                 minmax.P = minmax,
                 Choquet.A = Choquet.A))
}

