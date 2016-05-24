   
BAGIDIS.dist.BD= function(Details1, Breakpoints1,
                       Details2, Breakpoints2, 
                       p = 2, 
                       wk=NULL, 
                       Param=0.5){

#=======================================================================
# Selection des details et breakpoints  des deux series
#=======================================================================
if (p!= Inf){
det.S1  = ((1-Param)^(1/p)) * Details1
break.S1= (Param^(1/p)) * Breakpoints1

det.S2  =  ((1-Param)^(1/p)) * Details2
break.S2= (Param^(1/p)) * Breakpoints2 } else
{    # if p=Inf
det.S1  = (1-Param) * Details1
break.S1= Param * Breakpoints1

det.S2  =  (1-Param) * Details2
break.S2= Param * Breakpoints2 }


#=======================================================================
# Definition des poids  (si non specifie dans l'entete de la fonction)
#=======================================================================
if (is.null(wk)){
N =length(det.S1)
wk = log(N+1-(1:N))/log(N+1) }

#=======================================================================
# Cacul des dissimilarites partielles (en norme p)
#=======================================================================
if (p!= Inf){
dk =  ( (abs(break.S1-break.S2))^p + (abs(det.S1-det.S2))^p  )^(1/p)
} else {   #if p= Inf
dk = max( abs(break.S1-break.S2), abs(det.S1-det.S2) )    }
 
#print(dk)
#=======================================================================
# Cacul de la dissimilarite (semi-distance)
#=======================================================================

dissimilarity =sum(wk*dk)



return(dissimilarity) }
