#=============================================================================
# BUUHWT_2D:
#----------------------------------------------------------------------------
# Object: Compute The Bottom Up Unbalanced Haar Wavelet Transform of an image.
#----------------------------------------------------------------------------
# PRE:
# Data: a matrix of dimension N1*N2 encoding the image
#
# POST: 
# a list with
#'BUUHWT_2D_Basis': an array containing the successive basis matrices of the Unbalanced Haar Basis.
# BUUHWT_2D_Basis[,,i] returns the basis matrix obtained at step i, i.e. encoding rank N1*N2-i +1
#'BUUHWT_2D_Coeff': the vector of associated coefficients.
#'DataInit': the initial image = matrix Data
#------------------------------------------------------------------------------
# Side-effect:
# the successive basis matrix are plotted
#------------------------------------------------------------------------------
# Requires:
# Blue2Red.Col.r (library BAGIDIS)
#==============================================================================

SHAH_Stepwise=function(Data,Plot=TRUE){
BUUHWE_2D_Stepwise(Data,Plot)
}

BUUHWE_2D_Stepwise=function(Data,Plot=TRUE){

#===============================================================================
# Prepare the procedure
#===============================================================================
DataInit=Data
  
N1= nrow(Data)
N2= ncol(Data)
  
Nb_Link = (N1-1)*N2 + (N2-1)*N1
BUUHWT_2D_Basis =array(NA, dim= c( N1,N2,N1*N2-1))
BUUHWT_2D_Coeff=rep(NA, N1*N2-1)
BUUHWT_2D_Approx =array(NA, dim= c( N1,N2,N1*N2-1))
BUUHWT_2D_Coord =list() 


# Locations of the data points
Indices=  1:(N1*N2) # on se deplace dans la matrice colonne par colonne. 
Coord =   1:(N1*N2) # identification initiale des groupes


# weights
W = matrix(1, nrow=N1,ncol=N2)
# Ndata
ND = matrix(1, nrow=N1,ncol=N2)

NeighborV = (Indices+1)   #Down
NeighborV[seq(N1,N1*N2,by=N1)] =NA
NeighborH = Indices+N1  #Right
NeighborH[(N1*N2):(N1*N2-N1+1)] =NA
# Pas besoin de considerer Bottom et Left, car deja pris en compte via l'indice precedent. 

#===============================================================================
# Data reduction algorithm
#===============================================================================

if(Plot){
oldpar = par(no.readonly = TRUE)
par(mfrow=c(N1,N2), oma=c(0,0,0,0), mar=c(1,1,1,1))}


for (Step in 1: (N1*N2-1)){
#for (Step in 1:6){
#print(Step)
#print(matrix(Coord,nrow=N1))

BUUHWT_2D_Coord[[Step]] =matrix(Coord,nrow=N1)

### Selection de l'indice et de la direction pour lequel on fait un rassemblement 

V_Diff = abs( ND[Coord[Indices]] * sqrt(  ND[Coord[NeighborV]] / (ND[Coord[NeighborV]]*ND[Coord[Indices]] + ND[Coord[Indices]]^2  ) ) *Data[Coord[Indices]]
         -ND[Coord[NeighborV]] * sqrt(  ND[Coord[Indices]]  /  (ND[Coord[Indices]]*ND[Coord[NeighborV]] + ND[Coord[NeighborV]]^2) )*Data[Coord[NeighborV]]  )

H_Diff = abs(ND[Coord[Indices]] * sqrt(ND[Coord[NeighborH]]   /  (ND[Coord[NeighborH]]*ND[Coord[Indices]] + ND[Coord[Indices]]^2   ) ) *Data[Coord[Indices]]
       - ND[Coord[NeighborH]] * sqrt(  ND[Coord[Indices]] /   (ND[Coord[Indices]]*ND[Coord[NeighborH]] + ND[Coord[NeighborH]]^2)  ) *Data[Coord[NeighborH]] )



#V_Diff = ND[Coord[Indices]]*W[Coord[Indices]]^-1*Data[Coord[Indices]]- ND[Coord[NeighborV]]*W[Coord[NeighborV]]^1*Data[Coord[NeighborV]] 

#V_Diff= abs(V_Diff/ sqrt(  W[Coord[Indices]]^(-2) +W[Coord[NeighborV]]^(-2)  ))##

#H_Diff= ND[Coord[Indices]]*W[Coord[Indices]]^-1*Data[Coord[Indices]]-ND[Coord[NeighborH]]*W[Coord[NeighborH]]^-1*Data[Coord[NeighborH]] 
#H_Diff= abs( H_Diff/ sqrt(W[Coord[Indices]]^(-2)+W[Coord[NeighborH]]^(-2))    )
##renvoie NA aux bords   
# je passe par la matrice Coord car cela permettra de rassembler des points.

V_Diff=(round(V_Diff,7))
H_Diff=round(H_Diff,7)

V_Diff[ Coord[Indices]== Coord[NeighborV]] =Inf
H_Diff[ Coord[Indices]==Coord[NeighborH]] =Inf
## on ne peut pas rassembler des choses qui le sont deja...

#print(Step)
#print(V_Diff)
#print(H_Diff)

## on evite les erreurs numeriques...


V_min_index = min(Indices[V_Diff==min(V_Diff, na.rm=TRUE)] ,na.rm=TRUE)  

H_min_index = min(Indices[H_Diff== min(H_Diff, na.rm=TRUE)],na.rm=TRUE) ### indices dans les matrices 

min_vec =c(min(V_Diff, na.rm=TRUE), min(H_Diff, na.rm=TRUE))
HV_min_index =which( min_vec == min(min_vec))[1]   ## 1--> V, 2--> H.   
# si egalite, je prends le 1er element (choix arbitraire)

### Data reduction and BUUHWT construction

if (HV_min_index ==1) {
WorkIndex = V_min_index
NeighborIndex= WorkIndex+1

}else {
WorkIndex = H_min_index
NeighborIndex= WorkIndex+N1

}

WeightPlus= sqrt(  ND[Coord[NeighborIndex]] / (ND[Coord[NeighborIndex]]*ND[Coord[WorkIndex]] + ND[Coord[WorkIndex]]^2  ) )
WeightMinus= - sqrt(  ND[Coord[WorkIndex]]  /  (ND[Coord[WorkIndex]]*ND[Coord[NeighborIndex]] + ND[Coord[NeighborIndex]]^2) )

### Actualisation des poids
#DimGroup = sum( Coord==Coord[WorkIndex])
#W_old=W
#W[Coord==Coord[WorkIndex]] = sqrt(W_old[WorkIndex]^2 + W_old[NeighborIndex]^2)
#W[Coord==Coord[NeighborIndex]] = sqrt(W_old[WorkIndex]^2 + W_old[NeighborIndex]^2)

### Sauvegarde de la matrice de base
BasisMatrix =matrix(0, nrow=N1,ncol=N2)

BasisMatrix[Coord==Coord[NeighborIndex]] = WeightMinus
BasisMatrix[Coord==Coord[WorkIndex]] =WeightPlus

Basis_Matrix = BasisMatrix/sqrt(sum(BasisMatrix^2)) #normalize

BUUHWT_2D_Basis[,,Step] = Basis_Matrix

#print(Basis_Matrix)
### Actualisation du voisinage dans les Coord

Coord[Coord==Coord[NeighborIndex]] = Coord[WorkIndex]


   
BUUHWT_2D_Coeff[Step] =  sum(Data * BUUHWT_2D_Basis[,,Step] )
#print(sum(Data * BUUHWT_2D_Basis[,,Step] ))


### Actualisation des data

Data[Coord==Coord[WorkIndex]] = mean((Data[Coord==Coord[WorkIndex]]))
#BUUHWT_2D_Approx[,,Step] = Data*Basis_Matrix

### Actualisation des nombres de membres des groupes

DimGroup = sum( Coord==Coord[WorkIndex])
ND[Coord==Coord[WorkIndex]] = sum( Coord==Coord[WorkIndex])




if(Plot){
#--------------------------------------------------------------------------------
## Check ##
#print(Data)
#print(Coord)
#image(x=1:N1, y=1:N2, z=Data,zlim=c(-15,15))
image(x=1:N2, y=1:N1, z=t(Basis_Matrix),zlim=c(-1,1),col=Blue2Red.Col(51),xaxp=c(1,N1,N1-1),yaxp=c(1,N2,N2-1))
#--------------------------------------------------------------------------------
}


} ## end FOR

Coeff0 = mean(Data) * sqrt(N1*N2)

if(Plot){
par(oldpar)
}
return(list('basis'=BUUHWT_2D_Basis,
            'details'=BUUHWT_2D_Coeff,
		'd0' =Coeff0,
            'im'=DataInit,
	      'labels.hist' =BUUHWT_2D_Coord))

} ## end function


