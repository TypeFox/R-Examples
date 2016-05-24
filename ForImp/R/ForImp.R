ForImp <-
function(mat, p=2) 
{
# default value for the tolerance error in the ALS procedure of homals

eps<-0.000001

# matrix with missing values

Table1<-as.matrix(cbind(1:nrow(mat),mat))

# read the number of rows and columns

h<-nrow(Table1)
s<-ncol(Table1)

# one additional column for identifying rows
k<-s-1

TableVarNames<-1:s
MatDat<- as.matrix(Table1)
countmiss<-0
nummatric<-0

# VectNumMiss
# vector with the number of missing value for each unit (the number of cells equals the maximum number of missing values in a record)
# first element: number of complete units
# second element: number of units with one missing value
# third element: number of units with two missing values
#...
# last element: number of units with k missing values

VectNumMiss<-matrix(data=rep(0,k+1), nrow = k+1, ncol=1, byrow=TRUE)

# VectIndeces
# vector of zeros and ones with the corresponding indexes: the element (i,j)  equals 1 if unit j has i-1 missing values, otherwise equals 0

VectIndeces<-matrix(data=rep(0,(k+1)*h), nrow = k+1, ncol=h, byrow=TRUE)

for(i in 1:h) {
# recording the current row
TempVect<-matrix(data=MatDat[i,],nrow = 1, ncol=s, byrow=TRUE)
for(j in 2:s) {
if (is.na(TempVect[j])) {
# if missing value, increase the counter
countmiss<-countmiss+1
}
        }
for(l in 0:k) {
if (countmiss==l) {
# increases the number of records with a number of missing values equal to the index of the corresponding cell of the vector with the number of missing values
VectNumMiss[l+1]<-VectNumMiss[l+1]+1
}
}
# structure of VectIndeces: row=number of missing-1, columns=index of the record with that number of missing
VectIndeces[countmiss+1,i]<-1
countmiss<-0
}

### Setting up disjoint (sub-)matrices A_k - with exactly k missing values per unit - by splitting the original matrix

# structure of VectComb: the first two columns=number of missing, number of records with that number of missing;
# other columns as VectIndeces
VectComb<-cbind(0:k,VectNumMiss)
VectComb<-cbind(VectComb,VectIndeces)

# Creating (sub-)matrices A_k (see algorithm in the paper), "cleaning" VectComb
# x-1=MAX number of matrices to be created
x<-1
# there may be records with one or three missing values, and no one with two!
while((VectComb[x,2]!=0 & x!=k+2)) {
x<-x+1
}
x<-max(which(VectComb[,2]>0))+1
# deletes useless rows
VectComb<-head(VectComb,x-1)

Matrices <- as.list(1:(x-1))
Matrices1 <- as.list(1:(x-1))

# formatting matrices A_k with heading row "9999"
for (w in 1:(x-1)){
Matrices[[w]]<- matrix(data=rep(9999,(1*s)), nrow = 1, ncol=s, byrow=TRUE)
}

# Matrices contains the matrices to be written, with a heading row formatted with "9999"

for (w in 3:(h+2)){
TempVector<-VectComb[VectComb[,w]==1]
Matrices[[TempVector[1]+1]]<-rbind(Matrices[[TempVector[1]+1]],MatDat[w-2,])
}

# Matrices1 contains matrices to be written
for (w in 1:(x-1)){
Matrices1[[w]]<-as.data.frame(Matrices[[w]])
Matrices1[[w]]<-tail(Matrices1[[w]],dim(Matrices1[[w]])[1]-1)
}

Matrix0ID<- as.data.frame(Matrices1[[1]])
# A_0(0) submatrix with only observed and no missing values
Matrix0<- as.data.frame(Matrices1[[1]][,2:s])
# Matrixlw is used to compute loadings
Matrixlw<-Matrix0

###
# iterative algorithm
###

z<-1
while (z<=x-2){
variabili<-attributes(Matrix0)$names
###
# run NLPCA, one dimension, ordinal variables, rank-1 restriction, on A_0(z) (sub-matrix without missing values at step z)
###
NPCA <- homals(Matrix0, ndim=1, rank=1, level ="ordinal", eps)
# category quantifications
CatQuant<-NPCA$catscores
CatQuant1 <- as.list(1:(k))
CatQuant1<-CatQuant
Matrices2 <- as.list(1:(k))
for (w in 1:(k)){
livelli<-as.integer(levels(NPCA$dframe[,w]))
numlev<-as.integer(nlevels(NPCA$dframe[,w]))
Matrices2[[w]]<-data.frame(Variable=rep(TableVarNames[w],numlev),Category=livelli,Quantification=CatQuant1[[w]])
}
# variable loadings
CompLoad1<-unlist(NPCA$loadings, recursive = TRUE, use.names=FALSE)
CompLoadApp<-numeric()
for (r in 1:(k)){
CompLoadApp[r]<-CompLoad1[r]
}
CompLoad0<-data.frame(VARIABLE=variabili,VALUE=CompLoadApp)

First_Matrix<-Matrix0
First_MatrixID<-Matrix0ID
First_MatrixModified<-First_Matrix
First_MatrixModifiedID<-First_MatrixID
# check out if submatrix A_z+1 is empty
# if true, skip and go on
while(nrow(Matrices[[z+1]])<=1)
{
z<-z+1
}
Second_MatrixID<-Matrices1[[z+1]]
Second_Matrix<-Matrices1[[z+1]][,2:s]
Second_MatrixModified<-Second_Matrix
Second_MatrixModifiedID<-Second_MatrixID

# Replace original categories in the first matrix with new category quantifications from the NPCA weight'd with corresponding component loadings

h<-nrow(CompLoad0)
J<-nrow(First_Matrix)
COMP_LOAD_VARIABLE<-matrix(data=CompLoad0[,1], nrow = h, ncol=1, byrow=TRUE)
COMP_LOAD_VALUE<-matrix(data=CompLoad0[,2], nrow = h, ncol=1, byrow=TRUE)
for(i in 1:h) {
currentname<-COMP_LOAD_VARIABLE[i]
CURRENT_matrix<-Matrices2[[i]]
l<-nrow(CURRENT_matrix)
CURRENT_CATEGORY<-matrix(data=CURRENT_matrix$Category, nrow = l, ncol=1, byrow=TRUE)
CURRENT_QUANTIFICATION<-matrix(data=CURRENT_matrix$D1, nrow = l, ncol=1, byrow=TRUE)
for(j in 1:l) {
category<-CURRENT_CATEGORY[j]
quantification<-CURRENT_QUANTIFICATION[j]
# containing quantifications for the complete matrix
First_MatrixModified[First_MatrixModified[[currentname]] == category,][[currentname]]<-COMP_LOAD_VALUE[i]*rep(CURRENT_QUANTIFICATION[j],nrow(First_MatrixModified[First_MatrixModified[[currentname]] == category,]))
First_MatrixModifiedID[First_MatrixModifiedID[[currentname]] == category,][[currentname]]<-COMP_LOAD_VALUE[i]*rep(CURRENT_QUANTIFICATION[j],nrow(First_MatrixModifiedID[First_MatrixModifiedID[[currentname]] == category,]))
}
}

# Replace original categories in the second matrix with new category quantifications from the original matrix NPCA , weight'd with corresponding component loadings

for(i in 1:h) {
currentname<-COMP_LOAD_VARIABLE[i]
CURRENT_matrix<-Matrices2[[i]]
l<-nrow(CURRENT_matrix)
CURRENT_CATEGORY<-matrix(data=CURRENT_matrix$Category, nrow = l, ncol=1, byrow=TRUE)
CURRENT_QUANTIFICATION<-matrix(data=CURRENT_matrix$D1, nrow = l, ncol=1, byrow=TRUE)
# g is the number of units with z missing
g<-nrow(Second_MatrixModified)
for (t in 1:g) {
if (!is.na(Second_MatrixModified[[currentname]][t])){
for(j in 1:l) {
if (Second_MatrixModified[[currentname]][t]==CURRENT_CATEGORY[j]){
Second_MatrixModified[[currentname]][t]<-COMP_LOAD_VALUE[i]*CURRENT_QUANTIFICATION[j]
}
}
}
if (!is.na(Second_MatrixModifiedID[[currentname]][t])){
for(j in 1:l) {
if (Second_MatrixModifiedID[[currentname]][t]==CURRENT_CATEGORY[j]){
Second_MatrixModifiedID[[currentname]][t]<-COMP_LOAD_VALUE[i]*CURRENT_QUANTIFICATION[j]
}
}
}
}
}

# Count the number of missing data in the second matrix
# countmiss contains the col numbers of missing data
TempVect<-numeric()
TempVect<-as.matrix(Second_MatrixModified)
TempVectID<-as.matrix(Second_MatrixModifiedID)
# counters
countmiss<-numeric()
countelements<-0
countmissID<-numeric()
countelementsID<-0

for(d in 1:g) {

for(j in 1:h) {
if (is.na(TempVect[d,j])) {
countelements<-countelements +1
countmiss[countelements]<-j
        }

        }
for(b in 2:(h+1)) {
if (is.na(TempVectID[d,b])) {
countelementsID<-countelementsID +1
countmissID[countelementsID]<-b-1
                            }
                              }

numelem<-countelements
numelemID<-countelementsID

# Bind the two matrices WITHOUT columns containing missing data

GeneralMatrix<-rbind(First_MatrixModified,Second_MatrixModified[d,])
GeneralMatrixID<-rbind(First_MatrixModifiedID,Second_MatrixModifiedID[d,])
countdeleted<-0
countdeletedID<-0
for(q in 1:numelem)   {
GeneralMatrix<-GeneralMatrix[,-(countmiss[q]-countdeleted)]
countdeleted<-countdeleted +1
                      }
for(q in 1:numelemID) {
GeneralMatrixID<-GeneralMatrixID[,-(countmissID[q]-countdeletedID+1)]
countdeletedID<-countdeletedID +1
                      }

# Obtain distance matrix (do) and vector of Minkowsky distances between the last row (of A_z+1) and the previous rows (in A_0(z))
            # Implementation of nearest neighbour imputation method
nrowmat<-nrow(as.matrix(GeneralMatrix))
nrowmatID<-nrow(as.matrix(GeneralMatrixID))
ncolmatID<-ncol(GeneralMatrixID)
# Minkowsky distance with parameter p (by default, p=2->euclidean distance)
if(p>0 & p!=Inf)
{
do<-dist(GeneralMatrix, method = "minkowski", diag = FALSE, upper = FALSE, p=p)
doID<-dist(GeneralMatrixID[,2:ncolmatID], method = "minkowski", diag = FALSE, upper = FALSE, p = p)
}
else
# maximum distance (supremum norm)
{
do<-dist(GeneralMatrix, method = "maximum", diag = FALSE, upper = FALSE)
doID<-dist(GeneralMatrixID[,2:ncolmatID], method = "maximum", diag = FALSE, upper = FALSE)
}
attrn <- attr(do, "Size")
attrnID <- attr(doID, "Size")
# Vector of distances
distances<-numeric()
distancesID<-numeric()
for(y in 1:nrowmat-1)   {
distances[y]<-do[attrn*(y-1) - y*(y-1)/2 + attrn-y]
                        }
for(y in 1:nrowmatID-1) {
distancesID[y]<-doID[attrnID*(y-1) - y*(y-1)/2 + attrnID-y]
                        }
# find the minimum distance
mindist<-which.min(distances)
# find the unit with the minimum distance
mindistID<-which.min(distancesID)
vect<-numeric()
vectID<-numeric()
vect<-Matrix0[mindist,]
vectID<-Matrix0ID[mindistID,]

for(f in 1:countelements)   {
Second_Matrix[d,countmiss[f]]<-Matrix0[mindist,countmiss[f]]
                            }
for(f in 1:countelementsID) {
Second_MatrixID[d,countmissID[f]+1]<-Matrix0ID[mindistID,countmissID[f]+1]
                            }
###
# Update the complete/imputed matrix, appending the units with imputed values
###
Matrix0<-rbind(Matrix0,Second_Matrix[d,])
Matrix0ID<-rbind(Matrix0ID,Second_MatrixID[d,])
# reset the counters
countelements<-0
countelementsID<-0
}
# next iteration step
z<-z+1
}
# final imputed matrix
IM<-Matrix0ID[,2:s][order(Matrix0ID[,1]),]
IM<-as.matrix(IM)
rownames(IM)<-NULL
colnames(IM)<-NULL
IM
}

