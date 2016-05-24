#' Discrete sparse descriptor.
#'
#' \code{DiscreteSparse} returns the sum of sparse descriptors of amino acids in a protein
#' sequence.
#'
#' @param x A string
#' @return A 20 dimensional numeric vector
#'
#' @export ShiftedDiscreteSparse
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' ShiftedDiscreteSparse(x)

ShiftedDiscreteSparse<-function (x)
{
AAs = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
n=nchar(x)
x=strsplit(x, split = "")

    sparse = t(data.frame(
	A = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	R = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	N = c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	D = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	C = c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	E = c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
	Q = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
      G = c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
	H = c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
	I = c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
	L = c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
	K = c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
	M = c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
	F = c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
      P = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
	S = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
	T = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
	W = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
	Y = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
	V = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)))

O=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
vx=unlist(x)

for(i in 1:n)
{

if (vx[i]=="A")
O<-O+sparse[1,]

if (vx[i]=="R")
O<-O+sparse[2,]

if (vx[i]=="N")
O<-O+sparse[3,]

if (vx[i]=="D")
O<-O+sparse[4,]

if (vx[i]=="C")
O<-O+sparse[5,]

if (vx[i]=="E")
O<-O+sparse[6,]

if (vx[i]=="Q")
O<-O+sparse[7,]

if (vx[i]=="G")
O<-O+sparse[8,]

if (vx[i]=="H")
O<-O+sparse[9,]

if (vx[i]=="L")
O<-O+sparse[10,]

if (vx[i]=="I")
O<-O+sparse[11,]

if (vx[i]=="K")
O<-O+sparse[12,]

if (vx[i]=="M")
O<-O+sparse[13,]

if (vx[i]=="F")
O<-O+sparse[14,]

if (vx[i]=="P")
O<-O+sparse[15,]

if (vx[i]=="S")
O<-O+sparse[16,]

if (vx[i]=="T")
O<-O+sparse[17,]

if (vx[i]=="W")
O<-O+sparse[18,]

if (vx[i]=="Y")
O<-O+sparse[19,]

if (vx[i]=="V")
O<-O+sparse[20,]

O<-c(O,0)
sparse<-cbind(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), sparse)

}

result=O/n
return(result)

}


