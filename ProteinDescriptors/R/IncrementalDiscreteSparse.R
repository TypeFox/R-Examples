#' Incremental discrete sparse descriptor.
#'
#' \code{incrementalDiscreteSparse} returns the sum of incremented sparse
#' descriptors of amino acids in a protein sequence.
#'
#' @param x A string of amino acid letters
#' @return A 20 dimensional numeric vector
#'
#' @export IncrementalDiscreteSparse
#'
#' @examples
#' x = "LALHLLLLHMHMMDRSLLLH"
#' IncrementalDiscreteSparse(x)

IncrementalDiscreteSparse<-function (x)
{
    AAs= c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
n=nchar(x)
x=strsplit(x, split = "")
fact=summary(factor(x[[1]],levels=AAs),maxsum=21)/n

O=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
vx=unlist(x)

for(i in 1:n)
{
if (vx[i]=="A")
O<-O+c(i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="R")
O<-O+c(0,i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="N")
O<-O+c(0,0,i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="D")
O<-O+c(0,0,0,i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="C")
O<-O+c(0,0,0,0,i,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="E")
O<-O+c(0,0,0,0,0,i,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="Q")
O<-O+c(0,0,0,0,0,0,i,0,0,0,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="G")
O<-O+c(0,0,0,0,0,0,0,i,0,0,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="H")
O<-O+c(0,0,0,0,0,0,0,0,i,0,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="L")
O<-O+c(0,0,0,0,0,0,0,0,0,i,0,0,0,0,0,0,0,0,0,0)

if (vx[i]=="I")
O<-O+c(0,0,0,0,0,0,0,0,0,0,i,0,0,0,0,0,0,0,0,0)

if (vx[i]=="K")
O<-O+c(0,0,0,0,0,0,0,0,0,0,0,i,0,0,0,0,0,0,0,0)

if (vx[i]=="M")
O<-O+c(0,0,0,0,0,0,0,0,0,0,0,0,i,0,0,0,0,0,0,0)

if (vx[i]=="F")
O<-O+c(0,0,0,0,0,0,0,0,0,0,0,0,0,i,0,0,0,0,0,0)

if (vx[i]=="P")
O<-O+c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,i,0,0,0,0,0)

if (vx[i]=="S")
O<-O+c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,i,0,0,0,0)

if (vx[i]=="T")
O<-O+c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,i,0,0,0)

if (vx[i]=="W")
O<-O+c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,i,0,0)

if (vx[i]=="Y")
O<-O+c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,i,0)

if (vx[i]=="V")
O<-O+c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,i)

}

O=O/n
result=c(fact,O)
return(result)
}
