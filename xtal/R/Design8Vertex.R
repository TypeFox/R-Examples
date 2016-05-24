### subclass of 'Design'
### divide 96-well in four quardants as the last dimension
### for 3 dimension grid(6X4X4), 
### the portion of vertex/stock1 in each dimension input as list dim

setClass(Class='Design8Vertex',
#	representation=representation(dim='list'),
	contains='Design'
)

design8Vertex=function(volume,stock,dim){
	portion=array(dim=c(12,8,8),dimnames=list(col=c(1:12),row=LETTERS[1:8],vertex=c(1:8)))
	portion[1:6,1,1]=dim[[1]]
	portion[1,1:4,1]=dim[[2]]
	portion[1:6,1:4,1]=portion[1:6,1,1]%o%portion[1,1:4,1]

	portion[1:6,1,2]=1-dim[[1]]
	portion[6,1:4,2]=dim[[2]]
	portion[1:6,1:4,2]=portion[1:6,1,2]%o%portion[6,1:4,2]

	portion[1:6,4,3]=dim[[1]]
	portion[1,1:4,3]=1-dim[[2]]
	portion[1:6,1:4,3]=portion[1:6,4,3]%o%portion[1,1:4,3]

	portion[1:6,4,4]=1-dim[[1]]
	portion[6,1:4,4]=1-dim[[2]]
	portion[1:6,1:4,4]=portion[1:6,4,4]%o%portion[6,1:4,4]

	portion[7:12,1:4,1:4]=portion[1:6,1:4,1:4]*dim[[3]][2]
	portion[1:6,5:8,1:4]=portion[1:6,1:4,1:4]*dim[[3]][3]
	portion[7:12,5:8,1:4]=portion[1:6,1:4,1:4]*dim[[3]][4]

	portion[7:12,5:8,5:8]=portion[1:6,1:4,1:4]
	portion[1:6,1:4,5:8]=portion[7:12,5:8,5:8]*(1-dim[[3]][1])
	portion[7:12,1:4,5:8]=portion[7:12,5:8,5:8]*(1-dim[[3]][2])
	portion[1:6,5:8,5:8]=portion[7:12,5:8,5:8]*(1-dim[[3]][3])
	
	return(new(Class='Design8Vertex',volume=volume,stock=stock,portion=portion))

}

setMethod(
    f="writeTecan",
    signature=c(object="Design8Vertex"),
    def=function(object,fileName,source='Source1',destination='Destination',liquidType){
      tecan=object@portion*object@volume
      for (i in 1:8) { # i for vertex stock id
         for (k in 1:12) { # k for col
            for (j in 1:8) { #j  for row
                if (tecan[k,j,i]!=0) {
                  well=(k-1)*8+j
                  out=paste(source,i, destination,well,tecan[k,j,i],liquidType[i],sep=',')
                  write(out,file=fileName,append=TRUE)
                }
             }
         }
      }
    }
)

setMethod(
    f="writeTecan",
    signature=c(object="Design8Vertex",liquidType="missing"),
    def=function(object,fileName,source='Source1',destination='Destination',liquidType){
      liquidType=rep('B',length(object@stock))
      tecan=object@portion*object@volume
      for (i in 1:8) { # i for vertex stock id
         for (k in 1:12) { # k for col
            for (j in 1:8) { #j  for row
                if (tecan[k,j,i]!=0) {
                  well=(k-1)*8+j
                  out=paste(source,i, destination,well,tecan[k,j,i],liquidType[i],sep=',')
                  write(out,file=fileName,append=TRUE)
                }
             }
         }
      }
    }
)

