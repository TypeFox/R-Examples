gcomp <-
function(inf=NULL, perimeter='G', no.kprm=TRUE)
{
### DP 17.10.13
### this function compiles all cleaned Goldmann data 
### returns a matrix of individual results (areas, etc)
##
### inf == NULL  is interpreted as being in example mode, i.e.
### reading the demographic information from example files

if(!is.null(inf))
goldmann.demogr <- inf
else
goldmann.demogr<- Goldmann.demogr



names(goldmann.demogr)[1]<- 'ID'

b0<- (1:dim(goldmann.demogr)[1])[!is.na(goldmann.demogr$Eye)]
goldmann.demogr<- goldmann.demogr[b0,] 

list.of.files<- paste(goldmann.demogr$ID,
                      substring(goldmann.demogr$Eye,1,1),sep='')

num.files<- length(list.of.files)                             

### declare matrix for results     
areas.goldmann<- matrix(NA, nrow=num.files,ncol=8) ### 8 refers to 4 areas and 4 dists

                               
 
   
for( i in 1:num.files)
{
print(list.of.files[i])
areas.goldmann[i,]<- 
kFsubj(obj.name=list.of.files[i], 
       no.kprm=no.kprm, perimeter=perimeter,
no.cleaning=TRUE)
}


dimnames(areas.goldmann)<- list( list.of.files,
      c('III4e','I4e','I2e','blind spot',
        'superior nasal','superior temporal',
        'inferior temporal','inferior nasal'))

### separate ID and side.eye
side.eye<- as.factor(substring( dimnames(areas.goldmann)[[1]], 6,6))
attr(side.eye, 'levels')<- c('R','L')

id<- as.factor(substring( dimnames(areas.goldmann)[[1]], 1,5))
attr(id, 'levels')<- substring( dimnames(areas.goldmann)[[1]], 1,5)


areas.goldmann<- as.data.frame( areas.goldmann)
areas.goldmann$id<- id
areas.goldmann$side.eye<- side.eye
n2<- dim( areas.goldmann)[2]
areas.goldmann<- areas.goldmann[, c(n2-1, n2, 1:(n2-2))]

### write.csv(areas.goldmann, file="Goldmann results.csv", row.names=FALSE, na='')        
invisible(areas.goldmann)

}
