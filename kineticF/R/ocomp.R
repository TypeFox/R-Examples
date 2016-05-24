ocomp <-
function(inf=NULL, no.kprm=TRUE, perimeter='O')
{
### DP 17.10.13
### this function compiles all cleaned Octopus data 
### returns a matrix of individual results (areas, etc)
### and writes a *.csv file

if(!is.null(inf))
octopus.demogr <- inf
else
octopus.demogr<- Octopus.demogr

names(octopus.demogr)[1]<- 'ID'

b0<- (1:dim(octopus.demogr)[1])[!is.na(octopus.demogr$Eye)]
octopus.demogr<- octopus.demogr[b0,] 


list.of.files<- paste(octopus.demogr$ID,
                      substring(octopus.demogr$Eye,1,1),sep='')

num.files<- length(list.of.files)                             

### declare matrix for results     
areas.octopus<- matrix(NA, nrow=num.files,ncol=8) ### 8 refers to 4 areas and 4 dists

                               
        
for( i in 1:num.files)
{
print(list.of.files[i])
areas.octopus[i,]<- 
kFsubj(obj.name=list.of.files[i], 
       no.kprm=no.kprm, perimeter=perimeter,
no.cleaning=TRUE)
}


dimnames(areas.octopus)<- list( list.of.files,
      c('III4e','I4e','I2e','blind spot',
        'superior nasal','superior temporal',
        'inferior temporal','inferior nasal'))

### separate ID and side.eye
side.eye<- as.factor(substring( dimnames(areas.octopus)[[1]], 6,6))
attr(side.eye, 'levels')<- c('R','L')


id<- as.factor(substring( dimnames(areas.octopus)[[1]], 1,5))
attr(id, 'levels')<- substring( dimnames(areas.octopus)[[1]], 1,5)


areas.octopus<- as.data.frame( areas.octopus)
areas.octopus$id<- id
areas.octopus$side.eye<- side.eye
n2<- dim( areas.octopus)[2]
areas.octopus<- areas.octopus[, c(n2-1, n2, 1:(n2-2))]

                              

# write.csv(areas.octopus, file='Octopus results.csv', row.names=FALSE, na='')        

invisible(areas.octopus)

}
