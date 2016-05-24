kf.sector <-
function(file.name, is.octopus=FALSE)
{
### This function realigns individual points along 24 sectors

### file.name is the file name using format defined for studyid
### note that file.name only contains e.g. "M2005R"
### output is a list with named elements "III4e" "I4e"   "I2e"   "blind"
### note that exactly one of the first three elements must be NULL
### because... 
### the output list is split by sector and isopter and
### contains the distance of each point to the origin and frequencies in 
### each sector (must be 0 or 1)


obj<- ls(envir=.GlobalEnv)[grep(file.name, ls(envir=.GlobalEnv))]

isopter.names<-c('III4e','I4e','I2e')
III4e<- I4e<- I2e<-  NULL
isopter.obj<- substring( obj, 8, nchar(obj))
isopter.obj<- isopter.obj[isopter.obj!='blind']

file.name1<-  substring ( file.name, 1, 5)
which.files<- is.element( isopter.names, isopter.obj)


if(which.files[1]) III4e<- get(obj[grep('III4e', obj)])
if(which.files[2]) I4e<-get(obj[grep('I4e', obj)])
if(which.files[3]) I2e<- get(obj[grep('I2e', obj)])

file0<- file.name
file1<- paste(file0, 'GIII4e', sep='')
file2<- paste(file0, 'GI4e', sep='')
file3<- paste(file0, 'GI2e', sep='')



if(is.octopus)
{
file0<- file.name
file1<- paste(file0, 'OIII4eproc', sep='')
file2<- paste(file0, 'OI4eproc', sep='')
file3<- paste(file0, 'OI2eproc', sep='')
}


which.files<- rep(FALSE,3)
which.files[1] <-exists(file1)
which.files[2] <-exists(file2)
which.files[3] <-exists(file3)

isopter.names<-c('III4e','I4e','I2e')

if(which.files[1]) III4e<- get(file1)
if(which.files[2]) I4e<-   get(file2)
if(which.files[3]) I2e<-   get(file3)


side.eye<- substring(file.name, 6,6) ## must be L or R only

if(side.eye=='L')
{
if( which.files[1]) III4e[,1]<- -III4e[,1]
if( which.files[2]) I4e[,1]<- -I4e[,1]
if( which.files[3]) I2e[,1]<- -I2e[,1]

}


### compute distaces

dIII4e<- dI4e<- dI2e<-  NULL

if(  which.files[1]) dIII4e<- apply( III4e, 1, function(x){sqrt(sum(x^2))})
if(  which.files[2]) dI4e<- apply( I4e, 1, function(x){sqrt(sum(x^2))})
if(  which.files[3]) dI2e<- apply( I2e, 1, function(x){sqrt(sum(x^2))})


angIII4e<- angI4e<- angI2e<-  NULL

### transform coordinates to radians - use coord2rad from library(circular) - uses atan2
### and returns an object of class circular

if( which.files[1]) angIII4e<- coord2rad(III4e)
if( which.files[2]) angI4e<- coord2rad(I4e)
if( which.files[3]) angI2e<- coord2rad(I2e)


df.angles<- list(angIII4e, angI4e, angI2e)
df.dist<- list( dIII4e, dI4e, dI2e)

### counts by sector

arcsM<- rad(seq( 0,  345, by=15)) ### angles for centra lines of arcs
arcsL<- rad(c(352.5, seq(7.5, 337.5, by=15))) ### lower
arcsU<- rad( c( 7.5, seq(22.5, 352.5, by=15))) ### upper

positions<- dists<- freqs<- matrix(0, ncol=3, nrow=24) 
dimnames(positions)<- dimnames(dists)<- dimnames(freqs)<- 
        list( 1:24, c('III4e', 'I4e', 'I2e'))

length.dists<- unlist(lapply( df.dist, length))

for (i in 1:3) ### III4e, I4e, I2e
if( which.files[i] ) 
{
    pos.aux<- rep(0, length.dists[i]) ### initialises vector of positions for each observation
    for (j in 1:24)  ### iterate over sectors
           {   
                   if( j==1) 
                        b0<- df.angles[[i]] <= arcsU[j] | df.angles[[i]] > arcsL[j]
                     else
                        b0<- df.angles[[i]] <=arcsU[j] & df.angles[[i]] > arcsL[j] ### find the index of observation in sector j
                  if(sum(b0)>0) ### if the sector is not empty
                   {
                       if(sum(b0)>1) 
                           {
                              print(paste('there are ',sum(b0),' points in sector ',j,
                                          ' of isopter ',isopter.names[i],
                                          ' in subject ', file.name,sep=''),quote=FALSE)                                                                                                                                                                 
                             positions[j,i]<- dists[j,i]<- freqs[j,i]<-NA                                                                                  
                            }                                                     
                       else
                         {                                                                    
                             positions[j,i]<- (1:length.dists[i])[b0] ### observation in sector j 
                             dists[j,i]<- df.dist[[i]][b0] ### puts distance in obs determined by b0 in correct sector
                                if(j==1)
                                   freqs[j,i]<- freqs[j,i] + sum(df.angles[[i]]<= arcsU[j] | df.angles[[i]] > arcsL[j])
                                else
                                    freqs[j,i]<- freqs[j,i] + sum(df.angles[[i]]<= arcsU[j] & df.angles[[i]] > arcsL[j])
                         }
              }   
           } ### close loop over sectors
               
}


### output is a list with elements distance and frequencies
invisible(list( dists=dists, freqs=freqs))
}
