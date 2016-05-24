preprocess.octopus <-
function( octopus.file, octopus.demogr=Octopus.demogr)
{
## This function is used to clean Octopus data - one subject at a time
### octopus.file is the raw file in the working directory; 
###         octopus.file==NULL uses the example
### demogr.file is the name of the demographic file with columns
### ID, Eye tested, Sex, age, test quality
###
### output is matrix Mfile containing coordinates and isopter labels
### for the subject 
### Note that this function writes  a subject-specific set of *.csv files
### one per isopter
### If there is 1 or 0 points it crashes!

Mfile<- unlist(strsplit(octopus.file, split=';'))


### get the 1st field and assign it as the name of the data.frame
m1<-match( substring(octopus.file,1,5), octopus.demogr[,1])
id<- as.character(octopus.demogr[m1,1]) ### subject id
eye.side<-substring(as.character(octopus.demogr[ m1 , 2]),1,1)
### this gives the position of the eye label for octopus.file in octopus.demogr


Mfile<- Mfile[-(1:32)] ### exclude identifiers ### 
### 32 characters is standard length for identifiers



Mfile<- as.numeric(Mfile)

pos.outliers<- (1:length(Mfile))[ abs(Mfile) > 100]
num.outliers<- length(pos.outliers)
outliers<-Mfile[pos.outliers]

if(length(pos.outliers)>0) 
  {
   Mfile<- Mfile[-pos.outliers] ### Values can't be > 100 degrees
   print(paste(num.outliers,' outliers removed', outliers, collapse='',
   sep=' : '),quote=FALSE)
  }

Mfile<- matrix( Mfile, ncol=11, byrow=TRUE)


### variable description
###  1 - 2                            start points
###  3 - 4                            direction - NB (0,0) for centripetal presentation
###  5 - 6                            responses to be plotted
###  7 - 8                            undefined (0's)
###  9                                intensity (decibels)
###  10                               size (mm^2)
###  11                               speed (deg/sec)

Mfile<- as.data.frame(Mfile)
names(Mfile)<- c('start1','start2','direction1','direction2',
                 'X', 'Y', 'null1','null2',
                 'intensity','size','speed')

### define labels for isopters
labels<- c('III4e', 'I4e', 'I2e', 'blind')
### patterns: 00015 = I4e;   00035 = III4e; 001015 = I2e
### anything ending in "2" is blind spot data

### chop Mfile into isopters
n.rows<- dim( Mfile)[1]
patterns.isopters<- apply(Mfile[,7:11],1, function( x){ paste(x, collapse='')})
num.isopters<- length( unique( patterns.isopters))
s1<- 2:n.rows

pos.changes<- (1:n.rows)[ patterns.isopters[s1-1] != patterns.isopters[s1]] 
pos.changes<- c(0, pos.changes, n.rows)

labels.isopters<- ifelse(patterns.isopters=='00015', 'I4e',
                  ifelse(patterns.isopters=='00035', 'III4e',
                  ifelse(patterns.isopters=='001015', 'I2e', 'blind')))
                                                              
Mfile$labels.isopters<- labels.isopters

### the output of preprocess.octopus must be num.isopters .csv files                        
list.output<-  list()
list.output[[1]]<-  Mfile[,c('X','Y', 'labels.isopters')]                                     
names(list.output)[1]<- 'Subject'
### initialises the list

for( i in 1:num.isopters)
{
obj.name<- paste(id, eye.side, 'O',Mfile$labels.isopters[pos.changes[i]+1],sep='')
s1<- (pos.changes[i]+1) : (pos.changes[i+1])
list.output[[i+1]]<- Mfile[s1, 5:6]                   
names(list.output)[i+1]<- obj.name
} ### close loop over number of isopters populating list.output


### returns the whole data.frame for inspection 
### but only if preprocess.octopus is assigned to an object
invisible( list.output)
}
