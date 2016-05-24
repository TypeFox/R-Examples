kFsubj <-
function(obj.name,  perimeter='G', 
         no.cleaning=TRUE,
         no.kprm=TRUE, no.flip=TRUE)
{

### obj.name is the file name using format defined for studyid
### note that argument obj.name only contains e.g. "M2005R"
### perimeter must be either "G" (Goldmann) or "O" (Octopus)
### if no.cleaning==TRUE it assumes a pre-ordered dataset (Goldmann)
### if no.cleaning==FALSE it assumes no ordered (Octopus data)
###    and starts cleaning procedure
### if no.kprm==FALSE then there's no kinetic perimetry
###    reliability measure (KPRM)
### if no.flip==TRUE then display it as the original left/right eye
### if no.flip==FALSE then flip and display as right eye
### Assumes that objects obj.name***** exist in the working directory
### the format is, subject_eye_perimeter_stimulus, e.g.  M0001RGI4e

if(!no.cleaning) print("You're now starting data cleaning")

if(perimeter!='G') set.template(void=FALSE)


if(perimeter=='G')set.template()

dists<- areas<- rep(NA, 4)

if(perimeter!= 'O' & perimeter!= 'G') stop('perimeter value is invalid')

if(perimeter=='G')

{


obj.name<- paste(obj.name, as.character(perimeter), sep='')
if (obj.name=='M0001R')
{
data(M0001RGIII4e, envir = environment())
data(M0001RGI4e, envir = environment())
data(M0001RGblind, envir = environment())

obj<- ls('M0001RGblind', 'M0001RGI4e', 'M0001RGIII4e')

}
else obj<- ls(envir=.GlobalEnv)[grep(obj.name, ls(envir=.GlobalEnv))]


num.obj<- length( obj)

files<- files0<- c('III4e','I4e', 'I2e','blind','rmeas')

files.aux<- substring(obj, 8, nchar(obj))
m0<- match(files.aux, files)

files<- paste(obj.name, files, sep='')

which.files<- rep(FALSE, 5)
for(i in 1:5) which.files[i]<- length(grep(files[i], obj)) > 0

m1<- match( substring(obj, 8, nchar(obj)), files0)

files00<- files0
files0<-  paste( obj.name, files0, sep='')

files1<- files0[sort(m1)]
m2<- match( obj, files0)




III4e<- if(exists(files0[1]))assign(files0[1], get(obj[match(files0[1],obj)]))
I4e<- if(exists(files0[2]))assign(files0[2], get(obj[match(files0[2],obj)]))
I2e<- if(exists(files0[3]))assign(files0[3], get(obj[match(files0[3],obj)]))
blind<- if(exists(files0[4]))assign(files0[4], get(obj[match(files0[4],obj)]))
rmeas<- if(exists(files0[5]))assign(files0[5], get(obj[match(files0[5],obj)]))

}



if(perimeter=='O')

{
obj.name<- paste(obj.name, as.character(perimeter), sep='')
obj<- ls(envir=.GlobalEnv)[grep(obj.name, ls(envir=.GlobalEnv))]
num.obj<- length( obj)

files<- files0<- c('III4e','I4e', 'I2e','blind','rmeas')

if(perimeter=='O' & !no.cleaning) 
 {
    for (i in 1:5) files[i]<- paste(files0[i], 'raw',sep='')
which.files<- rep(FALSE, 5)
for(i in 1:5) which.files[i]<- length(grep(files[i], obj)) > 0
m1<- match( substring(obj, 8, nchar(obj)-3), files0)
files00<- files0
    files0<-  paste( obj.name, files0, sep='')
files2<- paste(files0, 'raw',sep='')

III4e<- if(exists(files2[1]))assign(files2[1], get(obj[match(files2[1],obj)]))
I4e<- if(exists(files2[2]))assign(files2[2], get(obj[match(files2[2],obj)]))
I2e<- if(exists(files2[3]))assign(files2[3], get(obj[match(files2[3],obj)]))
blind<- if(exists(files2[4]))assign(files2[4], get(obj[match(files2[4],obj)]))
rmeas<- if(exists(files2[5]))assign(files2[5], get(obj[match(files2[5],obj)]))

 }
if(perimeter=='O' & no.cleaning) 
 {
    for (i in 1:5) files[i]<- paste(files0[i], 'proc',sep='')
which.files<- rep(FALSE, 5)
for(i in 1:5) which.files[i]<- length(grep(files[i], obj)) > 0
m1<- match( substring(obj, 8, nchar(obj)-4), files0)
    files00<- files0
    files0<-  paste( obj.name, files0, sep='')
files2<- paste(files0, 'proc',sep='')

III4e<- if(exists(files2[1]))assign(files2[1], get(obj[match(files2[1],obj)]))
I4e<- if(exists(files2[2]))assign(files2[2], get(obj[match(files2[2],obj)]))
I2e<- if(exists(files2[3]))assign(files2[3], get(obj[match(files2[3],obj)]))
blind<- if(exists(files2[4]))assign(files2[4], get(obj[match(files2[4],obj)]))
rmeas<- if(exists(files2[5]))assign(files2[5], get(obj[match(files2[5],obj)]))

 }

}



### get objects' names associated with obj.name






#for(i in 1:num.obj) assign(files1[i], get(obj[i]))


if(which.files[1])
{
III4e <- rbind(III4e, III4e[1,])
III4e <- cbind(III4e[,1], III4e[,2])
}

if(which.files[2])
{
I4e <- rbind(I4e, I4e[1,])
I4e <- cbind(I4e[,1], I4e[,2])
}

if(which.files[3])
{  
I2e <- rbind(I2e, I2e[1,])
I2e <- cbind(I2e[,1], I2e[,2])
}

if(which.files[4])
{   
blind <- rbind(blind, blind[1,])
blind <- cbind(blind[,1], blind[,2])   
}  


side.eye<- substring(obj.name, 6,6) ## must be L or R only


if(no.cleaning) ### starts Goldmann
{


if(side.eye=='L' & !no.flip)
{
   if( which.files[1]) III4e[,1]<- -III4e[,1]
   if( which.files[2]) I4e[,1]<- -I4e[,1]
   if( which.files[3]) I2e[,1]<- -I2e[,1]
   if( which.files[4]) blind[,1]<- -blind[,1]
   if( which.files[5]) rmeas[,1] <- -rmeas[,1]
}


   if(which.files[1])
   {
   ##III4e 
   points(III4e, pch=4, col='black', lwd=2, cex=2/3)
   if(no.cleaning) lines(III4e, col= "black", lwd=2)
   areas[1]<- areapl(III4e)
   }

   if(which.files[2])
{
##I4e
points(I4e, pch=19, col='royalblue3')
if(no.cleaning) lines(I4e, col= "royalblue3", lwd=2)
areas[2]<- areapl(I4e)
}

if(which.files[3])
{
##I2e
points(I2e, pch=15, col='maroon', lwd=2, cex=2/3)
if(no.cleaning) lines(I2e, col="maroon", lwd=2)
areas[3]<- areapl(I2e)
}

if(which.files[4])
{
##Blind spot
points(blind, pch=4, col='maroon', lwd=2, cex=1/3)
if(no.cleaning) lines(blind, col="maroon", lwd=2)
areas[4]<- areapl(blind)
}

if(which.files[5] & !no.kprm)
{
## Reliability measure points
points(rmeas, pch=2, cex=2/3, col="red", lwd=2)
print(
cat('please identify four points on the plot','\n',
   'starting in the upper-left quadrant working clockwise','\n'),
quote=FALSE)

if(which.files[1])
{
qc.outer<- III4e[identify( III4e, n=4), ]
test<-dist( rbind( data.matrix(rmeas), qc.outer))
dists<-diag(dist2full( test)[1:4, 5:8]) ### four distances
}
else
{
qc.outer<- I4e[identify( I4e, n=4), ]
test<-dist( rbind( data.matrix(rmeas), qc.outer))
dists<-diag(dist2full( test)[1:4, 5:8]) ### four distances
}
if(side.eye=='L') dists[c(2,1, 4,3)]<- dists
}

names(areas)<- c('III4e','I4e','I2e','blind spot')
names(dists)<- c('superior nasal','superior temporal',
                 'inferior temporal','inferior nasal')
                      

mtext(side=3, line=2, paste('ID=',obj.name,sep=''), cex=1.5)
title(sub=
cat( paste(names(areas), round(areas,1), sep=' = '), sep=' '))
print(' ',quote=FALSE)
invisible(c(round(areas,1), round(dists,2)) )#### output - 4 areas and 4 dists
} ### stops Goldmann

else ### starts Octopus
{ 
output.list<- list(III4e, I4e, I2e, blind, rmeas)
names(output.list)<- c('III4e','I4e','I2e','blind', 'rmeas')
print(output.list)


for (i in 1:5) 
{
   if(!is.null(output.list[[i]])) 
   { 
     mat.aux<- output.list[[i]]
     points(mat.aux[,1], mat.aux[,2],col=i, cex=1, pch=19)
                          
               mat.aux2<- delete.points(mat.aux)
               print(mat.aux2)
               mat.aux3<- stop.identify(mat.aux2) ### ordered list of chosen points
               print(mat.aux3)
               mat.aux3<-as.data.frame(mat.aux2[mat.aux3, ]) ### 
               names(mat.aux3)<-c('X','Y') 
               output.list[[i]]<- mat.aux3
    }
}## closes the loop
invisible(output.list)
} ### closes if then else for no.cleaning

}
