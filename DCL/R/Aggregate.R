Aggregate<-function(triangle,freq=4)
{
  # check if there are complete years
  # if not then remove the last periods
  m<-nrow(triangle)/freq
  if (floor(m)-m!=0)   print('The given triangle does not consist of complete high level
  periods. The last periods have been removed')

  m<-floor(m)
  k<-m*freq
  my.triangle<-triangle
  triangle<-matrix(NA,k,k)
  for (i in 1:k)
  {
    for (j in 1:k)
    {  
      triangle[i,j]<-my.triangle[i,j] 
    }
  }
  
  my.triangle<-triangle
  triangle<-matrix(NA,m,m)
  for (rr in 1:m)
  {
    for (cc in 1:m)
    {
      dd<-rr+cc-1
      set.ind<-(1+(rr-1)*freq <=row(my.triangle) & row(my.triangle)<= freq+(rr-1)*freq) & 
        (col(my.triangle)<= freq+(cc-1)*freq) &
        ((1+(dd-1)*freq <= (row(my.triangle)+col(my.triangle)-1)) &
           ((row(my.triangle)+col(my.triangle)-1) <= freq+(dd-1)*freq))
      my.triangle.res<-as.numeric(my.triangle[set.ind])
      triangle[rr,cc]<-sum(my.triangle.res,na.rm=T)
    } 
  } 
  triangle[row(triangle)+col(triangle)>(m+1)]<-NA
  return(triangle)
  
}