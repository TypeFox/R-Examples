#------------------------------------------------------------#

### prediction in multicategory classification

pred.cz=function(f)
{
y=min(which(f==max(f)))
return(y)
} 