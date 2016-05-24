mmnorm <-
structure(function (data, minval = 0, maxval = 1) 
{
    d = dim(data)
    c = class(data)
    cnames = colnames(data)
    classes = data[, d[2]]
    data = data[, -d[2]]
    minvect = apply(data, 2, min)
    maxvect = apply(data, 2, max)
    rangevect = maxvect - minvect
    zdata = scale(data, center = minvect, scale = rangevect)
    newminvect = rep(minval, d[2] - 1)
    newmaxvect = rep(maxval, d[2] - 1)
    newrangevect = newmaxvect - newminvect
    zdata2 = scale(zdata, center = FALSE, scale = (1/newrangevect))
    zdata3 = zdata2 + newminvect
    zdata3 = cbind(zdata3, classes)
    if (c == "data.frame") 
        zdata3 = as.data.frame(zdata3)
    colnames(zdata3) = cnames
    return(zdata3)
}, source = c("function (data,minval=0,maxval=1) ", "{", "#This is a function to apply min-max normalization to a matrix or dataframe.", 
"#Min-max normalization subtracts the minimum of an attribute from each value", 
"#of the attribute and then divides the difference by the range of the attribute.", 
"#These new values are multiplied by the given range of the attribute", 
"#and finally added to the given minimum value of the attribute.", 
"#These operations transform the data into [minval,mxval].", 
"#Usually minval=0 and maxval=1.", "#Uses the scale function found in the R base package.", 
"#Input: data= The matrix or dataframe to be scaled", "", "", 
"#store all attributes of the original data", "d=dim(data)", 
"c=class(data)", "cnames=colnames(data)", "", "#remove classes from dataset", 
"classes=data[,d[2]]", "data=data[,-d[2]]", "", "minvect=apply(data,2,min)", 
"maxvect=apply(data,2,max)", "rangevect=maxvect-minvect", "zdata=scale(data,center=minvect,scale=rangevect)", 
"", "#remove attributes added by the function scale and turn resulting", 
"#vector back into a matrix with original dimensions", "#attributes(zdata)=NULL", 
"#zdata=matrix(zdata,dim(data)[1],dim(data)[2])", "", "newminvect=rep(minval,d[2]-1)", 
"newmaxvect=rep(maxval,d[2]-1)", "newrangevect=newmaxvect-newminvect", 
"zdata2=scale(zdata,center=FALSE,scale=(1/newrangevect))", "zdata3=zdata2+newminvect", 
"", "zdata3=cbind(zdata3,classes)", "", "if (c==\"data.frame\") zdata3=as.data.frame(zdata3)", 
"colnames(zdata3)=cnames", "return(zdata3)", "", "}"))
