relief <-
function (data, nosample, threshold,repet=1) 
{#Modified by Edgar Acuna in october 2015
#repet greater than 1 must be  used only when nosample is less than
#the total number of instances. Otherwise, it is good enough to use 
#a 10% of the total number of instances with 5 to 10 repetitions
#It is not neccesary to imput the nominal features, the program
#identifies them.        
    if (sum(is.na(data))> 0) 
        stop("This dataset has missing values, impute them before running this function.\n",call.=FALSE)
    p=dim(data)[2]
    n=dim(data)[1]
    f=p-1
    vnom=rep(0,f)
    for(i in 1:f)
        vnom[i]=is.factor(data[,i]) 
    fc=(1:f)[vnom==0]
    fn=(1:f)[vnom==1]
    if (length(fn) > 0) {
        reliefcat(data, nosample, threshold, fn,repet)
    }
    else {
        reliefcont(data, nosample, threshold,repet)
    }
}
