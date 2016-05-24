context("LIBLINEAR")	



test_that("LiblineaR.ACF works as expected (similar to LIBLINEAR)", {

    RSquared<-function(pred,true){
        # Pierre Gramme
        # Jan 2015
        # Compute the R-squared coefficient of determination
        
        ss_tot = sum((true-mean(true)) ^ 2)
        ss_res = sum((true-pred) ^ 2)
        return(1 - ss_res/ss_tot)
    }
    
    set.seed(1)
    N=10
    # Note that 1st record is negative
    df=data.frame(x1 = (1:N)/N*10 + 2*rnorm(N), x2 = (1:N)/N*10 + 2*rnorm(N), x3 = (1:N)/N*10 + 2*rnorm(N))
    df$y.regr = apply(as.matrix(df),1,mean) + 2*rnorm(N)
    df$y.logical = df$y.regr > 5.5
    df$y.int = ifelse(df$y.logical, 1L, -1L)
    df$y.double = as.double(df$y.int)
    df$y.char = as.character(df$y.int)
    df$y.factor = factor(df$y.int)
    df$y.factorRev = factor(df$y.int, levels=rev(levels(df$y.factor)))
    df$y.factorExtra = factor(df$y.int, levels=c(-1,1,99), labels=c("no","yes","maybe"))
    df$y.multiclass = cut(df$y.regr, breaks=c(-99,4,7,99))
    # Expected: good classification performance, negative bias and positive sum of weights 

    regrTargets = "y.regr"
    classifTargets = setdiff(grep("^y",colnames(df), value=TRUE), regrTargets)
    binTargets = setdiff(classifTargets,"y.multiclass")

    testClassif = function(rev,yy,weighted,tt) {
        #cat("Testing",rev,yy,weighted,tt,"\n")
        if(rev)
            is=1:nrow(df)
        else
            is=nrow(df):1
        nis = which(!df[is,"y.logical"])
        
        if(weighted) # class weights giving more importance to negative classes
            wi=c("1"=2,"TRUE"=2,"yes"=2, "(7,99]"=1, 
                    "(4,7]"=50,
                    "-1"=100,"FALSE"=100,"no"=100,"(-99,4]"=150)
        else
            wi=NULL
        
        y = df[is,yy]
        x = df[is,1:3]
        m = LiblineaR.ACF(x, y, type = tt, wi=wi)
        p = predict(m, newx = x)
        res=c(
            type=tt,
            target=yy,
            weighted=weighted,
            y1 = as.character(y[1]),
            perf=(mean(as.character(y)==as.character(p$predictions))), # >= 0.75
            perfNeg=(mean(as.character(y[nis])==as.character(p$predictions[nis]))), # >= 0.75
            dimW=paste(dim(m$W), collapse = " "),
            sumW=sum(m$W[,1:3]),
            biasW=m$W[1,][["Bias"]],
            classNames=paste(m$ClassNames, collapse = " "),
            yLev=paste(levels(y), collapse=" "),
            predLev=paste(levels(p$predictions), collapse=" "),
            yClass=class(y),
            predClass=class(p$predictions),
            weights=paste(colnames(m$W),"=",round(m$W[1,],3),collapse = " ; ")
        )
        
        #table(true=as.character(y),pred=as.character(p$predictions))
        
        return(res)
    }



    allRes=NULL
    for(tt in c(1,3,4)) {
        for(weighted in c(FALSE,TRUE)) {
            for(rev in c(FALSE,TRUE)) {
                for (yy in classifTargets) {
                    res = testClassif(rev,yy,weighted,tt)
                    allRes=rbind(allRes,res)
                }
            }
        }
    }

    allRes = as.data.frame(allRes,stringsAsFactors = F)
    # Simple tests
    allRes$dimOK=(allRes$type=="4" | allRes$target=="y.multiclass" | allRes$dimW=="1 4")
    allRes$perfOK=(allRes$target=="y.multiclass" & allRes$perf>=.6 | 
                                    ifelse(allRes$weighted, allRes$perfNeg>=.9, allRes$perf>=.75))
    allRes$sumOK=(!allRes$target%in%c("y.int","y.double") | allRes$type=="4" | as.numeric(allRes$sumW)>0)
    allRes$biasOK=(!allRes$target%in%c("y.int","y.double") | allRes$type=="4" | as.numeric(allRes$biasW)<0)
    allRes$levelsOK=(allRes$target%in%c("y.char","y.double") | (allRes$yLev==allRes$predLev & allRes$yClass==allRes$predClass))

    expect_equal(sum(!allRes$dimOK), 0)
    expect_equal(sum(!allRes$perfOK), 0)
    expect_equal(sum(!allRes$sumOK), 0)
    expect_equal(sum(!allRes$biasOK), 0)
    expect_equal(sum(!allRes$levelsOK), 0)
})
