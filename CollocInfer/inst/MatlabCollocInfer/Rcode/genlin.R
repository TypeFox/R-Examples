

make.genlin <- function()
{
genlin<-list()
genlin$fun.ode <- function(t,y,p,more)
{
    n = ncol(y)

    more = checkmore.genlin(more,n)

    pmat = more$mat

    pmat[more$sub[,1:2,drop=F]] = pmat[more$sub[,1:2,drop=F],drop=F] + p[more$sub[,3,drop=F]]

    r = pmat%*%y

    if(!is$null(more$force)){

        fs = matrix(0,length(t),length(more$force))

        for(i in 1:length(more$force)){
            if(is.fd(more$force[[i]])){
                fs[,i] = eval.fd(t,more$force[[i]])
            }
            else{
                fs[,i] = more$force[[i]](t,more$force.input)
            }
        }

        b = more$force.mat;
        b[more$force.sub[,1:2,drop=F]] = p[more$forc.sub[,3,drop=F]];

        r = r + b%*%t(fs);
    }
    return(r)
}






genlin$fn <- function(t,y,p,more)   # We get z = Ax for pmat = A
{

    n = ncol(y)

    more = checkmore.genlin(more,n)
    pmat = more$mat


    pmat[more$sub[,1:2,drop=F]] = pmat[more$sub[,1:2,drop=F],drop=F] + p[more$sub[,3,drop=F]]

    r1 = y%*%t(pmat)


    if( !is.null(more$force) ){
    
        fs = matrix(0,length(t),length(more$force))
        
        for(i in 1:length(more$force)){
            if(is.fd(more$force[[i]])){
                fs[,i] = eval.fd(t,more$force[[i]])
            }
            else{
                fs[,i] = more$force[[i]](t,more$force_input)
            }
        }

        b = more$force.mat
        b[more$force.sub[,1:2,drop=F]] = p[more$force.sub[,3,drop=F]]

        r1 = r1 + fs%*%t(b)
    }
    return(r1)
}






genlin$dfdx <- function(t,y,p,more)
{
    n = ncol(y)

    more = checkmore.genlin(more,n)
    pmat = more$mat
    pmat[more$sub[,1:2,drop=F]] = pmat[more$sub[,1:2,drop=F],drop=F] + p[more$sub[,3,drop=F]];

    r = array(0,c(dim(y)[1],nrow(pmat),n))

    for(i in 1:nrow(pmat)){
        for(j in 1:n){
            r[,i,j] = pmat[i,j]
        }
    }
    return(r)
}







genlin$dfdp <- function(t,y,p,more)
{
    n = ncol(y)

    more = checkmore.genlin(more,n)


    r = array(0,c(dim(y)[1],nrow(more$mat),length(p)))

    if(nrow(more$sub)>0){
        for(i in 1:nrow(more$sub)){
            r[,more$sub[i,1],more$sub[i,3]] = y[,more$sub[i,2]]
        }
    }

    if(!is.null(more$force)){
    

        fs = matrix(0,length(t),length(more$force))

        for(i in 1:length(more$force)){
            if(is.fd(more$forc[[i]])){
                fs[,i] = eval.fd(t,more$force[[i]])
            }
            else{
                fs[,i] = more$force[[i]](t,more$force.input)
            }
        }

        for(i in 1:nrow(more$force.sub)){
           r[,more$force.sub[i,1],more$force.sub[i,3]] =  r[,more$force.sub[i,1],more$force.sub[i,3]] +
               fs[,more$force.sub[i,2]]
        }
    }
    return(r)
}







genlin$d2fdxdp <- function(t,y,p,more)
{
    n = ncol(y)

    more = checkmore.genlin(more,n)

    r = array(0,c(dim(y)[1],nrow(more$mat),n,length(p)))

    if(nrow(more$sub)>0){
        for(i in 1:nrow(more$sub)){
            r[,more$sub[i,1],more$sub[i,2],more$sub[i,3]] = 1
        }
    }
    return(r)
}

genlin$d2fdx2 <- function(t,y,p,more)
{
    more = checkmore.genlin(more,n)
    n = ncol(y)
    r = array(0,c(dim(y)[1],nrow(more$mat), n, n))
}

checkmore.genlin <- function(more,n)    # checks additional arguments to genlin
{
    if(is.null(more)){
        more$sub = cbind(kronecker(matrix(1:n,n,1),matrix(1,n,1)),
            kronecker(matrix(1,n,1),matrix(1:n,n,1)),1:n^2)  
        more$mat = matrix(0,n,n)      
        }
    else{
        if(is.null(more$sub)){
            more$sub = cbind(kronecker(matrix(1:n,n,1),matrix(1,n,1)),
                kronecker(matrix(1,n,1),matrix(1:n,n,1)),1:n^2)
        }
        else{
            if(ncol(more$sub)==1){                      # can just specify which params to use
                more$sub = cbind(kronecker(matrix(1:n,n,1),matrix(1,n,1)),
                    kronecker(matrix(1,n,1),matrix(1:n,n,1)),more$sub)
            }
        }
        if(is.null(more$mat)){ more$mat = matrix(0,n,n) }
    
        if(!is.null(more$force)){
            m = length(more$force)
            if(is.null(more$force.mat)){ more$force.mat = matrix(0,n,m) }
            if(is.null(more$force.sub)){
                more$force.sub = cbind(kronecker(matrix(1:n,m,1),matrix(1,m,1)),
                    kronecker(matrix(1,n,1),matrix(1:m,n,1)),1:(n*m))
            }
            if(ncol(more$force.sub)==1 | is.vector(more$force.sub)){
                more$force.sub = 
                    cbind(more$force.sub,more$force.sub,nrow(more$sub)+length(more$force.sub))
            }
            else{                               # one force per collumn, and specify which params
                if(ncol(more$force.sub)==2){
                    more$force.sub = cbind(more$force.sub[,1],more$force.sub[,1],more$force.sub[,2])
                }
            }
        }
    }   
    return(more)
}



    return(genlin)
}


