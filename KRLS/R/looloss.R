looloss <-
function(y=NULL,Eigenobject=NULL,lambda=NULL,eigtrunc=NULL){
            return(solveforc(y=y,Eigenobject=Eigenobject,lambda=lambda,eigtrunc=eigtrunc)$Le)
        }

