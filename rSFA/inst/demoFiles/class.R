# % class_demo1: demonstrate how classification can be done with SFA
# %
# % see [Berkes05] Pietro Berkes: Pattern recognition with Slow Feature
# % Analysis. Cognitive Sciences EPrint Archive (CogPrint) 4104, http://cogprints.org/4104/ (2005)
# %
# % Instead of handwritten digits we generate artificial (synthetic) classification data 
# %
# % INPUT
# %   opts        struct with several parameter settings for CLASSIFY. If omitted on
# %               input, meaningful defaults will be set (see source code)
# %   pars        struct with several parameter settings for training and test data generation. 
# %               If omitted on input, meaningful defaults will be set (see source code)
# % OUTPUT
# %   res         struct with fields
# %       errtrn  1 x 2 matrix: error rate with / w/o SFA on training set 
# %       errtst  1 x 2 matrix: error rate with / w/o SFA on test set 
# %
# % Author: Wolfgang Konen, Sep 2009
# %
# % See also: CLASS_DEMO2, CLASSIFY, CLASS_SWEEP
class1 <- function (opts=NULL,pars=NULL,path=NULL){ #TODO umbenennen
	if(is.null(opts)){#TODO: check for unset default values even if not null
	    opts$idim=5;         #% number of input dimensions
        opts$nclass=4;       #% number of classes
        opts$ALIGNED=1;      #% =1: align the Gaussian classifiers with coordinate axes 
                             #% (recommended, this is more robust if the data distribution deviates
                             #% largely from multivariate Gaussian). 
                             #% =0: don't align, can be more precise if data are Gaussian. 
        opts$dographics=2;   #% =0: no plots, =1: some plots, =2: more plots
        opts$epsC=0;
	}
    if(is.null(pars)){
        pars$npatt=100;      #% number of patterns per class
        pars$sigma=0.08;     #% sigma of Gaussian noise added to prototypes
        pars$generateP=1;   #% =1: generate synthetic class patterns x (& save on path/classPat.dat)
                             #% =0: load synthetic class patterns x from path/classPat.dat  
                             #% =5: generate totally random synthetic data
        pars$MIXING=1;       #% =1: nonlinear mixing of input variables, same for all classes
                             #% =2: nonlinear mixing in x1 and x2, different for all classes
                             #% other: no mixing
        pars$tstfrac=0.40;   #% which fraction of the patterns becomes test patterns
        pars$tstvarseed=1;   #% =1: vary seed, =0: always seed 0, for random number generator
                             #% in case of random divide
		opts$doPB=1;         #% =1/0: do/don't automatic parametric bootstrap, if training set is too small
	}

    
    if (pars$generateP==5){
        #% totaly random data (points taken from the opts$idim-dimensional
        #% [0,1]-cube) --> nothing can be learned
        pars$npatt=30;
        opts$idim=20;
        opts$nclass=5;
        opts$ppRange=opts$idim;    #% number of preprocessing dimensions     
        #%opts$ppRange=5;
        realclass = as.vector(customRepmat (1:opts$nclass, pars$npatt,1)) # matlab::reshape(customRepmat (1:opts$nclass, pars$npatt,1),1,opts$nclass*pars$npatt);   #% [111..222..333....]
        x = matrix(runif(pars$npatt*opts$nclass*opts$idim),pars$npatt*opts$nclass,opts$idim);
    }

    idim    =opts$idim;
    nclass  =opts$nclass;
    npatt   =pars$npatt;

    if (pars$generateP==1){

        #% for the very first time one pass with generateP==1 is necessary
        #%
        prototype = 0.2+matrix(runif(nclass*idim),nclass,idim);  #% prototype(n,:) is prototype vector for class n
        x=NULL
        for (n in 1:nclass){
			x = rbind(x, customRepmat (t(prototype[n,]),npatt,1)+matrix(rnorm(npatt*idim),npatt,idim)*pars$sigma );
        }
		tmp=customSize(x[,5:idim])
        x[,5:idim] = matrix(runif(tmp[1]*tmp[2]),tmp[1],tmp[2]);    #%distractive dimension(s)
        realclass = as.vector(customRepmat (t(1:nclass), npatt,1))#matlab::reshape(customRepmat (1:nclass, npatt,1),1,nclass*npatt);   #% [111..222..333....]   

        if (pars$MIXING==1){
            #% make the class separation a little harder in input space:
            if (idim>2){
                x[,1]=x[,1]/x[,3];
                x[,2]=x[,2]/x[,3];
            }
            if (idim>3){
                x[,2]=x[,2]/x[,4];
                x[,3]=x[,3]/x[,4];
            }
		}
        else if (pars$MIXING==2){
            #% make the class separation considerably harder in input space:
            for (n in 1:nclass){
                m=n%%idim+1;
                n2=(n+2)%%idim+1;
                m2=n2%%idim+1;
                ind=which(realclass==n);
				cclass=length(ind)
                x[ind,n]=runif(cclass);
                x[ind,m]=runif(1)+x[ind,n]^2*runif(1)+runif(cclass)*pars$sigma; 
                x[ind,n2]=runif(cclass);
                x[ind,m2]=runif(1)+x[ind,n2]^2*runif(1)+runif(cclass)*pars$sigma;
            }
        }# % pars$MIXING
		if(!is.null(path)){
			classPat=list(x=x,realclass=realclass)
			save(classPat,file=paste(path,"/classPat.dat",sep=""));
			classPat<-NULL
		}
    }
	else if (pars$generateP==0){
        load(file=paste(path,"/classPat.dat",sep=""))
		x<-classPat$x
		realclass<-classPat$realclass
		classPat<-NULL
    }
		
	tmp<-divideRand(x,realclass,pars);#[x,xtst,realclass,realcTst] = divide_rand(x,realclass,pars);

    opts$gaussdim=opts$nclass-1;        #% number of SFA-outputs used in Gauss classifier    
    opts$classes = unique(as.vector(realclass))   #% labels of the classes (we assume that each class appears
                                        #% at least once in the training data)
    if(is.null(opts$ppRange))opts$ppRange=opts$idim;        #% dims after preproc = number of input dims
    resClass = sfaClassify(tmp$x,tmp$realclass,tmp$xtst,tmp$realcTst ,opts);
	return(resClass)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % helper function for class_*demo*.m:
# % Divide data x into training set x and test set xtst.
# %
divideRand <- function(x,realclass,pars){
	x<-as.matrix(x)
    #% divide x and realclass in test and training patterns by random split
    if (pars$tstvarseed==0){              
        if (is.null(pars$seed)) pars$seed=4; 
        set.seed(pars$seed);       		 #% each time the same random numbers
    }else set.seed(proc.time()[3]*100);  # % each time different random numbers
    sz=customSize(x,1);
    tstind=matrix(0,1,sz)
    tstind[1,1:(trunc(sz*pars$tstfrac)+1)] = 1;
    tstind = tstind[sample(sz,replace=F)]; #random permutation
    res=list()
    res$xtst = x[tstind==1,];
    res$x    = x[tstind==0,];
    res$realcTst = realclass[tstind==1];
    res$realclass = realclass[tstind==0];
	return(res)
}  
