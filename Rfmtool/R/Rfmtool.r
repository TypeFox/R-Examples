# Rfmtool Package v1.1


fm <- function()
{
    # This function outputs a list of all functions included in this toolbox.
    
	print("The list of functions in Rfmtool Tool Box:")
	
	print("fm.Init([number of variables])")
	print("fm.Banzhaf([general fuzzy measure],[environment])")
	print("fm.BanzhafMob([mobius fuzzy measure],[environment])")
	print("fm.Choquet([input criteria], [fuzzy measure (general represenation)],[environment])")
	print("fm.ChoquetMob([input criteria], [fuzzy measure (mobius represenation),[environment]])")
	print("fm.ConstructLambdaMeasure([singletons (array of size n of the values of fuzzy measure at singletons)],[environment])")
	print("fm.ConstructLambdaMeasureMob([singletons (array of size n of the values of fuzzy measure at singletons)],[environment])")
	print("fm.dualm([fuzzy measure (general represenation)],[environment])")
	print("fm.dualmMob([fuzzy measure (mobius represenation)],[environment])")
	print("fm.EntropyChoquet([fuzzy measure (general represenation)],[environment])")
	print("fm.EntropyChoquetMob([fuzzy measure (mobius represenation)])")
	print("fm.fitting <- function([emprical data], [k-additive]) - returns Mobius fuzzy measure")
	print("fm.fittingMob <- function([emprical data], [k-additive])")
	print("fm.FuzzyMeasureFitLP <- function([emprical data], [k-additive], [other options]) - returns standard fuzzy measure")
	print("fm.FuzzyMeasureFitLPMob <- function([emprical data], [k-additive], [other options]) - returns Mobius fuzzy measure")
	print("fm.Interaction([standard fuzzy measure],[environment])")
	print("fm.InteractionMob([mobius fuzzy measure],[environment])")
	print("fm.InteractionB([standard fuzzy measure],[environment])")
	print("fm.InteractionBMob([mobius fuzzy measure],[environment])")
	print("fm.IsMeasureAdditive([fuzzy measure (general represenation)],[environment])")
	print("fm.IsMeasureAdditiveMob([fuzzy measure (mobius represenation),[environment]])")
	print("fm.IsMeasureBalanced([fuzzy measure (general represenation),[environment]])")
	print("fm.IsMeasureBalancedMob([fuzzy measure (mobius represenation)],[environment])")
	print("fm.IsMeasureSelfdual([fuzzy measure (general represenation)],[environment])")
	print("fm.IsMeasureSelfdualMob([fuzzy measure (mobius represenation)],[environment])")
	print("fm.IsMeasureSubadditive([fuzzy measure (general represenation)],[environment])")
	print("fm.IsMeasureSubadditiveMob([fuzzy measure (mobius represenation)],[environment])")
	print("fm.IsMeasureSubmodular([fuzzy measure (general represenation)],[environment])")
	print("fm.IsMeasureSubmodularMob([fuzzy measure (mobius represenation)],[environment])")
	print("fm.IsMeasureSuperadditive([fuzzy measure (general represenation)],[environment])")
	print("fm.IsMeasureSuperadditiveMob([fuzzy measure (mobius represenation)],[environment])")
	print("fm.IsMeasureSupermodular([fuzzy measure (general represenation)],[environment])")
	print("fm.IsMeasureSupermodularMob([fuzzy measure (mobius represenation)],[environment])")
	print("fm.IsMeasureSymmetric([fuzzy measure (general represenation)],[environment])")
	print("fm.IsMeasureSymmetricMob([fuzzy measure (mobius represenation)],[environment])")
	print("fm.Mobius([general fuzzy measure],[environment])")
	print("fm.OrnessChoquet([fuzzy measure (standard represenation)],[environment])")
	print("fm.OrnessChoquetMob([fuzzy measure (mobius represenation)],[environment])")
	print("fm.Shapley([general fuzzy measure],[environment])")
	print("fm.ShapleyMob([mobius fuzzy measure],[environment])")
	print("fm.Sugeno([input criteria], [fuzzy measure (general represenation)],[environment])")
	print("fm.SugenoMob([input criteria], [fuzzy measure (mobius represenation)],[environment])")
	print("fm.test()")
	print("fm.Zeta([mobius fuzzy measure],[environment])")	
}


fm.Init <- function(n1)
{
    # Calculates an array of Banzhaf indices

    n<-as.integer(n1);
    m1 <-2^n1;

    out<-.C("Preparations_FMCall",n=as.integer(n), m=as.integer(m1), card=as.integer(1:m1),cardpos=as.integer(1:(n+1)),
	bit2card=as.integer(1:m1),card2bit=as.integer(1:m1), factorials=as.double(1:(n+1))
 );
 
					
    return (out);
}

fm.errorcheck <- function(env=NULL)
{
	if(is.null(env)) return(TRUE);
	if(env$n!=log2(env$m)) return(TRUE);
	if((env$n+1) !=length(env$cardpos)) return(TRUE);
	if(env$m!=length(env$card)) return(TRUE);
	if(env$m!=length(env$card2bit)) return(TRUE);
	if(env$m!=length(env$bit2card)) return(TRUE);

	if((env$n+1)!=length(env$factorials)) return(TRUE);

	return(FALSE);
}

fm.Banzhaf <- function(v,env=NULL)
{
    # Calculates an array of Banzhaf indices
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}

    BanzhafVal <- array(0,log2(length(v)));

	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}


	BanzhafValue <- .C("BanzhafCall", as.numeric(v), 
        out = as.numeric(BanzhafVal),
        as.integer(log2(length(v))), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
					
    return (round(BanzhafValue$out, digits=4));
}


fm.BanzhafMob <- function(Mob,env=NULL)
{
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    # Calculates an array of Banzhaf indices for Mobius fuzzy measure
    v = fm.Zeta(Mob,env);

    BanzhafMobVal <- array(0,log2(length(v)));
	BanzhafMobValue <- .C("BanzhafCall", as.numeric(v), 
        out = as.numeric(BanzhafMobVal),
        as.integer(log2(length(v))), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
					
    return (round(BanzhafMobValue$out, digits=4));
}


fm.Choquet <- function(x, v,env=NULL)
{
    # Calculates the value of a discrete Choquet integral of the input x, with provided fuzzy measure v (in general representation)
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}	
	if(env$m!=length(v)|| env$n!=length(x)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    ChoquetVal <- -1; #this is just a initial value.
    ChoquetValue <- .C("ChoquetCall", as.numeric(x),
        as.numeric(v),
        as.integer(length(x)),
        out = as.numeric(ChoquetVal), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
    return (ChoquetValue$out);
}


fm.ChoquetMob <- function(x, Mob,env=NULL)
{
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)||env$n!=length(x)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

	# This is an alternative calculation of the Choquet integral from the Mobius fuzzy measure.
    v = fm.Zeta(Mob,env);

    ChoquetVal <- -1; #this is just a initial value.
	ChoquetMobValue <- .C("ChoquetMobCall", as.numeric(x),
					        as.numeric(Mob),
				 	        as.integer(length(x)),
			                  out = as.numeric(ChoquetVal), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
	return (ChoquetMobValue$out);
}


fm.ConstructLambdaMeasure <- function(singletons,env=NULL)
{
	# Finds the value of lambda and calculates the rest of the values of the fuzzy measure,
    # given its values at singletons. singletons is an array of size n.
    # lambda and v are the outputs, v is in standard representation and binary ordering 
    # (array v of size m should be allocated by the calling routine).
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}	
	if(env$n!=length(singletons)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    lambda <- array(-1, 1);      # initial value of lambda: array of length 1 with value -1
    v <- array(0, 2^length(singletons));   # array of m zeros
	ConstructLambdaMeasureValue <- .C("ConstructLambdaMeasureCall", 
        as.numeric(singletons),
        out1 = as.numeric(lambda),
        out2 = as.numeric(v),
        as.integer(length(singletons)), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );

    return(list(lambda=ConstructLambdaMeasureValue$out1, measure=ConstructLambdaMeasureValue$out2));
}


fm.ConstructLambdaMeasureMob <- function(singletons,env=NULL)
{
	# Finds the value of lambda and calculates the rest of the values of the fuzzy measure,
    # given its values at singletons. singletons is an array of size n.
    # lambda and Mob are the outputs, Mob is in standard representation and binary ordering 
    # (array Mob of size m should be allocated by the calling routine).
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}	
	if(env$n!=length(singletons)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

	out<-fm.ConstructLambdaMeasure(singletons,env);

	out$measure<-fm.Mobius(out$measure,env);

    return(out);
}


fm.dualm <- function(v,env=NULL)
{
    # Calculates the dual of fuzzy measure v, returns it as value of the function (array of size m).
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    dualmVal <- array(0,length(v));  # array of m zeros
    dualmValue <- .C("dualmCall", 
        as.numeric(v),
        out = as.numeric(dualmVal),
        #as.integer(log2(length(v))),
        as.integer(length(v)), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
    return (dualmValue$out);
}


fm.dualmMob <- function(Mob,env=NULL)
{
    # Calculates the dual of fuzzy measure Mob in Mobius representation, returns it as value of the function (array of size m).
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}

	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    v = fm.Zeta(Mob,env);

    dualmVal <- array(0,length(v));  # array of m zeros
    dualmValue <- .C("dualmCall", 
        as.numeric(v),
        out = as.numeric(dualmVal),
        #as.integer(log2(length(v))),
        as.integer(length(v)), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
    return (fm.Mobius(dualmValue$out,env));
}


fm.EntropyChoquet <- function(v,env=NULL)
{
    # Calculates entropy value of the Choquet integral of fuzzy measure in general representation. 
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}	
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}


    EntropyChoquetVal <- -1;  # this is just a initial value.
    EntropyChoquetValue <- .C("EntropyChoquetCall", 
        as.numeric(v),
        as.integer(log2(length(v))),
        out = as.numeric(EntropyChoquetVal), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
    return (EntropyChoquetValue$out);
}


fm.EntropyChoquetMob <- function(Mob,env=NULL)
{
    # Calculates entropy value of the Choquet integral of fuzzy measure in general representation.
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	} 

	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    v = fm.Zeta(Mob,env);
	
    EntropyChoquetVal <- -1;  # this is just a initial value.
    EntropyChoquetValue <- .C("EntropyChoquetCall", 
        as.numeric(v),
        as.integer(log2(length(v))),
        out = as.numeric(EntropyChoquetVal), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
    return (EntropyChoquetValue$out);
}


fm.fitting<- function(data, env=NULL, kadd="NA")
{
	# This function estimates the values of a k-additive standard fuzzy measure based on empirical data. 
	# The result is an array containing the values of the fuzzy measure in Mobius, ordered according to set cardinalities.
	# kadd define the complexity of fuzzy measure. if kadd is not provided, its default value is equal to the number of inputs.

	size <- dim(as.matrix(data));
	n <- size[2] - 1;
	datanum <- size[1];
	m = 2^n;
	MobiusVal <- array(0,m);


	if (kadd == "NA") 
	{
		kadd = n;
    }
  	
    MobiusValue <- .C("fittingCall", as.integer(n),
        as.integer(datanum),
        as.integer(kadd),
        out = as.numeric(MobiusVal),
        as.numeric(t(data))
    );
	MobiusValue$out<-fm.Zeta(MobiusValue$out,env)				
	return (MobiusValue$out);
}


fm.fittingMob<- function(data, env=NULL, kadd="NA")
{
	# This function estimates the values of a k-additive Mobius fuzzy measure based on empirical data. 
	# The result is an array containing the values of the fuzzy measure in Mobius, ordered according to set cardinalities.
	# kadd define the complexity of fuzzy measure. if kadd is not provided, its default value is equal to the number of inputs.

	size <- dim(as.matrix(data));
	n <- size[2] - 1;
	datanum <- size[1];
	m = 2^n;
	MobiusVal <- array(0,m);
	
	if (kadd == "NA") 
	{
		kadd = n;
    }
  	
    MobiusValue <- .C("fittingCall", as.integer(n),
        as.integer(datanum),
        as.integer(kadd),
        out = as.numeric(MobiusVal),
        as.numeric(t(data))
    );
					
	return (MobiusValue$out);
}


fm.FuzzyMeasureFitLP <- function(data, env=NULL, kadd="NA", 
        options=0, indexlow=(NULL), indexhigh=(NULL) , option1=0, orness=(NULL))
{
	# This function estimates the values of a k-additive fuzzy measure based on empirical data. 
	# The result is an array containing the values of a standard fuzzy measure, ordered according to set cardinalities.
	# kadd define the complexity of fuzzy measure. if kadd is not provided, its default value is equal to the number of inputs.
    # int FuzzyMeasureFitLP(int n, int m, int K, int Kadd, double *v, double* XYData, int options=0, 
    #    double* indexlow=(NULL), double* indexhigh=(NULL) , int option1=0, double* orness=(NULL));
    # Input parameters: 
    # n - the dimension of inputs, m = 2^n - the number of fuzzy measure values
    # K - the number of empirical data
    # Kadd - k in k-additive f. measures, 1 < Kadd < n+1. Kdd=n - f.m. is unrestricted
    # XYData - an array of size K x (n+1), where each row is the pair (x,y), K data altogether
    # options (default value is 0)
    #    1 - lower bounds on Shapley values supplied in indexlow
    #    2 - upper bounds on Shapley values supplied in indexhigh
    #    3 - lower and upper bounds on Shapley values supplied in indexlow and indexhigh
    #    4 - lower bounds on all interaction indices supplied in indexlow
    #    5 - upper bounds on all interaction indices supplied in indexhigh
    #    6 - lower and upper bounds on all interaction indices supplied inindexlow and indexhigh
    #    all these value will be treated as additional constraints in the LP
    # indexlow, indexhigh - array of size n (options =1,2,3) or m (options=4,5,6)
    # containing the lower and upper bounds on the Shapley values or interaction indices
    # Example of orness in C:
	# double orness[2];
	# orness[0]=0; 
	# orness[1]=1;

	size <- dim(as.matrix(data));
	n <- size[2] - 1;
	datanum <- size[1];
	m = 2^n;
	MobiusVal <- array(0,m);
	
	if (kadd == "NA") 
	{
		kadd = n;
    }
  	
    MobiusValue <- .C("FuzzyMeasureFitLPCall", as.integer(n),
        as.integer(datanum),
        as.integer(kadd),
        out = as.numeric(MobiusVal),
        as.numeric(t(data)),
        as.integer(options), 
        as.numeric(indexlow), 
        as.numeric(indexhigh), 
        as.integer(option1), 
        as.numeric(orness)
    );
#	print(	MobiusValue );
			
	return (fm.Zeta(MobiusValue$out,env));
}


fm.FuzzyMeasureFitLPMob <- function(data, env=NULL, kadd="NA", 
        options=0, indexlow=(NULL), indexhigh=(NULL) , option1=0, orness=(NULL))
{
	# This function estimates the values of a k-additive fuzzy measure based on empirical data. 
	# The result is an array containing the values of the fuzzy measure in Mobius, ordered according to set cardinalities.
	# kadd define the complexity of fuzzy measure. if kadd is not provided, its default value is equal to the number of inputs.
    # int FuzzyMeasureFitLP(int n, int m, int K, int Kadd, double *v, double* XYData, int options=0, 
    #    double* indexlow=(NULL), double* indexhigh=(NULL) , int option1=0, double* orness=(NULL));
    # Input parameters: 
    # n - the dimension of inputs, m = 2^n - the number of fuzzy measure values
    # K - the number of empirical data
    # Kadd - k in k-additive f. measures, 1 < Kadd < n+1. Kdd=n - f.m. is unrestricted
    # XYData - an array of size K x (n+1), where each row is the pair (x,y), K data altogether
    # options (default value is 0)
    #    1 - lower bounds on Shapley values supplied in indexlow
    #    2 - upper bounds on Shapley values supplied in indexhigh
    #    3 - lower and upper bounds on Shapley values supplied in indexlow and indexhigh
    #    4 - lower bounds on all interaction indices supplied in indexlow
    #    5 - upper bounds on all interaction indices supplied in indexhigh
    #    6 - lower and upper bounds on all interaction indices supplied inindexlow and indexhigh
    #    all these value will be treated as additional constraints in the LP
    # indexlow, indexhigh - array of size n (options =1,2,3) or m (options=4,5,6)
    # containing the lower and upper bounds on the Shapley values or interaction indices
    # Example of orness in C:
	# double orness[2];
	# orness[0]=0; 
	# orness[1]=1;

	size <- dim(as.matrix(data));
	n <- size[2] - 1;
	datanum <- size[1];
	m = 2^n;
	MobiusVal <- array(0,m);
	
	if (kadd == "NA") 
	{
		kadd = n;
    }
  	
    MobiusValue <- .C("FuzzyMeasureFitLPCall", as.integer(n),
        as.integer(datanum),
        as.integer(kadd),
        out = as.numeric(MobiusVal),
        as.numeric(t(data)),
        as.integer(options), 
        as.numeric(indexlow), 
        as.numeric(indexhigh), 
        as.integer(option1), 
        as.numeric(orness)
    );
					
	return (MobiusValue$out);
}


fm.Interaction <- function(v,env=NULL)
{
	# calculates all interaction indices 
	# result is a matrix, whose first column is the interaction index
	# and second column is the index of coliation.
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}

	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    Mob <- fm.Mobius(v,env);

    coliation <- array(0,length(Mob));
    InteractionVal <- array(0,length(Mob));
    InteractionValue <- .C("InteractionCall", as.numeric(Mob), 
 		                   inter = as.numeric(InteractionVal),
						   as.integer(log2(length(Mob))),
					   colia = as.integer(coliation), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
    inteIndex <- as.matrix(InteractionValue$inter);
    coliaIndex <- as.matrix(InteractionValue$colia);
    index <- cbind(inteIndex,coliaIndex); 				
    return (round(index, digits=4));
}


fm.InteractionMob <- function(Mob,env=NULL)
{
	# calculates all interaction indices 
	# result is a matrix, whose first column is the interaction index
	# and second column is the index of coliation.
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

	 coalition <- array(0,length(Mob));
	 InteractionVal <- array(0,length(Mob));
	 InteractionValue <- .C("InteractionCall", as.numeric(Mob), 
			                   inter = as.numeric(InteractionVal),
						   as.integer(log2(length(Mob))),
					   colia = as.integer(coalition), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
	inteIndex <- as.matrix(InteractionValue$inter);
	coliaIndex <- as.matrix(InteractionValue$colia);
	index <- cbind(inteIndex,coliaIndex); 				
	return (round(index, digits=4));
}


fm.InteractionB <- function(v,env=NULL)
{
	# calculates all InteractionB indices 
	# result is a matrix, whose first column is the InteractionB index
	# and second column is the index of coliation.
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    Mob <- fm.Mobius(v,env);

    coalition <- array(0,length(Mob));
    InteractionBVal <- array(0,length(Mob));
    InteractionBValue <- .C("InteractionBCall", as.numeric(Mob), 
		                   inter = as.numeric(InteractionBVal),
						   as.integer(log2(length(Mob))),
                           colia = as.integer(coalition), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
	inteIndex <- as.matrix(InteractionBValue$inter);
	coliaIndex <- as.matrix(InteractionBValue$colia);
	index <- cbind(inteIndex,coliaIndex); 				
	return (round(index, digits=4));
}


fm.InteractionBMob <- function(Mob,env=NULL)
{
	# calculates all InteractionB indices 
	# result is a matrix, whose first column is the InteractionB index
	# and second column is the index of coliation.
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    coalition <- array(0,length(Mob));
    InteractionBVal <- array(0,length(Mob));
    InteractionBValue <- .C("InteractionBCall", as.numeric(Mob), 
		                   inter = as.numeric(InteractionBVal),
						   as.integer(log2(length(Mob))),
	                       colia = as.integer(coalition), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
	inteIndex <- as.matrix(InteractionBValue$inter);
	coliaIndex <- as.matrix(InteractionBValue$colia);
	index <- cbind(inteIndex,coliaIndex); 				
	return (round(index, digits=4));
}


fm.IsMeasureAdditive <- function(v,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    # v is a fuzzy measure in standard representation.
	result <- 1;
    res <- .C("IsMeasureAdditiveCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureAdditiveMob <- function(Mob,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}

	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    v = fm.Zeta(Mob,env);
	result <- 1;

    # v is a fuzzy measure in standard representation.
    res <- .C("IsMeasureAdditiveCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureBalanced <- function(v,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    # v is a fuzzy measure in standard representation.
	result <- 1;
    res <- .C("IsMeasureBalancedCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureBalancedMob <- function(Mob,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    # Mob is a fuzzy measure in Mobius representation.
    v = fm.Zeta(Mob,env);

	result <- 1;
    res <- .C("IsMeasureBalancedCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSelfdual <- function(v,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    # v is a fuzzy measure in standard representation.
	result <- 1;
    res <- .C("IsMeasureSelfdualCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSelfdualMob <- function(Mob,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    # Mob is a fuzzy measure in Mobius representation.
    v = fm.Zeta(Mob,env);

	result <- 1;
    res <- .C("IsMeasureSelfdualCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSubadditive <- function(v,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    # v is a fuzzy measure in standard representation.
	result <- 1;
    res <- .C("IsMeasureSubadditiveCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSubadditiveMob <- function(Mob,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    # v is a fuzzy measure in Mobius representation.
    v = fm.Zeta(Mob,env);

	result <- 1;
    res <- .C("IsMeasureSubadditiveCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSubmodular <- function(v,env=NULL)
{
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
	# Returns 1 if yes, 0 if no;
    # v is a fuzzy measure in standard representation.
	result <- 1;
    res <- .C("IsMeasureSubmodularCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSubmodularMob <- function(Mob,env=NULL)
{
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
	# Returns 1 if yes, 0 if no;
    # v is a fuzzy measure in standard representation.
    v = fm.Zeta(Mob,env);

	result <- 1;
    res <- .C("IsMeasureSubmodularCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSuperadditive <- function(v,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    # v is a fuzzy measure in standard representation.
	result <- 1;
    res <- .C("IsMeasureSuperadditiveCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSuperadditiveMob <- function(Mob,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    # v is a fuzzy measure in standard representation.
    v = fm.Zeta(Mob,env);

	result <- 1;
    res <- .C("IsMeasureSuperadditiveCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSupermodular <- function(v,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    # v is a fuzzy measure in standard representation.
	result <- 1;
    res <- .C("IsMeasureSupermodularCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSupermodularMob <- function(Mob,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    # Mob is a fuzzy measure in Mobius representation.
	result <- 1;
    v = fm.Zeta(Mob,env);

    res <- .C("IsMeasureSupermodularCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.IsMeasureSymmetric <- function(v, env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    # v is a fuzzy measure in standard representation.
	result <- 1;
    res <- .C("IsMeasureSymmetricCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

       );

	return (as.logical(res$result));
}


fm.IsMeasureSymmetricMob <- function(Mob,env=NULL)
{
	# Returns 1 if yes, 0 if no;
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    # Mob is a fuzzy measure in Mobius representation.
    v = fm.Zeta(Mob,env);
	result <- 1;
    res <- .C("IsMeasureSymmetricCall", 
        as.numeric(v), 
        as.integer(length(v)), result=as.integer(result), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (as.logical(res$result));
}


fm.Mobius <- function(v,env=NULL)
{
    # Calculates Mobius representation of the general fuzzy measure v
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    MobiusVal <-  array(0,length(v));
    MobiusValue <- .C("MobiusCall", as.numeric(v), 
        out = as.numeric(MobiusVal),
        as.integer(env$n), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );			
    return (MobiusValue$out);
}


fm.OrnessChoquet <- function(v,env=NULL)
{
	# Calculates the orness value of the Choquet integral for a standard fuzzy measure.
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
	Mob = fm.Mobius(v,env);

    OrnessChoquetMobVal <- -1;  # this is just a initial value.
	OrnessChoquetMobValue <- .C("OrnessChoquetMobCall", 
        as.numeric(Mob),
        as.integer(log2(length(Mob))),
        out = as.numeric(OrnessChoquetMobVal), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (OrnessChoquetMobValue$out);
}


fm.OrnessChoquetMob <- function(Mob,env=NULL)
{
	# Calculates the orness value of the Choquet integral for the Mobius fuzzy measure.
		if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    OrnessChoquetMobVal <- -1;  # this is just a initial value.
	OrnessChoquetMobValue <- .C("OrnessChoquetMobCall", 
        as.numeric(Mob),
        as.integer(log2(length(Mob))),
        out = as.numeric(OrnessChoquetMobVal), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
	return (OrnessChoquetMobValue$out);
}


fm.Shapley<- function(v,env=NULL)
{
    # Calculates an array of Shapley values.
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    ShapleyVal <- array(0,log2(length(v)));
    ShapleyValue <- .C("ShapleyCall", as.numeric(v), 
        out = as.numeric(ShapleyVal),
        as.integer(log2(length(v))), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
					
    return (round(ShapleyValue$out, digits=4));
}


fm.ShapleyMob<- function(Mob,env=NULL)
{
    # Calculates an array of Shapley values.
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
    v = fm.Zeta(Mob,env);

    ShapleyVal <- array(0,log2(length(v)));
    ShapleyValue <- .C("ShapleyCall", as.numeric(v), 
        out = as.numeric(ShapleyVal),
        as.integer(log2(length(v))), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)

    );
					
    return (round(ShapleyValue$out, digits=4));
}


fm.Sugeno <- function(x, v,env=NULL)
{
    # Calculates the value of a Sugeno integral of an input x, with provided fuzzy measure v (in general representation). 
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(v)||env$n!=length(x)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}
	
    SugenoVal <- -1;  # this is just a initial value.
    SugenoValue <- .C("SugenoCall", as.numeric(x),
        as.numeric(v),
        as.integer(length(x)),
        out = as.numeric(SugenoVal), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
    return (SugenoValue$out);
}


fm.SugenoMob <- function(x, Mob,env=NULL)
{
    # Calculates the value of a Sugeno integral of an input x, with provided fuzzy measure v (in Mobius representation). 
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)||env$n!=length(x)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    v = fm.Zeta(Mob,env);
	
    SugenoVal <- -1;  # this is just a initial value.
    SugenoValue <- .C("SugenoCall", as.numeric(x),
        as.numeric(v),
        as.integer(length(x)),
        out = as.numeric(SugenoVal), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
    return (SugenoValue$out);
}


fm.test <- function ()
{
	# Checking that the toolbox has been installed succeffully and the functions find correct values. 
	
	print("initialisation for n=3 env<-fm.Init(3)")
	env<-fm.Init(3); print(env);

	print("Banzhaf indices fm.Banzhaf(c(0, 0.3, 0.5, 0.6, 0.4, 0.8, 0.7, 1),env)")
	print(fm.Banzhaf(c(0, 0.3, 0.5, 0.6, 0.4, 0.8, 0.7, 1),env))	

	print("Banzhaf indices for Mobius fuzzy measure fm.BanzhafMob(c(0.0, 0.3, 0.5, -0.2, 0.4, 0.1, -0.2, 0.1),env)")
	print(fm.BanzhafMob(c(0.0, 0.3, 0.5, -0.2, 0.4, 0.1, -0.2, 0.1),env))	

	print("Choquet integral for a general fuzzy measure")
	print(fm.Choquet(c(0.6, 0.3, 0.8), c(0, 0.3, 0.5, 0.6, 0.4, 0.8, 0.7, 1),env))
	
    print("Choquet integral for a Mobius fuzzy measure")
	print(fm.ChoquetMob(c(0.6, 0.3, 0.8), c(0.0, 0.3, 0.5, -0.2, 0.4, 0.1, -0.2, 0.1),env))
	
    print("ConstructLambdaMeasure in standard representation")
	print(mea<-fm.ConstructLambdaMeasure(c(0, 0.3, 0.5),env))

    print("ConstructLambdaMeasure in Mobius representation")
	print(meamob<-fm.ConstructLambdaMeasureMob(c(0, 0.3, 0.5),env))

    print("Dual measure for a fuzzy measure in standard representation")
	print(fm.dualm(mea$measure,env))

    print("Dual measure for a fuzzy measure in Mobius representation")
	print(fm.dualmMob(meamob$measure,env))

	print("Entropy value of the Choquet integral for a general fuzzy measure")
	print(fm.EntropyChoquet(c(0, 0.3, 0.5, 0.6, 0.4, 0.8, 0.7, 1),env))

    print("Entropy value of the Choquet integral for a Mobius fuzzy measure")
	print(fm.EntropyChoquetMob(c(0.0, 0.3, 0.5, -0.2, 0.4, 0.1, -0.2, 0.1),env))
	
    print("Fitting a standard fuzzy measure to data")
    d <-  matrix( c( 0.00125122, 0.563568, 0.193298, 0.164338, 
            0.808716, 0.584991, 0.479858, 0.544309, 
            0.350281, 0.895935, 0.822815, 0.625868, 
            0.746582, 0.174103, 0.858917, 0.480347, 
            0.71048, 0.513519, 0.303986, 0.387631, 
            0.0149841, 0.0914001, 0.364441, 0.134229, 
            0.147308, 0.165894, 0.988495, 0.388044, 
            0.445679, 0.11908, 0.00466919, 0.0897714, 
            0.00891113, 0.377869, 0.531647, 0.258585, 
            0.571167, 0.601746, 0.607147, 0.589803, 
            0.166229, 0.663025, 0.450775, 0.357412, 
            0.352112, 0.0570374, 0.607666, 0.270228, 
            0.783295, 0.802582, 0.519867, 0.583348, 
            0.301941, 0.875946, 0.726654, 0.562174, 
            0.955872, 0.92569, 0.539337, 0.633631, 
            0.142334, 0.462067, 0.235321, 0.228419, 
            0.862213, 0.209595, 0.779633, 0.498077, 
            0.843628, 0.996765, 0.999664, 0.930197, 
            0.611481, 0.92426, 0.266205, 0.334666, 
            0.297272, 0.840118, 0.0237427, 0.168081), 
       nrow=20, 
       ncol=4);
    print(mea1<-fm.fitting(d,env))
#    mea1<-fm.fitting(d,env)

    print("Fitting a Mobius fuzzy measure to data")
    print(fm.fittingMob(d,env))
    fm.fitting(d,env)

    print("Transform to the Mobius representation")
    print(mea1mob<-fm.Mobius(mea1,env))

	

 
	print("Interaction index for standard fuzzy measure")
	print(fm.Interaction(mea1,env))

	print("Interaction index for Mobius fuzzy measure")
	print(fm.InteractionMob(mea1mob,env))

	print("InteractionB index for a standard fuzzy measure")
	print(fm.InteractionB(mea1,env))

	print("InteractionB index for a Mobius fuzzy measure")
	print(fm.InteractionBMob(mea1mob,env))

	print("Is a standard measure additive?")
	print(fm.IsMeasureAdditive(mea1,env))

	print("Is a Mobius measure additive?")
	print(fm.IsMeasureAdditiveMob(mea1mob,env))

	print("Is a standard measure balanced?")
	print(fm.IsMeasureBalanced(mea1,env))

	print("Is a Mobius measure balanced?")
	print(fm.IsMeasureBalancedMob(mea1mob,env)) 

	print("Is a standard measure selfdual?")
	print(fm.IsMeasureSelfdual(mea1,env))

	print("Is a Mobius measure selfdual?")
	print(fm.IsMeasureSelfdualMob(mea1mob,env))

	print("Is a standard measure subadditive?")
	print(fm.IsMeasureSubadditive(mea1,env))

	print("Is a Mobius measure subadditive?")
	print(fm.IsMeasureSubadditiveMob(mea1mob,env))

	print("Is a standard measure submodular?")
	print(fm.IsMeasureSubmodular(mea1,env))

	print("Is a Mobius measure submodular?")
	print(fm.IsMeasureSubmodularMob(mea1mob,env))

	print("Is a standard measure superadditive?")
	print(fm.IsMeasureSuperadditive(mea1,env))

	print("Is a Mobius measure superadditive?")
	print(fm.IsMeasureSuperadditiveMob(mea1mob,env))

	print("Is a standard measure supermodular?")
	print(fm.IsMeasureSupermodular(mea1,env))

	print("Is a Mobius measure supermodular?")
	print(fm.IsMeasureSupermodularMob(mea1mob,env))

	print("Is a standard measure symmetric?")
	print(fm.IsMeasureSymmetric(mea1,env))

	print("Is a Mobius measure symmetric?")
	print(fm.IsMeasureSymmetricMob(mea1mob,env))

	
    print("Orness value of the Choquet integral for a standard fuzzy measure")
	print(fm.OrnessChoquet(mea1,env))

    print("Orness value of the Choquet integral for a Mobius fuzzy measure")
	print(fm.OrnessChoquetMob(mea1mob,env))

	print("Shapley value for standard representation")
	print(fm.Shapley(mea1,env))

	print("Shapley value for Mobius representation")
	print(fm.ShapleyMob(mea1mob,env))

	print("Sugeno integral from general fuzzy measure")
	print(fm.Sugeno(c(0.6, 0.3, 0.8), c(0, 0.3, 0.5, 0.6, 0.4, 0.8, 0.7, 1),env))

	print("Sugeno integral from Mobius fuzzy measure")
	print(fm.SugenoMob(c(0.6, 0.3, 0.8), c(0.0, 0.3, 0.5, -0.2, 0.4, 0.1, -0.2, 0.1),env)) 

	print("Zeta transform")
	print(fm.Zeta(mea1mob,env))
}


fm.Zeta<- function(Mob,env)
{
	if(fm.errorcheck(env)) {
		print("Incorrect environment specified, call env<-fm.Init(n) first.");
		return (NULL);
	}
	if(env$m!=length(Mob)) {
		print("The environment mismatches the dimension to the fuzzy measure.");
		return (NULL);
	}

    # Calculates the general fuzzy measure from its Mobius representation.
    ZetaVal <- array(0,length(Mob));
    ZetaValue <- .C("ZetaCall", as.numeric(Mob), 
        out = as.numeric(ZetaVal),
        as.integer(env$n), 
	 as.integer(env$m), as.integer(env$card), as.integer(env$cardpos),as.integer(env$bit2card),as.integer(env$card2bit),as.double(env$factorials)
);
					
    return (ZetaValue$out);
}


fm.fittingWAM<- function(data, env=NULL)
{
	# This function estimates the values of a k-additive standard fuzzy measure based on empirical data. 
	# The result is an array containing the values of the fuzzy measure in Mobius, ordered according to set cardinalities.
	# kadd define the complexity of fuzzy measure. if kadd is not provided, its default value is equal to the number of inputs.

	size <- dim(as.matrix(data));
	n <- size[2] - 1;
	datanum <- size[1];
	WeightVal <- array(0,n);
 
	WeightValue <- .C("fittingWAMCall", as.integer(n),
	  				    as.integer(datanum),
			              out = as.numeric(WeightVal),
			                    as.numeric(t(data)));
					
	return (WeightValue$out);

}
fm.fittingOWA<- function(data, env=NULL)
{
	# This function estimates the values of a k-additive standard fuzzy measure based on empirical data. 
	# The result is an array containing the values of the fuzzy measure in Mobius, ordered according to set cardinalities.
	# kadd define the complexity of fuzzy measure. if kadd is not provided, its default value is equal to the number of inputs.

	size <- dim(as.matrix(data));
	n <- size[2] - 1;
	datanum <- size[1];
	WeightVal <- array(0,n);

	WeightValue <- .C("fittingOWACall", as.integer(n),
	  				    as.integer(datanum),
			              out = as.numeric(WeightVal),
			                    as.numeric(t(data)));
 					
	return (WeightValue$out);

}
