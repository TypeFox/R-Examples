check<-
  function(test=NULL,x=NULL,B=NULL,alpha=NULL,prop=NULL,S=NULL,bit=NULL,k=NULL,m=NULL,n=NULL,lambda=NULL,
           KS=TRUE,CSQ=TRUE,AD=TRUE,JB=TRUE,test.k=TRUE,test.g=TRUE,mu=NULL,sd=NULL,
           Excursion=TRUE,Expansion=TRUE,Height=TRUE,critical.value=NULL,num.class=NULL){
    # test shows inputs of the test will be checked
    # 1:Achi; 2:BDS; 3:BS; 4:GCD; 5:RW; 6:TBT
    # k will be num.class for BDS test and k BS test
    # B is the bit-length
    # alpha is the significance level
    # bit is logical shows the type of input
    
    if (is.numeric(x)==FALSE){
      stop("Number sequence must contain numeric entries!")
    }        
    
    if ((alpha<=0) | (alpha>=1)){
      stop("The significance level alpha must be entered between 0 and 1!")
    }
    if (test==1){
      if ((prop<=0) | (prop>=1)){
        stop("The proportion of test and training data sets in adaptive chi-square test must be entered between 0 and 1!")
      }
      if ((abs(S - round(S)) < (.Machine$double.eps^0.5))==FALSE){
        stop("An integer value should be entered for the number of sub-alphabets in adaptive chi-square test!")
      }
    }    
    if (test!=2){
      if (is.null(bit)==FALSE){
        if (is.logical(bit)==FALSE){
          stop("A logical value must be entered for bit!")
        }
        if ((bit==TRUE) & (all(x %in% 0:1)==FALSE)){ #Are all elements of x either 0 or 1?
          stop("Input sequence must contain binary entries when bit is set to TRUE!")
        }
        if ((bit==FALSE) & (is.numeric(x)==FALSE)){ #Are all elements of x numeric?
          stop("Input sequence must contain numeric entries when bit is set to FALSE!")
        }
      }            
      if (is.numeric(B)==FALSE){
        stop("Bit length must be a numeric entry!")
      }
      if ((abs(B - round(B)) < (.Machine$double.eps^0.5))==FALSE){ #is B a wholenumber
        stop("Bit length must be an integer!")
      }
      
    }   
    if (test==2){
      if ((abs(num.class - round(num.class)) < (.Machine$double.eps^0.5))==FALSE){
        stop("An integer value should be entered for the number of classes in birthday spacings test!")
      }
      if ((abs(m - round(m)) < (.Machine$double.eps^0.5))==FALSE){
        stop("An integer value should be entered for the number of birthdays in birthday spacings test!")
      }
      if (lambda<0){
        stop("A positive value should be entered for mean rate of theoretical Poisson distribution in birthday spacings test!")
      }
    }    
    if (test==3){ 
      if ((abs(k - round(k)) < (.Machine$double.eps^0.5))==FALSE){
        stop("An integer value should be entered for k in book stack test!")
      }
    }
    if (test==4){
      if ((is.logical(KS)==FALSE)|(is.logical(CSQ)==FALSE)|(is.logical(AD)==FALSE)|(is.logical(JB)==FALSE)|
            (is.logical(test.k)==FALSE)|(is.logical(test.g)==FALSE)){
        stop("At least one of the logical arguments of GCD test is entered unlogical!")
      }
      if (sd<0){
        stop("A positive value should be entered for standard deviation of theoretical normal distribution in GCD test!")
      }
    }    
    if (test==5){
      if ((is.logical(Excursion)==FALSE)|(is.logical(Expansion)==FALSE)|(is.logical(Height)==FALSE)){
        stop("At least one of the logical arguments of random walk tests is entered unlogical!")
      }
    }
    if (test==6){
      if (critical.value<=0){
        stop("Critical value of topological binary test must be positive!")
      }
    }
    
  }