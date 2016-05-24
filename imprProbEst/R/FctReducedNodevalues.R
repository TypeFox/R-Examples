`FctReducedNodevalues` <-
function(n,m,s,x.nodevalues,omega.nodevalues){

  # Load the R package 'inline':

  library(inline)

  # The variables which are passed to the C++ programm:

  Variables <- signature(n="integer", m="integer", s="integer",
    vecofnodevalues="numeric",
    xrem="integer", omegarem="integer")

  # The C - source code:

  Code <- " int i,j,l; int compare;
    int zrem[*n+*m];  /* This is a vector that will contain xrem and omegarem */
    zrem[0]=1;
    int ss=*s;  /* getting rid of those pointers */
    /* The following for-loop checks if some of the observations in x cannot be distinguished from any other
       observation. Those observations which cannot be distinguished are superfluous and may be deleted.
    */
    for(i=1; i<*n; i++){
        zrem[i]=1;  /* At the beginning, every component of zrem should be equal to 1 */
        for(j=1; j<=i; j++){
            if(zrem[i-j]!=0){ /* if omega_(i-j) is relevant (i.e. not assigned for being deleted) ... */
                for(l=0; l<ss; l++){
                    if( fabs(vecofnodevalues[i*ss+l]-vecofnodevalues[(i-j)*ss+l])>1E-7 ){
                        compare=1;
                        l=ss;
                    }
                    else {
                        compare=0;
                    }
                }
                /* 'compare==0' means: omega_i is superfluous */
                if(compare==0){
                    zrem[i]=0;  /* So, here it is assigned for being deleted.  */
                    zrem[i-j]=zrem[i-j]+1;  /*  Observation x_{i-j} is observed once more now.  */
                    j=i+1;   /* omega_i will already be deleted. Hence, we may stop the j-loop.  */
                }
            }
        }
    }
    /*  The following for-loop checks if some of the nodes are superfluous. A node omega_1 is superflouus if all
        of the evaluations of the functions f (l-loop) in omega_1 are not smaller than the evaluations of the
        functions f in some omega_2.
    */
    for(i=*n; i<*n+*m; i++){
        zrem[i]=1;  /* At the beginning, every component of zrem should be equal to 1 */
        /* The following j-loop checks for all j<=i
           if (compared with each other) omega_(i-j) or omega_i is superfluous.
        */
        for(j=1; j<=i; j++){
            if(zrem[i-j]!=0){ /* if omega_(i-j) is relevant (i.e. not assigned for being deleted) ... */
                for(l=0; l<ss; l++){
                    if(vecofnodevalues[i*ss+l] < vecofnodevalues[(i-j)*ss+l]){
                        compare=1;
                        l=ss;
                    }
                    else {
                        compare=0;
                    }
                }
                if(compare==0){  /* if omega_{i} should be assigned for being deleted ... */
                    zrem[i]=0;  /* ... it is assigned for being deleted ... */
                    j=i+1;       /* and we may stop the j-loop */
                }
                /* If omega_{i} is superfluous, it is already assigned for being deleted.
                   However, we can go further on: Perhaps omega_{i} makes some omega_{i-j} superfluous so that
                   we can assigne some omega_{i-j} for being deleted. */
                else {
                    if(i-j>=*n){ /* CAUTION: omega_{i} cannot make an observation x superfluous; thats why this
                                             if-codition is neccessary  */
                        for(l=0; l<ss; l++){
                            if(vecofnodevalues[i*ss+l] > vecofnodevalues[(i-j)*ss+l]){
                                compare=1;
                                l=ss;
                            }
                            else {
                                compare=0;
                            }
                        }
                        if(compare==0){
                           zrem[i-j]=0;
                        }
                    }
                }
            }
        }
    }
    /* Finally, zrem is plit into xrem and omegarem again:  */
    for(i=0; i<*n; i++){xrem[i]=zrem[i];}
    for(i=*n; i<*n+*m; i++){omegarem[i-*n]=zrem[i];}
  "


  # Use the R package 'inline':

  FctDelNotes <- cfunction(Variables, Code, language="C", convention=".C")

  # Preparing some variables:

  nodevalues <- cbind(x.nodevalues,omega.nodevalues)
  vecofnodevalues <- matrix(nodevalues, ncol=1)
  xrem <- matrix(1, ncol=1, nrow=n)
  omegarem <- matrix(1, ncol=1, nrow=m)

  # Use of the function 'FctDelNotes':

  DelNotes <- FctDelNotes(n,m,s,vecofnodevalues,xrem,omegarem)
  # DelNotes <- list(n,m,s,vecofnodevalues,xrem,omegarem)    Hiermit laest sich die C-Funktion ausschalten

  x.rem <- matrix(DelNotes[[5]],ncol=1)   # this vector stores, which observations have occured how often
  omega.rem <- matrix(DelNotes[[6]],ncol=1)  # this is the vector which stores the remaining
                                             # additional supporting nodes
  if (sum(omega.rem)==0) {omega.rem[1]=1} # This line prevents 'omega.nodevalues' from becoming an empty
                                          # vector which may lead to a programming error further on.
  mr <- sum(omega.rem)                       # number of remaining additional supporting nodes

  x.nodevalues <- x.nodevalues[,x.rem!=0]   # the values of the functions f at the remaining observed nodes of x
  x.frequency <- x.rem[x.rem!=0]             # the (absolute) frequency of the observations at the remaining
                                            # observed nodes of x

  omega.nodevalues <- omega.nodevalues[,omega.rem!=0]    # remaining additional supporting nodes omega

  # The return of the function 'FctReducedNodevalues'

  list(x.nodevalues,x.frequency,omega.nodevalues,mr)
}

