# Code for the monotonic detection function fitting!
# 2 functions here:
#  first for creating the reference points
#  second is the constraint function to pass to solnp
# Credit for this goes to Lorenzo Milazzo, any bugs are
# Dave Miller's fault.


# function to evaluate a set of `reference points' (distances) 
# used to define the inequality constraints
getRefPoints<- function(no_d, int.range){

  # previous versions just used width and assumed that left truncation
  # was at zero, now using int.range to tell us the interval

  # note this currently doesn't work with multiple integration ranges
  # i.e. when int.range is a matrix w. > 1 row. This should have been caught
  # and an error throw before now though.

  xlat <- (int.range[2]-int.range[1])/(no_d^1.5)
  ref_points <- double(no_d)
  for(i in 1:no_d){
    ref_points[i] <- (i^1.5) * xlat
  }

  ref_points <- ref_points+int.range[1]

  return(ref_points)
}


#
# set of equations associated with the Inequality Constraints
#
# input:
#  pars           - parameters
#  ddfobj         - ddf object with almost everything we need in it
#  misc.options   - everything else...
flnl.constr<- function(pars, ddfobj, misc.options,...){

  if(is.null(ddfobj$adjustment)){
    # this never gets called from ddf()
    ineq_constr <- rep(10,2*misc.options$mono.points)
  }else{
    ddfobj <- assign.par(ddfobj,pars)

    # apply the constraints?
    constr <- misc.options$mono
    # apply strict monotonicity?
    strict <- misc.options$mono.strict

    ### Constraint stuff here:
    # number of points (distances)
    # at which the DF is evaluated
    no_d <- misc.options$mono.points
    # reference points (distances)
    ref_p <- getRefPoints(no_d, misc.options$int.range)
    # to get detfct to play nice need to mudge ddfobj a bit...
    if(!is.null(ddfobj$scale)){
      ddfobj$scale$dm <- rep(1,no_d)
    }
    if(!is.null(ddfobj$shape)){
      ddfobj$shape$dm <- rep(1,no_d)
    }

    # evaluate the detection function at the reference points
    # note that we must standardize so 0<=g(x)<=1
    df_v_rp <- as.vector(detfct(ref_p,ddfobj,width=misc.options$width,
                               standardize=TRUE))

    # reference point associated with distance=0
    ref_p0 <- 0
    # again, to get detfct to play nice need to mudge ddfobj a bit...
    if(!is.null(ddfobj$scale)){
      ddfobj$scale$dm <- 1
    }
    if(!is.null(ddfobj$shape)){
      ddfobj$shape$dm <- 1
    }

    # evaluate the detection function at 0
    df_v_rp0 <- as.vector(detfct(ref_p0,ddfobj,width=misc.options$width,
                                standardize=TRUE))

    # inequality constraints ensuring the
    # (weak or strict) monotonicity
    ic_m <- NULL
    if(constr){
      # set the reference point to be the detection function
      # value at 0
      df_v_rp_p <- df_v_rp0
      ic_m <- double(no_d)
      for(i in 1:no_d){
        ic_m[i] <- (df_v_rp_p - df_v_rp[i])
        if(strict){
          # if we have strict monotonicity, then change the ref
          # point to be the last point
          df_v_rp_p <- df_v_rp[i]
        }
      }
    }
    # inequality constraints ensuring that
    # the detection function is always >=0
    ic_p <- double(no_d)
    ic_p <- df_v_rp

    #  set of inequality constraints
    ineq_constr <- c(ic_m, ic_p)
  }
  return(ineq_constr)
}

