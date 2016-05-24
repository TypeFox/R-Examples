sample.from.polytope <-
function# CORE function to draw the samples from the polytope
### internal function for sisus
(Ab
### internal variable
, M
### internal variable
, skip
### internal variable
, burnin
### internal variable
, warning.sw
### internal variable
, i.samples.isotope.mvn = 1
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  # attached via DESCRIPTION # library("polyapost");   # for drawing samples uniformly from solution polytope

  sam = NULL;

  # determine if running in windows or unix environment
  OS = .Platform$OS.type;
#  if (OS == "unix") {
    # attached via DESCRIPTION # library("rcdd");      # for determining vertices of solution polytope
    # Using RCDD as formal check of whether a solution exists
    # get verticies

      # p.o = paste("Determine vertices of solution polytope", "\n"); write.out(p.o);
      p.o = paste("(v)"); write.out(p.o);
    # Polytope Vertices (cdd), add to top of sample and remove same number of vertices from bottom of sample
      # format is column 1: 1=equality, 0=inequality
      #           column 2: b vector for Ax<=b, therefore -b for Ax>=b
      #  remaining columns: -A for Ax<=b, therefore A for Ax>=b
      # from polytope.constraints():  1:  Ax=b;  2:  Ax<=b;  3:  Ax>=b;
    # define halfspace H-representation using the linear constraints
    H.rep = rbind(cbind(rep(1,length(Ab$b1)),  Ab$b1, -Ab$A1),
                  cbind(rep(0,length(Ab$b2)),  Ab$b2, -Ab$A2),
                  cbind(rep(0,length(Ab$b3)), -Ab$b3,  Ab$A3));
    H.rep <- d2q(H.rep) # convert decimal to rational 3/5/2014
    valid.cdd.check = validcdd(H.rep, representation="H");  # check that representation is valid
    if (valid.cdd.check == 0) {
      p.o = paste("\n"); write.out(p.o);
      p.o = paste("           WARNING: CDD not valid H-representation of solution polytope.", "\n"); write.out(p.o);
      warning(p.o);
    }; # end if valid.cdd.check

    V.rep = scdd(H.rep, representation = "H");    # vertex V-representation
    V.sam = V.rep$output[,3:(2+dim(Ab$A2)[2])];   # extract the vertices
    V.sam = matrix(V.sam,ncol=dim(Ab$A3)[2]);     # make a matrix if it is a vector
    n.vertices = dim(V.sam)[1];                   # number of vertices
    if (n.vertices == 0) { # no solution
      V.sam = NULL;
      warning.sw = 1;
      sol.feasible = 0;   # indicate an unfeasible solution
      p.o = paste("\n"); write.out(p.o);
      p.o = paste("           WARNING: No vertices in solution polytope, Unfeasible Solution -- check Convex Hull plots.", "\n"); write.out(p.o);

      ############### return and quit if no vertices to solution polytope
      SAMPLE = new.env();  # create a list to return sample (sam) with whether it was feasible (sol.feasible)
      SAMPLE$sam          = sam         ;
      SAMPLE$sol.feasible = sol.feasible;
      SAMPLE$warning.sw   = warning.sw  ;
      SAMPLE$n.vertices   = n.vertices  ;
      ##SAMPLE$V.sam        = V.sam       ;
      return( as.list(SAMPLE) );
    }; # end if n.vertices == 0

    if (n.vertices == 1) { # unique solution
      M=1;
      sol.feasible = 1;   # indicate a feasible solution
      sam = V.sam;
    }; # end if n.vertices == 1

    #if (i.samples.isotope.mvn == 1) { # only get vertices for first sample
    #}; # end if i.samples.isotope.mvn
#  }; # end if OS unix

#  if (OS == "windows") {
#    n.vertices = 0;
#    V.sam = NULL  ;
#
#    #### RECODE THIS SECTION WITH THIS LOGIC 2/28/2008 8:32AM
#    # Try feasible with eps>0, if solution then M != 0 and sol.feasible = 1
#    #   else, try feasible with eps=0, if solution then M = 0 and sol.feasible = 1
#    #     else, sol.feasible = 0.
#
#    # unique solution or no solution under windows
#    if (dim(as.matrix(Ab$A1))[1] == dim(as.matrix(Ab$A1))[2]) {  # rows equal is unique solution
#      # unique solution
#      ### n=15;k=7;eps=0; initsol = feasible(Ab$A1[1:(k+1),1:n],Ab$A2[1:n,1:n],Ab$A3[1:n,1:n],Ab$b1[1:(k+1)],Ab$b2[1:n],Ab$b3[1:n],eps);initsol
#      M=1; eps=0; initsol = feasible(Ab$A1,Ab$A2,Ab$A3,Ab$b1,Ab$b2,Ab$b3,eps);
#      if (initsol[1] < 0 ) {
#        sol.feasible = 0;   # indicate an unfeasible solution
#        sam = NULL;
#        if (warning.sw == 0) {
#          warning.sw = 1;
#          p.o = paste("\n"); write.out(p.o);
#          p.o = paste("           WARNING: Encountered an Unfeasible Solution -- check Convex Hull plots.", "\n"); write.out(p.o);
#          warning(p.o);
#        }
#      } else {
#        sol.feasible = 1;   # indicate a feasible solution
#        sam = initsol;
#      };
#    };  # end M=1 unique solution
#  }; # end if OS windows

  ########################################
  ## Perform uniform sampling of solution polytope
  if (M != 1) { # nonunique solution
    eps = 0; #1e-2; # for finding an initial solution on boundary
    ## Ab <- d2q(Ab) # convert decimal to rational 3/5/2014
    ## Ab$A1 <- d2q(Ab$A1) # convert decimal to rational 3/5/2014
    ## Ab$A2 <- d2q(Ab$A2) # convert decimal to rational 3/5/2014
    ## Ab$A3 <- d2q(Ab$A3) # convert decimal to rational 3/5/2014
    ## Ab$b1 <- d2q(Ab$b1) # convert decimal to rational 3/5/2014
    ## Ab$b2 <- d2q(Ab$b2) # convert decimal to rational 3/5/2014
    ## Ab$b3 <- d2q(Ab$b3) # convert decimal to rational 3/5/2014
    ## eps   <- d2q(eps  ) # convert decimal to rational 3/5/2014
    initsol = feasible(Ab$A1,Ab$A2,Ab$A3,Ab$b1,Ab$b2,Ab$b3,eps);        # a starting solution within the solution polytope
    #if (initsol[1] < 0 ) {
    if ((initsol[1] < 0 ) & !all.equal(initsol[1], 0)) {  # negative, but not basically 0
      sol.feasible = 0;   # indicate an unfeasible solution
      sam = NULL;
      if (warning.sw == 0) {
        warning.sw = 1;
        p.o = paste("\n"); write.out(p.o);
        p.o = paste("           WARNING: Encountered an Unfeasible Solution -- check Convex Hull plots.", "\n"); write.out(p.o);
        warning(p.o);
      }
    } else {
      sol.feasible = 1;   # indicate a feasible solution
      # burnin
      burnin.count = 0; burnin.max.tries = 1000;
      burninsam = matrix(c(initsol,initsol),nrow=2,byrow=TRUE);
      while ( identical(burninsam[dim(burninsam)[1],], initsol) ) { # loop until off boundary
        burninsam = constrppprob(Ab$A1,Ab$A2,Ab$A3,Ab$b1,Ab$b2,Ab$b3,burninsam[dim(burninsam)[1],],skip,burnin); # draw BURNIN uniform samples within the solution polytope
        burnin.count = burnin.count+1;
        if (burnin.count > burnin.max.tries) {
          p.o = paste("WARNING: Burnin failed to get into solution polytope interior after", burnin*burnin.max.tries, "attempts.", "\n"); write.out(p.o);
          p.o = paste("         Extend burnin period, then contact SISUS author to attempt resolution.", "\n"); write.out(p.o);
          break;
        }
      }
      initsol = burninsam[dim(burninsam)[1],]; # assign initsol a point inside the polytope
      # sample
      sam = constrppprob(Ab$A1,Ab$A2,Ab$A3,Ab$b1,Ab$b2,Ab$b3,initsol,skip,M); # draw uniform samples within the solution polytope

    }; # if initsol

  }; # end M!=1
              # convert rational to decimal 3/5/2014
  sam = rbind(q2d(V.sam), q2d(sam));                      # append vertices to beginning of sample

  SAMPLE = new.env();  # create a list to return sample (sam) with whether it was feasible (sol.feasible)
  SAMPLE$sam          = sam         ;
  SAMPLE$sol.feasible = sol.feasible;
  SAMPLE$warning.sw   = warning.sw  ;
  SAMPLE$n.vertices   = n.vertices  ;
  ##SAMPLE$V.sam        = V.sam       ;

  return( as.list(SAMPLE) );

  ### internal variable
}
