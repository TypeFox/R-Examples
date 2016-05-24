# Class for spline basis objects
# see SplinParam.R for methods

setClass(".SplineBasis",
         representation(knots="numeric",
                        degree="integer",
                        nbases="integer",
                        log="logical",
                        "VIRTUAL"))
setClass("BSplineBasis",
         representation(Matrices="array"),
         contains=".SplineBasis")

# slot SplineBasis of class SplineBasis of package "orthogonalsplinebasis" 
setClass("MSplineBasis",
         representation(min="numeric",
                        max="numeric",
                        SplineBasis="SplineBasis"),
         contains="BSplineBasis")

setClass("TPSplineBasis",
         # knots are interior knots
         representation(min="numeric",
                        max="numeric",
                        coef="numeric",
                        degrees="integer",
                        type="character"),
         contains=".SplineBasis")


# idem but first bases are (x-ref)^i
setClass("TPRSplineBasis",
         representation(ref="numeric"),
         contains="TPSplineBasis",
         prototype=prototype(ref=0))

setClass("C0BSplineBasis", representation("BSplineBasis",
                                          ref="numeric"))
setClass("C0TPSplineBasis", representation("TPSplineBasis",
                                          ref="numeric"))

setClassUnion("AnySplineBasis", c("BSplineBasis", "MSplineBasis", "TPSplineBasis"))

##########################################################################################################
## Class DesignMatrix*
# see methods in DesignMatrix.R


setClass("DesignMatrix",
         representation(DM="array",
                        nObs="integer"
                        ),
         prototype=prototype(
           DM=NULL,
           nObs=0L)
)

setClass("DesignMatrixNPH",
         representation("DesignMatrix",
                        nX="integer",
                        TSplineBasis="AnySplineBasis",
                        nTbasis="integer",
                        intercept="logical",
                        names="character"
                        ),
         prototype=prototype(
           nX=0L,
           TSplineBasis=NULL,
           nTbasis=0L,
           intercept=TRUE,
           names=NULL
           )
         )

  # splinebasis : list of splines parameter of Z (one element per Zi)
  # signature : design matrix for spline parameters parameters :
  #                (DM %*% diag(param) %*%signature is the matrix of alpha(Z_i))
  #               alpha[index[i,1]:index[i,2]] is the vector parameter for Z_i
  # names are the names of (Z_i) 
setClass("DesignMatrixNPHNLL",
         representation("DesignMatrix",
                        nZ="integer",
                        nparam="integer",
                        signature="array",
                        index="array",
                        TSplineBasis="AnySplineBasis",
                        listSplineBasis="list",
                        nTbasis="integer",
                        names="character"
                        ),
         prototype=prototype(
           nZ=0L,
           nparam=0L,
           signature=NULL,
           index=NULL,
           TSplineBasis=NULL,
           listSplineBasis=NULL,
           nTbasis=0L,
           names=NULL)
         )

		 
		 setClass("NCStepParam",
         representation(step="numeric",
                        min="numeric",
                        max="numeric")    # the uniq step
         )

		 
		 
##########################################################################################################
## Class *StepParam*
# see methods in StepParam.R
		 
setClass("GLMStepParam",
         representation(
                        nbands="integer", # number og bands
                        ncuts="integer", # number og cuts = nbands+1
                        cuts="numeric",  # the cuts of the bands =c(min, ..., max)
                        steps="numeric",   # the steps c(b2-min, b3-b2, ..., max-b_(n-1)
                        points="numeric",  # the evaluation points c((b2+min)/2, (b3+b2)/2, ..., (max-b+(n-1))/
                        min="numeric",
                        max="numeric")
         )

setClass("NCAdaptedStepParam",
         representation(Nstep="integer",        # the number of steps
               theSteps="numeric",
               from="numeric",
               to="numeric"),    # the vector of steps
         contains="NCStepParam" )

#setClassUnion("StepParam", 
#         c("NCStepParam", "NCAdaptedStepParam", "GLMStepParam")
#         )

