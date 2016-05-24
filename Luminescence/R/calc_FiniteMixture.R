#' Apply the finite mixture model (FMM) after Galbraith (2005) to a given De
#' distribution
#'
#' This function fits a k-component mixture to a De distribution with differing
#' known standard errors. Parameters (doses and mixing proportions) are
#' estimated by maximum likelihood assuming that the log dose estimates are
#' from a mixture of normal distributions.
#'
#' This model uses the maximum likelihood and Bayesian Information Criterion
#' (BIC) approaches. \cr\cr Indications of overfitting are: \cr\cr - increasing
#' BIC \cr - repeated dose estimates \cr - covariance matrix not positive
#' definite \cr - covariance matrix produces NaNs\cr - convergence problems
#' \cr\cr \bold{Plot} \cr\cr If a vector (\code{c(k.min:k.max)}) is provided
#' for \code{n.components} a plot is generated showing the the k components
#' equivalent doses as normal distributions. By default \code{pdf.weight} is
#' set to \code{FALSE}, so that the area under each normal distribution is
#' always 1. If \code{TRUE}, the probability density functions are weighted by
#' the components proportion for each iteration of k components, so the sum of
#' areas of each component equals 1. While the density values are on the same
#' scale when no weights are used, the y-axis are individually scaled if the
#' probability density are weighted by the components proportion. \cr The
#' standard deviation (sigma) of the normal distributions is by default
#' determined by a common \code{sigmab} (see \code{pdf.sigma}). For
#' \code{pdf.sigma = "se"} the standard error of each component is taken
#' instead.\cr The stacked barplot shows the proportion of each component (in
#' per cent) calculated by the FFM. The last plot shows the achieved BIC scores
#' and maximum log-likelihood estimates for each iteration of k.
#'
#' @param data \code{\linkS4class{RLum.Results}} or \link{data.frame}
#' (\bold{required}): for \code{data.frame}: two columns with De
#' \code{(data[,1])} and De error \code{(values[,2])}
#' @param sigmab \code{\link{numeric}} (\bold{required}): spread in De values
#' given as a fraction (e.g. 0.2). This value represents the expected
#' overdispersion in the data should the sample be well-bleached (Cunningham &
#' Wallinga 2012, p. 100).
#' @param n.components \code{\link{numeric}} (\bold{required}): number of
#' components to be fitted. If a vector is provided (e.g. \code{c(2:8)}) the
#' finite mixtures for 2, 3 ... 8 components are calculated and a plot and a
#' statistical evaluation of the model performance (BIC score and maximum
#' log-likelihood) is provided.
#' @param grain.probability \code{\link{logical}} (with default): prints the
#' estimated probabilities of which component each grain is in
#' @param dose.scale \code{\link{numeric}}: manually set the scaling of the
#' y-axis of the first plot with a vector in the form of \code{c(min,max)}
#' @param pdf.weight \code{\link{logical}} (with default): weight the
#' probability density functions by the components proportion (applies only
#' when a vector is provided for \code{n.components})
#' @param pdf.sigma \code{\link{character}} (with default): if \code{"sigmab"}
#' the components normal distributions are plotted with a common standard
#' deviation (i.e. \code{sigmab}) as assumed by the FFM. Alternatively,
#' \code{"se"} takes the standard error of each component for the sigma
#' parameter of the normal distribution
#' @param pdf.colors \code{\link{character}} (with default): color coding of
#' the components in the the plot. Possible options are "gray", "colors" and
#' "none"
#' @param pdf.scale \code{\link{numeric}}: manually set the max density value
#' for proper scaling of the x-axis of the first plot
#' @param plot.proportions \code{\link{logical}} (with default): plot barplot
#' showing the proportions of components
#' @param plot \code{\link{logical}} (with default): plot output
#' @param \dots further arguments to pass.  See details for their usage.
#' @return Returns a plot (optional) and terminal output. In addition an
#' \code{\linkS4class{RLum.Results}} object is returned containing the
#' following elements:
#'
#' \item{summary}{\link{data.frame} summary of all relevant model results.}
#' \item{data}{\link{data.frame} original input data} \item{args}{\link{list}
#' used arguments} \item{call}{\link{call} the function call} \item{mle}{
#' covariance matrices of the log likelhoods} \item{BIC}{ BIC score}
#' \item{llik}{ maximum log likelihood} \item{grain.probability}{ probabilities
#' of a grain belonging to a component} \item{components}{\link{matrix}
#' estimates of the de, de error and proportion for each component}
#' \item{single.comp}{\link{data.frame} single componente FFM estimate}
#'
#' If a vector for \code{n.components} is provided (e.g.  \code{c(2:8)}),
#' \code{mle} and \code{grain.probability} are lists containing matrices of the
#' results for each iteration of the model.
#'
#' The output should be accessed using the function
#' \code{\link{get_RLum}}
#' @section Function version: 0.4
#' @author Christoph Burow, University of Cologne (Germany) \cr Based on a
#' rewritten S script of Rex Galbraith, 2006. \cr
#' @seealso \code{\link{calc_CentralDose}}, \code{\link{calc_CommonDose}},
#' \code{\link{calc_FuchsLang2001}}, \code{\link{calc_MinDose}}
#' @references Galbraith, R.F. & Green, P.F., 1990. Estimating the component
#' ages in a finite mixture. Nuclear Tracks and Radiation Measurements 17,
#' 197-206. \cr\cr Galbraith, R.F. & Laslett, G.M., 1993. Statistical models
#' for mixed fission track ages. Nuclear Tracks Radiation Measurements 4,
#' 459-470.\cr\cr Galbraith, R.F. & Roberts, R.G., 2012. Statistical aspects of
#' equivalent dose and error calculation and display in OSL dating: An overview
#' and some recommendations. Quaternary Geochronology 11, 1-27.\cr\cr Roberts,
#' R.G., Galbraith, R.F., Yoshida, H., Laslett, G.M. & Olley, J.M., 2000.
#' Distinguishing dose populations in sediment mixtures: a test of single-grain
#' optical dating procedures using mixtures of laboratory-dosed quartz.
#' Radiation Measurements 32, 459-465.\cr\cr Galbraith, R.F., 2005. Statistics
#' for Fission Track Analysis, Chapman & Hall/CRC, Boca Raton.\cr\cr
#' \bold{Further reading}\cr\cr Arnold, L.J. & Roberts, R.G., 2009. Stochastic
#' modelling of multi-grain equivalent dose (De) distributions: Implications
#' for OSL dating of sediment mixtures. Quaternary Geochronology 4,
#' 204-230.\cr\cr Cunningham, A.C. & Wallinga, J., 2012. Realizing the
#' potential of fluvial archives using robust OSL chronologies. Quaternary
#' Geochronology 12, 98-106.\cr\cr Rodnight, H., Duller, G.A.T., Wintle, A.G. &
#' Tooth, S., 2006. Assessing the reproducibility and accuracy of optical
#' dating of fluvial deposits.  Quaternary Geochronology 1, 109-120.\cr\cr
#' Rodnight, H. 2008. How many equivalent dose values are needed to obtain a
#' reproducible distribution?. Ancient TL 26, 3-10.
#' @examples
#'
#' ## load example data
#' data(ExampleData.DeValues, envir = environment())
#'
#' ## (1) apply the finite mixture model
#' ## NOTE: the data set is not suitable for the finite mixture model,
#' ## which is why a very small sigmab is necessary
#' calc_FiniteMixture(ExampleData.DeValues$CA1,
#'                    sigmab = 0.2, n.components = 2,
#'                    grain.probability = TRUE)
#'
#' ## (2) repeat the finite mixture model for 2, 3 and 4 maximum number of fitted
#' ## components and save results
#' ## NOTE: The following example is computationally intensive. Please un-comment
#' ## the following lines to make the example work.
#' FMM<- calc_FiniteMixture(ExampleData.DeValues$CA1,
#'                          sigmab = 0.2, n.components = c(2:4),
#'                          pdf.weight = TRUE, dose.scale = c(0, 100))
#'
#' ## show structure of the results
#' FMM
#'
#' ## show the results on equivalent dose, standard error and proportion of
#' ## fitted components
#' get_RLum(object = FMM, data.object = "components")
#'
#' @export
calc_FiniteMixture <- function(
  data,
  sigmab,
  n.components,
  grain.probability = FALSE,
  dose.scale,
  pdf.weight = TRUE,
  pdf.sigma = "sigmab",
  pdf.colors = "gray",
  pdf.scale,
  plot.proportions = TRUE,
  plot=TRUE,
  ...
){

  ##============================================================================##
  ## CONSISTENCY CHECK OF INPUT DATA
  ##============================================================================##

  if(missing(data)==FALSE){

    if(is(data, "data.frame") == FALSE & is(data,"RLum.Results") == FALSE){
      stop("[calc_FiniteMixture] Error: 'data' object has to be of type
           'data.frame' or 'RLum.Results'!")
    } else {
      if(is(data, "RLum.Results") == TRUE){
        data <- get_RLum(data, signature(object = "De.values"))
      }
    }
  }
  try(colnames(data)<- c("ED","ED_Error"),silent=TRUE)
  if(colnames(data[1])!="ED"||colnames(data[2])!="ED_Error") {
    cat(paste("Columns must be named 'ED' and 'ED_Error'"), fill = FALSE)
    stop(domain=NA)
  }
  if(sigmab <0 | sigmab >1) {
    cat(paste("sigmab needs to be given as a fraction between",
              "0 and 1 (e.g. 0.2)"), fill = FALSE)
    stop(domain=NA)
  }
  if(any(n.components<2) == TRUE) {
    cat(paste("Atleast two components need to be fitted"), fill = FALSE)
    stop(domain=NA)
  }
  if(pdf.sigma!="se" ) {
    if(pdf.sigma!="sigmab") {
      cat(paste("Only 'se' or 'sigmab' allowed for the pdf.sigma argument"),
          fill = FALSE)
      stop(domain=NA)
    }
  }

  ##============================================================================##
  ## ... ARGUMENTS
  ##============================================================================##

  extraArgs <- list(...)

  ## console output
  if("verbose" %in% names(extraArgs)) {
    verbose<- extraArgs$verbose
  } else {
    verbose<- TRUE
  }
  # trace calculations
  if("trace" %in% names(extraArgs)) {
    trace<- extraArgs$trace
  } else {
    trace<- FALSE
  }
  # plot title
  if("main" %in% names(extraArgs)) {
    main<- extraArgs$main
  } else {
    main<- "Finite Mixture Model"
  }

  ##============================================================================##
  ## CALCULATIONS
  ##============================================================================##

  ## create storage variables if more than one k is provided
  if(length(n.components)>1) {

    # counter needed for various purposes
    cnt<- 1

    # create summary matrix containing DE, standard error (se) and proportion
    # for each component
    comp.n<- matrix(data = NA, ncol = length(n.components),
                    nrow = n.components[length(n.components)] * 3,
                    byrow = TRUE)

    # create empty vector as storage for BIC and LLIK scores
    BIC.n<- vector(mode = "double")
    LLIK.n<- vector(mode = "double")

    # create empty vectors of type "lists" as storage for mle matrices and
    # grain probabilities
    vmat.n<- vector(mode = "list", length = length(n.components))
    grain.probability.n<- vector(mode = "list", length = length(n.components))

  }

  ## start actual calculation (loop) for each provided maximum components to
  ## be fitted.
  for(i in 1:length(n.components)) {

    k<- n.components[i]

    # calculate yu = log(ED),  su = se(logED),  n = number of grains
    yu<- log(data$ED)
    su<- data$ED_Error/data$ED
    n<- length(yu)

    # compute starting values
    fui<- matrix(0,n,k)
    pui<- matrix(0,n,k)
    nui<- matrix(0,n,k)
    pii<- rep(1/k,k)
    mu<- min(yu) + (max(yu)-min(yu))*(1:k)/(k+1)

    # remove the # in the line below to get alternative starting values
    # (useful to check that the algorithm converges to the same values)
    #	mu<- quantile(yu,(1:k)/(k+1))

    # compute maximum log likelihood estimates
    nit<- 499L
    wu<- 1/(sigmab^2 + su^2)
    rwu<- sqrt(wu)

    for(j in 1:nit){
      for(i in 1:k)
      {
        fui[,i]<-  rwu*exp(-0.5*wu*(yu-mu[i])^2)
        nui[,i]<-  pii[i]*fui[,i]
      }
      pui<- nui/apply(nui,1,sum)
      mu<- apply(wu*yu*pui,2,sum)/apply(wu*pui,2,sum)
      pii<- apply(pui,2,mean)
    }

    # calculate the log likelihood and BIC
    llik<- sum( log( (1/sqrt(2*pi))*apply(nui,1,sum) ))
    bic<- -2*llik + (2*k - 1)*log(n)

    # calculate the covariance matrix and standard errors of the estimates
    # i.e., the dose estimtes in Gy and relative standard errors, and
    # the mixing proportions and standard errors.
    aui<- matrix(0,n,k)
    bui<- matrix(0,n,k)
    for(i in 1:k)
    {
      aui[,i]<- wu*(yu-mu[i])
      bui[,i]<- -wu + (wu*(yu-mu[i]))^2
    }
    delta<- diag(rep(1,k))

    Au<- matrix(0,k-1,k-1)
    Bu<- matrix(0,k-1,k)
    Cu<- matrix(0,k,k)

    for(i in 1:(k-1)){ for(j in 1:(k-1)){
      Au[i,j]<- sum( (pui[,i]/pii[i] - pui[,k]/pii[k])*(pui[,j]/pii[j] -
                                                          pui[,k]/pii[k]) )}}

    for(i in 1:(k-1)){ for(j in 1:k){
      Bu[i,j]<- sum( pui[,j]*aui[,j]*(pui[,i]/pii[i] - pui[,k]/pii[k] -
                                        delta[i,j]/pii[i] + delta[k,j]/pii[k] ) )}}

    for(i in 1:k){ for(j in 1:k){
      Cu[i,j]<- sum( pui[,i]*pui[,j]*aui[,i]*aui[,j] - delta[i,j]*bui[,i]*
                       pui[,i] ) }}

    invvmat<- rbind(cbind(Au,Bu),cbind(t(Bu),Cu))
    vmat<- solve(invvmat, tol=.Machine$double.xmin)
    rek<- sqrt(sum(vmat[1:(k-1),1:(k-1)]))


    # calculate DE, relative standard error, standard error
    dose<- exp(mu)
    re<- sqrt(diag(vmat))[-c(1:(k-1))]
    sed<- dose*re
    estd<- rbind(dose,re,sed)

    # rename proportion
    prop<- pii

    # this calculates the proportional standard error of the proportion of grains
    # in the fitted components. However, the calculation is most likely erroneous.
    # sep<-  c(sqrt(diag(vmat))[c(1:(k-1))],rek)

    # rename proportion
    estp<- prop

    # merge results to a data frame
    blk<- rep("    ",k)
    comp<- rbind(blk,round(estd,4),blk,round(estp,4))
    comp<- data.frame(comp,row.names=c("","dose (Gy)    ","rse(dose)    ",
                                       "se(dose)(Gy)"," ","proportion   "))

    # label results data frame
    cp<- rep("comp",k)
    cn<- c(1:k)
    names(comp)<- paste(cp,cn,sep="")

    # calculate the log likelihood and BIC for a single component -- can
    # be useful to see if there is evidence of more than one component
    mu0<- sum(wu*yu)/sum(wu)
    fu0<-  rwu*exp(-0.5*wu*(yu-mu0)^2)
    L0<- sum( log((1/sqrt(2*pi))*fu0 ) )
    bic0<- -2*L0 + log(n)
    comp0<- round(c(exp(mu0),sigmab,L0,bic0),4)


    ## save results for k components in storage variables
    if(length(n.components)>1) {

      # vector of indices needed for finding the dose rows of the summary
      # matrix - position 1,4,7...n
      pos.n<- seq(from = 1, to = n.components[cnt]*3, by = 3)

      # save results of each iteration to summary matrix
      for(i in 1:n.components[cnt]) {
        comp.n[pos.n[i], cnt]<- round(dose[i], 2) #De
        comp.n[pos.n[i]+1, cnt]<- round(sed[i], 2) #SE
        comp.n[pos.n[i]+2, cnt]<- round(estp[i], 2) #Proportion
      }

      # save BIC and llik of each iteration to corresponding vector
      BIC.n[cnt]<- bic
      LLIK.n[cnt]<- llik

      # merge BIC and llik scores to a single data frame
      results.n<- rbind(BIC = round(BIC.n, 3),
                        llik = round(LLIK.n, 3))

      # save mle matrix and grain probabilities to corresponding vector
      vmat.n[[cnt]]<- vmat
      grain.probability.n[[cnt]]<- as.data.frame(pui)

      # increase counter by one for next iteration
      cnt<- cnt+1
    }#EndOf::save intermediate results
  }##EndOf::calculation loop

  ##============================================================================##
  ## STATISTICAL CHECK
  ##============================================================================##

  if(length(n.components)>1) {

    ## Evaluate maximum log likelihood estimates
    LLIK.significant<- vector(mode = "logical")

    # check if llik is at least three times greater when adding a further
    # component
    for(i in 1:c(length(LLIK.n)-1)) {
      LLIK.significant[i]<- (LLIK.n[i+1]/LLIK.n[i])>3
    }

    ## Find lowest BIC score
    BIC.lowest<- n.components[which.min(BIC.n)]
  }

  ##============================================================================##
  ## OUTPUT
  ##============================================================================##

  if(verbose==TRUE) {

    ## HEADER (always printed)
    cat("\n [calc_FiniteMixture]")

    ##----------------------------------------------------------------------------
    ## OUTPUT WHEN ONLY ONE VALUE FOR n.components IS PROVIDED

    if(length(n.components) == 1) {

      # covariance matrix
      cat(paste("\n\n--- covariance matrix of mle's ---\n\n"))
      print(round(vmat,6))

      # general information on sample and model performance
      cat(paste("\n----------- meta data ------------"))
      cat(paste("\n n:                    ",n))
      cat(paste("\n sigmab:               ",sigmab))
      cat(paste("\n number of components: ",k))
      cat(paste("\n llik:                 ",round(llik,4)))
      cat(paste("\n BIC:                   ",round(bic,3)))

      # fitted components
      cat(paste("\n\n----------- components -----------\n\n"))
      print(comp)


      # print (to 2 decimal places) the estimated probabilities of which component
      # each grain is in -- sometimes useful for diagnostic purposes
      if(grain.probability==TRUE) {
        cat(paste("\n-------- grain probability -------\n\n"))
        print(round(pui,2))
      }

      # output for single component
      cat(paste("\n-------- single component --------"))
      cat(paste("\n mu:                    ", comp0[1]))
      cat(paste("\n sigmab:                ", comp0[2]))
      cat(paste("\n llik:                  ", comp0[3]))
      cat(paste("\n BIC:                   ", comp0[4]))
      cat(paste("\n----------------------------------\n\n"))

    }#EndOf::Output for length(n.components) == 1

    ##----------------------------------------------------------------------------
    ## OUTPUT WHEN ONLY >1 VALUE FOR n.components IS PROVIDED
    if(length(n.components) > 1) {

      ## final labeling of component and BIC/llik matrices
      # create labels
      dose.lab<- paste("c", 1:n.components[length(n.components)],"_dose", sep="")
      se.lab<- paste("c", 1:n.components[length(n.components)],"_se", sep="")
      prop.lab<- paste("c", 1:n.components[length(n.components)],"_prop", sep="")

      # empty vector which stores the labeles in correct order (dose, se, prop)
      n.lab<- vector(mode = "expression",
                     n.components[length(n.components)]*3)

      # loop to store the labels in correct order (dose, se, prop)
      cnt<- 1
      for(i in pos.n) {
        n.lab[i]<- dose.lab[cnt]
        n.lab[i+1]<- se.lab[cnt]
        n.lab[i+2]<- prop.lab[cnt]
        cnt<- cnt+1
      }

      # label columns and rows of summary matrix and BIC/LLIK data frame
      colnames(comp.n)<- n.components[1]:n.components[length(n.components)]
      rownames(comp.n)<- n.lab
      colnames(results.n)<- n.components[1]:n.components[length(n.components)]

      ## CONSOLE OUTPUT
      # general information on sample and model performance
      cat(paste("\n\n----------- meta data ------------"))
      cat(paste("\n n:                    ",n))
      cat(paste("\n sigmab:               ",sigmab))
      cat(paste("\n number of components:  ",n.components[1],"-",
                n.components[length(n.components)], sep=""))

      # output for single component
      cat(paste("\n\n-------- single component --------"))
      cat(paste("\n mu:                    ", comp0[1]))
      cat(paste("\n sigmab:                ", comp0[2]))
      cat(paste("\n llik:                  ", comp0[3]))
      cat(paste("\n BIC:                   ", comp0[4]))

      # print component matrix
      cat(paste("\n\n----------- k components -----------\n"))
      print(comp.n, na.print="")

      # print BIC scores and LLIK estimates
      cat(paste("\n----------- statistical criteria -----------\n"))
      print(results.n)

      ## print evaluation of statistical criteria
      # lowest BIC score
      cat(paste("\n Lowest BIC score for k =", BIC.lowest))

      # first significant increase in LLIK estimates
      if(any(LLIK.significant)!=TRUE) {
        cat(paste("\n No significant increase in maximum log",
                  "likelihood estimates. \n"))
      } else {
        cat(paste("\n First significant increase in maximum log likelihood for",
                  "k =", which(LLIK.significant==TRUE)[1], "\n\n"))
      }

      cat(paste("\n"))
    }#EndOf::Output for length(n.components) > 1
  }

  ##============================================================================##
  ## RETURN VALUES
  ##============================================================================##

  # .@data$meta
  BIC<- data.frame(n.components=k, BIC=bic)
  llik<- data.frame(n.components=k, llik=llik)

  if(length(n.components)>1) {
    BIC.n<- data.frame(n.components=n.components, BIC=BIC.n)
    llik.n<- data.frame(n.components=n.components, llik=LLIK.n)
  }

  # .@data$single.comp
  single.comp<- data.frame(mu=comp0[1],sigmab=comp0[2],
                           llik=comp0[3],BIC=comp0[4])

  # .@data$components
  comp.re<- t(rbind(round(estd,4),round(estp,4)))
  colnames(comp.re)<- c("de","rel_de_err","de_err","proportion")
  comp.re<- comp.re[,-2] # remove the relative error column

  # .@data$grain.probability
  grain.probability<- round(pui, 2)

  summary<- data.frame(comp.re)
  call<- sys.call()
  args<- list(sigmab = sigmab, n.components = n.components)

  # create S4 object
  newRLumResults.calc_FiniteMixture <- set_RLum(
    class = "RLum.Results",
    data = list(
      summary=summary,
      data=data,
      args=args,
      call=call,
      mle=if(length(n.components)==1){vmat}else{vmat.n},
      BIC=if(length(n.components)==1){BIC}else{BIC.n},
      llik=if(length(n.components)==1){llik}else{llik.n},
      grain.probability=if(length(n.components)==1){grain.probability}else{grain.probability.n},
      components=if(length(n.components)==1){comp.re}else{comp.n},
      single.comp=single.comp))

  ##=========##
  ## PLOTTING
  if(plot==TRUE) {
    try(plot_RLum.Results(newRLumResults.calc_FiniteMixture, ...))
  }#endif::plot

  # Return values
  invisible(newRLumResults.calc_FiniteMixture)

}
