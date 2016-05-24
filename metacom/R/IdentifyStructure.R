#' Classifies metacommunities
#'
#' Identifies structure (or quasi-structure) and outputs a classification.
#'
#'
#' @param metacom.obj The result of the `Metacommunity` function, containing a
#' list of 4 elements; the empirical matrix being tested, and results for
#' coherence, turnover, and boundary clumping.
#' @return Ouputs a classification of the metacommunity.
#' @note Quasi structures, as well as 'random' and 'Gleasonian' structures, may
#' not strictly be discernable through the EMS approach, as they rely on
#' inferring a result from a non-significant test ('accepting the null'), which
#' is typically a bad idea.
#' @author John Lefcheck and Tad Dallas
#' @export
#' @examples
#'
#' # Define an interaction matrix
#' data(TestMatrices)
#' pres3c=TestMatrices[[6]]
#'
#' # Analyze metacommunity
#' pres3c.metacom=Metacommunity(pres3c, sims=10, method='r1')
#'
#' # Identify the structure
#' IdentifyStructure(pres3c.metacom)
#'
#'
IdentifyStructure=function(metacom.obj) {
  #Coherence
  if(as.numeric(t(metacom.obj$Coherence)[,3]) >= 0.05) "Random" else
    if(as.numeric(t(metacom.obj$Coherence)[,1]) < as.numeric(t(metacom.obj$Coherence)[,4]) &
         as.numeric(t(metacom.obj$Coherence)[,3]) < 0.05) "Checkerboard (negative coherence)" else
           if(as.numeric(t(metacom.obj$Coherence)[,1]) >= as.numeric(t(metacom.obj$Coherence)[,4]) &
                as.numeric(t(metacom.obj$Coherence)[,3]) < 0.05) {
             print("Positive coherence...")
             #Significant positive turnover
             if(as.numeric(t(metacom.obj$Turnover)[,1]) >= as.numeric(t(metacom.obj$Turnover)[,4]) &
                  as.numeric(t(metacom.obj$Turnover)[,3]) <= 0.05 &
                  metacom.obj$Boundary[,1] >=0 & metacom.obj$Boundary[,2] < 0.05) "Clementsian" else
                    if(as.numeric(t(metacom.obj$Turnover)[,1]) >= as.numeric(t(metacom.obj$Turnover)[,4]) &
                         as.numeric(t(metacom.obj$Turnover)[,3]) <= 0.05 &
                         metacom.obj$Boundary[,1] < 0 & metacom.obj$Boundary[,2] >= 0.05) "Gleasonian" else
                           if(as.numeric(t(metacom.obj$Turnover)[,1]) >= as.numeric(t(metacom.obj$Turnover)[,4]) &
                                as.numeric(t(metacom.obj$Turnover)[,3]) <= 0.05 &
                                metacom.obj$Boundary[,1] < 0 & metacom.obj$Boundary[,2] < 0.05) "Evenly spaced" else
                                  #Significant negative turnover
                                  if(as.numeric(t(metacom.obj$Turnover)[,1]) < as.numeric(t(metacom.obj$Turnover)[,4]) &
                                       as.numeric(t(metacom.obj$Turnover)[,3]) <= 0.05 &
                                       metacom.obj$Boundary[,1] >=0 & metacom.obj$Boundary[,2] < 0.05) "Nested (clumped)" else
                                         if(as.numeric(t(metacom.obj$Turnover)[,1]) < as.numeric(t(metacom.obj$Turnover)[,4]) &
                                              as.numeric(t(metacom.obj$Turnover)[,3]) <= 0.05 &
                                              metacom.obj$Boundary[,1] < 0 & metacom.obj$Boundary[,2] >= 0.05) "Nested (random)" else
                                                if(as.numeric(t(metacom.obj$Turnover)[,1]) < as.numeric(t(metacom.obj$Turnover)[,4]) &
                                                     as.numeric(t(metacom.obj$Turnover)[,3]) <= 0.05 &
                                                     metacom.obj$Boundary[,1] < 0 & metacom.obj$Boundary[,2] < 0.05) "Nested (hyperdispersed" else
                                                       #Non-significant positive turnover
                                                       if(as.numeric(t(metacom.obj$Turnover)[,1]) >= as.numeric(t(metacom.obj$Turnover)[,4]) &
                                                            as.numeric(t(metacom.obj$Turnover)[,3]) > 0.05 &
                                                            metacom.obj$Boundary[,1] >=0 & metacom.obj$Boundary[,2] < 0.05) "Quasi-clementsian" else
                                                              if(as.numeric(t(metacom.obj$Turnover)[,1]) >= as.numeric(t(metacom.obj$Turnover)[,4]) &
                                                                   as.numeric(t(metacom.obj$Turnover)[,3]) > 0.05 &
                                                                   metacom.obj$Boundary[,1] < 0 & metacom.obj$Boundary[,2] >= 0.05) "Quasi-gleasonian" else
                                                                     if(as.numeric(t(metacom.obj$Turnover)[,1]) >= as.numeric(t(metacom.obj$Turnover)[,4]) &
                                                                          as.numeric(t(metacom.obj$Turnover)[,3]) > 0.05 &
                                                                          metacom.obj$Boundary[,1] < 0 & metacom.obj$Boundary[,2] < 0.05) "Quasi-evenly spaced" else
                                                                            #Non-significant negative turnover
                                                                            if(as.numeric(t(metacom.obj$Turnover)[,1]) < as.numeric(t(metacom.obj$Turnover)[,4]) &
                                                                                 as.numeric(t(metacom.obj$Turnover)[,3]) > 0.05 &
                                                                                 metacom.obj$Boundary[,1] >=0 & metacom.obj$Boundary[,2] < 0.05) "Quasi-nested (clumped)" else
                                                                                   if(as.numeric(t(metacom.obj$Turnover)[,1]) < as.numeric(t(metacom.obj$Turnover)[,4]) &
                                                                                        as.numeric(t(metacom.obj$Turnover)[,3]) > 0.05 &
                                                                                        metacom.obj$Boundary[,1] < 0 & metacom.obj$Boundary[,2] >= 0.05) "Quasi-nested (random)" else
                                                                                          if(as.numeric(t(metacom.obj$Turnover)[,1]) < as.numeric(t(metacom.obj$Turnover)[,4]) &
                                                                                               as.numeric(t(metacom.obj$Turnover)[,3]) > 0.05 &
                                                                                               metacom.obj$Boundary[,1] < 0 & metacom.obj$Boundary[,2] < 0.05) "Quasi-nested (hyperdispersed)" } }
