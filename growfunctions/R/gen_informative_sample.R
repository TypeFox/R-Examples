#' Generate a finite population and take an informative single or two-stage sample.
#'
#' Used to compare performance of sample design-weighted and unweighted estimation procedures.
#'
#' @param clustering Boolean input on whether want population generated from clusters of covariance 
#'        parameters.  Defaults to \code{clustering = FALSE}
#' @param two_stage Boolean input on whether want two stage sampling, with first stage defining set
#'        of \code{L} blocks, where membership in blocks determined by quantiles of observation unit
#'        variance functions.  (They are structured like strata, though they are sub-sampled).
#' @param theta A numeric vector of global covariance parameters in the case of \code{clustering = FALSE}.
#'        The length, \code{P}, of \code{theta} must be consistent with the selected \code{gp_type}.  
#'        Defaults to \code{theta = c(0.30.7,1.0)} in the case of \code{clustering = FALSE}.
#' @param M Scalar input denoting number of clusters to employ if \code{clustering = TRUE}. Defaults to
#'        \code{M = 3}
#' @param theta_star An \emph{P x M} matrix of cluster location values associated with the choice of 
#'        \code{M} and the selected \code{gp_type}. Defaults to 
#'        \code{matrix(c(0.3,0.3,0.3,0.31,0.72,2.04,0.58,0.83,1.00),3,3,byrow=TRUE))}.
#' @param gp_type Input of choice for covariance matrix formulation to be used to generate the functions
#'        for the \code{N} population units.  Choices are \code{c("se","rq")}, where \code{"se"} denotes
#'        the squared exponential covariance function and \code{"rq"} denotes the rational quadratic.
#'        Defaults to \code{gp_type = "se"} 
#' @param N A scalar input denoting the number of population units (or establishments).
#' @param T A scalar input denoting the number of time points in each of \code{N}, \emph{T x 1} functions
#'        that contribute to the \emph{N x T} population data matrix, \code{y}.  Defaults to \code{T = 15}.
#' @param L A scalar input that denotes the number of blocks in which to assign the population
#'        units to be sub-sampled in the first stage of sampling.  
#'        Defaults to \code{L = 10}.
#' @param R  A scalar input that denotes the number of blocks to sample from \code{L  = 10} with
#'        probability proportional to the average variance of member functions in each block.
#' @param I A scalar input denoting the number of strata to form within each block.  Population units
#'        are divided into equally-sized strata based on variance quantiles. Defaults to \code{I  = 4}.
#' @param n Sample size to be generated.  Both an informative sample under either single
#'        (\code{two_stage = FALSE}) or 2-stage (\code{two_stage = TRUE}) sample is taken, along with
#'        a non-informative, \emph{iid} sample of the same size (\code{n}) from the finite population
#'        (generated with (\code{clustering = TRUE}) or without clustering). Defaults to \code{n = 770}.
#' @param incl_gradient A character input on whether stratum probabilities from lowest-to-highest
#'        is to \code{"high"}, in which case they are proportional to the exponential of the
#'        cluster number.  If set to \code{"medium"} , the inclusion probabilities are proportional
#'        to the square of the cluster number.  Note that population units are assigned to each
#'        stratum proportional to a progressively increasing quantile variance.  The 
#'        \code{incl_gradient} setting is used for both \code{two_stage = TRUE}, in which
#'        case it is applied to strata within block, as well as \code{two_stage = FALSE},
#'        in which case a simple stratified random sample is conducted.  Defaults to 
#'        \code{incl_gradient = "medium"}   
#' @param noise_to_signal A numeric input in the interval, \code{(0,1)}, denoting the ratio of noise
#'        variance to the average variance of the generated functions, \code{bb_i}.  Defaults to 
#'        \code{noise_to_signal = 0.05}               
#' @return A list object named \code{dat_sim} containing objects related to the generated sample
#'         finite population, the informative sample and the non-informative, \emph{iid}, sample. 
#'         Some important objects, include:
#'     \item{H}{A vector of length \code{N}, the population size, with cluster assignments
#'             for each establishment (unit) in \code{1,..M} clusters.}
#'     \item{map.tot}{A \code{data.frame} object including unit label identifiers
#'                   (under \code{establishment}),
#'                    the cluster assignment (if \code{clustering = TRUE}), 
#'                    the block (if\code{two_stage = TRUE}) and stratum assignments 
#'                    and the sample inclusion probabilities.}
#'     \item{map.obs}{A \code{data.frame} object configured the same as \code{map.tot}, only
#'                  confined to those establishments/units selected into the \emph{informative}
#'                  sample of size \code{n}.}
#'     \item{map.iid}{A \code{data.frame} object configured the same as \code{map.tot}, only
#'                  confined to those establishments/units selected into the \emph{non-informative},
#'                  iid sample of size \code{n}.}
#'     \item{(y,bb)}{\emph{N x T} \code{matrix} objects containing data responses and de-noised '
#'                  functions, respectively, for each of the \code{N} population units. The order
#'                  of the \code{N} units is consistent with \code{map}.}
#'     \item{(y_obs,bb_obs)}{\emph{N x T} \code{matrix} objects containing observed responses and de-noised '
#'                  functions, respectively, for each of the \code{n} units sampled under an
#'                  informative sampling design. The order of the \code{n} units is consistent
#'                  with \code{map_obs}.}
#'     \item{(y_iid,bb_iid)}{\emph{N x T} \code{matrix} objects containing observed responses and de-noised '
#'                  functions, respectively, for each of the \code{n} units sampled under a
#'                  non-informative / iid sampling design. The order of the \code{n} units is consistent
#'                  with \code{map_iid}.}
#' @seealso \code{\link{gpdpgrow}}, \code{\link{gmrfdpgrow}}
#' @examples 
#' \dontrun{
#' library(growfunctions)
#' ## use gen_informative_sample() to generate an 
#' ## N X T population drawn from a dependent GP
#' ## By default, 3 clusters are used to generate 
#' ## the population.
#' ## A single stage stratified random sample of size n 
#' ## is drawn from the population using I = 4 strata. 
#' ## The resulting sample is informative in that the 
#' ## distribution for this sample is
#' ## different from the population from which 
#' ## it was drawn because the strata inclusion
#' ## probabilities are proportional to a feature 
#' ## of the response, y (in the case, the variance.
#' ## The stratified random sample over-samples 
#' ## large variance strata).
#' ## (The user may also select a 2-stage 
#' ## sample with the first stage
#' ## sampling "blocks" of the population and 
#' ## the second stage sampling strata within blocks). 
#' dat_sim        <- gen_informative_sample(N = 10000, 
#'                                 n = 500, T = 10,
#'                                 noise_to_signal = 0.1)
#'
#' ## extract n x T observed sample under informative
#' ## stratified sampling design.
#' y_obs                       <- dat_sim$y_obs
#' T                           <- ncol(y_obs)
#' }
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases gen_informative_sample
#' @export
gen_informative_sample        <- function(clustering = TRUE, two_stage = FALSE, theta = c(0.2,.7,1.0), 
                                          M = 3,
                                          theta_star = matrix(c(0.3,0.3,0.3,0.31,0.72,2.04,0.58,
                                                                0.83,1.00),3,3,byrow=TRUE),
                                          gp_type = "rq", N = 10000, T = 15, L = 10, R = 8, I = 4, n = 750, 
                                          noise_to_signal = 0.05, incl_gradient = "medium")
{
     ###############
     ## some quick input checks
     ###############
     ## number of covariance parameters
     if( (gp_type == "rq")  )
     {
          if(clustering == FALSE)
          {
               if( length(theta) != 3 ){stop("\nMust specify P = 3 parameters for rational quadratic option.\n")}
          }else{
               if( nrow(theta_star) != 3 ){stop("\nMust specify P = 3 parameters for rational 
                                                quadratic option.\n")}
          }
          
     }else{ ## "se"
          if(clustering == FALSE)
          {
               if( length(theta) != 2 ){stop("\nMust specify P = 2 parameters for 
                                                  squared exponential option.\n")}
          }else{
               if( nrow(theta) != 2 ){stop("\nMust specify P = 2 parameters for squared 
                                                  exponential option.\n")}
          }
     } ## end conditional check on number of GP covariance parameters, P
     
     ## number of clusters
     if( clustering == TRUE )
     {
          if( ncol(theta_star) != M ){stop("\nMust specify M columns for theta_star.\n")}
     }
     
     #############
     ## generate finite population of size N
     ############# 
     ## Create T x T GP distance matrix, Omega
     Omega     <- sapply(1:T,function(i){
          sapply(1:T,function(j){
               d = (i-j)^2
          })
     })
     zeros.T          <- rep(0,T)
     
     if(clustering == TRUE)
     {
          ## randomly assign establishments to clusters
          cluster                       <- vector("list",M)
          units_left                    <- 1:N
          n_sample                      <- rep(floor(N/M),M)
          n_left                        <- N - sum(n_sample)
          if( n_left > 0 )
          {
               pos_m               <- sample(1:M,n_left,replace = TRUE)
               n_sample[pos_m]     <- n_sample[pos_m] + 1
          } ## end conditional statement on whether any unallocated units to clusters
          
          H         <- vector("numeric",N)
          for( m in 1:M )
          {
               cluster[[m]]             <- sample(units_left,n_sample[m],replace=FALSE)
               H[cluster[[m]]]          <- m
               units_left               <- setdiff(units_left,cluster[[m]])    
          } ## end loop m over clusters to sample units to each cluster
          
          if(gp_type == "rq")
          {
               G         <- lapply(1:M,function(m){
                              1/theta_star[1,m] * 
                              (1+Omega/(theta_star[2,m]*theta_star[3,m]))^
                              (-theta_star[3,m])
                         })
               ## Generate set of N x T functions based on hyperparameters associated 
               ## with cluster memberships 
               bb        <- t(sapply(1:N,function(i){     
                              ut	<- rmvnorm( 1,zeros.T, G[[H[i]]]  )
                              ut}))    
          }else{ ## gp_type == "se"
               G         <- lapply(1:M,function(m){
                              1/theta_star[1,m] * 
                                   exp(-Omega/theta_star[2,m])
               })
               ## Generate set of N x T functions based on hyperparameters associated 
               ## with cluster memberships 
               bb        <- t(sapply(1:N,function(i){     
                              ut     <- rmvnorm( 1,zeros.T, G[[H[i]]]  )
                              ut}))    
          } ## end generation of N X T functions, bb, for clustering
     }else{ ## clustering == FALSE
          if( gp_type == "rq" )
          {
               bb        <- t(sapply(1:N,function(i){     
                              G    <- 1/theta[1] * 
                                   (1+Omega/(theta[2]*theta[3]))^
                                   (-theta[3])
                              ut	<- rmvnorm( 1,zeros.T, G  )
                              ut}))
          }else{ ## se
               bb        <- t(sapply(1:N,function(i){     
                              G    <- 1/theta[1] * 
                                   exp(-Omega/theta[2])
                              ut	<- rmvnorm( 1,zeros.T, G  )
                              ut}))
               
          } ## end generation of N X T functions, bb, for non-clustering 
     } ##  end generation of N x T functions, bb
     
     ## generate observed N x T data, y
     ## add some noise
     bb_mvar             <- mean(apply(bb,1,var))
     sig_e               <- sqrt(bb_mvar) * noise_to_signal
     y                   <- bb + matrix(rnorm(N*T,0,sig_e),N,T)
     
     ## plots of population level functions
     ## plot emphasizing type of curves in each cluster
     
     ## true functions, f
     res_pattern     <- plot_cluster(bb,H,y.axis.label = bquote("Population true "~Delta~" Counts, f "),
                                     sample_rate = 0.05)
     res_shape       <- plot_cluster(bb,H,y.axis.label = bquote("Population true "~Delta~" Counts, f "),
                                     smoother = FALSE, fade = 1, sample_rate = .001)
     ## noisy data values, y
     res_y          <- plot_cluster(y,H,y.axis.label = bquote("Population noisy "~Delta~" Counts, y "),
                                         sample_rate = 0.05)
     
     ## assign block (2-stage) labels to population
     ## based on variance quantiles of the T x 1, {y_{i}}
     ## also assign block inclusion probabilities to be proportional to the average block variance
     
     ## create mapping of N establishments to cluster (clustering, they are real; no clustering, just dummy)
     ## and to variance of each T x 1, y_{i}
     ## dummy cluster assignment that puts all units in a single cluster for 1-stage
     if(clustering == FALSE){H  <- rep(1,N)}
     map_tot             <- data.frame(1:N,H)
     names(map_tot)      <- c("establishment","cluster")
     map_tot             <- map_tot[order(map_tot$establishment),]
     map_tot$y_var       <- apply(y,1,var)
     
     if( two_stage == TRUE)
     {
          ####
          ## split finite population into L disjoint blocks
          ####     
          ## create L representative sub-populations from which to draw blocks
          est_clu <- map_clu <- map_sel <- y_obs <- bb_obs <- res_yobs_shape <- vector("list",L)
          boundaries     <- incl_prob      <- n_ind <- N_ind <- vector("list",L)
          ## disjointly (randomly) divide finite population of size N into L equally-sized sub-populations 
          cut_cats                      <- rep(1/L,L)
          cut_probs                     <- cumsum( cut_cats )[1:(L-1)]
          cut_points                    <- quantile(map_tot$y_var,probs = cut_probs)
          est_clu[[1]]                  <- map_tot$establishment[ map_tot$y_var <= cut_points[1] ] ## low var
          for( l in 2:(L-1) )
          {
               est_clu[[l]]                  <- map_tot$establishment[ (map_tot$y_var > cut_points[l-1]) &
                                                                      (map_tot$y_var <= cut_points[l]) ]
          }
          est_clu[[L]]                  <- map_tot$establishment[ map_tot$y_var > cut_points[(L-1)] ] ## hi var
               
          ## use est_clu to sub-divide map_tot data.frame, by block
          for( l in 1:L )
          {
               map_clu[[l]]            <- subset( map_tot, establishment %in% est_clu[[l]] )
          }
                   
          ## set selection probabilities for blocks
          block_probs      <- sapply(1:L,function(l){
                    x = mean( map_clu[[l]]$y_var )
                    x})
          block_probs      <- block_probs / sum( block_probs )
     } ## end assignment of blocks for 2-stage sampling
     
     ## Assign strata (within block for 2-stage) and conduct stratified sample (within
     ## each of the L = 10 blocks for 2-stage)
     
     ## define quantile probabilities for I strata (in each of L clusters for 2-stage)
     strat_cats                         <- rep(1/I,I)
     strat_probs                        <- cumsum( strat_cats )[1:(I-1)]
     
     if( two_stage == TRUE )
     {  
          for( l in 1:L )
          {
               ## assign I Industry strata to establishments in each cluster based on variance quantiles    
               boundaries[[l]]                    <- quantile(map_clu[[l]]$y_var,probs = strat_probs)
               map_clu[[l]]$Industry              <- 1 ## lowest variance stratum
               map_clu[[l]]$Industry[map_clu[[l]]$y_var > boundaries[[l]][(I-1)]]  <- I ## highest var cluster
               for( i in 2:(I-1) )
               {
                    map_clu[[l]]$Industry[(map_clu[[l]]$y_var > boundaries[[l]][i-1]) 
                                          & (map_clu[[l]]$y_var <= boundaries[[l]][i])]  <- i
               }
               
               ## draw informative sample from map_clu[[l]] within each of I strata
               ## strata population and sample sizes
               n_max_l                  <- floor(nrow(map_clu[[l]]))
               n_l                      <- min( n_max_l, (750/R) ) ## will only select R of L clusters
               ## more samples from high variance strata
               strat_labs               <- 1:I
               ## inclusion proportions to allocate sample of size n_l
               if( incl_gradient == "high" )
               {
                    incl_prop_l              <- exp(strat_labs) ## sampling more of higher var strata
               }else{ ## medium rate of increased sample rate across strata
                    incl_prop_l              <- strat_labs^2
               } ## end setting sample allocation proportions across I strata
               incl_prop_l              <- incl_prop_l/(incl_prop_l)
               n_ind_l                  <- floor( n_l*incl_prop_l )
               if( sum(n_ind_l) < n_l ){n_ind_l[I] <- (n_ind_l[I] + n_l - sum(n_ind_l))}
               N_ind_l                  <- sapply(1:I,function(i){
                    x = length(map_clu[[l]]$Industry[map_clu[[l]]$Industry == i])      
               })
               ## compose inclusion probabilities from product of 1st stage and 2nd stage weights
               incl_prob[[l]]              <- (n_ind_l / N_ind_l) * block_probs[l]
               n_ind[[l]]                  <- n_ind_l
               N_ind[[l]]                  <- N_ind_l
               map_clu[[l]]$incl_prob      <- 0
               obs_units_l                 <- vector("list",I) ## by stratum
               ## draw sample of n_ind_l[i] units from each stratum
               for( i in 1:I )
               {
                    map_clu[[l]]$incl_prob[map_clu[[l]]$Industry == i] <- incl_prob[[l]][i]
                    map_clu_li               <- subset(map_clu[[l]], Industry == i)
                    obs_units_l[[i]]         <- sample(map_clu_li$establishment,n_ind_l[i],replace=FALSE)
               }
               
               ## draw observational units and write inclusion indicator to map_ind
               obs_units_l              <- unlist(obs_units_l)
               map_clu[[l]]$incl_ind    <- 0
               map_clu[[l]]$incl_ind[map_clu[[l]]$establishment %in% obs_units_l] <- 1
               ## pull sub-set of map_clu for observed units
               map_sel[[l]]             <- subset(map_clu[[l]], incl_ind == 1)
               ## must sort map_sel by establishment in order to ensure that order in (bb_obs,y_obs)
               ## corresponds to establishment order
               map_sel[[l]]                  <- map_sel[[l]][order(map_sel[[l]]$establishment),]  
               bb_obs[[l]]                   <- bb[map_sel[[l]]$establishment,] 
               y_obs[[l]]                    <- y[map_sel[[l]]$establishment,]   
          } ## end loop l over clusters
          
          ###
          ## sample cluster 
          ###
          sel_blocks           <- sort( sample(1:L,R,prob=block_probs) )
          map_obs              <- lapply(1:R,function(r){
               x <- map_sel[[ sel_blocks[r] ]]
               x})
          map_obs             <- do.call("rbind", map_obs) ## retains data.frame 
                                                        ##structure with names from map_clu
          
          bb_sel <- y_sel <- vector("list",R)
          for( r in 1:R )
          {
               bb_sel[[r]]    <- bb_obs[[ sel_blocks[r] ]]
               y_sel[[r]]   <- y_obs[[ sel_blocks[r] ]]
          }
          bb_obs         <- do.call("rbind", bb_sel)
          y_obs          <- do.call("rbind", y_sel)
     }else{ ## 1-stage
          ## map contains c("establishment","cluster","y_var") labels and is sorted by establishment
          map_tot                 <- map_tot[order(map_tot$establishment),]
          
          ## assign strata to population units
          boundaries                         <- quantile(map_tot$y_var,probs = strat_probs)
          map_tot$Industry                   <- 1 ## lowest variance stratum
          map_tot$Industry[map_tot$y_var > boundaries[(I-1)]]  <- I ## highest var cluster
          for( i in 2:(I-1) )
          {
               map_tot$Industry[(map_tot$y_var > boundaries[i-1]) 
                                     & (map_tot$y_var <= boundaries[i])]  <- i
          } ## end loop to assign strata to N units
          
          ## set stratum sample sizes
          ## incl_prop            <- c(0.1,0.15,0.3,0.45) 
          strat_labels         <- 1:I
          if(incl_gradient == "high")
          {
               incl_prop            <- exp(strat_labels)  
          }else{  ## use ^2 instead of exp()
               incl_prop            <- strat_labels^2
          }
          incl_prop            <- incl_prop / sum(incl_prop)
          n_ind                <- floor( n*incl_prop )
          if( sum(n_ind) < n ){n_ind[I] <- (n_ind[I] + n - sum(n_ind))}
          map_tot$incl_prob    <- 0
          ## draw a stratified random sample
          obs_units            <- vector("list",I)
          N_ind                <- vector("numeric",I)
          for( i in 1:I )
          {
               est_cand_i       <- map_tot$establishment[map_tot$Industry == i]
               N_ind[i]         <- length(est_cand_i)
               obs_units[[i]]   <- sample(est_cand_i,n_ind[i],replace = FALSE)
               ## write inclusion probabilities to population units
               map_tot$incl_prob[map_tot$Industry == i] <- ( n_ind[i] / N_ind[i] )
          } ## end loop over strata for drawing sample
          obs_units           <- unlist(obs_units)
          
          ## write inclusion indicator to map_tot
          map_tot$incl_ind    <- 0
          map_tot$incl_ind[map_tot$establishment %in% obs_units] <- 1
          ## pull sub-set of map_tot for observed units
          map_obs             <- subset(map_tot, incl_ind == 1)
          ## must sort map_obs by establishment in order to ensure that order in (bb_obs,y_obs)
          ## corresponds to establishment order
          map_obs             <- map_obs[order(map_obs$establishment),]  
          bb_obs              <- bb[map_obs$establishment,] 
          y_obs               <- y[map_obs$establishment,]
          
          ## plot observed functions by cluster - pattern attenuated, but remains
          res_y_obs           <- plot_cluster(y_obs,map_obs$cluster,
                                              y.axis.label = bquote("1-Stage Sample "~Delta~" Counts, y "),
                                              sample_rate = 0.5)
          
          ## plot emphasizes shape of functions generated in each cluster
          res_yobs_shape      <- plot_cluster(y_obs,map_obs$cluster,
                                               y.axis.label = bquote("1-Stage Sample "~Delta~" Counts, b "),
                                               smoother = FALSE, fade = 1, sample_rate = .1)
          
     } ## end assignment of strata labels and drawing of informative sample (deposited into map_obs object)
     
     ## adjust incl_prob for observed units to sum to n
     map_obs$incl_prob_orig        <- map_obs$incl_prob
     map_obs$incl_prob             <- map_obs$incl_prob  * (sum(1/map_obs$incl_prob)/nrow(map_obs))
     
     ## draw iid sample
     iid_units            <- sample(map_tot$establishment,n,replace=FALSE)
     map_tot$iid_ind      <- 0
     map_tot$iid_ind[map_tot$establishment %in% iid_units] <- 1
     ## pull sub-set of map_ind for observed units
     map_iid             <- subset(map_tot, iid_ind == 1)
     ## must sort map_obs by establishment in order to ensure that order in (bb_obs,y_obs)
     ## corresponds to establishment order
     map_iid             <- map_iid[order(map_iid$establishment),]  
     bb_iid              <- bb[map_iid$establishment,] 
     y_iid               <- y[map_iid$establishment,]
     ## set indicator in map_obs (informative sample) where overlaps with iid sample.
     map_overlap         <- subset(map_iid, establishment %in% map_obs$establishment)
     map_obs$overlap     <- 0
     map_obs$overlap[map_obs$establishment %in% 
                          map_overlap$establishment] = 1
     
     ## define dat_sim return object
     dat_sim   <- list(bb = bb, y = y, H = H, map_tot = map_tot, map_obs = map_obs, 
                         bb_obs = bb_obs, y_obs = y_obs, n_ind = n_ind, 
                         map_iid = map_iid, bb_iid = bb_iid, y_iid = y_iid,
                         sig_e = sig_e, noise_to_signal = noise_to_signal)
     ## memo: for 2-stage, n_ind is an L length list with I elements each (ordered by variance quantile).
     ## for 1-stage it is an I length vector, ordered by variance quantile
     if(two_stage == TRUE)
     {
          dat_sim$N_ind = N_ind; dat_sim$sel_blocks = sel_blocks; 
          dat_sim$block_probs = block_probs
     }
     
     if(clustering == TRUE)
     {
          dat_sim$theta_star = theta_star
     }else{ ## no clustering
          dat_sim$theta <- theta
     } ## end writing out dat_sim return object
     return(dat_sim)  
     
     establishment <- Industry <- incl_ind <- iid_ind <- NULL
     
} ## end function gen_informative_sample