#' Extract Metapop management action details
#' 
#' Extract management action details from RAMAS Metapop .mp files.
#' 
#' @param mp A character string containing the path to a RAMAS Metapop .mp file.
#' @return A \code{data.frame} containing one row per management action, with 
#'   columns: \item{do.action}{Logical. Will the action be performed 
#'   (\code{TRUE}) or ignored (\code{FALSE}).} \item{action}{Factor. The type of
#'   action to be performed.} \item{sourcepop}{The identity of the source 
#'   population.} \item{targetpop}{The identity of the target population.} 
#'   \item{start}{The timestep at which the action will commence.} 
#'   \item{end}{The timestep at which the action will end.} \item{freq}{The 
#'   frequency of the action, in timestep units.} 
#'   \item{after.dispersal}{Logical. Whether the action will be performed after 
#'   (\code{TRUE}) or before (\code{FALSE}) dispersal has taken place.} 
#'   \item{quantity}{Factor. Whether the action affects an absolute 
#'   \code{number} of individuals, or a \code{proportion} of the source 
#'   population.} \item{number}{The absolute number of individuals involved in 
#'   the action.} \item{proportion}{The proportion of the source population 
#'   involved in the action.} \item{fromstage}{The lowest stage involved in the 
#'   action.} \item{tostage}{The highest stage involved in the action.} 
#'   \item{condition}{The condition under which the action will be performed.} 
#'   \item{thr1}{If \code{condition} is either \code{N<thr1} or 
#'   \code{N<thr1_and_N>thr2}, this is the abundance threshold \code{thr1}.} 
#'   \item{thr2}{If \code{condition} is either \code{N>thr2} or 
#'   \code{N<thr1_and_N>thr2}, this is the abundance threshold \code{thr2}.} 
#'   \item{unknownfield}{Unknown.} \item{linear.to}{If \code{condition} is 
#'   \code{linear}, this is the upper quantity (absolute number, or proportion, 
#'   depending on \code{quantity}) towards which linear change will move.} 
#'   \item{linear.lowerN}{If \code{condition} is \code{linear}, this is the 
#'   abundance at which the quantity affected is equal to \code{number} or 
#'   \code{proportion}, depending on the value of \code{quantity}.} 
#'   \item{linear.upperN}{If \code{condition} is \code{linear}, this is the 
#'   abundance at which the quantity affected is equal to \code{linear.to}.} 
#'   \item{N.comprises.stages}{Factor. The stages included in the definition of 
#'   N, when calculating \code{thr1}, \code{thr2}, \code{linear.lowerN} and 
#'   \code{linear.upperN}.} \item{N.comprises.pops}{Factor. The populations
#'   included in the definition of N, when calculating \code{thr1}, \code{thr2},
#'   \code{linear.lowerN} and \code{linear.upperN}.}
#' @export
actions <- function(mp) {
  message("Extracting population management action info from file:\n", mp)
  metapop <- check_mp(mp)
  metapop <- metapop[-(1:6)]
  mgmt.start <- grep('pop mgmnt', metapop)
  n.actions <- as.numeric(gsub('\\D', '', metapop[mgmt.start]))
  if(n.actions==0) stop(sprintf('No management actions in %s.', mp))
  metapop <- metapop[(mgmt.start + 1):(mgmt.start + n.actions)]
  metapop <- gsub('\\s+', ' ', metapop)
  metapop <- as.data.frame(apply(do.call(rbind, strsplit(metapop, ' ')), 
                                 2, as.numeric))
  colnames(metapop) <- c(
    'do.action', 'action', 'sourcepop', 'targetpop', 'start', 'end', 'freq', 
    'after.dispersal', 'quantity', 'number', 'proportion', 'fromstage', 
    'tostage', 'condition', 'thr1', 'thr2', 'unknownfield', 'linear.to', 
    'linear.lowerN', 'linear.upperN', 'N.comprises.stages', 'N.comprises.pops')
  metapop$do.action <- metapop$do.action == 1
  metapop$action <- factor(metapop$action, 0:2, c('harvest', 'intro', 'translo'))
  metapop$after.dispersal <- metapop$after.dispersal == 1 
  metapop$quantity <- factor(metapop$quantity, 0:1, c('number', 'proportion'))
  metapop$number <- ifelse(metapop$quantity == 'proportion', 
                           NA, metapop$quantity)
  metapop$proportion <- ifelse(metapop$quantity == 'proportion', 
                               metapop$proportion, NA)
  metapop$condition <- factor(
    metapop$condition, 0:4, 
    c('none', 'N<thr1', 'N>thr2', 'N<thr1_and_N>thr2', 'linear'))
  metapop$thr1 <- ifelse(metapop$condition %in% 
                           c('none', 'N>thr2', 'linear'), NA, metapop$thr1)
  metapop$thr2 <- ifelse(metapop$condition %in% 
                           c('none', 'N<thr1', 'linear'), NA, metapop$thr2)
  metapop$linear.to <- ifelse(metapop$condition != 'linear', 
                              NA, metapop$linear.to)
  metapop$linear.lowerN <- ifelse(metapop$condition != 'linear', 
                                  NA, metapop$linear.lowerN)
  metapop$linear.upperN <- ifelse(metapop$condition != 'linear', 
                                  NA, metapop$linear.upperN)  
  metapop$N.comprises.stages <- ifelse(metapop$condition == 'none',
                                       NA, metapop$N.comprises.stages)  
  metapop$N.comprises.pops <- ifelse(metapop$condition == 'none', 
                                     NA, metapop$N.comprises.pops)  
  metapop$N.comprises.stages <- factor(
    metapop$N.comprises.stages, 0:2, 
    c('each stage', 'all selected stages', 'all stages'))
  metapop$N.comprises.pops <- 
    factor(metapop$N.comprises.stages, c(0, 2), c('each pop', 'all pops'))
  
  metapop
}
