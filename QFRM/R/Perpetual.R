#'@title        Perpetual option valuation via Black-Scholes (BS) model
#'@description  An exotic option is an option which has features making it more complex than commonly traded options.
#'              A perpetual option is non-standard financial option with no fixed maturity and no exercise limit. 
#'              While the life of a standard option can vary from a few days to several years, 
#'              a perpetual option (XPO) can be exercised at any time. 
#'              Perpetual options are considered an American option. European options can be exercised only on the option's maturity date. 
#'@author       Kim Raath, Department of Statistics, Rice University, Spring 2015.  
#'
#'@param        o AN object of class \code{OptPx}
#'@return       A list of class \code{Perpetual.BS} consisting of the input object \code{OptPx}
#'
#'@references   Chi-Guhn Lee, \emph{The Black-Scholes Formula}, Courses, Notes, Note2, Sec 1.5 and 1.6
#'              \url{http://www.mie.utoronto.ca/courses/mie566f/materials/note2.pdf} 
#'@examples     
#'
#' #Perpetual American Call and Put
#' #Verify pricing with \url{http://www.coggit.com/freetools}
#' (o <- PerpetualBS())$PxBS # Approximately valued at $8.54 
#'              
#' #This example should produce approximately $33.66
#' o = Opt(Style="Perpetual", Right='Put', S0=50, K=55)
#' o = OptPx(o, r = .03, q = 0.1, vol = .4)
#' (o = PerpetualBS(o))$PxBS
#'             
#' #This example should produce approximately $10.87
#' o = Opt(Style="Perpetual", Right='Call', S0=50, K=55)
#' o = OptPx(o, r = .03, q = 0.1, vol = .4)
#' (o <- PerpetualBS(o))$PxBS
#'              
#' @export
#' 
PerpetualBS = function(o=OptPx(Opt(Style='Perpetual'), q=0.1)){   #q cannot be zero for perpetual options, specify the style of this option
  stopifnot(is.OptPx(o), o$Style$Perpetual); #Do not complete code if the style is not perpetual
  #Perpetual option valuation using Black-Scholes (BS) model - see references above for more detail 
  sigma_sqr =with(o, vol^2)   #volitility (std Deviation) is sigma, variance is sigma squared
  
  h1  =with(o, 0.5 - ( (r-q)/sigma_sqr ))    #adding the first two values in the formula together
  
  #Value of $1 received when S first reaches H from below = (S/H)^h2 where h2 is the value below
  h2 = with(o, h1 + sqrt( ( ((r-q)/sigma_sqr-0.5)^2 ) + 2*r/sigma_sqr ))   
  
  #Value of $1 received when S reaches H from above  = (S/H)^h3 where h3 is the value below
  h3 =with(o, h1 - sqrt( ( ((r-q)/sigma_sqr-0.5)^2 ) + 2*r/sigma_sqr ))
  
  c = with(o, (K/(h2-1))*( ( ((h2-1)/h2)*(S0/K) )^h2 ))   #Price of perpetual call
  p = with(o, (K/(h3-1))*( ( ((h3-1)/h3)*(S0/K) )^h3 ))   #Price of perpetual put
  
  if (o$Right$Call) o$PxBS = c          #if Call give price for Perpetual call option
  if (o$Right$Put) o$PxBS = p           #if Put give price for Perpetual put option
  
  return(o)       #return object o
};
