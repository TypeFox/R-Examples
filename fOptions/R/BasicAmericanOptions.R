
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA


################################################################################
# FUNCTION:                  DESCRIPTION:
#  RollGeskeWhaleyOption      Roll-Geske-Whaley Calls on Dividend Paying Stocks
#  BAWAmericanApproxOption    Barone-Adesi and Whaley Approximation
#  BSAmericanApproxOption     Bjerksund and Stensland Approximation
################################################################################


RollGeskeWhaleyOption = 
function(S, X, time1, Time2, r, D, sigma, title = NULL, description = NULL) 
{   # A function implemented by Diethelm Wuertz
 
    # Description:
    #   Calculates the option price of an American call on a stock
    #   paying a single dividend with specified time to divident
    #   payout. The option valuation formula derived by Roll, Geske 
    #   and Whaley is used.
    
    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas
    
    # FUNCTION:
    
    # Settings:
    big = 100000000
    eps = 1.0e-5
    t1 = time1
    T2 = Time2
    
    # Compute:
    Sx = S - D * exp(-r * t1)
    if(D <= X * (1 - exp(-r*(T2-t1)))) {         
        result = GBSOption("c", Sx, X, T2, r, b=r, sigma)@price
        cat("\nWarning: Not optimal to exercise\n")
        return(result) }
    ci = GBSOption("c", S, X, T2-t1, r, b=r, sigma)@price
    HighS = S
    while ( ci-HighS-D+X > 0 && HighS < big ) {
        HighS = HighS * 2
        ci = GBSOption("c", HighS, X, T2-t1, r, b=r, sigma)@price }
    if(HighS > big) {
        result = GBSOption("c", Sx, X, T2, r, b=r, sigma)@price
        stop()}
    LowS = 0
    I = HighS * 0.5
    ci = GBSOption("c", I, X, T2-t1, r, b=r, sigma)@price 
    # Search algorithm to find the critical stock price I
    while ( abs(ci-I-D+X) > eps && HighS - LowS > eps ) {
         if(ci-I-D+X < 0 ) { HighS = I }
        else { LowS = I }
        I = (HighS + LowS) / 2
        ci = GBSOption("c", I, X, T2-t1, r, b=r, sigma)@price }
    a1 = (log(Sx/X) + (r+sigma^2/2)*T2) / (sigma*sqrt(T2))
    a2 = a1 - sigma*sqrt(T2)
    b1 = (log(Sx/I) + (r+sigma^2/2)*t1) / (sigma*sqrt(t1))
    b2 = b1 - sigma*sqrt(t1)
    result = Sx*CND(b1) + Sx*CBND(a1,-b1,-sqrt(t1/T2)) -
        X*exp(-r*T2)*CBND(a2,-b2,-sqrt(t1/T2)) - 
            (X-D)*exp(-r*t1)*CND(b2)
    
    # Parameters:
    # S, X, time1, Time2, r, D, sigma
    param = list()
    param$S = S
    param$X = X
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$D = D
    param$sigma = sigma
    
    # Add title and description:
    if(is.null(title)) title = "Roll Geske Whaley Option"
    if(is.null(description)) description = as.character(date())
    
    # Return Value:
    new("fOPTION", 
        call = match.call(),
        parameters = param,
        price = result, 
        title = title,
        description = description
        )      
}


# ******************************************************************************


BAWAmericanApproxOption = 
function(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma, title = NULL, 
description = NULL)
{   # A function implemented by Diethelm Wuertz
 
    # Description:
    #   Calculates the option price of an American call or put
    #   option on an underlying asset for a given cost-of-carry rate.
    #   The quadratic approximation method by Barone-Adesi and
    #   Whaley is used.

    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas
    
    # FUNCTION:
    
    # Settings:
    TypeFlag = TypeFlag[1]
    
    # Compute:
    if(TypeFlag == "c") {
        result = .BAWAmCallApproxOption(S, X, Time, r, b, sigma) }
    if(TypeFlag == "p") {      
        result = .BAWAmPutApproxOption(S, X, Time, r, b, sigma) }
       
    # Parameters:
    # TypeFlag = c("c", "p"), S, X, Time, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma
    
    # Add title and description:
    if(is.null(title)) title = "BAW American Approximated Option"
    if(is.null(description)) description = as.character(date())
    
    # Return Value:
    new("fOPTION", 
        call = match.call(),
        parameters = param,
        price = result, 
        title = title,
        description = description
        )      
}


.BAWAmCallApproxOption <- 
function(S, X, Time, r, b, sigma) 
{
    # Internal Function - The Call:
        
    # Compute:
    if(b >= r) {
        result = GBSOption("c", S, X, Time, r, b, sigma)@price }
    else {
        Sk = .bawKc(X, Time, r, b, sigma)
        n = 2*b/sigma^2
        K = 2*r/(sigma^2*(1-exp(-r*Time)))
        d1 = (log(Sk/X)+(b+sigma^2/2)*Time)/(sigma*sqrt(Time))
        Q2 = (-(n-1)+sqrt((n-1)^2+4*K))/2
        a2 = (Sk/Q2)*(1-exp((b-r)*Time)*CND(d1))
        if(S < Sk) {
            result = GBSOption("c", S, X, Time, r, b, sigma)@price +
                a2*(S/Sk)^Q2 
        } else {
            result = S-X 
        } 
    }
    
    # Return Value:
    result 
}


.bawKc <- 
function(X, Time, r, b, sigma) 
{   
    # Newton Raphson algorithm to solve for the critical commodity 
    # price for a Call.
    # Calculation of seed value, Si
    n = 2*b/sigma^2
    m = 2*r/sigma^2
    q2u = (-(n-1)+sqrt((n-1)^2+4*m))/2
    Su = X/(1-1/q2u)
    h2 = -(b*Time+2*sigma*sqrt(Time))*X/(Su-X)
    Si = X+(Su-X)*(1-exp(h2))
    K = 2*r/(sigma^2*(1-exp(-r*Time)))
    d1 = (log(Si/X)+(b+sigma^2/2)*Time)/(sigma*sqrt(Time))
    Q2 = (-(n-1)+sqrt((n-1)^2+4*K))/2
    LHS = Si-X
    RHS = GBSOption("c", Si, X, Time, r, b, sigma)@price + 
        (1-exp((b-r)*Time)*CND(d1))*Si/Q2
    bi = exp((b-r)*Time)*CND(d1)*(1-1/Q2) +
        (1-exp((b-r)*Time)*CND(d1)/(sigma*sqrt(Time)))/Q2
    E = 0.000001
    
    # Newton Raphson algorithm for finding critical price Si
    while (abs(LHS-RHS)/X > E) {
        Si = (X+RHS-bi*Si)/(1-bi)
        d1 = (log(Si/X)+(b+sigma^2/2)*Time)/(sigma*sqrt(Time))
        LHS = Si-X
        RHS = GBSOption("c", Si, X, Time, r, b, sigma)@price + 
            (1-exp((b-r)*Time)*CND(d1))*Si/Q2
        bi = exp((b-r)*Time)*CND(d1)*(1-1/Q2) + 
        (   1-exp((b-r)*Time)*CND(d1)/(sigma*sqrt(Time)))/Q2 }
    
    # Return Value:
    Si
}


.BAWAmPutApproxOption <- 
function(S, X, Time, r, b, sigma) 
{
    # Internal Function - The Put:
    
    # Compute:
    Sk = .bawKp(X, Time, r, b, sigma)
    n = 2*b/sigma^2
    K = 2*r/(sigma^2*(1-exp(-r*Time)))
    d1 = (log(Sk/X)+(b+sigma^2/2)*Time)/(sigma*sqrt(Time))
    Q1 = (-(n-1)-sqrt((n-1)^2+4*K))/2
    a1 = -(Sk/Q1)*(1-exp((b-r)*Time)*CND(-d1))
    if(S > Sk) {
        result = GBSOption("p", S, X, Time, r, b, sigma)@price + a1*(S/Sk)^Q1 
    } else {
        result = X-S 
    }  
    
    # Return Value:
    result
}


.bawKp <- 
function(X, Time, r, b, sigma) 
{   
    # Internal Function - used for the Put:
    
    # Newton Raphson algorithm to solve for the critical commodity 
    # price for a Put.
    # Calculation of seed value, Si
    n = 2*b/sigma^2
    m = 2*r/sigma^2
    q1u = (-(n-1)-sqrt((n-1)^2+4*m))/2
    Su = X/(1-1/q1u)
    h1 = (b*Time-2*sigma*sqrt(Time))*X/(X-Su)
    Si = Su+(X-Su)*exp(h1) 
    K = 2*r/(sigma^2*(1-exp(-r*Time)))
    d1 = (log(Si/X)+(b+sigma^2/2)*Time)/(sigma*sqrt(Time))
    Q1 = (-(n-1)-sqrt((n-1)^2+4*K))/2
    LHS = X-Si
    RHS = GBSOption("p", Si, X, Time, r, b, sigma)@price -
        (1-exp((b-r)*Time)*CND(-d1))*Si/Q1
    bi = -exp((b-r)*Time)*CND(-d1)*(1-1/Q1) -
        (1+exp((b-r)*Time)*CND(-d1)/(sigma*sqrt(Time)))/Q1
    E = 0.000001
    # Newton Raphson algorithm for finding critical price Si
    while (abs(LHS-RHS)/X > E ) {
        Si = (X-RHS+bi*Si)/(1+bi)
        d1 = (log(Si/X)+(b+sigma^2/2)*Time)/(sigma*sqrt(Time))
        LHS = X-Si
        RHS = GBSOption("p", Si, X, Time, r, b, sigma)@price -
            (1-exp((b-r)*Time)*CND(-d1))*Si/Q1
        bi = -exp((b-r)*Time)*CND(-d1)*(1-1/Q1) -
            (1+exp((b-r)*Time)*CND(-d1)/(sigma*sqrt(Time)))/Q1 }
    # Return Value:
    Si
}


# ------------------------------------------------------------------------------


BSAmericanApproxOption = 
function(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma, title = NULL, 
description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the option price of an American call or 
    #   put option stocks, futures, and currencies. The 
    #   approximation method by Bjerksund and Stensland is used.
    
    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas
    
    # FUNCTION:
    
    # Settings:
    TypeFlag = TypeFlag[1]
    
    # The Bjerksund and Stensland (1993) American approximation:
    if(TypeFlag == "c") {
      result = .BSAmericanCallApprox(S, X, Time, r, b, sigma) }
    if(TypeFlag == "p") {
      # Use the Bjerksund and Stensland put-call transformation
      result = .BSAmericanCallApprox(X, S, Time, r - b, -b, sigma) }
    
    # Parameters:
    # TypeFlag = c("c", "p"), S, X, Time, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma
    if(!is.na(result$TriggerPrice)) param$TrigerPrice = result$TriggerPrice 
    
    # Add title and description:
    if(is.null(title)) title = "BS American Approximated Option"
    if(is.null(description)) description = as.character(date())
    
    # Return Value:
    new("fOPTION", 
        call = match.call(),
        parameters = param,
        price = result$Premium, 
        title = title,
        description = description
        )      
}


.BSAmericanCallApprox <- 
function(S, X, Time, r, b, sigma) 
{ 
    # Call Approximation:
    
    if(b >= r) { 
        # Never optimal to exersice before maturity
        result = list(
            Premium = GBSOption("c", S, X, Time, r, b, sigma)@price,
            TriggerPrice = NA)
    } else {
    Beta = (1/2 - b/sigma^2) + sqrt((b/sigma^2 - 1/2)^2 + 2*r/sigma^2)
    BInfinity = Beta/(Beta-1) * X
    B0 = max(X, r/(r-b) * X)
    ht = -(b*Time + 2*sigma*sqrt(Time)) * B0/(BInfinity-B0)
    # Trigger Price I:
    I = B0 + (BInfinity-B0) * (1 - exp(ht))
    alpha = (I-X) * I^(-Beta)
    if(S >= I) { 
        result = list(
            Premium = S-X, 
            TriggerPrice = I) }
    else {
        result = list(
            Premium = alpha*S^Beta - alpha*.bsPhi(S,Time,Beta,I,I,r,b,sigma) + 
            .bsPhi(S,Time,1,I,I,r,b,sigma) - .bsPhi(S,Time,1,X,I,r,b,sigma) - 
            X*.bsPhi(S,Time,0,I,I,r,b,sigma) + X*.bsPhi(S,Time,0,X,I,r,b,sigma), 
            TriggerPrice = I) } }
    result}
      

.bsPhi <- 
function(S, Time, gamma, H, I, r, b, sigma) 
{
    # Utility function phi:

    lambda = (-r + gamma*b + 0.5*gamma * (gamma-1)*sigma^2) * Time
    d = -(log(S/H) + (b + (gamma-0.5)*sigma^2)*Time) / 
        (sigma*sqrt(Time))
    kappa = 2 * b / (sigma^2) + (2*gamma - 1)
    result = exp(lambda)*S^gamma * 
    (CND(d)-(I/S)^kappa*CND(d-2*log(I/S)/(sigma*sqrt(Time))))
    
    # Return Value:
    result 
}


################################################################################

