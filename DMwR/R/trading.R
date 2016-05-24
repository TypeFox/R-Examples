# ===============================================
# TRADING SIMULATOR
# -----------------------------------------------
# This function simulates the trading given the market
# prices for a future period and a set of model signals
# The function calls a user-suplied function (policy.func)
# that makes the decisions at the end of each trading day
# -----------------------------------------------
# L. Torgo, Nov 2009
# -----------------------------------------------
trading.simulator <- function(market,
                              signals,
                              policy.func,
                              policy.pars=list(),
                              trans.cost=5,
                              init.cap=1000000)
  {

    #-----------------
    # Aux Functions
    # ----------------
    
    # ---- Opening a position
    open.position <- function(type,prc,quant) {
      #quant <- round(bet*(trading[d,'Money']-trans.cost)/prc,0)
      .currPos <<- .currPos + 1
      if (.currPos > .maxPos) {
        n <- .maxPos %/% 4
        .maxPos <<- .maxPos + n
        positions <<- rbind(positions,
                            matrix(NA,nrow=n,ncol=7,
                                   dimnames=list(.currPos:.maxPos,
                                     c('pos.type','N.stocks','Odate',
                                       'Oprice','Cdate','Cprice','result')))
                            )
      }
      positions[.currPos,] <<- c(pos.type=type,
                                 N.stocks = quant,
                                 Odate=d,Oprice=prc,
                                 Cdate=NA,Cprice=NA,
                                 result=NA)
      trading[d,'Order'] <<- type
      trading[d,'N.Stocks'] <<- trading[d,'N.Stocks']+type*quant
      trading[d,'Money'] <<- trading[d,'Money'] -
                                   type*quant*prc -
                                   trans.cost
      if (trading[d,'Money'] < 0)
        cat('Borrowing money (',abs(trading[d,'Money']),
            ') for opening a long position (PosID=',.currPos,')\n')
      open.positions <<- c(open.positions,.currPos)
      return(.currPos)
    }
    # ---- Closing a position
    close.position <- function(ID,prc) {
      #browser()
      quant <- positions[ID,'N.stocks']
      value <- positions[ID,'pos.type']*quant*prc
               - trans.cost
      trading[d,'Order'] <<- -positions[ID,'pos.type']
      trading[d,'Money'] <<- trading[d,'Money'] +
                             positions[ID,'pos.type']*quant*prc -
                             trans.cost
      if (trading[d,'Money'] < 0)
        cat('Borrowing money (',abs(trading[d,'Money']),
            ') for closing a short position (PosID=',ID,')\n')
      trading[d,'N.Stocks'] <<- trading[d,'N.Stocks']-
          positions[ID,'pos.type']*quant
      positions[ID,'Cdate'] <<- d
      positions[ID,'Cprice'] <<- prc
      init <- if (positions[ID,'pos.type'] == 1) positions[ID,'Oprice'] else positions[ID,'Cprice']
      fin <- if (positions[ID,'pos.type'] == 1) positions[ID,'Cprice'] else positions[ID,'Oprice']
      positions[ID,'result'] <<- 100*(fin/init - 1)
      open.positions <<- open.positions[-which(open.positions == ID)]

    }
    
    # ----------------
    # Initialization stuff
    # ----------------
    dates <- index(market)
    market <- as.data.frame(market)
    N.days <- nrow(market)

    res <- list()
    
    trading <- matrix(0,nrow=N.days,ncol=5)
    colnames(trading) <- c('Close','Order','Money','N.Stocks','Equity')
    trading[,'Close'] <- market$Close
    trading[1,'Money'] <- init.cap

    .maxPos <- N.days %/% 2
    positions <- matrix(NA,nrow=.maxPos,ncol=7,
                        dimnames=list(1:.maxPos,
                          c('pos.type','N.stocks','Odate',
                            'Oprice','Cdate','Cprice','result')))
    .currPos <- 0
    open.positions <- c() # positions currently opened


    pending.orders <- NULL
    .orderID <- 1
    
    # -----------------
    # Main Loop through all days
    # -----------------

    for(d in 1:N.days) {
      
      # ------------------------
      # Pre-Open market actions
      #if (d>=43) browser()
      # update money and n.stocks
      if (d > 1) {
        trading[d,'Money'] <- trading[d-1,'Money']
        trading[d,'N.Stocks'] <- trading[d-1,'N.Stocks']
      }
      

      # ------------------------
      # During the day actions
      
      # check pending orders
      if (NROW(pending.orders)) {
        #browser()
        # * Market orders
        mkts <- which(pending.orders$order.type == 1) # market orders are first
        closed <- c()
        #if (d >= 488) browser()
        for(i in mkts) {
          # - opening a position
          if (pending.orders[i,'action'] == 'open') {
            idP <- open.position(pending.orders[i,'order'],
                                market[d,'Open'],
                                pending.orders[i,'val'])
            # associate the respective limit and stop with this opened position
            pending.orders[pending.orders$ID==pending.orders[i,'ID'],
                           'posID'] <- idP
            
          # - closing a position  
          } else { 
            close.position(pending.orders[i,'posID'],market[d,'Open'])
            # closed positions
            closed <- c(closed,pending.orders[i,'posID'])
          }
        }
        if (length(mkts)) {
          pending.orders <- pending.orders[-mkts,]
          toRem <- which(pending.orders$posID %in% closed)
          if (length(toRem)) pending.orders <- pending.orders[-toRem,]
        }
        
        if (NROW(pending.orders)) {
          done <- c()
          for(i in 1:NROW(pending.orders)) {
            if (! i %in% done) {
              # * Limit orders (always for closing positions)
              if (pending.orders[i,'order.type'] == 2) {

                if (pending.orders[i,'order'] == 1) { # it's a buy to close a short
                  if (market[d,'Low'] < pending.orders[i,'val']) {
                    close.position(pending.orders[i,'posID'],
                                   pending.orders[i,'val'])
                    done <- c(done,which(pending.orders$posID == pending.orders[i,'posID']))
                  }
                  
                } else { # it's a sell to close a long
                  if (market[d,'High'] > pending.orders[i,'val']) {
                    close.position(pending.orders[i,'posID'],
                                   pending.orders[i,'val'])
                    done <- c(done,which(pending.orders$posID == pending.orders[i,'posID']))
                  }
                }
                
              # * Stop orders (always for closing positions)
              } else if (pending.orders[i,'order.type'] == 3) {

                if (pending.orders[i,'order'] == 1) { # it's a buy to close a short
                  if (market[d,'High'] > pending.orders[i,'val']) {
                    close.position(pending.orders[i,'posID'],
                                   pending.orders[i,'val'])
                    done <- c(done,which(pending.orders$posID == pending.orders[i,'posID']))
                  }
                } else { # it's a sell to close a long
                  if (market[d,'Low'] < pending.orders[i,'val']) {
                    close.position(pending.orders[i,'posID'],
                                   pending.orders[i,'val'])
                    done <- c(done,which(pending.orders$posID == pending.orders[i,'posID']))
                  }
                }

              }
            }
          }
          
          if (length(done)) pending.orders <- pending.orders[-done,]
        }
        
      }

      # ------------------------
      # Post-Close market actions

      # call user-define decision policy function
      # orders is a data.frame with (order,type,val,action,posID)

      orders <- do.call(policy.func,
                        c(list(signals[1:d],
                               market[1:d,],
                               if (.currPos) positions[open.positions,1:4,drop=F] else NULL,
                               trading[d,'Money']),
                          policy.pars
                          ))
      if (no <- NROW(orders)) orders <- cbind(ID=rep(.orderID,no),orders)
      .orderID <- .orderID + NROW(orders)
      pending.orders <- rbind(pending.orders,orders)
      
      #if (NROW(pending.orders)) {
      #  cat('\nEnd of day ',d,'\n')
      #  print(pending.orders)
      #  print(tail(positions))
      #}

      # Update our current Equity
      trading[d,'Equity'] <- trading[d,'Money'] +
        trading[d,'N.Stocks']*market[d,'Close']
      
    }


#    trading <- xts(trading,dates)
    trading <- zoo(trading,dates)
    tradeRecord(trading,if (.currPos) positions[1:.currPos,,drop=FALSE] else as.matrix(vector()),trans.cost,init.cap,policy.func,policy.pars)
  }




# =====================================================
# Function that obtains the precision and recall of the
# trading signals of a model, given the "true" signals.
# It produces a matrix with three rows ('s','b','s+b'),
# and two columns ('precision', 'recall'), with the
# respective scores.
# =====================================================
# Luis Torgo, Nov 2009
# =====================================================
sigs.PR <- function(preds,trues) {
  confM <- table(preds,trues)
  sigs <- c('s','b')
  r <- matrix(NA,ncol=2,nrow=3,dimnames=list(c(sigs,'s+b'),c('precision','recall')))
  for(i in seq(along=sigs)) {
    r[i,'precision'] <- confM[sigs[i],sigs[i]]/sum(confM[sigs[i],])
    r[i,'recall'] <- confM[sigs[i],sigs[i]]/sum(confM[,sigs[i]])
  }
  r['s+b','precision'] <- (confM[1,1]+confM[3,3])/sum(confM[c(1,3),])
  r['s+b','recall'] <- (confM[1,1]+confM[3,3])/sum(confM[,c(1,3)])
  r
}


# =====================================================
# Function that obtains the trading signals corresponding
# to a set of numeric values.
# The user specifies the buy and sell trheshold that guide
# the discretization
# =====================================================
# Luis Torgo, Nov 2009
# =====================================================
trading.signals <- function(vs,b.t,s.t) 
  factor(ifelse(vs < s.t,'s',ifelse(vs > b.t,'b','h')), levels=c('s','h','b'))



# =====================================================
# Function that obtains a vector of econometric evaluation
# statistics of a trading record.
# The trading record of class tradeRecord that is generated
# by the trading.simulator function
# =====================================================
# Luis Torgo, Nov 2009
# =====================================================
tradingEvaluation <- function(t) {
  if (!inherits(t,'tradeRecord')) stop(t,' is not of class "tradeRecord".\n')
  s <- c(NTrades=0,NProf=0,PercProf=0,
         PL=0,Ret=0,RetOverBH=NA,
         MaxDD=NA,SharpeRatio=NA,
         AvgProf=NA,AvgLoss=NA,AvgPL=NA,MaxProf=NA,MaxLoss=NA)

  s['NTrades'] <- NROW(t@positions)
  if (s['NTrades']) {

    s['NProf'] <- length(which(t@positions[,'result'] > 0))
  
    s['PercProf'] <- 100*s['NProf']/s['NTrades']

    s['AvgProf'] <- sum(t@positions[t@positions[,'result'] > 0,'result'],
                        na.rm=T)/s['NProf']
    s['AvgLoss'] <- sum(t@positions[t@positions[,'result'] < 0,'result'],
                        na.rm=T)/(s['NTrades']-s['NProf'])
    s['AvgPL'] <- sum(t@positions[,'result'],na.rm=T)/s['NTrades']
    
    s['MaxProf'] <- max(0,max(t@positions[,'result'],na.rm=T))
    s['MaxLoss'] <- min(0,min(t@positions[,'result'],na.rm=T))
  }
  
  s['PL'] <- t@trading[nrow(t@trading),'Equity']-t@init.cap
  s['Ret'] <- 100*(t@trading[nrow(t@trading),'Equity']-t@init.cap)/t@init.cap
  s['RetOverBH'] <- s['Ret'] - 100*(coredata(t@trading[nrow(t@trading),'Close'])-t@trading[1,'Close'])/t@trading[1,'Close']
 
  s['MaxDD'] <- max(cummax(t@trading[,'Equity'])-t@trading[,'Equity'])
  
    
  s['SharpeRatio'] <- mean(100*Delt(t@trading[,'Equity']),na.rm=T)/
                       sd(100*Delt(t@trading[,'Equity']),na.rm=T)

  round(s,2)
}


.Eq <- function(p) p[,'Equity']
.St <- function(p) p[,'N.Stocks']
.addEq <- newTA(FUN = .Eq, col = 'red', legend = "Equity")
.addSt <- newTA(FUN = .St, col = 'green', legend = "N.Stocks")
