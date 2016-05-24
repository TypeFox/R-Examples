if (getRversion() >= "2.15.1") utils::globalVariables(c("s0","strike","u","d", "q", "maturity",  "Rate", "nsteps", "exercisetype","opttype"))


stockoptiontreegui <-
  function(){
    
    my.draw <- function(panel){
      
      s0=as.numeric(panel$s0)
      u=as.numeric(panel$u)
      d=as.numeric(panel$d)
      q=as.numeric(panel$q)
      maturity=as.numeric(panel$nsteps)
      Rate=as.numeric(panel$Rate)
      strike=as.numeric(panel$strike)
      opttype=panel$opttype
      exercisetype=panel$exercisetype
      nsteps=as.numeric(panel$nsteps)
      
      
      
      stocktree <- function(s0,u,d,q,nsteps,optional){
        
        if(missing(optional)) {
          timepoints = nsteps + 1
          gap = 1 # gap between rows
          margin = 0
          nrows = (gap+1) * nsteps + 1 + 2 * margin
          ncols = 2 * (nsteps) + 1 + 2 * margin
          dt = 1
          startrow = margin + 1
          startcol = margin + 1
          lastrow = nrows - margin
        }
        
        #optional=c(timepoints,gap,margin,nrows,ncols,dt,startrow,startcol,lastrow)
        
        else{
          timepoints=optional[1]
          gap=optional[2]
          margin=optional[3]
          nrows=optional[4]
          ncols=optional[5]
          dt=optional[6]
          startrow=optional[7]
          startcol=optional[8]
          lastrow=optional[9]
        }
        
        S = matrix(data = NA, nrow = nrows, ncol = ncols)
        S[startrow, nsteps + margin + 1] = s0
        
        thisstep= 0
        
        for (row in seq(from = startrow + gap+1, to = nrows - margin, by = gap+1)) {
          thisstep = thisstep+1
          for (col in seq(from = ncols-nsteps-margin-thisstep, to=ncols-nsteps-margin+thisstep, by = 2)) {
            if (col <= ceiling(ncols/2)) {
              S[row, col] = S[row - (gap+1), col + 1] * d
            }
            else {
              S[row, col] = S[row - (gap+1), col - 1] * u
            }
          }
        }
        #S=round(S,2)
        S
      }
      
      #stocktree(s0=100,u=1.25,d=0.8,q=0.6,nsteps=3)
      
      
      
      
      
      
      stockoptiontree <- function(strike,opttype,exercisetype,s0,u,d,q,maturity,Rate){
        
        
        nsteps = maturity
        timepoints = nsteps + 1
        gap = 1 # gap between rows
        margin = 0
        nrows = (gap+1) * nsteps + 1 + 2 * margin
        ncols = 2 * (nsteps) + 1 + 2 * margin
        dt = 1
        startrow = margin + 1
        startcol = margin + 1
        lastrow = nrows - margin
        
        optional=c(timepoints,gap,margin,nrows,ncols,dt,startrow,startcol,lastrow)
        
        
        P <- stocktree(s0,u,d,q,nsteps,optional)
        
        
        
        O = matrix(data = NA, nrow = nrows, ncol = ncols)
        
        # final nodes
        thisstep=nsteps   
        #for (col in seq(from = (startcol + timepoints - lastrow/2), by = 2, length.out = lastrow/2)) {
        for (col in seq(from = ncols-nsteps-margin-thisstep, to=ncols-nsteps-margin+thisstep, by = 2)) {
          if (opttype == "Call"){
            O[lastrow, col] = pmax(P[lastrow, col] - strike, 0)
          }
          else{
            O[lastrow, col] = pmax(strike- P[lastrow, col], 0)
          }
          
        }
        # intermediate nodes
        thisstep=nsteps-1       
        for (row in seq(from = (lastrow - (gap+1)), to = startrow, by = -(gap+1))) {
          #for (col in seq(from = (startcol + timepoints - row/2), by = 2, length.out = row/2)) {
          for (col in seq(from = ncols-nsteps-margin-thisstep, to=ncols-nsteps-margin+thisstep, by = 2)) {
            if (opttype == "Call"){
              O[row, col] =  (O[row + gap+1, col + 1] * q + O[row + gap+1, col - 1] * (1 - q))/(1+Rate/100)
            }
            else{
              O[row, col] =  (O[row + gap+1, col + 1] * q + O[row + gap+1, col - 1] * (1 - q))/(1+Rate/100)
            }
          }
          thisstep=thisstep-1
        }
        
        
        thisstep=nsteps-1          
        if (exercisetype == "American") {
          for (row in seq(from = (lastrow - (gap+1)), to = startrow, by = -(gap+1))) {
            #for (col in seq(from = (startcol + timepoints -  row/2), by = 2, length.out = row/2)) {
            for (col in seq(from = ncols-nsteps-margin-thisstep, to=ncols-nsteps-margin+thisstep, by = 2)) {
              if (opttype == "Call"){
                O[row, col] = pmax(pmax(P[row, col] - strike,  0), O[row, col])
              }
              else {
                O[row, col] = pmax(pmax(strike - P[row, col],  0), O[row, col])
              }
            }
            thisstep=thisstep-1
          }
        }
        
        else { # European
        }
        
        #O=round(O,2)
        O
      }
      
      
      
      #stockoptiontree(strike=110,opttype="Call",exercisetype="European",s0=100,u=1.25,d=0.8,q=0.6,maturity=3,Rate=5.0)
      
      
      
      
      
      
      
      S=stocktree(s0,u,d,q,nsteps)
      C=stockoptiontree(strike,opttype,exercisetype,s0,u,d,q,maturity,Rate)
      
      # set graphs options
      if (nsteps>= 2){
        cex=0.9
      }
      else{
        cex=1
      }
      
      
      nrows = dim(C)[1]
      ncols = dim(C)[2]
      
      if (length(dev.list()) == 0) 
        dev.new()
      plot(1:nrows, 1:ncols, type="n",ylab="",xlab="", 
           axes=FALSE, frame = FALSE)
      
      for (i in 1:nrows){
        for (j in 1:ncols){
          if (panel$plot == "Stock Tree"){
            topaste = "Stock"  
            text(i, j, round(S[i,j],2),cex=cex) # ,col="red") 
          }
          else{
            topaste = paste(panel$exercisetype, panel$opttype)
            text(i, j, round(C[i,j],2),cex=cex) # ,col="red")  
          }
        }
      }
      title(main = paste(nsteps,"Step ", topaste, " Tree"))
      panel
    }   
    
    
    
    
    my.panel <- rp.control(title = "Stock Option Tree")
    
    rp.radiogroup(panel = my.panel, variable= opttype,
                  vals = c("Call", "Put"), 
                  action = my.draw, title = "Type of Option")
    
    rp.radiogroup(panel = my.panel, variable= exercisetype,
                  vals = c("European", "American"), 
                  action = my.draw, title = "Exercise style")
    
    rp.textentry(panel=my.panel,variable=s0,action=my.draw,title="Stock price     ",initval=100)
    rp.textentry(panel=my.panel,variable=strike,action=my.draw,title="Strike price     ",initval=110)
    #rp.textentry(panel=my.panel,variable=Time,action=my.draw,title="Time               ",initval=0.25)
    #rp.textentry(panel=my.panel,variable=sigma,action=my.draw,title="Volatility         ",initval=0.25)
    rp.textentry(panel=my.panel,variable=Rate,action=my.draw,title="Risk free rate  ",initval=4)
    rp.textentry(panel=my.panel,variable=u,action=my.draw,title="u                   ",initval=1.1)
    rp.textentry(panel=my.panel,variable=d,action=my.draw,title="d                   ",initval=0.9)
    rp.textentry(panel=my.panel,variable=q,action=my.draw,title="q                   ",initval=0.5)
    #rp.textentry(panel=my.panel,variable=Div,action=my.draw,title="Dividend rate ",initval=0)
    rp.doublebutton(panel = my.panel, showvalue=TRUE, variable= nsteps, step = 1, range = c(1, 15),initval=3,
                    title = "No. of Steps", action = my.draw)
    rp.radiogroup(panel = my.panel, variable= plot,
                  vals = c("Stock Tree", "Option Tree"), 
                  action = my.draw, title = "Plot Type")
    rp.do(my.panel, my.draw)
  }


