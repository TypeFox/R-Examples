if (getRversion() >= "2.15.1") utils::globalVariables(c("FV","strike", "r0", "u","d", "q", "optmaturity","coupon", "nsteps", "ratesteps"))

floortreegui <-
  function(){
    
    my.draw <- function(panel){
      
      FV=as.numeric(panel$FV)
      strike=as.numeric(panel$strike)
      r0=as.numeric(panel$r0)
      u=as.numeric(panel$u)
      d=as.numeric(panel$d)
      q=as.numeric(panel$q)
      optmaturity=as.numeric(panel$optmaturity)
      nsteps=as.numeric(panel$optmaturity)
      ratesteps=as.numeric(panel$optmaturity)
      
      ratetree <- function(r0,u,d,q,ratesteps,optional){
        
        if(missing(optional)) {
          timepoints = ratesteps + 1
          gap = 1 # gap between rows
          margin = 0
          nrows = (gap+1) * ratesteps + 1 + 2 * margin
          ncols = 2 * (ratesteps) + 1 + 2 * margin
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
        
        Rates = matrix(data = NA, nrow = nrows, ncol = ncols)
        Rates[startrow, ratesteps + margin + 1] = r0
        
        thisstep= 0
        
        for (row in seq(from = startrow + gap+1, to = nrows - margin, by = gap+1)) {
          thisstep = thisstep+1
          for (col in seq(from = ncols-ratesteps-margin-thisstep, to=ncols-ratesteps-margin+thisstep, by = 2)) {
            if (col <= ceiling(ncols/2)) {
              Rates[row, col] = Rates[row - (gap+1), col + 1] * d
            }
            else {
              Rates[row, col] = Rates[row - (gap+1), col - 1] * u
            }
          }
        }
        #Rates=round(Rates,2)
        Rates
      }
      #ratetree(r0=6,u=1.25,d=0.9,q=0.5,ratesteps=4)
      
      
      
      
      R = ratetree(r0,u,d,q,ratesteps)
      
      floortree <- function(FV,strike,optmaturity,coupon,r0,u,d,q){
        
        nsteps = optmaturity-1
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
        ratesteps=optmaturity-1
        
        Rates <- ratetree(r0,u,d,q,ratesteps,optional)
        
        
        
        O = matrix(data = NA, nrow = nrows, ncol = ncols)
        
        # final nodes
        thisstep=nsteps   
        
        for (col in seq(from = ncols-nsteps-margin-thisstep, to=ncols-nsteps-margin+thisstep, by = 2)) {
          
          O[lastrow, col] = max(strike-Rates[lastrow, col],0)/(1+Rates[lastrow, col]/100)
          
        }
        
        # intermediate nodes
        thisstep=nsteps-1       
        for (row in seq(from = (lastrow - (gap+1)), to = startrow, by = -(gap+1))) {
          
          for (col in seq(from = ncols-nsteps-margin-thisstep, to=ncols-nsteps-margin+thisstep, by = 2)) {
            if(thisstep!=0){
              O[row, col] = (max(strike-Rates[row, col] ,0) + O[row + gap+1, col + 1] * q + O[row + gap+1, col - 1] * (1 - q)) /(1+Rates[row, col]/100)
            }
            else{
              O[row, col] = ( O[row + gap+1, col + 1] * q + O[row + gap+1, col - 1] * (1 - q)) /(1+Rates[row, col]/100)
            }
            
          }
          thisstep=thisstep-1
        }
        
        
        
        floorvalue= O[startrow,ceiling(ncols/2)]*FV
        floorvalue=round(floorvalue,2)
        
        ans = list(floorvalue,O)
        
        ans
        
      }
      
      
      
      
      C = floortree(FV,strike,optmaturity,coupon,r0,u,d,q)[[2]]
      val = floortree(FV,strike,optmaturity,coupon,r0,u,d,q)[[1]]
      
      
      
      
      
      # set graphs options
      if (nsteps>= 2){
        cex=0.9
      }
      else{
        cex=1
      }
      
      
      
      if (panel$plot == "Floor Tree"){
        topaste = "Floor" 
        M = C 
        nrows = dim(C)[1]
        ncols = dim(C)[2]
      }
      else{
        topaste = "Rate"
        M = R  
        nrows = dim(R)[1]
        ncols = dim(R)[2]
      }
      
      
      if (length(dev.list()) == 0) 
        dev.new()
      plot(1:nrows, 1:ncols, type="n",ylab="",xlab="", 
           axes=FALSE, frame = FALSE)
      
      for (i in 1:nrows){
        for (j in 1:ncols){
          text(i, j, round(M[i,j],3),cex=cex) # ,col="red")  
        }
      }
      title(main = paste(floor(nrows/2),"Step ", topaste, " Tree. Floor value = ", val))
      panel
    }   
    
    
    
    
    my.panel <- rp.control(title = "Floor Tree")
    
    rp.textentry(panel=my.panel,variable=FV,action=my.draw,title="Face value    ",initval=100)
    rp.textentry(panel=my.panel,variable=strike,action=my.draw,title="strike            ",initval=8.0)
    rp.textentry(panel=my.panel,variable=r0,action=my.draw,title="Rate              ",initval=6.0)
    rp.textentry(panel=my.panel,variable=u,action=my.draw,title="u                   ",initval=1.25)
    rp.textentry(panel=my.panel,variable=d,action=my.draw,title="d                   ",initval=0.9)
    rp.textentry(panel=my.panel,variable=q,action=my.draw,title="q                   ",initval=0.5)
    rp.textentry(panel=my.panel,variable=coupon,action=my.draw,title="Coupon        ",initval=5.0)
    rp.doublebutton(panel = my.panel, showvalue=TRUE, variable= optmaturity, step = 1, range = c(1, 15),initval=6,
                    title = "Floor Maturity", action = my.draw)
    rp.radiogroup(panel = my.panel, variable= plot,
                  vals = c("Floor Tree", "Rate Tree"), 
                  action = my.draw, title = "Plot Type")
    rp.do(my.panel, my.draw)
  }

