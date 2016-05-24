if (getRversion() >= "2.15.1") utils::globalVariables(c("r0", "u","d", "q", "nsteps", "ratesteps"))



ratetreegui <-
  function(){
    
    my.draw <- function(panel){
      
      
      r0=as.numeric(panel$r0)
      u=as.numeric(panel$u)
      d=as.numeric(panel$d)
      q=as.numeric(panel$q)
      nsteps=as.numeric(panel$ratesteps)
      ratesteps=as.numeric(panel$ratesteps)
      
      
      
      
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
      
      
      
      # set graphs options
      if (nsteps>= 2){
        cex=0.9
      }
      else{
        cex=1
      }
      
      
      topaste = "Rate"
      M = R  
      nrows = dim(R)[1]
      ncols = dim(R)[2]
      
      
      if (length(dev.list()) == 0) 
        dev.new()
      plot(1:nrows, 1:ncols, type="n",ylab="",xlab="", 
           axes=FALSE, frame = FALSE)
      
      for (i in 1:nrows){
        for (j in 1:ncols){
          text(i, j, round(M[i,j],2),cex=cex) # ,col="red")  
        }
      }
      title(main = paste(floor(nrows/2),"Step ", topaste, " Tree"))
      panel
    }   
    
    
    
    
    my.panel <- rp.control(title = "Rate Tree")
    
    rp.textentry(panel=my.panel,variable=r0,action=my.draw,title="Rate              ",initval=5.0)
    rp.textentry(panel=my.panel,variable=u,action=my.draw,title="u                   ",initval=1.1)
    rp.textentry(panel=my.panel,variable=d,action=my.draw,title="d                   ",initval=0.9)
    rp.textentry(panel=my.panel,variable=q,action=my.draw,title="q                   ",initval=0.5)
    rp.doublebutton(panel = my.panel, showvalue=TRUE, variable= ratesteps, step = 1, range = c(1, 16),initval=10,
                    title = "Rate tree steps", action = my.draw)
    rp.do(my.panel, my.draw)
  }






