## this file contains the function originally from chapter 7. simulations.
## the function generate.sims( 1, blah, blah, blah ) returns a data.frame
## that contains various metrics associated with forest projections.
## volumes are computed using the method described simulations.rnw (section

## this file contains functions for
## the following functions are documented 
## generate.sims.Rd


## this function returns a data.frame object that represents
## the results of combined simulations.


generate.sims <- function( s,
                          ch.site=125.0, ch.pnba=1.0,
                          ##t.treats=c(200),
                          s.treats=c(0,0.25,0.50,0.75,1.0),
                          ##n.rep=1,
                          age.start=10, age.end=50)
{
  
  res <- NULL
  rcon.res <- NULL
  chambers.res <- NULL
  rcch.res <- NULL ## combined rconifers->25->chambers
  
  results <- data.frame()
  
    ## for each level of shrub removal
    for( j in 1:(length(s.treats)) ) {
      
      shrub.rem <- s.treats[j]
      
      ##        msg <- sprintf( "running t.treats[%d]=%f, s.treats[%d]=%f rep %d", i, t.treats[i], j, s.treats[j], k )
      msg <- sprintf( "running s.treats[%d]=%f", j, s.treats[j] )
      print( msg )
      
################################################################################
      ## 1) run the rconifers model first
      print( "running rconifers" )
      s3 <- s
      
      ## remove the competing vegetation
      s3$plants[s3$plants$sp.code == "CV",]$expf <-
        s3$plants[s3$plants$sp.code == "CV",]$expf * ( 1.0 - shrub.rem )        
      ## temporarily assign a sample.data object for annual simulation
      ## the stand is currently 3 years old
      s0 <- s3
      
      ## project the stand to age 153
      for( m in 1:200 ) {
        ## project sample.3 in one year intervals
        s1 <- project( s0, 1, control=list(rand.err=0,
                                rand.seed=0,
                                endemic.mort=1,
                                sdi.mort=1) )                  
        ## compute summaries for the projection and add it to the results data.frame
        ##rep.res <- data.frame( veg.ctrl=i, pct.crtl=j, rep=k, age=s1$age, sp.sums.2( s1 )["DF",], col="grey" )
        ##rep.res <- data.frame( veg.ctrl=i, pct.ctrl=j, rep=k, age=s1$age, sp.sums.2( s1 )["DF",], col="rconifers" )
        ##rep.res <- data.frame( veg.ctrl=s.treats[j], pct.ctrl=t.treats[i], rep=k, age=s1$age, sp.sums.2( s1 )["DF",], col="rconifers" )
        rep.res <- data.frame( veg.ctrl=s.treats[j], age=s1$age, sp.sums.2( s1 )["DF",], col="rconifers" )
        print( rep.res )
        rcon.res <- rbind( rcon.res, rep.res )       
        s0 <- s1      
      }
      
      print( "done with the rconifers runs" )
        
################################################################################
      ## 2) run the chambers model second
      ## the projected stand is 30 years old and should have the nescessary 
      ## data to use to feed into the chambers-1980 model
      ##kings.site <- s1$site.index ## this needs to be converted to kings 1966 from smc @ 30        
      ## for each of the ages in the projection, you need to perform the following tasks:
      ## 1) construct a dummy sample.data object 
      ## 2) compute the volumes using sp.sums.2
      ## 3) append the resulting data.frame onto the results (res)
      print( "running chambers" )
      osr <- chambers.1980( ages=1:200, site=ch.site, pnba=ch.pnba )
      for( m in 10:200 ) {          
        s2.stem <- data.frame( plot=1,
                              sp.code="DF", 
                              d6=NA,
                              dbh=osr[m,]$qmd, 
                              tht=osr[m,]$tht,
                              cr=0.6,
                              n.stems=1,
                              expf=osr[m,]$expf,
                              crown.width=NA )
        s2.plot <- data.frame( plot=1,elevation=NA,slope=NA,
                              aspect=NA,whc=NA,map=NA,si30=ch.site )
        s2 <- list( plots=s2.plot, plants=s2.stem, age=m, x0=s$x0 )
        class( s2 ) <- "sample.data"          
        ## you need the volume here.
        ##rep.res <- data.frame( veg.ctrl=i, pct.crtl=j, rep=k, age=m,
        ##                      sp.sums.2( s2 )["DF",], col="black" )
        ##rep.res <- data.frame( veg.ctrl=s.treats[j], pct.ctrl=t.treats[i], rep=k, age=s2$age, sp.sums.2( s2 )["DF",], col="chambers" )
        ##rep.res <- data.frame( veg.ctrl=i, pct.ctrl=j, rep=k, age=s2$age, sp.sums.2( s2 )["DF",], col="chambers" )
        ## rbind the results to the larger results table.
        rep.res <- data.frame( veg.ctrl=s.treats[j], age=s2$age, sp.sums.2( s2 )["DF",], col="chambers" )
        print( rep.res )
        chambers.res <- rbind( chambers.res, rep.res )                 
      }
      print( "running chambers complete" )
      
      
      

################################################################################
        ## 3) run the rconifers model, and feed the results at age 25 into chambers
        ##s3 <- sample.3       
##       print( "running chambers from rconifers @ age 25" )
##         s3 <- s
##         ## remove the competing vegetation
##         s3$plants[s3$plants$sp.code == "CV",]$expf <-
##           s3$plants[s3$plants$sp.code == "CV",]$expf * ( 1.0 - shrub.rem )        
##         ## temporarily assign a sample.data object for annual simulation
##         ## the stand is currently 3 years old
##         s0 <- s3        
##         ## project the stand to age 25 (3+22)
##         for( m in 1:22 ) {          
##           ## project sample.3 in one year intervals
##           s1 <- project( s0, 1, control=list(rand.err=0,
##                                   rand.seed=0,
##                                   endemic.mort=1,
##                                   sdi.mort=1) )
##           ## compute summaries for the projection and add it to the results data.frame
##           ##rep.res <- data.frame( veg.ctrl=i, pct.crtl=j, rep=k, age=s1$age, sp.sums.2( s1 )["DF",], col="grey" )
##           ##rep.res <- data.frame( veg.ctrl=s.treats[j], pct.ctrl=t.treats[i], rep=k, age=s1$age, sp.sums.2( s1 )["DF",], col="rconifers->chambers" )
##           ##rep.res <- data.frame( veg.ctrl=s.treats[j], pct.ctrl=t.treats[i], rep=k, age=s1$age, sp.sums.2( s1 )["DF",], col="rconifers->chambers" )
##           ##rep.res <- data.frame( veg.ctrl=i, pct.ctrl=j, rep=k, age=s1$age, sp.sums.2( s1 )["DF",], col="rconifers->chambers" )
##           rep.res <- data.frame( veg.ctrl=s.treats[j], age=s1$age, sp.sums.2( s1 )["DF",], col="rconifers->chambers" )
##           print( rep.res )
##           rcch.res <- rbind( rcch.res, rep.res )       
##           s0 <- s1      
##         }
##         ## compute the percent normal basal area to hand-off to the chambers model
##         rcch.pnba <- sp.sums.2( s1 )["DF",]$ba / chambers.1980.nba( ch.site, s1$age )
##         osr <- chambers.1980( ages=1:200, site=ch.site, pnba=rcch.pnba )
##         for( m in 26:200 ) {
          
##           s2.stem <- data.frame( plot=1,
##                                 sp.code="DF", 
##                                 d6=NA,
##                                 dbh=osr[m,]$qmd, 
##                                 tht=osr[m,]$tht,
##                                 cr=0.6,
##                                 n.stems=1,
##                                 expf=osr[m,]$expf,
##                                 crown.width=NA )
##           s2.plot <- data.frame( plot=1,elevation=NA,slope=NA,
##                                 aspect=NA,whc=NA,map=NA,si30=site )
##           s2 <- list( plots=s2.plot, plants=s2.stem, age=m, x0=s1$x0 )
##           class( s2 ) <- "sample.data"
          
##           ## compute the summaries for the sample.data object.
##           ##rep.res <- data.frame( veg.ctrl=i, pct.crtl=j, rep=k, age=m,
##           ##                      sp.sums.2( s2 )["DF",], col="black" )
## ##           rep.res <- data.frame( veg.ctrl=s.treats[j], pct.ctrl=t.treats[i], rep=k, age=s2$age,
## ##                                 sp.sums.2( s2 )["DF",], col="rconifers->chambers" )
##           ##rep.res <- data.frame( veg.ctrl=i, pct.ctrl=j, rep=k, age=s2$age, sp.sums.2( s2 )["DF",], col="rconifers->chambers" )
##           rep.res <- data.frame( veg.ctrl=s.treats[j], age=s2$age,
##                                 sp.sums.2( s2 )["DF",], col="rconifers->chambers" )
##           print( rep.res )

##           ## rbind the results to the larger results table.
##           rcch.res <- rbind( rcch.res, rep.res )                 
##         }
        
##      } ## done with the replications        
      ## for the current level of shrub removal    
    } ## done with the different shrub removal levels
##  } ## done with the different pct thinning levels
  
  ## you now have chambers results

  rcon.res$mai <- rcon.res$sm.vol / rcon.res$age
  chambers.res$mai <- chambers.res$sm.vol / chambers.res$age
##  rcch.res$mai <- rcch.res$sm.vol / rcch.res$age
  
  ## only return the data.frames for which you request data
  res$rconifers <- rcon.res[rcon.res$age >= age.start & rcon.res$age <= age.end,]
  res$chambers <- chambers.res[chambers.res$age >= age.start & chambers.res$age <= age.end,]
  ##res$rcch <- rcch.res[rcch.res$age >= age.start & rcch.res$age <= age.end,]

  ## only return the results where the two are the same length.
  ## for each of the runs, go back and 
  
  res
}




