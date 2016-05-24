ibdhap.compare.segs <- function( calls, true, data.type, seg.cutoff, pos = NA){
### Compare called states and true states by segment of IBD

	## A little function for the mode
	Mode = function(x){
		ux = unique(x)
		ux[which.max(tabulate(match(x,ux)))]	
	}			

    ## Check that data type has been entered 
    if(length(data.type)>1){ stop("data.type improperly indicated")
    }else if( is.element("h", data.type)){ no.ibd.ind = 15
    }else if( is.element("g", data.type)){no.ibd.ind = 9
    }else if( is.element("r", data.type)){no.ibd.ind = 4
    }else{ stop("data.type improperly indicated")}

    ## Set up parameters
    num.sets <- dim( calls )[2]             
    num.snps <- dim( calls )[1]

    if( is.na(pos)) pos <- 1:num.snps
    seg.info <- NULL
    
    for( j in 1:num.sets ){ ## for each set of ibd states
        
        callsj <- calls[,j]
        truej <- true[,j]
        startpos <- c( 0, which(diff(truej)!=0) )+1  ## starts of each segment
        endpos <- c( tail( startpos,-1)-1, num.snps ) ## ends of each segment
        
        for( i in 1:length( startpos )){

            ## true value in segment
            true.state <- truej[ startpos[i] ] 
            ## called segment value
            mode.call <- Mode( callsj[ startpos[i]:endpos[i] ] )
            prop.corr <- mean( mode.call == callsj[ startpos[i]:endpos[i] ] )
            if( max( prop.corr, 1-prop.corr) > seg.cutoff ){
                seg.call <- mode.call 
            }else{
                seg.call <- 0 
            }
            ## seg length 
            seg.length <- pos[ endpos[i] ] - pos[ startpos[i] ]
            
            ## save
            seg.info <- rbind( seg.info,
                              c(seg.length, true.state, seg.call, mode.call, prop.corr ))
        }
    }
    
    ## Stats for all the segments
    colnames( seg.info ) <- c("seg.length","true.state","seg.call","mode.call","prop.corr" ) 
    
    seg.stats <- matrix( NA, ncol = 1, nrow = 6 )
    rownames( seg.stats ) <- c("Number segments",
                               "Number of IBD segments",
                               "IBD Segs called correctly",
                               "IBD Segs called no-IBD",
                               "IBD Segs called wrong IBD",
                               "IBD Segs with no call")
    
    
    ibd <- seg.info[,2]!=no.ibd.ind  ## Subset IBD segments
    seg.t <- seg.info[ibd,2]
    seg.c <- seg.info[ibd,3]
    
    seg.stats[1,] <- dim( seg.info )[1]
    seg.stats[2,] <- length( seg.t ) 
    seg.stats[3,] <- mean( seg.t == seg.c )*100
    seg.stats[4,] <- mean( seg.t != seg.c & seg.c == no.ibd.ind )*100
    seg.stats[5,] <- mean( seg.t != seg.c & seg.c < no.ibd.ind & seg.c > 0 )*100
    seg.stats[6,] <- mean( seg.t != seg.c & seg.c == 0 )*100
   
    return( list( seg.stats = seg.stats, seg.info = seg.info ) )
}
