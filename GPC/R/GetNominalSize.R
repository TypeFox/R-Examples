GetNominalSize <- function(L,Growth,Rule){
    # L         = quadrature level
    # Growth    = 'Lin'
    #             'LinOdd'
    #             'LindEven'
    #             'ExpOdd'
    #             'ExpEven'
    # Rule      = 'Gauss
    #             'Fejer'
    #             'ClenshawCurtis'
    #             'GaussKP' or 'LegendreKP'
    #             'HertotalteKP' 
    #             'HertotalteKP126' or 'HertotalteKP128'
    unique <- 0
    uniquesymmetric <- 0
    total <- 1:L
    # Determine the size, cf Table 47 Ko 2009
    if (Growth == 'Lin'){
        total <- 1:L
    } else if (Growth == 'LinOdd'){
        total <- total*2-1        
    } else if (Growth == 'LinEven') {
        total <- total*2        
    } else if (Growth == 'ExpOdd') {
        if (Rule == 'Fejer'){
            total <- 2^total - 1 
            total[1] <- 1            
        } else {
            total <- 2^(total-1) + 1
            total[1] <- 1            
        }
    } else if(Growth == 'ExpEven'){
        total <- 2^(total)
    }
    # Determine the total of number of nominal quadrature formula
    if(Rule=='HertotalteKP126'){
        total <- (1:L)*2-1        
    } else if(Rule == 'HertotalteKP128'){
        total <- (1:L)*2-1        
    } else if(Rule == 'HertotalteKP'){
        total <- (1:L)*2-1        
    }
    # Determine the total of unique nominal quadrature formula
    if((Rule=='ClenshawCurtis') | grepl('Fejer',Rule) | grepl('KP',Rule)){
        unique <- c(1,diff(total))
        uniquesymmetric <- unique/2 
        uniquesymmetric[1] <- 1            
    } else {
        if(Growth=='Lin'){
            # unique <- total - mod(total,2) 
            unique <- total - (total %% 2) 
            unique[1] <- 1
            uniquesymmetric <- unique/2 
            uniquesymmetric[1] <- 1            
        } else if (grepl('Odd',Growth)){
            unique <- total - 1 
            unique[1] <- 1
            uniquesymmetric <- unique/2 
            uniquesymmetric[1] <- 1            
        } else if (grepl('Even',Growth)){
            unique <- total
            uniquesymmetric <- unique/2            
        }
    }
    if(Rule=='GaussLaguerre'){
        total <- 1:L        
        if(Growth == 'Lin'){
            # total <- 1:L            
        } else if (Growth == 'LinOdd'){
            total <- total*2
        } else if (Growth == 'LinEven'){
            total <- total*2            
        } else if(Growth == 'ExpOdd'){
            total <- 2^(total-1) + 1 
            total[1] <- 1   
        } else if(Growth == 'ExpEven'){
            total <- 2^(total)            
        }
        unique <- total
        uniquesymmetric <- rep(0,length(total))
    }
    size <- list(Total = total, Unique = unique, UniqueSymmetric = uniquesymmetric)
    return(size)
}


