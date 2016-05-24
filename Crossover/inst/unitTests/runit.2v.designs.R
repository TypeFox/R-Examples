checkResult <- function (result, refDesign, model) {
    allDesigns <- result@misc$designs
    df <- data.frame(run=numeric(0), step=numeric(0), crit=numeric(0))
    for (i in 1:length(allDesigns)) {
        designs <- allDesigns[[i]]
        for (j in 1:length(designs)) {
            design <- designs[[j]]
            crit <- sum(general.carryover(design, model=model)$Var.trt.pair)/2
            df <- rbind(df, data.frame(run=i, step=j, crit=crit))
        }
    }
    plot <- ggplot(df, aes(x=step, y=crit)) + geom_point() + facet_wrap( ~ run)
    print(plot)
    return(df)
}

test.2v.designs <- function() {    
    if (!"extended" %in% strsplit(Sys.getenv("CROSSOVER_UNIT_TESTS"),",")[[1]]) {
        cat("Skipping design tests for v=2.\n")
        return()
    }    
    
    design1 <- t(rbind(c(1,1,2,2),
                       c(2,2,1,1),
                       c(1,1,2,2),
                       c(2,2,1,1),
                       c(1,2,2,1),
                       c(2,1,1,2)))
    design2 <- t(rbind(c(1,1,2,1),
                       c(2,2,1,2),
                       c(1,1,2,1),
                       c(2,2,1,2),
                       c(1,2,2,1),
                       c(2,1,1,2)))
    design3 <- t(rbind(c(1,1,2,2),
                       c(2,2,1,1),
                       c(1,1,2,2),
                       c(2,2,1,1),
                       c(1,1,2,2),
                       c(2,2,1,1)))
    design4 <- t(rbind(c(1,1,2,2),
                       c(2,2,1,1),
                       c(1,1,2,2),
                       c(2,2,1,1),
                       c(1,1,2,2),
                       c(2,2,1,1))) # Why is there no design 4 in the book?
    design5 <- t(rbind(c(1,2,2,2),
                       c(2,1,1,1),
                       c(1,2,2,2),
                       c(2,1,1,1),
                       c(1,1,2,2),
                       c(2,2,1,1)))
    design6 <- t(rbind(c(1,2,1,2),
                       c(2,1,2,1),
                       c(1,1,2,1),
                       c(2,2,1,2),
                       c(1,2,1,2),
                       c(2,1,2,1)))
    design7 <- t(rbind(c(1,1,2,1),
                       c(2,2,1,2),
                       c(1,1,1,2),
                       c(2,2,2,1),
                       c(1,2,1,1),
                       c(2,1,2,2)))
    design8 <- t(rbind(c(1,2,2,2),
                       c(2,1,1,1),
                       c(1,1,2,2),
                       c(2,2,1,1),
                       c(1,2,1,2),
                       c(2,1,2,1)))
    
    designs <- list(design1, design2, design3, design4, design5, design6, design7, design8)
    
    models <- c("Standard additive model",
                "Self-adjacency model",
                "Proportionality model",
                "Placebo model",
                "No carry-over into self model",
                "Treatment decay model",
                "Full set of interactions",
                "Second-order carry-over effects")
    
    s <- 6; p <- 4; v <- 2
    
    results <- list()
    
    for (i in c(1,2,3,4,5,6)) {  
        model <- models[i]
        cat("======= ", model, " =======","\n")
        
        result <- searchCrossOverDesign(s=s, p=p, v=v, model=model, start.designs=list(designs[[i]]))
        print(result)
        print(plot(result))
        
        print(general.carryover(designs[[i]], model=i))  
        
        checkResult(result, designs[[i]], model=i)
    }
    
    estimable(design3, v, model=3)
    
    searchCrossOverDesign(s=s, p=p, v=v, model=3)
    searchCrossOverDesign(s=s, p=p, v=v, model=3, start.designs=list(design3))
    
}

