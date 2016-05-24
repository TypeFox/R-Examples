getCutpoint <-
function(dataset, base.form, cov, median = TRUE){
    if(!median){
        mod.form <- as.formula(paste(as.character(base.form[2]),
                                     as.character(base.form[1]),
                                     cov))
        
        if(length(unique(dataset[[as.character(base.form[2])]])) == 2){
            base.mod <- glm(mod.form, data = dataset, family = 'binomial')
        } else{
            base.mod <- lm(mod.form, data = dataset)
        }
        
        seg.reg <- segmented(base.mod, seg.Z=mod.form[c(1,3)], psi = median(dataset[[cov]]), control = seg.control(it.max = 10000))
        cutpoint <- seg.reg$psi[2]
        print(cov)
        print(cutpoint)
    }
    else{
        cutpoint <- median(dataset[[cov]])
    }
    return(cutpoint)
}
