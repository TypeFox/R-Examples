posterior_probability <-
function(prevalence,lrpos=-1, lrneg=-1 ,posterior_odds=FALSE)
{
ORPPV <- lrpos *  (prevalence/(1-prevalence))
PPV <-  ORPPV / (1+ORPPV)
ORNPV <- lrneg *  (prevalence/(1-prevalence))
NPV <- 1/(ORNPV+1)
values1 <- list(posterior_probability_neg=NPV)
values2 <- list(posterior_probability_pos=PPV)
values3 <- list(posterior_probability_neg=NPV, posterior_odds_neg=ORNPV)
values4 <- list(posterior_probability_pos=PPV, posterior_odds_pos=ORPPV)
if (posterior_odds == FALSE & lrpos==-1 & lrneg >= 0) {
        return(values1)
    }
if (posterior_odds == FALSE & lrneg==-1 & lrpos >= 0) {
        return(values2)
    }
if (posterior_odds == TRUE & lrpos==-1 & lrneg >= 0) {
        return(values3)
    }
if (posterior_odds ==  TRUE & lrneg==-1 & lrpos >= 0) {
        return(values4)
    }
}
