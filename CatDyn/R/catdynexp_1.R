catdynexp <-
function(x, p, par, dates, distr)
    {
     fleet.name <- x$Properties$Fleets$Fleet
     distr.set  <- c("poisson","negbin","normal","apnormal","lognormal","aplnormal","gamma")
     if(class(x) != "CatDynData") {stop("Pass an object of class CatDynData as first argument - see the help for as.CatDynData")}
     if(any(is.na(p))) {stop("NAs are not allowed in the perturbations per fleet vector 'p'")}
     if(length(p) < 1 || length(p) > 2) {stop("The integer vector 'p' determines the number of perturbations per fleet; its length must be 1 for single fleet, 2 for two fleets")}
     if(any(p < -20) || any(p > 100)) {stop("The number of perturbations per fleet shall not be more than 100, nor more than 20 for transit fisheries")}
     if(class(p) != "numeric") {stop("'p' must be a numeric vector or a scalar")}
     if(any(is.na(par))) {stop("NAs are not allowed in the parameter vector 'par'")}
     if(length(par) < 4 || length(par) > 208) {stop("For any of the model versions the number of parameters is > 3 and < 208")}
     if(class(par) != "numeric") {stop("'par' must be a numeric vector")}
     if(any(is.na(dates))){stop("NAs are not allowed in the dates vector")}
     if(dates[1] > dates[2]) {stop("Initial date must not be less than next date")}
     if(tail(dates,1) < tail(dates,2)[1]) {stop("Final date must not be less than previous date")}
     options(warn=-1)
     if(length(fleet.name) == 1)
       {
        if(p >= 0)
          {
           if(sum(sort(p) == p)<length(p)) {stop("Number of perturbations per fleet must not be arranged in descending order")}
           if(length(dates) != (sum(p)+2)) {stop("The dates vector must contain initial time step, time steps of all perturbations, and final time step")}
           if(sum(sort(dates) == dates)<length(dates)){stop("Dates must be arranged in ascending order")}
           #if(length(fleet.name) == 1)
           #  {
           #   if(sum(sort(dates) == dates)<length(dates)){stop("Dates must be arranged in ascending order")}
           #  }
           #else if(length(fleet.name) == 2)
           #  {
           #   if(sum(sort(dates[1:(p[1]+1)]) == dates[1:(p[1]+1)]) < length(dates[1:(p[1]+1)]) || sum(sort(dates[(p[1]+2):(p[1]+p[2]+2)]) == dates[(p[1]+2):(p[1]+p[2]+2)]) < length(dates[(p[1]+2):(p[1]+p[2]+2)]))
           #     {stop("Perturbation dates for each fleet must be arranged in ascending order")}
           #   if(length(unique(dates[1:(p[1]+1)])) != length(dates[1:(p[1]+1)]) || length(unique(dates[(p[1]+2):(p[1]+p[2]+2)])) != length(dates[(p[1]+2):(p[1]+p[2]+2)]))
           #     {stop("Perturbation dates inside each fleet should be all distinct and distinct from initial and final date")}
           #  }
          }
        else
          {
          if(length(dates) != (2*abs(p)+2)) {stop("The dates vector must contain initial time step, time step of entry of each perturbation immediately followed by time step of its exit, for all perturbations ordered by time of entry, and final time step")}
          }           
        #if(class(dates) != "integer") {stop("'dates' must be a vector of integers")}
        if(!distr%in%distr.set)
          {stop("distr must be a length-1 character vector with either 'poisson', 'negbin', 'normal', 'apnormal', 'lognormal', 'aplnormal', or 'gamma' ")}
        parlist <- list(par=par, 
                        dates=dates, 
                        obseff1=x$Data[[fleet.name]][,2], 
                        obscat1=x$Data[[fleet.name]][,5], 
                        obsmbm1=x$Data[[fleet.name]][,4], 
                        distr=distr, 
                        properties=x$Properties,
                        output="predict");
        if(p==0)
          {
          if(length(dates) != 2)
            {
             stop("For a 1-fleet 0 perturbation (simple depletion model) dates must be a vector with the following time step marks: initial, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 5)
               {
                stop("For a 1-fleet 0 perturbation (simple depletion model) fit with poisson or adjusted profile likelihood par must be a vector of length 5")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 6)
               {
                stop("For a 1-fleet 0 perturbation (simple depletion model) fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 6")
               }
            }
          results <- do.call(.CDMN0P, parlist);
          }
        else if(p==1)
          {
          if(length(dates) != 3)
            {
             stop("For a 1-fleet 1-perturbation model dates must be a vector with the following time step marks: initial, perturbation 1, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 6)
               {
                stop("For a 1-fleet 1 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 6")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 7)
               {
                stop("For a 1-fleet 1 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 7")
               }
            }
          results <- do.call(.CDMN1P, parlist);
          }
        else if(p==2)
          {
          if(length(dates) != 4)
            {
             stop("For a 1-fleet 2-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, perturbation 2, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 7)
               {
                stop("For a 1-fleet 2 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 7")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 8)
               {
                stop("For a 1-fleet 2 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 8")
               }
            }
          results <- do.call(.CDMN2P, parlist);
          }
        else if(p==3)
          {
          if(length(dates) != 5)
            {
             stop("For a 1-fleet 3-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, perturbation 2, perturbation 3, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 8)
               {
                stop("For a 1-fleet 3 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 8")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 9)
               {
                stop("For a 1-fleet 3 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 9")
               }
            }
          results <- do.call(.CDMN3P, parlist);
          }
        else if(p==4)
          {
          if(length(dates) != 6)
            {
             stop("For a 1-fleet 4-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 4, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 9)
               {
                stop("For a 1-fleet 4 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 9")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 10)
               {
                stop("For a 1-fleet 4 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 10")
               }
            }
          results <- do.call(.CDMN4P, parlist);
          }
        else if(p==5)
          {
          if(length(dates) != 7)
            {
             stop("For a 1-fleet 5-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 5, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 10)
               {
                stop("For a 1-fleet 5 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 10")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 11)
               {
                stop("For a 1-fleet 5 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 11")
               }
            }
          results <- do.call(.CDMN5P, parlist);
          }
        else if(p==6)
          {
          if(length(dates) != 8)
            {
             stop("For a 1-fleet 6-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 6, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 11)
               {
                stop("For a 1-fleet 6 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 11")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 12)
               {
                stop("For a 1-fleet 6 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 12")
               }
            }
          results <- do.call(.CDMN6P, parlist);
          }
        else if(p==7)
          {
          if(length(dates) != 9)
            {
             stop("For a 1-fleet 7-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 7, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 12)
               {
                stop("For a 1-fleet 7 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 12")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 13)
               {
                stop("For a 1-fleet 7 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 13")
               }
            }
          results <- do.call(.CDMN7P, parlist);
          }
        else if(p==8)
          {
          if(length(dates) != 10)
            {
             stop("For a 1-fleet 8-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 8, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 13)
               {
                stop("For a 1-fleet 8 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 13")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 14)
               {
                stop("For a 1-fleet 8 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 14")
               }
            }
          results <- do.call(.CDMN8P, parlist);
          }
        else if(p==9)
          {
          if(length(dates) != 11)
            {
             stop("For a 1-fleet 9-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 9, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 14)
               {
                stop("For a 1-fleet 9 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 14")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 15)
               {
                stop("For a 1-fleet 9 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 15")
               }
            }
          results <- do.call(.CDMN9P, parlist);
          }
        else if(p==10)
          {
          if(length(dates) != 12)
            {
             stop("For a 1-fleet 10-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 10, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 15)
               {
                stop("For a 1-fleet 10 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 15")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 16)
               {
                stop("For a 1-fleet 10 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 16")
               }
            }
          results <- do.call(.CDMN10P, parlist);
          }
        else if(p==11)
          {
          if(length(dates) != 13)
            {
             stop("For a 1-fleet 11-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 11, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 16)
               {
                stop("For a 1-fleet 11 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 16")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 17)
               {
                stop("For a 1-fleet 11 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 17")
               }
            }
          results <- do.call(.CDMN11P, parlist);
          }
        else if(p==12)
          {
          if(length(dates) != 14)
            {
             stop("For a 1-fleet 12-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 12, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 17)
               {
                stop("For a 1-fleet 12 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 17")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 18)
               {
                stop("For a 1-fleet 12 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 18")
               }
            }
          results <- do.call(.CDMN12P, parlist);
          }
        else if(p==13)
          {
          if(length(dates) != 15)
            {
             stop("For a 1-fleet 13-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 13, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 18)
               {
                stop("For a 1-fleet 13 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 18")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 19)
               {
                stop("For a 1-fleet 13 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 19")
               }
            }
          results <- do.call(.CDMN13P, parlist);
          }
        else if(p==14)
          {
          if(length(dates) != 16)
            {
             stop("For a 1-fleet 14-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 14, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 19)
               {
                stop("For a 1-fleet 14 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 19")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 20)
               {
                stop("For a 1-fleet 14 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 20")
               }
            }
          results <- do.call(.CDMN14P, parlist);
          }
        else if(p==15)
          {
          if(length(dates) != 17)
            {
             stop("For a 1-fleet 15-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 15, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 20)
               {
                stop("For a 1-fleet 15 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 20")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 21)
               {
                stop("For a 1-fleet 15 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 21")
               }
            }
          results <- do.call(.CDMN15P, parlist);
          }
        else if(p==16)
          {
          if(length(dates) != 18)
            {
             stop("For a 1-fleet 16-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 16, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 21)
               {
                stop("For a 1-fleet 16 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 21")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 22)
               {
                stop("For a 1-fleet 16 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 22")
               }
            }
          results <- do.call(.CDMN16P, parlist);
          }
        else if(p==17)
          {
          if(length(dates) != 19)
            {
             stop("For a 1-fleet 17-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 17, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 22)
               {
                stop("For a 1-fleet 17 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 22")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 23)
               {
                stop("For a 1-fleet 17 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 23")
               }
            }
          results <- do.call(.CDMN17P, parlist);
          }
        else if(p==18)
          {
          if(length(dates) != 20)
            {
             stop("For a 1-fleet 18-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 18, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 23)
               {
                stop("For a 1-fleet 18 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 23")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 24)
               {
                stop("For a 1-fleet 18 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 24")
               }
            }
          results <- do.call(.CDMN18P, parlist);
          }
        else if(p==19)
          {
          if(length(dates) != 21)
            {
             stop("For a 1-fleet 19-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 19, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 24)
               {
                stop("For a 1-fleet 19 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 24")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 25)
               {
                stop("For a 1-fleet 19 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 25")
               }
            }
          results <- do.call(.CDMN19P, parlist);
          }
        else if(p==20)
          {
          if(length(dates) != 22)
            {
             stop("For a 1-fleet 20-perturbations model dates must be a vector with the following time step marks: initial, perturbation 1, ..., perturbation 20, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 25)
               {
                stop("For a 1-fleet 20 perturbation fit with poisson or adjusted profile likelihood par must be a vector of length 25")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 26)
               {
                stop("For a 1-fleet 20 perturbation fit with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 26")
               }
            }
          results <- do.call(.CDMN20P, parlist);
          }
        #Transit models
        else if(p==-1)
          {
          if(length(dates) != 4)
            {
             stop("For a 1-fleet 1-perturbation transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 6)
               {
                stop("For a 1-fleet 1 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 6")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 7)
               {
                stop("For a 1-fleet 1 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 7")
               }
            }
          results <- do.call(.CDMNT1P, parlist);
          }
        else if(p==-2)
          {
          if(length(dates) != 6)
            {
             stop("For a 1-fleet 2-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, perturbation 2, exit 2, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 7)
               {
                stop("For a 1-fleet 2 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 7")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 8)
               {
                stop("For a 1-fleet 2 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 8")
               }
            }
          results <- do.call(.CDMNT2P, parlist);
          }
        else if(p==-3)
          {
          if(length(dates) != 8)
            {
             stop("For a 1-fleet 3-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, perturbation 2, exit 2, perturbation 3, exit 3, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 8)
               {
                stop("For a 1-fleet 3 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 8")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 9)
               {
                stop("For a 1-fleet 3 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 9")
               }
            }
          results <- do.call(.CDMNT3P, parlist);
          }
        else if(p==-4)
          {
          if(length(dates) != 10)
            {
             stop("For a 1-fleet 4-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 4, exit 4, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 9)
               {
                stop("For a 1-fleet 4 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 9")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 10)
               {
                stop("For a 1-fleet 4 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 10")
               }
            }
          results <- do.call(.CDMNT4P, parlist);
          }
        else if(p==-5)
          {
          if(length(dates) != 12)
            {
             stop("For a 1-fleet 5-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 5, exit 6, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 10)
               {
                stop("For a 1-fleet 5 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 10")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 11)
               {
                stop("For a 1-fleet 5 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 11")
               }
            }
          results <- do.call(.CDMNT5P, parlist);
          }
        else if(p==-6)
          {
          if(length(dates) != 14)
            {
             stop("For a 1-fleet 6-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 6, exit 6, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 11)
               {
                stop("For a 1-fleet 6 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 11")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 12)
               {
                stop("For a 1-fleet 6 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 12")
               }
            }
          results <- do.call(.CDMNT6P, parlist);
          }
        else if(p==-7)
          {
          if(length(dates) != 16)
            {
             stop("For a 1-fleet 7-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 7, exit 7, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 12)
               {
                stop("For a 1-fleet 7 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 12")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 13)
               {
                stop("For a 1-fleet 7 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 13")
               }
            }
          results <- do.call(.CDMNT7P, parlist);
          }
        else if(p==-8)
          {
          if(length(dates) != 18)
            {
             stop("For a 1-fleet 8-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 8, exit 8, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 13)
               {
                stop("For a 1-fleet 8 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 13")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 14)
               {
                stop("For a 1-fleet 8 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 14")
               }
            }
          results <- do.call(.CDMNT8P, parlist);
          }
        else if(p==-9)
          {
          if(length(dates) != 20)
            {
             stop("For a 1-fleet 9-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 9, exit 9, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 14)
               {
                stop("For a 1-fleet 9 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 14")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 15)
               {
                stop("For a 1-fleet 9 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 15")
               }
            }
          results <- do.call(.CDMNT9P, parlist);
          }
        else if(p==-10)
          {
          if(length(dates) != 22)
            {
             stop("For a 1-fleet 10-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 10, exit 10, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 15)
               {
                stop("For a 1-fleet 10 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 15")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 16)
               {
                stop("For a 1-fleet 10 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 16")
               }
            }
          results <- do.call(.CDMNT10P, parlist);
          }
        else if(p==-11)
          {
          if(length(dates) != 24)
            {
             stop("For a 1-fleet 11-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 11, exit 11, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 16)
               {
                stop("For a 1-fleet 11 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 16")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 17)
               {
                stop("For a 1-fleet 11 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 17")
               }
            }
          results <- do.call(.CDMNT11P, parlist);
          }
        else if(p==-12)
          {
          if(length(dates) != 26)
            {
             stop("For a 1-fleet 12-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 12, exit 12, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 17)
               {
                stop("For a 1-fleet 12 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 17")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 18)
               {
                stop("For a 1-fleet 12 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 18")
               }
            }
          results <- do.call(.CDMNT12P, parlist);
          }
        else if(p==-13)
          {
          if(length(dates) != 28)
            {
             stop("For a 1-fleet 13-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 13, exit 13, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 18)
               {
                stop("For a 1-fleet 13 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 18")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 19)
               {
                stop("For a 1-fleet 13 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 19")
               }
            }
          results <- do.call(.CDMNT13P, parlist);
          }
        else if(p==-14)
          {
          if(length(dates) != 30)
            {
             stop("For a 1-fleet 14-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 14, exit 14, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 19)
               {
                stop("For a 1-fleet 14 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 19")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 20)
               {
                stop("For a 1-fleet 14 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 20")
               }
            }
          results <- do.call(.CDMNT14P, parlist);
          }
        else if(p==-15)
          {
          if(length(dates) != 32)
            {
             stop("For a 1-fleet 15-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 15, exit 15, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 20)
               {
                stop("For a 1-fleet 15 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 20")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 21)
               {
                stop("For a 1-fleet 15 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 21")
               }
            }
          results <- do.call(.CDMNT15P, parlist);
          }
        else if(p==-16)
          {
          if(length(dates) != 34)
            {
             stop("For a 1-fleet 16-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 16, exit 16, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 21)
               {
                stop("For a 1-fleet 16 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 21")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 22)
               {
                stop("For a 1-fleet 16 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 22")
               }
            }
          results <- do.call(.CDMNT16P, parlist);
          }
        else if(p==-17)
          {
          if(length(dates) != 36)
            {
             stop("For a 1-fleet 17-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 17, exit 17, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 22)
               {
                stop("For a 1-fleet 17 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 22")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 23)
               {
                stop("For a 1-fleet 17 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 23")
               }
            }
          results <- do.call(.CDMNT17P, parlist);
          }
        else if(p==-18)
          {
          if(length(dates) != 38)
            {
             stop("For a 1-fleet 18-perturbations model transit dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 18, exit 18, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 23)
               {
                stop("For a 1-fleet 18 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 23")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 24)
               {
                stop("For a 1-fleet 18 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 24")
               }
            }
          results <- do.call(.CDMNT18P, parlist);
          }
        else if(p==-19)
          {
          if(length(dates) != 40)
            {
             stop("For a 1-fleet 19-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 19, exit 18, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 24)
               {
                stop("For a 1-fleet 19 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 24")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 25)
               {
                stop("For a 1-fleet 19 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 25")
               }
            }
          results <- do.call(.CDMNT19P, parlist);
          }
        else if(p==-20)
          {
          if(length(dates) != 42)
            {
             stop("For a 1-fleet 20-perturbations transit model dates must be a vector with the following time step marks: initial, perturbation 1, exit 1, ..., perturbation 20, exit 18, final")
            }
          if(distr == "apnormal" | distr == "aplnormal" | distr == "poisson")
            {
             if(length(par) != 25)
               {
                stop("For a 1-fleet 20 perturbation transit model with poisson or adjusted profile likelihood par must be a vector of length 25")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "negbin" | distr == "gamma")
            {
             if(length(par) != 26)
               {
                stop("For a 1-fleet 20 perturbation transit model with negative binomial, normal, lognormal or gamma distribution, par must be a vector of length 26")
               }
            }
          results <- do.call(.CDMNT20P, parlist);
          }
       }
     #Two fleets
     else 
       {
        if(sum(distr%in%distr.set) != 2)
          {stop("'distr' must be a length-2 character vector with any of the 49 pairs that can be made with 'poisson','negbin','normal','apnormal','lognormal','aplnormal', and 'gamma')")}
        parlist <- list(par=par,
                        dates=dates, 
                        obseff1=x$Data[[fleet.name[1]]][,2], 
                        obscat1=x$Data[[fleet.name[1]]][,5], 
                        obsmbm1=x$Data[[fleet.name[1]]][,4], 
                        obseff2=x$Data[[fleet.name[2]]][,2], 
                        obscat2=x$Data[[fleet.name[2]]][,5], 
                        obsmbm2=x$Data[[fleet.name[2]]][,4],
                        distr=distr, 
                        properties=x$Properties,
                        output="predict");
        sdistr.set <- c("poisson","apnormal","aplnormal")
        if(sum(p == c(0,0)) == length(p))
          {
          if(length(dates) != 2)
            {
             stop("For a 2-fleet 0-perturbation 0-perturbation model dates must be a vector with the following time step marks: initial, and final time steps in the season")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 10)
               {
                stop("For a 2-fleet 0 perturbation model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 10 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 9)
               {
                stop("For a 2-fleet 0 perturbation model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 9 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 8)
               {
                stop("For a 2-fleet 0 perturbation model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 8")
               }
            }
          results <- do.call(.CDMN0P0P, parlist);
          }
        else if(sum(p == c(0,1)) == length(p))
          {
          if(length(dates) != 3)
            {
             stop("For a 2-fleet 0P 1P model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 11)
               {
                stop("For a 2-fleet 0P 1P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 11 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 10)
               {
                stop("For a 2-fleet 0P 1P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 10 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 9)
               {
                stop("For a 2-fleet 0P 1P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 9")
               }
            }
          results <- do.call(.CDMN0P1P, parlist);
          }
        else if(sum(p == c(0,2)) == length(p))
          {
          if(length(dates) != 4)
            {
             stop("For a 2-fleet 0-perturbation 2-perturbation model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, fleet 2 perturbation 2, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 12)
               {
                stop("For a 2-fleet 0P 2P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 12 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 11)
               {
                stop("For a 2-fleet 0P 2P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 11 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 10)
               {
                stop("For a 2-fleet 0P 2P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 10")
               }
            }
          results <- do.call(.CDMN0P2P, parlist);
          }
        else if(sum(p == c(0,3)) == length(p))
          {
          if(length(dates) != 5)
            {
             stop("For a 2-fleet 0-perturbation 3-perturbation model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, fleet 2 perturbation 2, fleet 2 perturbation 3, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 13)
               {
                stop("For a 2-fleet 0P 3P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 13 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 12)
               {
                stop("For a 2-fleet 0P 3P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 12 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 11)
               {
                stop("For a 2-fleet 0P 3P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 11")
               }
            }
          results <- do.call(.CDMN0P3P, parlist);
          }
        else if(sum(p == c(0,4)) == length(p))
          {
          if(length(dates) != 6)
            {
             stop("For a 2-fleet 0-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, ..., fleet 2 perturbation 4, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 14)
               {
                stop("For a 2-fleet 0P 4P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 14 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 13)
               {
                stop("For a 2-fleet 0P 4P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 13 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 12)
               {
                stop("For a 2-fleet 0P 4P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 12")
               }
            }
          results <- do.call(.CDMN0P4P, parlist);
          }
        else if(sum(p == c(0,5)) == length(p))
          {
          if(length(dates) != 7)
            {
             stop("For a 2-fleet 0-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 15)
               {
                stop("For a 2-fleet 0P 5P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 15 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 14)
               {
                stop("For a 2-fleet 0P 5P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 14 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 13)
               {
                stop("For a 2-fleet 0P 5P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 13")
               }
            }
          results <- do.call(.CDMN0P5P, parlist);
          }
        else if(sum(p == c(1,1)) == length(p))
          {
          if(length(dates) != 4)
            {
             stop("For a 2-fleet 1-perturbation 1-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 12)
               {
                stop("For a 2-fleet 1P 1P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 12 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 11)
               {
                stop("For a 2-fleet 1P 1P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 11 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 10)
               {
                stop("For a 2-fleet 1P 1P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 10")
               }
            }
          results <- do.call(.CDMN1P1P, parlist);
          }                                                                                                                                            
        else if(sum(p == c(1,2)) == length(p))
          {
          if(length(dates) != 5)
            {
             stop("For a 2-fleet 1-perturbation 2-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, fleet 2 perturbation 2, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 13)
               {
                stop("For a 2-fleet 1P 2P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 13 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 12)
               {
                stop("For a 2-fleet 1P 2P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 12 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 11)
               {
                stop("For a 2-fleet 1P 2P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 11")
               }
            }
          results <- do.call(.CDMN1P2P, parlist);
          }
        else if(sum(p == c(1,3)) == length(p))
          {
          if(length(dates) != 6)
            {
             stop("For a 2-fleet 1-perturbation 3-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, ..., fleet 2 perturbation 3, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 14)
               {
                stop("For a 2-fleet 1P 3P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 14 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 13)
               {
                stop("For a 2-fleet 1P 3P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 13 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 12)
               {
                stop("For a 2-fleet 1P 3P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 12")
               }
            }
          results <- do.call(.CDMN1P3P, parlist);
          }
        else if(sum(p == c(1,4)) == length(p))
          {
          if(length(dates) != 7)
            {
             stop("For a 2-fleet 1-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, ..., fleet 2 perturbation 4, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 15)
               {
                stop("For a 2-fleet 1P 4P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 15 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 14)
               {
                stop("For a 2-fleet 1P 4P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 14 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 13)
               {
                stop("For a 2-fleet 1P 4P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 13")
               }
            }
          results <- do.call(.CDMN1P4P, parlist);
          }
        else if(sum(p == c(1,5)) == length(p))
          {
          if(length(dates) != 8)
            {
             stop("For a 2-fleet 1-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 16)
               {
                stop("For a 2-fleet 1P 5P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 16 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 15)
               {
                stop("For a 2-fleet 1P 5P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 15 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 14)
               {
                stop("For a 2-fleet 1P 5P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 14")
               }
            }
          results <- do.call(.CDMN1P5P, parlist);
          }
        else if(sum(p == c(2,2)) == length(p))
          {
          if(length(dates) != 6)
            {
             stop("For a 2-fleet 2-perturbation 2-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 1 perturbation 2, fleet 2 perturbation 1, fleet 2 perturbation 2, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 14)
               {
                stop("For a 2-fleet 2P 2P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 14 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 13)
               {
                stop("For a 2-fleet 2P 2P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 13 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 12)
               {
                stop("For a 2-fleet 2P 2P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 12")
               }
            }
          results <- do.call(.CDMN2P2P, parlist);
          }
        else if(sum(p == c(2,3)) == length(p))
          {
          if(length(dates) != 7)
            {
             stop("For a 2-fleet 2-perturbation 3-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 1 perturbation 2, fleet 2 perturbation 1, ..., fleet 2 perturbation 3, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 15)
               {
                stop("For a 2-fleet 2P 3P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 15 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 14)
               {
                stop("For a 2-fleet 2P 3P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 14 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 13)
               {
                stop("For a 2-fleet 2P 3P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 13")
               }
            }
          results <- do.call(.CDMN2P3P, parlist);
          }
        else if(sum(p == c(2,4)) == length(p))
          {
          if(length(dates) != 8)
            {
             stop("For a 2-fleet 2-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 1 perturbation 2, fleet 2 perturbation 1, ..., fleet 2 perturbation 4, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 16)
               {
                stop("For a 2-fleet 2P 4P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 16 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 15)
               {
                stop("For a 2-fleet 2P 4P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 15 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 14)
               {
                stop("For a 2-fleet 2P 4P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 14")
               }
            }
          results <- do.call(.CDMN2P4P, parlist);
          }
        else if(sum(p == c(2,5)) == length(p))
          {
          if(length(dates) != 9)
            {
             stop("For a 2-fleet 2-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, fleet 1 perturbation 2, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 17)
               {
                stop("For a 2-fleet 2P 5P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 17 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 16)
               {
                stop("For a 2-fleet 2P 5P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 16 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 15)
               {
                stop("For a 2-fleet 2P 5P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 15")
               }
            }
          results <- do.call(.CDMN2P5P, parlist);
          }
        else if(sum(p == c(3,3)) == length(p))
          {
          if(length(dates) != 8)
            {
             stop("For a 2-fleet 3-perturbation 3-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 3, ..., fleet 2 perturbation 2, fleet 2 perturbation 3, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 16)
               {
                stop("For a 2-fleet 3P 3P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 16 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 15)
               {
                stop("For a 2-fleet 3P 3P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 15 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 14)
               {
                stop("For a 2-fleet 3P 3P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 14")
               }
            }
          results <- do.call(.CDMN3P3P, parlist);
          }
        else if(sum(p == c(3,4)) == length(p))
          {
          if(length(dates) != 9)
            {
             stop("For a 2-fleet 3-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 3, fleet 2 perturbation 1, ..., fleet 2 perturbation 4, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 17)
               {
                stop("For a 2-fleet 3P 4P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 17 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 16)
               {
                stop("For a 2-fleet 3P 4P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 16 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 15)
               {
                stop("For a 2-fleet 3P 4P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 15")
               }
            }
          results <- do.call(.CDMN3P4P, parlist);
          }
        else if(sum(p == c(3,5)) == length(p))
          {
          if(length(dates) != 10)
            {
             stop("For a 2-fleet 3-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 3, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 18)
               {
                stop("For a 2-fleet 3P 5P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 18 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 17)
               {
                stop("For a 2-fleet 3P 5P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 17 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 16)
               {
                stop("For a 2-fleet 3P 5P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 16")
               }
            }
          results <- do.call(.CDMN3P5P, parlist);
          }
        else if(sum(p == c(4,4)) == length(p))
          {
          if(length(dates) != 10)
            {
             stop("For a 2-fleet 4-perturbation 4-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 3, fleet 1 perturbation 4, fleet 2 perturbation 1, fleet 2 perturbation 2, fleet 2 perturbation 3, fleet 2 perturbation 4, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 18)
               {
                stop("For a 2-fleet 4P 4P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 18 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 17)
               {
                stop("For a 2-fleet 4P 4P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 17 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 16)
               {
                stop("For a 2-fleet 4P 4P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 16")
               }
            }
          results <- do.call(.CDMN4P4P, parlist);
          }
        else if(sum(p == c(4,5)) == length(p))
          {
          if(length(dates) != 11)
            {
             stop("For a 2-fleet 4-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 4, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 19)
               {
                stop("For a 2-fleet 4P 5P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 19 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 18)
               {
                stop("For a 2-fleet 4P 5P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 18 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 17)
               {
                stop("For a 2-fleet 4P 5P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 17")
               }
            }
          results <- do.call(.CDMN4P5P, parlist);
          }
        else if(sum(p == c(5,5)) == length(p))
          {
          if(length(dates) != 12)
            {
             stop("For a 2-fleet 5-perturbation 5-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 5, fleet 2 perturbation 1, ..., fleet 2 perturbation 5, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 20)
               {
                stop("For a 2-fleet 5P 5P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 20 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 19)
               {
                stop("For a 2-fleet 5P 5P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 19 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 18)
               {
                stop("For a 2-fleet 5P 5P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 18")
               }
            }
          results <- do.call(.CDMN5P5P, parlist);
          }
        else if(sum(p == c(6,6)) == length(p))
          {
          if(length(dates) != 14)
            {
             stop("For a 2-fleet 6-perturbation 6-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 6, fleet 2 perturbation 1, ..., fleet 2 perturbation 6, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 22)
               {
                stop("For a 2-fleet 6P 6P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 22 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 21)
               {
                stop("For a 2-fleet 6P 6P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 21 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 20)
               {
                stop("For a 2-fleet 6P 6P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 20")
               }
            }
          results <- do.call(.CDMN6P6P, parlist);
          }
        else if(sum(p == c(7,7)) == length(p))
          {
          if(length(dates) != 16)
            {
             stop("For a 2-fleet 7-perturbation 7-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 7, fleet 2 perturbation 1, ..., fleet 2 perturbation 7, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 24)
               {
                stop("For a 2-fleet 7P 7P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 24 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 23)
               {
                stop("For a 2-fleet 7P 7P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 23 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 22)
               {
                stop("For a 2-fleet 7P 7P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 22")
               }
            }
          results <- do.call(.CDMN7P7P, parlist);
          }
        else if(sum(p == c(8,8)) == length(p))
          {
          if(length(dates) != 18)
            {
             stop("For a 2-fleet 8-perturbation 8-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 8, fleet 2 perturbation 1, ..., fleet 2 perturbation 8, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 26)
               {
                stop("For a 2-fleet 8P 8P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 26 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 25)
               {
                stop("For a 2-fleet 8P 8P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 25 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 24)
               {
                stop("For a 2-fleet 8P 8P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 24")
               }
            }
          results <- do.call(.CDMN8P8P, parlist);
          }
        else if(sum(p == c(9,9)) == length(p))
          {
          if(length(dates) != 20)
            {
             stop("For a 2-fleet 9-perturbation 9-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 9, fleet 2 perturbation 1, ..., fleet 2 perturbation 9, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 28)
               {
                stop("For a 2-fleet 9P 9P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 28 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 27)
               {
                stop("For a 2-fleet 9P 9P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 27 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 26)
               {
                stop("For a 2-fleet 9P 9P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 26")
               }
            }
          results <- do.call(.CDMN9P9P, parlist);
          }
        else if(sum(p == c(10,10)) == length(p))
          {
          if(length(dates) != 22)
            {
             stop("For a 2-fleet 10-perturbation 10-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 10, fleet 2 perturbation 1, ..., fleet 2 perturbation 10, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 30)
               {
                stop("For a 2-fleet 10P 10P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 30 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 29)
               {
                stop("For a 2-fleet 10P 10P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 29 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 28)
               {
                stop("For a 2-fleet 10P 10P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 28")
               }
            }
          results <- do.call(.CDMN10P10P, parlist);
          }
        else if(sum(p == c(11,11)) == length(p))
          {
          if(length(dates) != 24)
            {
             stop("For a 2-fleet 11-perturbation 11-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 11, fleet 2 perturbation 1, ..., fleet 2 perturbation 11, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 32)
               {
                stop("For a 2-fleet 11P 11P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 32 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 31)
               {
                stop("For a 2-fleet 11P 11P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 31 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 30)
               {
                stop("For a 2-fleet 11P 11P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 30")
               }
            }
          results <- do.call(.CDMN11P11P, parlist);
          }
        else if(sum(p == c(12,12)) == length(p))
          {
          if(length(dates) != 26)
            {
             stop("For a 2-fleet 12-perturbation 12-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 12, fleet 2 perturbation 1, ..., fleet 2 perturbation 12, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 34)
               {
                stop("For a 2-fleet 12P 12P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 34 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 33)
               {
                stop("For a 2-fleet 12P 12P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 33 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 32)
               {
                stop("For a 2-fleet 12P 12P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 32")
               }
            }
          results <- do.call(.CDMN12P12P, parlist);
          }
        else if(sum(p == c(13,13)) == length(p))
          {
          if(length(dates) != 28)
            {
             stop("For a 2-fleet 13-perturbation 13-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 13, fleet 2 perturbation 1, ..., fleet 2 perturbation 13, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 36)
               {
                stop("For a 2-fleet 13P 13P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 36 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 35)
               {
                stop("For a 2-fleet 13P 13P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 35 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 34)
               {
                stop("For a 2-fleet 13P 13P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 34")
               }
            }
          results <- do.call(.CDMN13P13P, parlist);
          }
        else if(sum(p == c(14,14)) == length(p))
          {
          if(length(dates) != 30)
            {
             stop("For a 2-fleet 14-perturbation 14-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 14, fleet 2 perturbation 1, ..., fleet 2 perturbation 14, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 38)
               {
                stop("For a 2-fleet 14P 14P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 38 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 37)
               {
                stop("For a 2-fleet 14P 14P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 37 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 36)
               {
                stop("For a 2-fleet 14P 14P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 36")
               }
            }
          results <- do.call(.CDMN14P14P, parlist);
          }
        else if(sum(p == c(15,15)) == length(p))
          {
          if(length(dates) != 32)
            {
             stop("For a 2-fleet 15-perturbation 15-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 15, fleet 2 perturbation 1, ..., fleet 2 perturbation 15, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 40)
               {
                stop("For a 2-fleet 15P 15P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 40 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 39)
               {
                stop("For a 2-fleet 15P 15P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 39 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 38)
               {
                stop("For a 2-fleet 15P 15P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 38")
               }
            }
          results <- do.call(.CDMN15P15P, parlist);
          }
        else if(sum(p == c(16,16)) == length(p))
          {
          if(length(dates) != 34)
            {
             stop("For a 2-fleet 16-perturbation 16-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 16, fleet 2 perturbation 1, ..., fleet 2 perturbation 16, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 42)
               {
                stop("For a 2-fleet 16P 16P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 42 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 41)
               {
                stop("For a 2-fleet 16P 16P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 41 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 40)
               {
                stop("For a 2-fleet 16P 16P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 40")
               }
            }
          results <- do.call(.CDMN16P16P, parlist);
          }
        else if(sum(p == c(17,17)) == length(p))
          {
          if(length(dates) != 36)
            {
             stop("For a 2-fleet 17-perturbation 17-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 17, fleet 2 perturbation 1, ..., fleet 2 perturbation 17, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 44)
               {
                stop("For a 2-fleet 17P 17P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 44 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 43)
               {
                stop("For a 2-fleet 17P 17P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 43 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 42)
               {
                stop("For a 2-fleet 17P 17P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 42")
               }
            }
          results <- do.call(.CDMN17P17P, parlist);
          }
        else if(sum(p == c(18,18)) == length(p))
          {
          if(length(dates) != 38)
            {
             stop("For a 2-fleet 18-perturbation 18-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 18, fleet 2 perturbation 1, ..., fleet 2 perturbation 18, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 46)
               {
                stop("For a 2-fleet 18P 18P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 46 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 45)
               {
                stop("For a 2-fleet 18P 18P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 45 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 44)
               {
                stop("For a 2-fleet 18P 18P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 44")
               }
            }
          results <- do.call(.CDMN18P18P, parlist);
          }
        else if(sum(p == c(19,19)) == length(p))
          {
          if(length(dates) != 40)
            {
             stop("For a 2-fleet 19-perturbation 19-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 19, fleet 2 perturbation 1, ..., fleet 2 perturbation 19, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 48)
               {
                stop("For a 2-fleet 19P 19P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 48 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 47)
               {
                stop("For a 2-fleet 19P 19P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 47 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 46)
               {
                stop("For a 2-fleet 19P 19P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 46")
               }
            }
          results <- do.call(.CDMN19P19P, parlist);
          }
        else if(sum(p == c(20,20)) == length(p))
          {
          if(length(dates) != 42)
            {
             stop("For a 2-fleet 20-perturbation 20-perturbation model dates must be a vector with the following time step marks: initial, fleet 1 perturbation 1, ..., fleet 1 perturbation 20, fleet 2 perturbation 1, ..., fleet 2 perturbation 20, final")
            }
          if(sum(distr%in%sdistr.set) == 0)
            {
             if(length(par) != 50)
               {
                stop("For a 2-fleet 20P 20P model with either negative binomial, normal, lognormal or gamma distribution for both fleets, par must be a vector of length 50 \n Put initial values for dispersion parameters for fleet 1 and fleet 2 in the last two positions of the parameter vector, in their obvious order")
               }
            }
          if(sum(distr%in%sdistr.set) == 1)
            {
             if(length(par) != 49)
               {
                stop("For a 2-fleet 20P 20P model with either negative binomial, normal, lognormal or gamma distribution for one of the fleets, par must be a vector of length 49 \n Put initial value for dispersion parameter of that fleet in the last position of the parameter vector")
               }
            }
          if(sum(distr%in%sdistr.set) == 2)
            {
             if(length(par) != 48)
               {
                stop("For a 2-fleet 20P 20P model with either poisson, apnormal, or aplnormal distribution for both fleets, par must be a vector of length 48")
               }
            }
          results <- do.call(.CDMN20P20P, parlist);
          }
       }
     if(any(results$Model$Results[,4] < 0))
       {stop("Predicted catch at some time steps is negative; consider decreasing natural mortality or increasing abundance. \n In transit fisheries you may also increase vulnerable abundance by delaying exit times")}
     return(results);
    }
