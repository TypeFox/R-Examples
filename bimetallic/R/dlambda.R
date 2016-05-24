dlambda <-
function(R, gam_ca, gam_co, ppv, npv, homRR, N_co, maf, prev, model, diff="gam_ca"){
env =  do.call(chisq.setup, as.list(environment()))
arg = match.arg(diff, c("gam_ca", "gam_co"))


P = env$combGtFreq
if(arg == "gam_ca"){
d0 = ((gam_co+1)*(2*((P[1,3]*ppv+P[2,3]*(1-ppv))/(gam_ca+1)
                     -(gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])
                      /(gam_ca+1)^2)
                    *((gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1)
                     -(gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])/(gam_co+1))
                  /((gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])*R/(gam_ca+1)
                   +(gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1))
                  -((gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1)
                   -(gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])/(gam_co+1))^2
                   *(-(gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])*R
                    /(gam_ca+1)^2
                    -(gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1)^2
                    +(P[1,3]*ppv+P[2,3]*(1-ppv))/(gam_ca+1))
                   /((gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])*R/(gam_ca+1)
                    +(gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1))^2
                  +2*((P[1,2]*ppv+P[2,2]*(1-ppv))/(gam_ca+1)
                     -(gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])
                      /(gam_ca+1)^2)
                    *((gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1)
                     -(gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])/(gam_co+1))
                   /((gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])*R/(gam_ca+1)
                    +(gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1))
                  -((gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1)
                   -(gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])/(gam_co+1))^2
                   *(-(gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])*R
                    /(gam_ca+1)^2
                    -(gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1)^2
                    +(P[1,2]*ppv+P[2,2]*(1-ppv))/(gam_ca+1))
                   /((gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])*R/(gam_ca+1)
                    +(gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1))^2
                  +2*((P[1,1]*ppv+P[2,1]*(1-ppv))/(gam_ca+1)
                     -(gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])
                      /(gam_ca+1)^2)
                    *((gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1)
                     -(gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])/(gam_co+1))
                   /((gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])*R/(gam_ca+1)
                    +(gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1))
                  -((gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1)
                   -(gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])/(gam_co+1))^2
                   *(-(gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])*R
                    /(gam_ca+1)^2
                    -(gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1)^2
                    +(P[1,1]*ppv+P[2,1]*(1-ppv))/(gam_ca+1))
                   /((gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])*R/(gam_ca+1)
                    +(gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1))^2))
} else if(arg == "gam_co"){
 d0 = ((gam_co+1)*(2*((gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])/(gam_co+1)^2
                     -(P[2,3]*npv+P[1,3]*(1-npv))/(gam_co+1))
                    *((gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1)
                     -(gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])/(gam_co+1))
                  /((gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])*R/(gam_ca+1)
                   +(gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1))
                  -(P[2,3]*npv+P[1,3]*(1-npv))
                   *((gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1)
                    -(gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])/(gam_co+1))^2*R
                   /((gam_ca+1)*((gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])*R
                                /(gam_ca+1)
                                +(gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])
                                 /(gam_ca+1))^2)
                  +2*((gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])/(gam_co+1)^2
                     -(P[2,2]*npv+P[1,2]*(1-npv))/(gam_co+1))
                    *((gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1)
                     -(gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])/(gam_co+1))
                   /((gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])*R/(gam_ca+1)
                    +(gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1))
                  -(P[2,2]*npv+P[1,2]*(1-npv))
                   *((gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1)
                    -(gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])/(gam_co+1))^2*R
                   /((gam_ca+1)*((gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])*R
                                /(gam_ca+1)
                                +(gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])
                                 /(gam_ca+1))^2)
                  +2*((gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])/(gam_co+1)^2
                     -(P[2,1]*npv+P[1,1]*(1-npv))/(gam_co+1))
                    *((gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1)
                     -(gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])/(gam_co+1))
                   /((gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])*R/(gam_ca+1)
                    +(gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1))
                  -(P[2,1]*npv+P[1,1]*(1-npv))
                   *((gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1)
                    -(gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])/(gam_co+1))^2*R
                   /((gam_ca+1)*((gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])*R
                                /(gam_ca+1)
                                +(gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])
                                 /(gam_ca+1))^2))
        +((gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1)
         -(gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])/(gam_co+1))^2
         /((gam_co*(P[2,3]*npv+P[1,3]*(1-npv))+P[2,3])*R/(gam_ca+1)
          +(gam_ca*(P[1,3]*ppv+P[2,3]*(1-ppv))+P[1,3])/(gam_ca+1))
        +((gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1)
         -(gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])/(gam_co+1))^2
         /((gam_co*(P[2,2]*npv+P[1,2]*(1-npv))+P[2,2])*R/(gam_ca+1)
          +(gam_ca*(P[1,2]*ppv+P[2,2]*(1-ppv))+P[1,2])/(gam_ca+1))
        +((gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1)
         -(gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])/(gam_co+1))^2
         /((gam_co*(P[2,1]*npv+P[1,1]*(1-npv))+P[2,1])*R/(gam_ca+1)
          +(gam_ca*(P[1,1]*ppv+P[2,1]*(1-ppv))+P[1,1])/(gam_ca+1)))
}

d0
}

