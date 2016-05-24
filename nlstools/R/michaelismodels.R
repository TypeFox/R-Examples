"michaelis"   <- as.formula(v ~ S/(S+Km) * Vmax)

"compet_mich"   <- as.formula(v ~ S/(S + Km*(1+I/Ki) ) * Vmax)

"non_compet_mich"   <- as.formula(v ~ S/( (S + Km)*(1+I/Ki) ) * Vmax)
