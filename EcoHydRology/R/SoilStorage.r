SoilStorage <-function(S_avg,field_capacity, soil_water_content, porosity){
w2<-(log(1-(field_capacity/(1-0.4348*S_avg/(2.381*S_avg)))-field_capacity)-log(1-(porosity/(1-2.54/(2.381*S_avg)))-porosity
))/(porosity-field_capacity)##
w1<-log(1-field_capacity/(1-0.4348*S_avg/(2.381*S_avg))-field_capacity)+w2*field_capacity
return(2.381*S_avg*(1-(soil_water_content/(soil_water_content+exp(w1-w2*soil_water_content)))))
}

