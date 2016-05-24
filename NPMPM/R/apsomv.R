apsomv <-
function(pstep){

currentmatrix <- matrix[(matrix$id==pstep$matrix_id),]

# NA => aw & pH VALUE
if (is.na(currentmatrix$pH_min)) {currentmatrix$pH_min<-7}
if (is.na(currentmatrix$pH_max)) {currentmatrix$pH_max<-7}
if (is.na(currentmatrix$aw_min)) {currentmatrix$aw_min<-1}
if (is.na(currentmatrix$aw_max)) {currentmatrix$aw_max<-1}

mvalues <- data.frame(id=numeric() ,  
temperature_C=numeric() , 
pH=numeric(), 
aw=numeric()
)

# EXACT VALUES - RANGE OF PARAMETERS IN INTERVAL [MIN,MAX]
suppl <- raw_data[(raw_data$temperature_C<=pstep$temp_max) & (raw_data$temperature_C>=pstep$temp_min) & 
(raw_data$pH<=currentmatrix$pH_max) & (raw_data$pH>=currentmatrix$pH_min) & (raw_data$aw<=currentmatrix$aw_max) & 
(raw_data$aw>=currentmatrix$aw_min),]
mvalues <- rbind(mvalues,suppl)

#
# WIDEN THE RANGE IN WHICH IS SEARCHED (aw AND pH)
if (length(mvalues$id)<30) {
# aw +- 0.05
mvalues <- data.frame(id=numeric() ,  
temperature_C=numeric() , 
pH=numeric(), 
aw=numeric()
)
suppl <- raw_data[(raw_data$temperature_C<=pstep$temp_max) & (raw_data$temperature_C>=pstep$temp_min) & 
(raw_data$pH<=currentmatrix$pH_max) & (raw_data$pH>=currentmatrix$pH_min) & (raw_data$aw<=currentmatrix$aw_max+0.05) & 
(raw_data$aw>=currentmatrix$aw_min-0.05),]
mvalues <- rbind(mvalues,suppl)
#
if (length(mvalues$id)<30) {
# pH +-1 AND aw +- 0.05
mvalues <- data.frame(id=numeric() ,  
temperature_C=numeric() , 
pH=numeric(), 
aw=numeric()
)
suppl <- raw_data[(raw_data$temperature_C<=pstep$temp_max) & (raw_data$temperature_C>=pstep$temp_min) & 
(raw_data$pH<=currentmatrix$pH_max+1) & (raw_data$pH>=currentmatrix$pH_min-1) & (raw_data$aw<=currentmatrix$aw_max+0.05) & 
(raw_data$aw>=currentmatrix$aw_min-0.05),]
mvalues <- rbind(mvalues,suppl)
}} # END IF

#
# INCLUDE NA FOR aw AND pH 
#
if (length(mvalues$id)<30) {
# INCLUDE NA FOR aw 
suppl <- raw_data[(raw_data$temperature_C<=pstep$temp_max) & (raw_data$temperature_C>=pstep$temp_min) & 
(raw_data$pH<=currentmatrix$pH_max) & (raw_data$pH>=currentmatrix$pH_min) &
(is.na(raw_data$aw)) 
,]
mvalues <- rbind(mvalues,suppl)
#
if (length(mvalues$id)<30) {
# INCLUDE NA FOR pH 
suppl <- raw_data[(raw_data$temperature_C<=pstep$temp_max) & (raw_data$temperature_C>=pstep$temp_min) & 
(is.na(raw_data$pH))  & 
(raw_data$aw<=currentmatrix$aw_max) & (raw_data$aw>=currentmatrix$aw_min)
,]
mvalues <- rbind(mvalues,suppl)
#
if (length(mvalues$id)<30) {
# INCLUDE NA FOR BOTH aw AND pH 
suppl <- raw_data[(raw_data$temperature_C<=pstep$temp_max) & (raw_data$temperature_C>=pstep$temp_min) & 
(is.na(raw_data$pH))  & 
(is.na(raw_data$aw)) 
,]
mvalues <- rbind(mvalues,suppl)
}}} # END IF


# 
# TAKE ALL MEASUREMENTS IN THE SPECIFIED TEMPERATURE RANGE
if (length(mvalues$id)<30) {
mvalues <- data.frame(id=numeric() ,  
temperature_C=numeric() , 
pH=numeric(), 
aw=numeric()
)
suppl <- raw_data[(raw_data$temperature_C<=pstep$temp_max) & (raw_data$temperature_C>=pstep$temp_min) ,]
mvalues <- rbind(mvalues,suppl)
}

# MAKE SURE THAT mvalues IS NOT EMPTY, THERE SHOULD BE AT LEAST ONE SERIES OF MEASURED VALUES!!!
# EXPANSION OF THE TEMPERATURE RANGE PLUS/MINUS 1°C


NAwarning <<- "No"
if (length(mvalues$id)==0) {
ii <- 0
while (length(mvalues$id)==0) {
ii <- ii+1
suppl <- raw_data[(raw_data$temperature_C<=pstep$temp_max+ii) & (raw_data$temperature_C>=pstep$temp_min-ii) ,]
mvalues <- rbind(mvalues,suppl)
print("Warning: No appropriate series of measured values for this process step")
NAwarning <<- "Yes"
} # END WHILE
} # END IF

mvaluesid <- mvalues$id
return(mvaluesid)
} # END FUNCTION

