data(CEMSwide)

# Define variable names containing the outcome of the comparisons
columns <- colnames(CEMSwide)[grep("V",colnames(CEMSwide))]

# Define names of the objects that are compared (shoud be separated with space)
# "London Paris" means that the first object is "London" and the second "Paris".
comparison <- c("London Paris", "London Milano", "Paris Milano", 
                "London St.Gallen", "Paris St.Gallen", 
                "Milano St.Gallen", "London Barcelona",
                "Paris Barcelona", "Milano Barcelona", 
                "St.Gallen Barcelona", "London Stockholm", 
                "Paris Stockholm", "Milano Stockholm",
                "St.Gallen Stockholm", "Barcelona Stockholm")

# Transform to 'long' format where v.names="Y" sets the name of 
# our response variable (see ?reshape)
CEMSlong <- wide2long(data=CEMSwide, paircomp=columns, names=comparison, v.names="Y")

head(CEMSlong)
head(CEMSwide)
