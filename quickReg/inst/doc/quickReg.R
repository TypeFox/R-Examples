## ----set, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------
# Change the width of html file
options(width = 150)


## ----data-------------------------------------------------------------------------------------------------------------------------------------------

# If you haven't install the package, you can download it from cran

# install.packages("quickReg")

library(quickReg)

# Load the dataset

data(diabetes)

# Show the first 6 rows of the data

head(diabetes)


## ----display----------------------------------------------------------------------------------------------------------------------------------------

show_data<-display(diabetes)

# We can show the results with indices or just the name of variables

show_data[1:2]

show_data$BMI


## ----quickReg---------------------------------------------------------------------------------------------------------------------------------------

# Apply univariate regression models

reg_glm<-reg(data = diabetes, y = 5, factor = c(1, 3, 4), model = 'glm')

# reg_glm have two componets, the regression models in detail and a concentrated data frame

# We can show the detail information with: reg_glm$detail, detail(reg_glm)

reg_glm$detail$BMI

# To show the concentrated data frame: reg_glm$dataframe, dataframe(reg_glm)

dataframe(reg_glm)

# Linear model and cox regression model are also avaiable

reg_lm<-reg(data = diabetes, x = c(1:6,8:12), y = 7, factor = c(1, 3, 4), model = 'lm')

# Use varible names
reg_coxph<-reg(data = diabetes, y = "diabetes", time = "age", factor = c("sex", "smoking", "education"), model = 'coxph')


# Display could be used to a reg class to summarize univariate models

display(reg_glm)

display(reg_lm)

display(reg_coxph)


## ----plot,fig.width=8,fig.height=5------------------------------------------------------------------------------------------------------------------

# `quickReg` package provides forest plot for univariate regression models

plot(reg_glm)

# One OR value is larger than others, we can set the limits
plot(reg_glm,limits=c(NA,3))

plot(reg_glm,limits=c(1,2))

# Sort the variables according to alphabetical

plot(reg_glm,limits=c(NA,3), sort ="alphabetical")

# Similarly, we can plot lm and cox regression results

plot(reg_lm,limits=c(-2,5))

plot(reg_coxph,limits=c(0.5,2))

# Modify plot.reg like ggplot2, add themes from package `ggthemes` 
library(ggplot2);library(ggthemes)

plot(reg_coxph,limits=c(0.5,2))+
  labs(list(title = "Logistic Regression Model", x = "variables"))+
  theme_classic() %+replace% 
  theme(legend.position ="none",axis.text.x=element_text(angle=45,size=rel(1.5)))


