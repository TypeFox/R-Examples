# Rattle: A GUI for Data Mining in R
#
# AUDIT DATASET
#
# Time-stamp: <2014-09-05 21:28:19 gjw>
#
# Copyright (c) 2009-2014 Togaware Pty Ltd
#
# This file is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.
#
########################################################################
#
# Generate an audit dataset that is fictional but illustrates the typcial
# financial audit.

acquireAuditData <- function(write.to.file=FALSE)
{
  if (!file.exists('survey.csv'))
  {
    UCI <- "ftp://ftp.ics.uci.edu/pub"
    REPOS <- "machine-learning-databases"
    survey.url <- sprintf("%s/%s/adult/adult.data", UCI, REPOS)
    download.file(survey.url, "survey.data")
    survey <- read.csv("survey.data", header=F, strip.white=TRUE,
                       na.strings="?",
                       col.names=c("Age", "Workclass", "fnlwgt", 
                         "Education", "Education.Num", "Marital.Status", 
                         "Occupation", "Relationship", "Race", "Gender", 
                         "Capital.Gain", "Capital.Loss", 
                         "Hours.Per.Week", "Native.Country", 
                         "Salary.Group"))
    write.table(survey, "survey.csv", sep=",", row.names=F)
  }

  survey <- read.csv("survey.csv")

  audit <- survey[,c(1,2,4,6,7,8,10,12,13,14,11,15)]

  colnames(audit)[2] <- "Employment"
  colnames(audit)[4] <- "Marital"
  colnames(audit)[6] <- "Income"
  colnames(audit)[8] <- "Deductions"
  colnames(audit)[9] <- "Hours"
  colnames(audit)[10] <- "Accounts"
  colnames(audit)[11] <- "Adjustment"
  colnames(audit)[12] <- "Adjusted"
  
  audit$Adjusted <- as.integer(audit$Adjusted)-1

  # Make sure most productive cases have an adjustment

  adj <- audit[audit$Adjusted==0 & audit$Adjustment != 0, 'Adjustment']
  a <- length(adj)
  m <- length(audit[audit$Adjusted==1 & audit$Adjustment==0,'Adjusted'])
  r <- m%/%a*a

  set.seed(12345)
  audit[audit$Adjusted==1 & audit$Adjustment==0, 'Adjustment'][sample(m, r)] <-
    as.integer(adj*(rnorm(r) + 2))

  # Make sure no nonproductive case has an adjustment

  audit[audit$Adjusted==0 & audit$Adjustment!=0,'Adjustment'] <- 0

  # Tidyup ForeignAccounts

  levels(audit$Accounts)[6] <- "NewZealand"
  levels(audit$Accounts)[8] <- "Singapore"
  levels(audit$Accounts)[15] <- "Holand"
  levels(audit$Accounts)[28] <- "Fiji"
  levels(audit$Accounts)[33] <- "Malaysia"
  levels(audit$Accounts)[35] <- "Vietnam"
  levels(audit$Accounts)[38] <- "Indonesia"
  levels(audit$Accounts)[39] <- "UnitedStates"

  # Tidyup Employment

  levels(audit$Employment)[1] <- "PSFederal"
  levels(audit$Employment)[2] <- "PSLocal"
  levels(audit$Employment)[3] <- "Unemployed"
  levels(audit$Employment)[5] <- "SelfEmp"
  levels(audit$Employment)[6] <- "Consultant"
  levels(audit$Employment)[7] <- "PSState"
  levels(audit$Employment)[8] <- "Volunteer"
  
  # Tidyup Marital
  
  levels(audit$Marital)[2] <- "Married"
  levels(audit$Marital)[3] <- "Married"
  levels(audit$Marital)[4] <- "Absent"
  levels(audit$Marital)[5] <- "Unmarried"
  
  # Tidyup Occupation
  
  levels(audit$Occupation)[1] <- "Clerical"
  levels(audit$Occupation)[2] <- "Military"
  levels(audit$Occupation)[3] <- "Repair"
  levels(audit$Occupation)[4] <- "Executive"
  levels(audit$Occupation)[5] <- "Farming"
  levels(audit$Occupation)[6] <- "Cleaner"
  levels(audit$Occupation)[7] <- "Machinist"
  levels(audit$Occupation)[8] <- "Service"
  levels(audit$Occupation)[9] <- "Home"
  levels(audit$Occupation)[10] <- "Professional"
  levels(audit$Occupation)[11] <- "Protective"
  levels(audit$Occupation)[12] <- "Sales"
  levels(audit$Occupation)[13] <- "Support"
  levels(audit$Occupation)[14] <- "Transport"
  
  levels(audit$Education)[1] <- "Yr10"
  levels(audit$Education)[2] <- "Yr11"
  levels(audit$Education)[3] <- "Yr12"
  levels(audit$Education)[4] <- "Yr1t4"
  levels(audit$Education)[5] <- "Yr5t6"
  levels(audit$Education)[6] <- "Yr7t8"
  levels(audit$Education)[7] <- "Yr9"
  levels(audit$Education)[8] <- "Associate"
  levels(audit$Education)[9] <- "Vocational"
  levels(audit$Education)[10] <- "Bachelor"
  levels(audit$Education)[11] <- "Doctorate"
  levels(audit$Education)[12] <- "HSgrad"
  levels(audit$Education)[13] <- "Master"
  levels(audit$Education)[14] <- "Preschool"
  levels(audit$Education)[15] <- "Professional"
  levels(audit$Education)[16] <- "College"
  
  # Turn Relationship into Income
  
  set.seed(12345)
  audit$Income <- round(abs(as.numeric(audit$Income)*rnorm(length(audit$Income),
                                                           35000, 15000)), 2)

  # Make deductions look more 0 for the non-productive cases!

  audit[audit$Adjusted==0,'Deductions'] <-
    audit[audit$Adjusted==0,'Deductions']/1.5

  # Sample just 2000 cases and add an Identifier - always the same

  set.seed(12345)
  cases <- sample(nrow(audit), 2000)
  set.seed(12345)
  idents <- as.integer(sort(runif(2000, 1000000, 9999999)))
  audit <- cbind(ID=idents, audit[cases,])
  
  # Use standard prefixes
  
  colnames(audit)[11] <- "IGNORE_Accounts" # randomForest can't handle
  colnames(audit)[12] <- "RISK_Adjustment" 
  colnames(audit)[13] <- "TARGET_Adjusted"

  audit.orig <- audit
  
  # Write out the data
  
  if (write.to.file) 
  {
    audit <- read.csv("audit.csv")
    save(audit, file="audit.RData", compress=TRUE)
    write.table(audit, "audit.csv", sep=",", row.names=FALSE)
  
    arff <- audit
    arff$TARGET_Adjusted <- as.factor(arff$TARGET_Adjusted)
    if (write.to.file) foreign::write.arff(arff, "audit.arff")
  
    # Create a dataset with special variable names.
    # 080709 I now do this as default.
  
    # colnames(audit)[11] <- "IGNORE_Accounts"
    # colnames(audit)[12] <- "RISK_Adjustment"
    # write.table(audit, "audit_auto.csv", sep=",", row.names=FALSE)
  
    # Create a dataset with many more missing values.
  
    mr <- sample(1:nrow(audit), nrow(audit)/4, replace=TRUE)
    mc <- sample(2:(ncol(audit)-1), nrow(audit)/4, replace=TRUE)
  
    for (i in 1:(nrow(audit)/4))
    {
      audit[mr[i], mc[i]] <- NA
    }
    write.table(audit, "audit_missing.csv", sep=",", row.names=FALSE)
  }
  if (write.to.file)
    invisible(audit.orig)
  else
    return(audit.orig)
}


