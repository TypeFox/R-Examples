#' Sample of dams from the NID database
#' 
#' The data dictionary from NID provides details on the data attributes. 
#' Below is only a short description of the attributes. Note that not all 
#' the attributes specified in the data dictionary are available in the 
#' data from NID. It appears that NID revised the data but not the 
#' data dictionary.
#' 
#' 
#' Variables:
#' 
#' \itemize{
#'  \item Dam_Name - The official name of the dam
#'  \item Other_Dam_Name - Names other than the official name (i.e., 
#'  reservoir name) of the dam in common use
#'  \item State_ID - The Official State or Agency identification number
#'  \item NID_ID - The official NID identification number, known formerly as 
#'  the National ID
#'  \item Num_Separate_Struct - Number of separate structures 
#'  \item Other_Structure_ID - The identification number (S001, S002, etc.) 
#'  for the saddle dam or dike associated with the larger dam project
#'  \item Longitude - Longitude at dam centerline, in decimal degrees, NAD83
#'  \item Latitude - Latitude at dam centerline, in decimal degrees, NAD83
#'  \item Section - The information is in any form that is understandable and 
#'  that clearly designates the individual values, i.e. S21, 73N, R69W
#'  \item County - The name of the county in which the dam is located
#'  \item River - The River or Stream designation
#'  \item Owner_Name - Name(s) of the dam owner
#'  \item Owner_Type - Code to indicate the type of owner
#'  \item Private_Dam - Y or N indicating whether or not the dam is privately 
#'  owned or not
#'  \item Dam_Designer - Name of the principal firm(s) or agency accomplishing 
#'  design of dam and major appurtenant operating features
#'  \item Dam_Type - Codes to indicate the type of dam
#'  \item Core - Code to indicate the position, type of watertight member 
#'  and certainty
#'  \item Foundation - Code for the material upon which dam is founded, 
#'  and certainty
#'  \item Primary_Purpose - Code(s) to indicate the current purpose(s) 
#'  for which the reservoir is used
#'  \item All_Purposes - Code(s) to indicate the current purpose(s) 
#'  for which the reservoir is used
#'  \item Year_Completed - Year when the original main dam structure 
#'  was completed
#'  \item Year_Modified - Year when major modifications or rehabilitation 
#'  of dam or major control structures were completed
#'  \item Dam_Length - Length of the dam, in feet
#'  \item Dam_Height - Height of the dam, in feet 
#'  \item Structural_Height - Structural height of the dam, in feet 
#'  \item Hydraulic_Height - Hydraulic height of the dam, in feet
#'  \item NID_Height - Maximum value of dam height, structural height, 
#'  and hydraulic height. Accepted as the general height of the dam
#'  \item Max_Discharge - Maximum discharge, in cubic feet per second 
#'  \item Max_Storage - Maximum storage, in acre-feet
#'  \item Normal_Storage - Normal storage, in acre-feet
#'  \item NID_Storage - Maximum value of normal storage and maximum storage. 
#'  Accepted as the general storage of the dam
#'  \item Surface_Area - Surface area, in acres
#'  \item Drainage_Area - Drainage area of the dam, in square miles
#'  \item EAP - Code indicating whether this dam has an Emergency Action 
#'  Plan (EAP) developed by the dam owner
#'  \item Inspection_Date - Date of the most recent inspection of the dam 
#'  prior to the transmittal of the data by the submitting agency
#'  \item Inspection_Frequency - The scheduled frequency interval for 
#'  periodic inspections, in years
#'  \item Spillway_Type - Code that describes the type of spillway
#'  \item Spillway_Width - The width of the spillway, in feet
#'  \item Outlet_Gates - Code(s) that describe the type of (1) spillway and 
#'  (2) controlled outlet gates
#'  \item Volume - Total number of cubic yards occupied by the materials 
#'  used in the dam structure
#'  \item Num_Locks - Number of existing navigation locks for the project
#'  \item Length_Locks - Length of the primary navigation lock, in feet
#'  \item Width_Locks - Width of the primary navigation lock, in feet
#'  \item Permitting_Authority - Yes if the state regulatory organization 
#'  has the authority to review
#'  \item Inspection_Authority - Yes if the state regulatory organization 
#'  has the authority to require or perform the inspection
#'  \item Enforcement_Authority - Yes if the state regulatory organization 
#'  has the authority to issue notices
#'  \item Jurisdictional_Dam - Yes if this dam meets the state regulatory 
#'  organization's definition of a jurisdictional dam
#'  \item State_Reg_Dam - Calculated field based on Permitting Authority, 
#'  Inspection Authority and Enforcement Authority
#'  \item State_Reg_Agency - Name of the primary state agency with regulatory 
#'  or approval authority over the dam
#'  \item Fed_Funding - Code identifying which federal agency was involved 
#'  in funding of the dam
#'  \item Fed_Design - Code identifying which federal agency was involved 
#'  in the design of the dam
#'  \item Fed_Construction - Code identifying which federal agency was 
#'  involved in the construction of the dam
#'  \item Fed_Regulatory - Code identifying which federal agency is involved 
#'  in the regulation of the dam
#'  \item Fed_Inspection - Code identifying which federal agency is involved 
#'  in the inspection of the dam
#'  \item Fed_Operation - Code identifying which federal agency is involved 
#'  in the operation of the dam
#'  \item Fed_Owner - Code identifying which federal agency partly or wholly 
#'  owns the dam
#'  \item Fed_Other - Code identifying which federal agency is involved in 
#'  other aspects of the dam
#'  \item Source_Agency - Primary state or federal agency responsible for data
#'  \item State - State where dam is located
#'  \item Submit_Date - Date data was submitted to the US Army Corps of 
#'  Engineers for inclusion to the National Inventory of Dams
#'  \item Url_Address - Web Site for more information on particular dam
#'  \item Congress_Rep - Name of congressional representative for the 
#'  congressional district where dam is located
#'  \item Political_Party - Name of political party associated with the 
#'  congressional representative for the congressional district 
#'  where dam is located
#'  \item Congress_District - Congressional District where dam is located
#' }
#'
#'
#' @references NID: The National Inventory of Dams from the United States 
#' Army Corps of Engineers, http://nid.usace.army.mil, data extracted from 
#' NID's website in March 2014.
#' 
#' @docType data
#' @name nid_sample
#' @usage data(nid_sample)
#' @format Data frame with 64 columns and 100 rows
#' @keywords datasets
NULL

#' Dams in the United States from the National Inventory of Dams (NID).
#'
#' Data from NID was extracted manually from http://nid.usace.army.mil. 
#' Subsequently, the raw data was checked for potential errors and cleaned. 
#' dams package provides sample cleaned data from the NID and provides 
#' functionality to access the entire cleaned NID data.
#'
#' @import RCurl
#' @docType package
#' @name dams
NULL
