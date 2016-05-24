#' @name haz_fl_occ
#' @title Building occupancy classes, specific to flood.
#' 
#' @description Modified from Table 3.1 from HAZUS-MH MR4 Flood Model Technical 
#' Manual.
#' 
#' @details
#' Variables:
#' 
#' \itemize{
#'  \item Occupy_Class - Generalized occupancy class, e.g., RES1 and RES2 
#'  occupanices were assigned to the RES class. One of - AGR, AUTO, COM, EDU, 
#'  GOV, IND, Other_Occupy, REL and RES
#'  \item Occupancy - Subclasses within each Occupy_Class
#'  \item Occ_Desc1 - Description of Occupy_Class
#'  \item Occ_Desc2 - Description of Occupancy provided by HAZUS
#'  \item SIC_code	- SIC (Standard Industrial Classification) code specified 
#'  by HAZUS
#' }
#'
#' @references HAZUS-MH MR4 Flood Model User Manual and Technical Manual, 
#' August 2009, \url{http://www.fema.gov/protecting-our-communities/hazus/hazus-user-technical-manuals}.
#' 
#' @author Gopi Goteti
#' 
#' @docType data
#' 
#' @usage data(haz_fl_occ)
#' 
#' @format Data frame with 5 columns and 33 rows
#' 
#' @keywords datasets
NULL

#' @name haz_fl_dept
#' @title Depth-damage functions from HAZUS, specific to flood
#' 
#' @description Tables D.28, D.30, D.32 (pg. D-22, D-26 and D-28 of the 
#' User Manual) describe  the attributes of these tables. Data was obtained 
#' from the tables flBldgStructDmgFn, flBldgContDmgFn, flBldgInvDmgFn, 
#' flEssntStructDmgFn, flEssntContDmgFn, flUtilFltyDmgFn, flVehicleDmgFn, 
#' in the MS Access Database flDmRsFn found in the HAZUS software package. 
#' Data from the above tables was combined into a single data frame.
#' 
#' @details
#' Variables:
#' 
#' \itemize{
#'  \item Occupancy - Subclasses within each Occupy_Class 
#'  (\code{\link{haz_fl_occ}})
#'  \item DmgFnId	- Identifier used by HAZUS
#'  \item Source - Source identified by HAZUS
#'  \item Description - Description from HAZUS
#'  \item Comment	- Comments from HAZUS, usually blank
#'  \item Columns beginning with ft - Percent damage at specified flood depth
#'  \item Source_Table - HAZUS table name from which the data was obtained
#'  \item Occupy_Class - Generalized occupancy class (\code{\link{haz_fl_occ}}), 
#'  e.g., RES1 and RES2 occupanices were assigned to the RES class. 
#'  One of - AGR, AUTO, COM, EDU, GOV, IND, Other_Occupy, REL and RES
#'  \item Cover_Class - Coverage class, one of building (Bldg), contents (Cont),
#'  inventory (Inv), or other (Other_Cover)
#' }
#'
#' @references HAZUS-MH MR4 Flood Model User Manual and Technical Manual, 
#' August 2009, \url{http://www.fema.gov/protecting-our-communities/hazus/hazus-user-technical-manuals}.
#' 
#' @author Gopi Goteti
#' 
#' @docType data
#' 
#' @usage data(haz_fl_dept)
#' 
#' @format Data frame with 51 columns and 1260 rows
#' 
#' @keywords datasets
NULL


#' @name haz_fl_velo
#' @title Velocity-depth-damage functions from HAZUS, specific to flood
#' 
#' @description These functions are specified in HAZUS as empirical equations 
#' and indicate whether or not the structure will collapse under a given 
#' combination of flood velocity and flood depth. The empirical equations 
#' have the form specified below, where the coefficients and exponents are 
#' dependent on the velocity and depth thresholds for a given structure type
#' \deqn{coef1 * V ^ {expo1} + coef2 * V + coef3}
#' 
#' @details
#' Variables:
#' 
#' \itemize{
#'  \item num_story - Number of stories (1, 2 or 3)
#'  \item struc_type - Structure type (wood, masonry, concrete or steel)
#'  \item thresh_d - Depth threshold (ft) used in the above equation
#'  \item thresh_v - Velocity threshold (ft/s) used in the above equation
#'  \item coef1	- Coefficient in the above equation
#'  \item coef2 - Similar to coef1
#'  \item coef3 - Similar to coef1
#'  \item expo1 - Exponent in the above equation
#' }
#'
#' @references Obtained from Tables 5.5, 5.6 and 5.7 of the HAZUS-MH MR4 Flood 
#' Model Technical Manual, August 2009, \url{http://www.fema.gov/protecting-our-communities/hazus/hazus-user-technical-manuals}
#' 
#' @author Gopi Goteti
#' 
#' @docType data
#' 
#' @usage data(haz_fl_velo)
#' 
#' @format Data frame with 8 columns and 12 rows
#' 
#' @keywords datasets
NULL


#' @name haz_fl_agri
#' @title Agriculture damage functions from HAZUS, specific to flood
#' 
#' @description Table D.31 (pg. D-27 of the User Manual) describes the 
#' attributes of these damage functions. Data was obtained from the table 
#' flAgDmgFn in the MS Access Database flDmRsFn found in the HAZUS 
#' software package.
#' 
#' @details
#' Variables:
#' 
#' \itemize{
#'  \item Crop - Name or type of the crop (currently 20 possibilities)
#'  \item FunctionSource - Source of the data (either HAZUS default or from USACE)
#'  \item JulianDay - Day of year (1 to 365)
#'  \item PctCropLoss - Maximum potential percentage crop damage
#'  \item PctLossDuration0_d - 0-Day flood duration damage modifier
#'  \item PctLossDuration3_d - 3-Day flood duration damage modifier
#'  \item PctLossDuration7_days - 7-Day flood duration damage modifier
#'  \item PctLossDuration14_days - 14-Day flood duration damage modifier
#' }
#'
#' @references HAZUS-MH MR4 Flood Model User Manual and Technical Manual, 
#' August 2009, \url{http://www.fema.gov/protecting-our-communities/hazus/hazus-user-technical-manuals}
#' 
#' @author Gopi Goteti
#' 
#' @docType data
#' 
#' @usage data(haz_fl_agri)
#' 
#' @format Data frame with 6 columns and 7300 rows
#' 
#' @keywords datasets
NULL

#' @name haz_fl_bridge
#' @title Damage functions from HAZUS for highway, railway and light rail 
#' bridges, specific to flood
#' 
#' @description Table D.29 (pg. D-24 of the User Manual) describes the 
#' attributes of the table. Data was obtained from the table flBridgeDmgFn in 
#' the MS Access Database flDmRsFn found in the HAZUS software package.
#' 
#' @details
#' Variables:
#' 
#' \itemize{
#'  \item BridgeDmgFnId - Identifier used by HAZUS
#'  \item Occupancy - Bridge-specific occupancy
#'  \item Source - Source of the data (currently only HAZUS default)
#'  \item Description - Single span or continuous span
#'  \item RP - Percent damage for return period in years
#' }
#'
#' @references HAZUS-MH MR4 Flood Model User Manual and Technical Manual, 
#' August 2009, \url{http://www.fema.gov/protecting-our-communities/hazus/hazus-user-technical-manuals}
#' 
#' @author Gopi Goteti
#' 
#' @docType data
#' 
#' @usage data(haz_fl_bridge)
#' 
#' @format Data frame with 45 columns and 8 rows
#' 
#' @keywords datasets
NULL

#' @name haz_fl_depr
#' @title Depreciation functions from HAZUS, specfic to flood
#' 
#' @description Table D.5 (pg. D-9 of the User Manual) describes the attributes 
#' of the table. Data was obtained from the table flDepFunction in the MS 
#' Access Database flDmRsFn found in the HAZUS software package.
#' 
#' @details
#' Variables:
#' 
#' \itemize{
#'  \item Age - Average age of the structure in years (0 to 100 years)
#'  \item Next 35 columns - Depreciation by age for 35 occupancy classes, where 
#'  occupancy is defined by \code{\link{haz_fl_occ}}
#' }
#'
#' @references HAZUS-MH MR4 Flood Model User Manual and Technical Manual, 
#' August 2009, \url{http://www.fema.gov/protecting-our-communities/hazus/hazus-user-technical-manuals}
#' 
#' @author Gopi Goteti
#' 
#' @docType data
#' 
#' @usage data(haz_fl_depr)
#' 
#' @format Data frame with 36 columns and 101 rows
#' 
#' @keywords datasets
NULL


#' @name hazus
#' 
#' @title Damage functions and other utilities from FEMA's HAZUS software for 
#' use in modeling financial risk from natural disasters.
#'
#' @details Damage functions are useful in modeling financial risk from natural 
#' disasters. \pkg{hazus} provides the damage functions used by FEMA's HAZUS 
#' software. \pkg{hazus} currently provides functionality to extract and visualize
#' damage functions specific to flood hazard. Hurricane and earthquake damage
#' functions would be included in the future.
#'
#' @import reshape2
#' 
#' @docType package
#' 
#' @author Gopi Goteti
#' 
#' @references HAZUS, The Federal Emergency Management Agency's (FEMA's) 
#' Methodology for Estimating Potential Losses from Disasters, \url{http://www.fema.gov/hazus}
NULL
