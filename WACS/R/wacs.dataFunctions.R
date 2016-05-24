###################################################################
#
# This function is part of WACSgen V1.0
# Copyright © 2013,2014,2015, D. Allard, BioSP,
# and Ronan Trépos MIA-T, INRA
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. http://www.gnu.org
#
###################################################################
  
# Internal function that builds the mapping of names between original 
# data frame names and the expected names for WACS model
#
#  mapping_user The mapping of names given by the user
#  names_data   The vector of column names into original data
#  Trange      TRUE if one uses the Trange option
# 
#  A dataframe with two columns, the first are the original names,
#  the second is the names in the matrix used by WACS model
# 

wacs.buildMapping = function(mapping_user, names_data, Trange)
{
    mapCol1 = names_data;
    mapCol2 =  c();

    for (i in 1:length(names_data)) {
        nd = names_data[i];
        nde = nd;
        if (! is.null(mapping_user) && (nd %in% names(mapping_user))) {
            nde = mapping_user[[nd]];
        }
        mapCol2 = c(mapCol2, nde);
    }
    if (Trange) {
        mapCol1 = c(mapCol1,"null");
        mapCol2 = c(mapCol2,"trange");
    }
    checkVars = c("year", "month", "day", "rain", "tmin");
    if (Trange) {
        checkVars = c(checkVars, "trange");
    }
    for (cv in checkVars) {
        if (! cv %in% mapCol2) {
            stop(paste("[WACSdata] column '",sep="",cv),sep="", "' missing");
        }
    }
    return(data.frame(names_data=mapCol1, wacs_names=mapCol2));
}

# Internal function that gets the wacs variable name from a column name of 
# original data
# 
# orig_var   The original variable name to get
# mapping    The mapping a returned by wacs.buildMapping
# 
# The wacs variable corresponding to var
#
wacs.getVarWacs = function(orig_var, mapping)
{
    i = which(mapping[,1] == orig_var);
    return(as.character(mapping[i,2]));
}

# Internal function that gets the original variable name from a wacs 
# variable name
# 
# wacs_var   The wacs variable name to get
# mapping    The mapping a returned by wacs.buildMapping
# 
# The wacs variable corresponding to var
#
wacs.getVarOrig = function(wacs_var, mapping)
{
    i = which(mapping[,2] == wacs_var);
    return(as.character(mapping[i,1]));
}


# Internal function that renames the original data columns according
# the WACS model expectations.
#
# data    The original data given by the user
# mapping The mapping of var names built with wacs.buildMapping
# 
# The data with columns renamed
# 

wacs.renameData = function(data, mapping)
{
    nd = names(data);
    nd = sub(wacs.getVarOrig("year",mapping), "year", nd);
    nd = sub(wacs.getVarOrig("month",mapping), "month", nd);
    nd = sub(wacs.getVarOrig("day",mapping), "day", nd);
    nd = sub(wacs.getVarOrig("rain",mapping), "rain", nd);
    nd = sub(wacs.getVarOrig("tmin",mapping), "tmin", nd);
    nd = sub(wacs.getVarOrig("tmax",mapping), "tmax", nd);
    names(data) = nd;
    return(data);
}

# Internal function that reorders columns according WACS model expectations.
#
# data    The original data given by the user
# mapping The mapping of var names built with wacs.buildMapping
# Trange TRUE if one uses the Trange option
# 
# The data with columns renamed
# 

wacs.reorderData = function(data, mapping, Trange)
{
    if (Trange) {
        data$trange = data$tmax - data$tmin;
        vars = c("year", "month", "day", "season", "rain", "tmin", "trange");
        data=data[c(vars, setdiff(names(data), vars))];
        data$tmax=NULL;
    } else {
        vars = c("year", "month", "day", "season", "rain", "tmin", "tmax");
        data=data[c(vars, setdiff(names(data), vars))];
    }
    return (data);
}


# Internal function that computes begin and end index of table
# 
# subdata Data (only date variables) to check
# from    Begin date
# to      End date
# 
# vector of begin and end dates to keep into data
# 

wacs.getTemporal =function(subdata, from, to) # Never called
{
    index_begin = 1;
    index_end = nrow(subdata);
    if (!is.null(from) | !is.null(to)) {
        
    } else {
        y = NULL;
        d = NULL;
        for (i in 1:nrow(subdata)) {
            y = subdata$year[i];
            d = as.Date(ISOdate(y, subdata$month[i], subdata$day[i]));
            if (! is.null(from)) {
                if (as.Date(from) == d) {
                    index_begin = i;
                }
            }
            if (! is.null(to)) {
                if (as.Date(to) == d) {
                    index_end = i;
                }
            }
        }
    }
    if (index_end <= index_begin) {
        stop ("[WACSdata] stop param 'from' is greater than 'to' ");
    }
    return (c(index_begin,index_end));
}

# Internal function that checks or computes the bounds 
# 
# data       The dataframe to study
# bounds     The bounds given by the user 
# mapping    The name mapping for special vars
# The bounds or sends an error 
# 
wacs.getBounds = function(data, bounds, mapping)
{
    var = setdiff(names(data), c("year","month","day","season"))
    res = data.frame(matrix(rep(0, 2*length(var)),nrow=2));
    names(res) <- var;
    has_trange = FALSE;
    for (c in var) {
        done_min =FALSE;
        done_max =FALSE;
        cmap=c;
        if (c =="rain"){
            cmap=wacs.getVarOrig("rain", mapping);
            res[1,c] = 0
            done_min = TRUE
        }
        if (c =="tmin"){
            cmap=wacs.getVarOrig("tmin", mapping);
        }
        if (c =="tmax"){
            cmap=wacs.getVarOrig("tmax", mapping);
        }
        if (c =="trange"){
            res[1,c] = 0;
            done_min = TRUE;
            has_trange = TRUE;
        } 
        if (c =="RG"){
          res[1,c] = 0
          done_min = TRUE
        }
        if (c =="ETPP"){
          res[1,c] = 0
          done_min = TRUE
        }   
        if (c =="V"){
          res[1,c] = 0
          done_min = TRUE
        }
        else {
            if((! is.null(bounds)) && exists(cmap, bounds)) {
                if (exists("min", bounds[[cmap]])) {
                    minv = bounds[[cmap]]$min;
                    if ((c =="rain") && (minv <0)){
                        stop ("[WACSdata] found a min bound < 0 for 'rain'");
                    }
                    res[1,c] = minv;
                    done_min = TRUE;
                }
                if (exists("max", bounds[[cmap]])) {
                    res[2,c] = bounds[[cmap]]$max;
                    done_max = TRUE;
                }
            }
        }
        if (! done_min) {
            minv  = min(data[[c]])
            delta = sort(data[[c]])[10] - minv
            res[1,c] = minv - delta
        }
        if (! done_max) {
            maxv  = max(data[[c]])
            delta = maxv - sort(data[[c]], decreasing=TRUE)[10]
            res[2,c] = maxv + delta
        }        
    }
    ##add Tmax if trange
    if (has_trange){
        done_min =FALSE;
        done_max =FALSE;
        cmap=wacs.getVarOrig("tmax", mapping);
        if((! is.null(bounds)) && exists(cmap, bounds)) {
            if (exists("min", bounds[[cmap]])) {
                res[1,"tmax"] = bounds[[cmap]]$min;
                done_min = TRUE;
            }
            if (exists("max", bounds[[cmap]])) {
                res[2,"tmax"] = bounds[[cmap]]$max;
                done_max = TRUE;
            }
        }
        zz = data[["tmin"]] + data[["trange"]]
        if (! done_min) {
            minv  = min(zz)
            delta = sort(zz)[10] - minv
            res[1,"tmax"] = minv - delta
        }
        if (! done_max) {
            maxv = max(zz)
            delta = maxv - sort(zz,decreasing=TRUE)[10]
            res[2,"tmax"] = maxv + delta
        }
    }
    return (res);
}
