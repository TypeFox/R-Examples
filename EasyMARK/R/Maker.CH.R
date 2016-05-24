Maker.CH <-
function (dates = dates, id = id, date.format = "%Y_%m_%d") 
{
    ids = as.character(unique(id))
    dates = as.Date(dates, date.format)
    D = data.frame(id = id, dates = dates)
    Start = D[order(D["dates"]), ][1, "dates"]
    End = D[order(D["dates"]), ][nrow(D), "dates"]
    StudyLength = as.numeric(abs(difftime(End, Start, units = "days"))) + 
        1
    DataList = list()
    for (i in 1:length(ids)) {
        id_iter = ids[i]
        idD = subset(D, id == id_iter)
        idD = idD[!duplicated(idD[c("dates")]), ]
        idD = idD[order(idD["dates"]), ]
        DataList[[i]] = idD
    }
    StartData = list()
    for (i in 1:length(DataList)) {
        FirstDay = DataList[[i]][1, "dates"]
        StartDifference = as.numeric(abs(difftime(Start, FirstDay, 
            units = "days")))
        id = DataList[[i]][1, "id"]
        StartData[[i]] = data.frame(id = id, StartDifference = StartDifference, 
            FirstDay = FirstDay)
    }
    DifferenceList = list()
    for (id_iter in 1:length(DataList)) {
        idD = DataList[[id_iter]]
        if (nrow(idD) > 1) {
            Differences = c()
            for (i in 1:nrow(idD) - 1) {
                date1 = idD[i, "dates"]
                date2 = idD[i + 1, "dates"]
                Differences[i] = as.numeric(abs(difftime(date1, 
                  date2, units = "days")))
                id = as.character(idD[1, "id"])
            }
            DifferenceList[[id_iter]] = list(id = id, Differences = Differences)
        }
        else {
            id = as.character(idD[1, "id"])
            DifferenceList[[id_iter]] = list(id = id, Differences = 0)
        }
    }
    chList = list()
    for (id_iter in 1:length(StartData)) {
        idStart = StartData[[id_iter]]
        idDiff = DifferenceList[[id_iter]]
        zeros = list()
        for (i in 1:length(idDiff$Differences)) {
            Differences = idDiff$Differences
            if (!(0 %in% Differences)) 
                zeros[[i]] = c(rep(0, Differences[i] - 1), 1)
        }
        ch = c(1, do.call("c", zeros))
        FrontZeros = idStart[1, "StartDifference"]
        if (FrontZeros != 0) 
            ch = c(rep(0, FrontZeros), ch)
        BackZeros = StudyLength - length(ch)
        if (BackZeros != 0) 
            ch = c(ch, rep(0, BackZeros))
        chList[[id_iter]] = t(ch)
    }
    StudyDays = seq(Start, End, "days")
    chData = data.frame(do.call("rbind", chList))
    colnames(chData) = StudyDays
    rownames(chData) = ids
    return(chData)
}
