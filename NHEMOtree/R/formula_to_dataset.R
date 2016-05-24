formula_to_dataset <-
function(formula, data){

  # Find data columns according to formula
    Spalten<- all.vars(formula)

    if (Spalten[2]=="."){
       (temp<- which(names(data)==Spalten[1]))
        Spalten<- c(Spalten[1], names(data)[-temp])
    }

    if (length(intersect(names(data), Spalten))!=length(Spalten)){ 
        stop("Input in formula and data do not match!")
    }else{

  # Create working dataset
    Spalten_num<- c()
    for (i in 1: length(Spalten)) Spalten_num[i]<- which(names(data)==Spalten[i])
    data_work<- data[,Spalten_num]
    }

  # Output
    return(data_work)
}
