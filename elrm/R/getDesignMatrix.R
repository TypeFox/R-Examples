`getDesignMatrix` <-
function(formula, interest, dataset)
{
    ## options(contrasts=c("contr.treatment","contr.poly"));

    yCol = as.character(as.list(as.formula(formula))[[2]])[2];
    mCol = as.character(as.list(as.formula(formula))[[2]])[3];

    response =  as.character(as.vector(formula)[[2]]);

    if(!('/' %in% response))
    {
        stop("model formula must be specified in the form y/n~x");
    }

    if(!(yCol %in% names(dataset)))
    {
        stop("binomial response variable, '", yCol,"', was not found");
    }

    if(!(mCol %in% names(dataset)))
    {
        stop("weights variable, '", mCol,"', was not found");
    }

    design.matrix = model.matrix(as.formula(formula),data=dataset);

    design.names = attr(terms.formula(as.formula(formula)),"term.labels");

    z.names = attr(terms.formula(as.formula(interest)),"term.labels");

    match.names = match(z.names,design.names,nomatch=-1);

    if(mCol %in% design.names)
    {
        stop("weights variable, '",mCol,"', cannot appear as a covariate in the formula");
    }

    if(yCol %in% design.names)
    {
        stop("response variable, '",yCol,"', cannot appear as a covariate in the formula");
    }
    
    if(-1 %in% match.names)
    {
        stop("the 'term.labels' attribute of 'terms.formula(interest)' must match those found in 'terms.formula(formula)'");
    }

    design.frame = as.data.frame(design.matrix);
    
    zCols = (1:(ncol(design.frame)))[pmatch(as.character(attr(design.matrix,"assign")),as.character(match.names),duplicates.ok=T,nomatch=-1)!=-1];
    wCols = (1:(ncol(design.frame)))[-zCols];

    design.matrix = as.data.frame(cbind(dataset[yCol],dataset[mCol],design.frame));

    yCol = 1;
    mCol = 2;

    zCols = zCols+2;
    wCols = wCols+2;

    result = list(design.matrix,yCol,mCol,zCols,wCols);

    names(result) = c("design.matrix","yCol","mCol","zCols","wCols");

    return(result);
}

