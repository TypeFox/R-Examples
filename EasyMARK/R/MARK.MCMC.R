MARK.MCMC <-
function (ch, cov, n.iter = 100, burn.in = 50, number.of.models = 10, 
    n.chains = 2, add = TRUE, quad = TRUE, corr = TRUE) 
{
    # require(stringr)
    # require(rjags)
    # load.module("glm")
    # require(coda)
    # require(foreach)
    # require(doParallel)
    # require(random)
    # load.module("lecuyer")
    
	cores = detectCores()
    if (n.chains > cores) {
        print(paste("WARNING: number of chains exceeds the number of cores. Running", 
            cores, " chains instead"))
        if (cores != 1) 
            n.chains = cores - 1
        else n.chains = 1
    }
    if (!(is.matrix(ch))) {
        cat("WARNING: ch should be a matrix \n Attempting to convert to matrix...\n")
        if (is.data.frame(ch) && ncol(ch) == 1) {
            ch = try(t(sapply(strsplit(ch[, 1], ""), function(x) as.numeric(x))))
            cat("done\n")
        }
        else stop("WARNING: could not covert into matrix")
    }
    if (n.iter <= burn.in) 
        stop("WARNING: n.iter should be higher than the burn.in")
    if (nrow(ch) != nrow(cov)) 
        stop("WARNING: number of rows of cov and ch should be equal")
    classes = apply(cov, 2, class)
    if (!(any(classes == "numeric"))) 
        stop("WARNING: all columns in the cov data.frame should be numeric")
    if (add == FALSE && quad == FALSE && corr == FALSE) 
        stop("WARNING: sorry you cannot test nothing, choose one add, quad, or corr to be true")
    if (corr == TRUE && ncol(cov) < 2) 
        stop("WARNING: sorry you need at least 2 covariates terms for an interaction, set corr == FALSE")
    if (TRUE) {
        N <- dim(ch)[[1]]
        Years <- dim(ch)[[2]]
        First <- NULL
        for (i in 1:N) {
            temp <- 1:Years
            First <- c(First, min(temp[ch[i, ] == 1]))
        }
        Xinit <- matrix(NA, nrow = N, ncol = Years)
        for (i in 1:N) {
            for (j in 1:(Years)) {
                if (j > First[i]) 
                  Xinit[i, j] <- 1
            }
        }
    }
    datax <- list(N = N, Years = Years, dat = as.matrix(ch), 
        First = First, cov = cov)
    covariates = colnames(cov)
    Model.Creator = function(covariates) {
        if (TRUE) {
            model = "\n\tmodel\n\t{\n\tw[1] ~ dbern(0.5) w[2] ~ dbern(0.5) w[3] ~ dbern(0.5) w[4] ~ dbern(0.5) w[5] ~ dbern(0.5)\n\tfor (i in 1:N)  # for each individual\n\t{\n\tphi[i] <- 1/(1+exp(-(mu + w[1] * beta[1] * cov[i,1] + w[2] * beta[2] * cov[i,2] + w[3] * beta[3] * pow(cov[i,1],2) + w[4] * beta[4] * pow(cov[i,2],2) + w[5] * beta[5] * cov[i,1] * cov[i,2] + epsilon[i])))\n\tepsilon[i] <- xi*eta[i]\n\teta[i] ~ dnorm (0, tau.eta) # hierarchical model for theta\n\t# '1' means it is alive the first time it was seen\n\talive[i, First[i]] ~ dbern(1)   \n\t# for each year after the first\n\tfor (j in (First[i]+1):Years)\n\t{ \n\t# state equation\n\talivep[i,j] <- phi[i] * alive[i, j-1]\n\talive[i,j] ~ dbern(alivep[i,j])\n\t# observation equation\n\tsightp[i,j] <- p * alive[i, j] \n\tdat[i, j] ~ dbern(sightp[i,j])\n\t}\n\t}\n\t############# PRIORS\n\tmu ~ dnorm(0,1)\n\tfor (j in 1:5){beta[j] ~ dnorm(0,1)}\n\tp ~ dunif(0,1)\n\t############ half cauchy, voir Gelman 2004 dans Bayesian Analysis\n\tprior.scale <- 1\n\txi ~ dnorm (0, tau.xi)\n\ttau.xi <- pow(prior.scale, -2)\n\ttau.eta ~ dgamma (.5, .5) # chi^2 with 1 d.f.\n\tsigmaeps <- abs(xi)/sqrt(tau.eta) # cauchy = normal/sqrt(chi^2)\n\t############ half cauchy\n\t}\n\t"
        }
        model = gsub(pattern = "\n", replacement = "\n__", model)
        model = unlist(strsplit(model, "\n"))
        model = gsub(pattern = "__", replacement = "\n", model)
        if (length(covariates) > 1) {
            combinations = length(combn(1:length(covariates), 
                m = 2, simplify = FALSE))
            num.vector = c(length(covariates), length(covariates), 
                combinations)
            bool = c(add, quad, corr)
            num.elements = sum(num.vector[bool])
            counter = 1:num.elements
            form.vector = c()
            if (add == TRUE) {
                insert = 1:length(covariates)
                lin.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * cov[i,", insert, "]")
                lin.form = paste(lin.vector, collapse = " + ")
                form.vector[1] = lin.form
            }
            if (quad == TRUE && add == TRUE) {
                insert = (length(covariates) + 1):(length(covariates) * 
                  2)
                quad.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * pow(cov[i,", 1:length(covariates), 
                  "],2)")
                quad.form = paste(quad.vector, collapse = " + ")
                form.vector[2] = quad.form
            }
            if (quad == TRUE && add == FALSE) {
                insert = 1:length(covariates)
                quad.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * pow(cov[i,", 1:length(covariates), 
                  "],2)")
                quad.form = paste(quad.vector, collapse = " + ")
                form.vector[1] = quad.form
            }
            if (corr == TRUE && add == TRUE && quad == TRUE) {
                pairwise.matrix = apply(combn(1:length(covariates), 
                  2), 2, c)
                insert = (length(covariates) * 2 + 1):num.elements
                corr.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * cov[i,", pairwise.matrix[1, ], 
                  "] * cov[i,", pairwise.matrix[2, ], "]")
                corr.form = paste(corr.vector, collapse = " + ")
                form.vector[3] = corr.form
            }
            if (corr == TRUE && add == TRUE && quad == FALSE) {
                pairwise.matrix = apply(combn(1:length(covariates), 
                  2), 2, c)
                insert = (length(covariates) + 1):num.elements
                corr.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * cov[i,", pairwise.matrix[1, ], 
                  "] * cov[i,", pairwise.matrix[2, ], "]")
                corr.form = paste(corr.vector, collapse = " + ")
                form.vector[2] = corr.form
            }
            if (corr == TRUE && add == FALSE && quad == FALSE) {
                pairwise.matrix = apply(combn(1:length(covariates), 
                  2), 2, c)
                insert = 1:num.elements
                corr.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * cov[i,", pairwise.matrix[1, ], 
                  "] * cov[i,", pairwise.matrix[2, ], "]")
                corr.form = paste(corr.vector, collapse = " + ")
                form.vector[1] = corr.form
            }
            if (corr == TRUE && add == FALSE && quad == TRUE) {
                pairwise.matrix = apply(combn(1:length(covariates), 
                  2), 2, c)
                insert = (length(covariates) + 1):num.elements
                corr.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * cov[i,", pairwise.matrix[1, ], 
                  "] * cov[i,", pairwise.matrix[2, ], "]")
                corr.form = paste(corr.vector, collapse = " + ")
                form.vector[2] = corr.form
            }
            full.form = paste("\nphi[i] <- 1/(1+exp(-(mu", paste(form.vector, 
                collapse = " + "), "epsilon[i])))", sep = " + ")
            model[7] = full.form
            bool_counter = unlist(counter[bool])
            w.vector = paste("w[", 1:num.elements, "]", "~ dbern(0.5)")
            model[4] = paste("\n", paste(w.vector, collapse = " "))
            beta.vector = paste("beta[", 1:num.elements, "] ~ dnorm(0,1)")
            model[25] = paste("\n", paste(beta.vector, collapse = " "))
        }
        else {
            num.vector = c(length(covariates), length(covariates))
            bool = c(add, quad)
            num.elements = sum(num.vector[bool])
            form.vector = c()
            if (add == TRUE) {
                insert = 1:length(covariates)
                lin.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * cov[i,", insert, "]")
                lin.form = paste(lin.vector, collapse = " + ")
                form.vector[1] = lin.form
            }
            if (quad == TRUE && add == TRUE) {
                insert = (length(covariates) + 1):(length(covariates) * 
                  2)
                quad.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * pow(cov[i,", 1:length(covariates), 
                  "],2)")
                quad.form = paste(quad.vector, collapse = " + ")
                form.vector[2] = quad.form
            }
            if (quad == TRUE && add == FALSE) {
                insert = 1:length(covariates)
                quad.vector = paste("w[", insert, "] * beta[", 
                  insert, "] * pow(cov[i,", 1:length(covariates), 
                  "],2)")
                quad.form = paste(quad.vector, collapse = " + ")
                form.vector[1] = quad.form
            }
            full.form = paste("\nphi[i] <- 1/(1+exp(-(mu", paste(form.vector, 
                collapse = " + "), "epsilon[i])))", sep = " + ")
            model[7] = full.form
            w.vector = paste("w[", 1:num.elements, "]", "~ dbern(0.5)")
            model[4] = paste("\n", paste(w.vector, collapse = " "))
            beta.vector = paste("beta[", 1:num.elements, "] ~ dnorm(0,1)")
            model[25] = paste("\n", paste(beta.vector, collapse = " "))
        }
        return(list(model = as.character(paste(model, collapse = "")), 
            num.elements = num.elements))
    }
    model_list = Model.Creator(covariates)
    model = model_list$model
    num.elements = model_list$num.elements
    if (TRUE) {
        deb = Sys.time()
        inits <- function() {
            list(p = 0.2, alive = as.matrix(Xinit), mu = 0.5, 
                w = rep(0, num.elements), .RNG.name = "lecuyer::RngStream", 
                .RNG.seed = as.numeric(randomNumbers(n = 1, min = 1, 
                  max = 1e+06, col = 1)))
        }
        cl <- makePSOCKcluster(n.chains)
        clusterSetRNGStream(cl)
        registerDoParallel(cl)
        mcmc <- foreach(i = 1:n.chains, .packages = c("rjags", 
            "random", "coda"), .multicombine = TRUE) %dopar% 
            {
                load.module("lecuyer")
                model.jags <- jags.model(file = textConnection(model), 
                  datax, inits = inits(), n.chains = 1, n.adapt = n.iter, 
                  quiet = TRUE)
                result <- coda.samples(model.jags, c("p", "sigmaeps", 
                  "mu", "beta", "w"), n.iter, thin = 1)[[1]]
                return(result)
            }
        stopCluster(cl)
        fin = Sys.time()
        duration = (fin - deb)
        print(paste("time taken in minutes: ", duration))
        mcmc = lapply(mcmc, function(x) x[burn.in:n.iter, ])
        obj = mcmc
        obj = lapply(obj, function(x) apply(x, 2, mcmc))
        obj = lapply(obj, function(x) mcmc(x))
        obj = mcmc.list(obj)
        gelman = try(gelman.diag(obj))
        mcmc.list = obj
        mcmc = do.call(rbind, mcmc)
    }
    gradients = apply(mcmc, 2, mean)[1:num.elements]
    p = mean(mcmc[, "p"])
    formula.Maker = function(covariates) {
        lin.form = covariates
        quad.form = paste(covariates, "^2", sep = "")
        if (length(covariates) > 1) {
            corr.form = apply(format(combn(covariates, 2)), 2, 
                paste, collapse = "*")
            bool = c(add, quad, corr)
            formula.elements = list(lin.form, quad.form, corr.form)
            formula.elements = unlist(formula.elements[bool])
        }
        else {
            formula.elements = list(lin.form, quad.form)
            bool = c(add, quad)
            formula.elements = unlist(formula.elements[bool])
        }
        return(formula.elements)
    }
    formula.elements = formula.Maker(covariates)
    estimates = data.frame(covariate = formula.elements, gradients = gradients)
    Model.Selection = function() {
        w = mcmc[, (num.elements + 4):ncol(mcmc)]
        if (num.elements > 1) {
            w.string = apply(w, 1, paste, collapse = "")
        }
        else {
            w.string = as.character(w)
        }
        w.models = unique(w.string)
        pp = c()
        total = length(w.string)
        for (i in 1:length(w.models)) {
            pp[i] = length(which(w.models[i] == w.string))/total
        }
        form = formula.elements
        idx = str_locate_all(string = w.models, pattern = "1")
        ff.list = list()
        for (i in 1:length(idx)) {
            focal = idx[[i]]
            f.list = c()
            for (j in 1:nrow(focal)) {
                foc_idx = focal[j]
                f.list[j] = form[foc_idx]
            }
            f.list = paste(f.list, collapse = "+")
            ff.list[i] = f.list
        }
        formulas = unlist(ff.list)
        pp.results = data.frame(model = formulas, probability = pp)
        pp.results = pp.results[order(pp.results$probability, 
            decreasing = TRUE), ]
        pp.results = pp.results[1:number.of.models, ]
        return(pp.results)
    }
    pp.results = Model.Selection()
    pp.results = as.data.frame(lapply(pp.results, function(x) if (is.character(x) | 
        is.factor(x)) 
        gsub("NA", "intercept", x)
    else x))
    output = list(mcmc = mcmc, mcmc.list = mcmc.list, pp = pp.results, 
        estimates = estimates, p = p, gelman = gelman)
    return(output)
}
