report = function(x, adj.method = "none", cutoff = ifelse(adj.method == "none", 0.01, 0.05), 
          template.file = paste(system.file(package = "CePa"), "/extdata/cepa.template", sep=""),
          only.sig = TRUE, dir.path = NULL, ...) {
    
    if(is.null(dir.path)) {
        dir.path = paste("cepa.report", gsub("[ :]", "-", as.character(Sys.time())), sep=".")
    }
    dir.create(dir.path, showWarnings=FALSE)
    
    file.copy(paste(system.file(package = "CePa"), "/extdata/js", sep=""), dir.path, recursive = TRUE)
    file.copy(paste(system.file(package = "CePa"), "/extdata/swf", sep=""), dir.path, recursive = TRUE)
    
    f1.path = paste(dir.path, "/f1.png", sep="")
    png(f1.path, width=800, height=200)
    plot(x, adj.method=adj.method, only.sig=FALSE, cutoff=cutoff)
    dev.off()
    
    p.value = p.table(x)
    np = ncol(p.value)
    
    p2 = p.value
    p2 = apply(p.value, 2, p.adjust, adj.method)
    f2.path = paste(dir.path, "/f2.png", sep="")
    pathway.name = rownames(p2)
    max.length = max(nchar(pathway.name))
    png(f2.path, width=800, height=400)
    plot(x, adj.method=adj.method, only.sig=TRUE, cutoff=cutoff)
    dev.off()

    # generate images for each pathway
    dir.create(paste(dir.path, "/image", sep=""), showWarnings=FALSE)
    dir.create(paste(dir.path, "/xml", sep=""), showWarnings=FALSE)
    cen = names(x[[1]])
    pathway.name = names(x)
    for(i in 1:length(pathway.name)) {
        if(only.sig){
            if(any(p2[i, ] < cutoff)) {
                cat("  generate images for", pathway.name[i], "...\n")
                for(ce in cen) {
                    image.path = paste(dir.path, "/image/", pathway.name[i], "-", ce, "-graph.png", sep="")
                    png(image.path, width=800, height=750)
                    gg = plot(x, pathway.name[i], ce, ...)
                    dev.off()
                    
                    image.path = paste(dir.path, "/image/", pathway.name[i], "-", ce, "-null.png", sep="")
                    png(image.path, width=800, height=500)
                    plot(x, pathway.name[i], ce, type="null", ...)
                    dev.off()
                    
                    write.graph(gg, file = paste(dir.path, "/xml/", pathway.name[i], "-", ce, ".xml", sep=""), format = "graphml")
                }
            }
        }
        else {
            cat("  generate images for", pathway.name[i], "...\n")
            for(ce in cen) {
                image.path = paste(dir.path, "/image/", pathway.name[i], "-", ce, "-graph.png", sep="")
                png(image.path, width=1200, height=800)
                gg = plot(x, pathway.name[i], ce, ...)
                dev.off()
                
                image.path = paste(dir.path, "/image/", pathway.name[i], "-", ce, "-null.png", sep="")
                png(image.path, width=1000, height=500)
                plot(x, pathway.name[i], ce, type="null", ...)
                dev.off()
                
                write.graph(gg, file = paste(dir.path, "/xml/", pathway.name[i], "-", ce, ".xml", sep=""), format = "graphml")
            }
        }
    }
    
    cat("\n")
    cat("  generate summary in HTML ...\n")
    
    write.table(p2, file = paste(dir.path, "/significance-by-", adj.method, ".txt", sep=""), quote=FALSE, sep = "\t")
    
    replacement = list(x = x,
                       f1.path = "f1.png",
                       f2.path = "f2.png",
                       sig.path = paste("significance-by-", adj.method, ".txt", sep=""),
                       adj.method = adj.method,
                       cutoff = cutoff,
                       only.sig = only.sig,
                       procedure = ifelse(is.ora(x), "Over-representation Analysis", "Gene-Set Analysis"))
    
    tt = readLines(template.file, n = -1)
    html = sapply(tt, function(x) template(x, replacement))
    
    output.file = "index.html"
    output.file = gsub("[:<>|*/\\?]", "", output.file)
    output.file = gsub(" +", "-", output.file)
    writeLines(html, paste(dir.path, "/", output.file, sep=""))
    
    cat(paste("\nVisit here: ", getwd(), "/", dir.path, "/", output.file, "\n\n", sep=""))
}

table.head = function(x, ...) {
    if(is.ora(x)) {
        h = c("Pathway", "Differential nodes", "All nodes", "Differential genes", "All genes", names(x[[1]]))
    } else {
        h = c("Pathway", names(x[[1]]))
    }
    return(h)
}

table.content = function(x, adj.method="none", cutoff = ifelse(adj.method == "none", 0.01, 0.05), only.sig = TRUE, ...) {
    p.value = apply(p.table(x), 2, p.adjust, adj.method)
    
    if(is.ora(x)) {
        count = t(sapply(x, function(x) x[[1]]$count))
    }
    
    pathway.name = names(x)
    cen = names(x[[1]])
    np = nrow(p.value)
    l = apply(p.value, 1, function(x) sum(x < cutoff) > 0)

    oa = 1:np
    i1 = oa[l]
    i2 = oa[!l]
    o1 = order(apply(-log(p.value[l,,drop=FALSE] + 1e-8), 1, mean), decreasing = TRUE)
    o2 = order(apply(-log(p.value[!l,,drop=FALSE] + 1e-8), 1, mean), decreasing = TRUE)
    o = c(i1[o1], i2[o2])
    
    pathway.name = pathway.name[o]
    if(is.ora(x)) {
        count = count[o, ]
    }
    l = l[o]
    p.value = p.value[o,,drop=FALSE]
    p.value = apply(p.value, 2, round, 4)
    p.text = p.value
    
    if(is.ora(x)) {
        for(i in 1:np) {
            for(j in 1:dim(count)[2]) {
                count[i, j] = paste("<td>", count[i, j], "</td>", sep = "")
            }
        }
    }
    for(i in 1:np) {
        for(j in 1:length(cen)) {
        
            if(p.value[i, j] < cutoff) {
                p.text[i, j] = paste("<strong>", p.text[i, j], "</strong>", sep = "")
            }
            
            if(l[i] || only.sig == FALSE) {
                p.text[i, j] = paste("<td><a href='#' onclick=\"getGraph('", pathway.name[i], "', '", cen[j], "');return false;\" >", p.text[i, j], "</a></td>", sep = "")
            }
            else {
                p.text[i, j] = paste("<td>", p.text[i, j], "</td>", sep = "")
            }
        }
    }
    p.text.row = apply(p.text, 1, function(x) paste(x, sep = "", collapse = ""))
    if(is.ora(x)) {
        count.row = apply(count, 1, function(x) paste(x, sep = "", collapse = ""))
    }
    tr.row = sapply(l, function(x) ifelse(x, "<tr class='significant-pathway'>", "<tr>"))
    if(is.ora(x)) {
        return(paste(tr.row, "<td>", pathway.name, "</td>", count.row, p.text.row, "</tr>", sep = ""))
    } else {
        return(paste(tr.row, "<td>", pathway.name, "</td>", p.text.row, "</tr>", sep = ""))
    }
}


template = function(text, replacement, code.pattern = NULL) {

    if (is.null(code.pattern)) {
        code.pattern = "\\@\\{CODE\\}"
    }
    if(length(text) != 1) {
        stop("length of the text should be 1.")
    }

    lines = strsplit(text, "\n")[[1]]
    if(length(lines) == 0) {
        lines = ""
    }
    newlines = character(length(lines))

    # import variables in replacement
    attach(replacement, warn.conflicts = FALSE)

    for (i in 1:length(lines)) {

        # check wether there are code replacements
        code = find_code(code.pattern, lines[i])
        code.template = code[[1]]
        code.variable = code[[2]]

        if(length(code.template)) {

            # if there is code replacement
            # replace the code with its value
            code.result = lapply(code.variable, function(code) eval(parse(text = code)))

            # length of the return value
            v.lines = sapply(code.result, function(x) length(x))

            if(max(v.lines) > 1) {
                current.line = rep(lines[i], max(v.lines))
                for(ai in 1:max(v.lines)) {
                    for(iv in 1:length(code.template)) {
                        current.line[ai] = gsub(code.template[iv],
                        code.result[[iv]][(ai-1) %% length(code.result[[iv]]) + 1],
                        current.line[ai], fixed = TRUE)
                    }
                }
                newlines[i] = paste(current.line, collapse = "\n")
            }
            else if(max(v.lines == 1)) {
                current.line = lines[i]
                for(iv in 1:length(code.template)) {
                    current.line = gsub(code.template[iv], code.result[[iv]],
                                   current.line, fixed = TRUE)
                }
                newlines[i] = current.line
            }
            else {
                newlines[i] = ""
            }
        }
        else {
            newlines[i] = lines[i]
        }
    }

    detach(replacement)

    return(paste(newlines, collapse="\n"))
}

find_code = function(m, text) {

    if(length(text) != 1) {
        stop("text must be length of 1.")
    }

    m2 = gsub("CODE", ".+?", m)

    reg = gregexpr(m2, text, perl = TRUE)[[1]]
    v1 = character(0)
    if(reg[1] > -1) {
        v1 = sapply(1:length(reg), function(i)substr(text, as.numeric(reg)[i], as.numeric(reg)[i]+ attr(reg, "match.length")[i] - 1))
    }
    edge = strsplit(m, "CODE")[[1]]
    v2 = gsub(paste("^", edge[1], "|", edge[2], "$", sep=""), "", v1)
    return(list(template=v1, variable=v2))
}
