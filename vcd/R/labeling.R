################################################################
## labeling

pexpand <- function(par, len, default_value, default_names, choices = NULL) {
  if (is.null(par))
      par <- default_value
  nam <- names(par)
  if (!is.null(choices))
      par <- sapply(par, match.arg, choices)
  if (is.null(nam)) {
      default_value <- par
      par <- rep(par, length.out = len)
      nam <- names(par) <- default_names
  } else if (length(nam[nam == ""])) {
      default_value <- par[nam == ""]
      nam <- nam[nam != ""]
  }
  ret <- rep(default_value, length.out = len)
  if (!is.null(nam)) {
      names(ret) <- default_names
      ret[nam] <- par[nam]
  }
  ret
}

labeling_list <- function(gp_text = gpar(),
                          just = "left",
                          pos = "left",
                          lsep = ": ", sep = " ",
                          offset = unit(c(2, 2), "lines"),
                          varnames = TRUE,
                          cols = 2,
                          ...) {
  function(d, split_vertical, condvars, prefix = "") {
    if (is.table(d))
      d <- dimnames(d)
    ld <- length(d)
    labeling_border(labels = FALSE, varnames = varnames)(d, split_vertical, condvars, prefix)
    seekViewport(paste(prefix, "margin_bottom", sep = ""))
    pos <- unit(switch(pos, left = 0, center = 0.5, 1) / cols, "npc")
    ind <- split(seq(ld), rep.int(seq(cols), ceiling(ld / cols))[seq(ld)])

    for (i in seq_along(ind))
      grid.text(x = offset[1] + pos + unit((i - 1) / cols, "npc"),
                y = unit(1, "npc") - offset[2],
                paste(names(d[ind[[i]]]),
                      sapply(d[ind[[i]]], paste, collapse = sep),
                      sep = lsep,
                      collapse = "\n"
                      ),
                just = c(just, "top"),
                gp = gp_text
                )
  }
}
class(labeling_list) <- "grapcon_generator"

labeling_conditional <- function(...) {
  function (d, split_vertical, condvars, prefix = "") {
    if (is.table(d))
      d <- dimnames(d)
    v <- rep.int(TRUE, length(d))
    v[seq(condvars)] <- FALSE
    labeling_border(labels = !v, ...)(d, split_vertical, condvars, prefix)
    labeling_cells(labels = v, ...)(d, split_vertical, condvars, prefix)
  }
}
class(labeling_conditional) <- "grapcon_generator"

labeling_cells <- function(labels = TRUE, varnames = TRUE,
                         abbreviate_labels = FALSE, abbreviate_varnames = FALSE,
                         gp_text = gpar(), lsep = ": ", lcollapse = "\n",
                         just = "center", pos = "center", rot = 0,
                         margin = unit(0.5, "lines"), clip_cells = TRUE,
                         text = NULL, ...) {
  function(d, split_vertical, condvars, prefix = "") {
    if (is.table(d))
      d <- dimnames(d)
    dn <- names(d)
    ld <- length(d)

    ## expand parameters
    if (length(pos) < 2) pos <- c(pos, pos)
    labels <- pexpand(labels, ld, TRUE, dn)
    varnames <- pexpand(varnames, ld, TRUE, dn)
    abbreviate_labels <- pexpand(abbreviate_labels, ld, FALSE, dn)
    abbreviate_varnames <- pexpand(abbreviate_varnames, ld, FALSE, dn)

    ## margin
    if (!is.unit(margin))
      margin <- unit(margin, "lines")

    prvars <- ifelse(abbreviate_varnames,
                     sapply(seq_along(dn),
                            function(i) abbreviate(dn[i], abbreviate_varnames[i])),
                     dn)
    prvars <- ifelse(varnames, paste(prvars, lsep, sep = ""), "")

    ## draw labels
    split <- function(vind = 1, labs = c()) {
      n <- d[[vind]]
      for (labind in seq_along(n)) {
        lab <- c(labs, n[labind])
        names(lab) <- names(d)[1:vind]
        mlab <- paste(prefix, "cell:", paste(dn[1:vind], lab, sep = "=", collapse = ","),
                      sep = "")

        if (vind < ld)
          split(vind + 1, lab)
        else {
          seekViewport(mlab)
          pushViewport(viewport(width = max(unit(0, "npc"), unit(1, "npc") - 2 * margin),
                                height = unit(1, "npc") - 2 * margin,
                                clip = clip_cells))
          txt <- if (!is.null(text)) {
            lab <- lab[names(dimnames(text))]
            do.call("[", c(list(text), as.list(lab)))
          } else {
            prlab <- ifelse(abbreviate_labels,
                            sapply(seq_along(lab),
                                   function(i) abbreviate(lab[i], abbreviate_labels[i])),
                            lab)
            prlab <- prlab[labels[1:ld]]
            paste(prvars[labels[1:ld]], prlab, sep = "", collapse = lcollapse)
          }

          grid.text(if(!is.na(txt)) txt,
                    x = switch(pos[1], left =, top = 0, center = 0.5, 1),
                    y = switch(pos[2], left =, top = 1, center = 0.5, 0),
                    gp = gp_text, just = just, rot = rot)
          popViewport()
        }
      }
    }
    split()
    seekViewport(paste(prefix, "base", sep = ""))
    upViewport(1)
  }
}
class(labeling_cells) <- "grapcon_generator"

labeling_border <- function(labels = TRUE, varnames = labels,
                            set_labels = NULL, set_varnames = NULL,
                            tl_labels = NULL, alternate_labels = FALSE,
                            tl_varnames = NULL,
                            gp_labels = gpar(fontsize = 12),
                            gp_varnames = gpar(fontsize = 12, fontface = 2),
                            rot_labels = c(0, 90, 0, 90),
                            rot_varnames = c(0, 90, 0, 90),
                            pos_labels = "center", pos_varnames = "center",
                            just_labels = "center", just_varnames = pos_varnames,
                            boxes = FALSE, fill_boxes = FALSE,
                            offset_labels = c(0, 0, 0, 0),
                            offset_varnames = offset_labels,

                            labbl_varnames = NULL,
                            labels_varnames = FALSE, sep = ": ",

                            abbreviate_labs = FALSE, rep = TRUE,
                            clip = FALSE, ...
                            ) {
    ## expand parameters that apply to the four table margins
    pos_labels <- pexpand(pos_labels, 4, "center", c("top", "right", "bottom", "left"),
                          c("left", "center", "right"))
    just_labels <- pexpand(just_labels, 4, "center", c("top", "right", "bottom", "left"),
                           c("left", "center", "right"))
    offset_varnames <- if (!is.unit(offset_varnames))
      unit(pexpand(offset_varnames, 4,
                   rep.int(0, 4), c("top","right","bottom","left")), "lines")
    else
      rep(offset_varnames, length.out = 4)

    offset_labels <- if (!is.unit(offset_labels))
      unit(pexpand(offset_labels, 4,
                   rep.int(0, 4), c("top","right","bottom","left")), "lines")
    else
      rep(offset_labels, length.out = 4)

    rot_labels <- pexpand(rot_labels, 4, c(0, 90, 0, 90),
                          c("top", "right", "bottom", "left"))

    if (inherits(gp_varnames, "gpar"))
      gp_varnames <- list(gp_varnames)
    gp_varnames <- pexpand(gp_varnames, 4, list(gpar(fontsize = 12, fontface = 2)),
                           c("top", "right", "bottom", "left"))

    rot_varnames <- pexpand(rot_varnames, 4, c(0, 90, 0, 90),
                          c("top", "right", "bottom", "left"))

    pos_varnames <- pexpand(pos_varnames, 4, "center",
                           c("top", "right", "bottom", "left"),
                            c("left", "center", "right"))

    just_varnames <- pexpand(just_varnames, 4, pos_varnames,
                             c("top", "right", "bottom", "left"),
                             c("left", "center", "right"))

  function(d, split_vertical, condvars, prefix = "") {
    if (is.table(d))
      d <- dimnames(d)
    dn <- names(d)
    ld <- length(d)

    ## expand table- (i.e., dimensionality)-dependent parameters
    clip <- pexpand(clip, ld, TRUE, dn)
    labels <- pexpand(labels, ld, TRUE, dn)
    labels_varnames <- pexpand(labels_varnames, ld, FALSE, dn)

    ## tl_labels
    def <- logical()
    def[split_vertical] <- rep(c(TRUE, FALSE), length.out = sum(split_vertical))
    def[!split_vertical] <- rep(c(TRUE, FALSE), length.out = sum(!split_vertical))
    tl_labels <- if (is.null(tl_labels))
      def
    else
      pexpand(tl_labels, ld, def, dn)

    ## rep labels
    rep <- pexpand(rep, ld, TRUE, dn)
    printed <- lapply(d, function(i) rep.int(FALSE, length(i)))

    ## alternate labels
    alternate_labels <- pexpand(alternate_labels, ld, FALSE, dn)

    ## abbreviate
    abbreviate_labs <- pexpand(abbreviate_labs, ld, FALSE, dn)
    labs <- d
    for (i in seq_along(d))
      if (abbreviate_labs[i])
        labs[[i]] <- abbreviate(labs[[i]], abbreviate_labs[i])

    ## gp_labels
    if (inherits(gp_labels, "gpar"))
      gp_labels <- list(gp_labels)
    gp_labels <- pexpand(gp_labels, ld, list(gpar(fontsize = 12)), dn)

    ## varnames
    varnames <- pexpand(varnames, ld, labels, dn)

    ## tl_varnames
    if (is.null(tl_varnames) && is.null(labbl_varnames))
      tl_varnames <- tl_labels
    tl_varnames <- pexpand(tl_varnames, ld, tl_labels, dn)

    ## labbl_varnames
    if (!is.null(labbl_varnames))
      labbl_varnames <- pexpand(labbl_varnames, ld, TRUE, dn)

    ## boxes
    boxes <- pexpand(boxes, ld, FALSE, dn)

    ## fill_boxes
    dnl <- sapply(d, length)
    fill_boxes <- if (is.atomic(fill_boxes)) {
      fill_boxes <- if (is.logical(fill_boxes))
        ifelse(pexpand(fill_boxes, ld, FALSE, dn), "grey", NA)
      else
        pexpand(fill_boxes, ld, "grey", dn)
      col <- rgb2hsv(col2rgb(fill_boxes))
      lapply(seq(along.with = dnl),
             function(i) if (is.na(fill_boxes[i])) "white" else
                         hsv(h = col["h",i],
                             s = col["s",i],
                             v = seq(from = col["v",i],
                                     to = 0.5 * col["v",i],
                                     length = dnl[i])
                             )
             )
    } else {
      fill_boxes <- pexpand(fill_boxes, ld, "white", dn)
      lapply(seq(ld),
             function(i) pexpand(fill_boxes[[i]], dnl[i], "white", d[[i]])
             )
    }


    ## precompute spaces
    lsp <- tsp <- bsp <- rsp <- 0
    labsp <- rep.int(0, ld)
    for (i in seq_along(dn)[tl_labels & labels])
      labsp[i] <- if (split_vertical[i]) {
        if (alternate_labels[i]) bsp <- bsp - 1
        tsp <- tsp + 1
      } else {
        if (alternate_labels[i]) rsp <- rsp + 1
        lsp <- lsp - 1
      }
    for (i in rev(seq_along(dn)[!tl_labels & labels]))
      labsp[i] <- if (split_vertical[i]) {
        if (alternate_labels[i]) tsp <- tsp + 1
        bsp <- bsp - 1
      } else {
        if (alternate_labels[i]) lsp <- lsp - 1
        rsp <- rsp + 1
      }

    if(is.null(labbl_varnames)) {
    ## varnames in the outer margin
      ## compute axis names
      tt <- bt <- lt <- rt <- ""
      for (i in seq_along(dn)) {
        var <- if (!is.null(set_varnames) && !is.na(set_varnames[dn[i]]))
          set_varnames[dn[i]]
        else
          dn[i]
        if (varnames[i]) {
          if (split_vertical[i]) {
            if (tl_varnames[i])
              tt <- paste(tt, var, sep = if (tt == "") "" else " / ")
            else
              bt <- paste(bt, var, sep = if (bt == "") "" else " / ")
          } else {
            if (tl_varnames[i])
              lt <- paste(lt, var, sep = if (lt == "") "" else " / ")
            else
              rt <- paste(rt, var, sep = if (rt == "") "" else " / ")
          }
        }
      }
      ## draw axis names
      if (tt != "")
        grid.text(tt, y = unit(1, "npc") + unit(tsp + 1, "lines") + offset_varnames[1],
                  x = switch(pos_varnames[1], left =, bottom = 0, center =, centre = 0.5, 1),
                  rot = rot_varnames[1], just = just_varnames[1], gp = gp_varnames[[1]])
      if (bt != "")
        grid.text(bt, y = unit(bsp - 1, "lines") + -1 * offset_varnames[3],
                  x = switch(pos_varnames[3], left =, bottom = 0, center =, centre = 0.5, 1),
                  rot = rot_varnames[3], just = just_varnames[3], gp = gp_varnames[[3]])
      if (lt != "")
        grid.text(lt, x = unit(lsp - 1, "lines") + -1 * offset_varnames[4],
                  y = switch(pos_varnames[4], left =, bottom = 0, center =, centre = 0.5, 1),
                  rot = rot_varnames[4], just = just_varnames[4], gp = gp_varnames[[4]])
      if (rt != "")
        grid.text(rt, x = unit(1, "npc") + unit(rsp + 1, "lines") + offset_varnames[2],
                  y = switch(pos_varnames[2], left =, bottom = 0, center =, centre = 0.5, 1),
                  rot = rot_varnames[2], just = just_varnames[2], gp = gp_varnames[[2]])
    } else {
    ## varnames beneath labels
      for (i in seq_along(dn)) {
        var <- if (!is.null(set_varnames) && !is.na(set_varnames[dn[i]]))
          set_varnames[dn[i]]
        else
          dn[i]
        if (varnames[i]) {
          if (split_vertical[i]) {
            if (tl_labels[i]) {
              if (labbl_varnames[i]) {
                grid.text(var,
                          y = unit(1, "npc") + unit(1 + tsp - labsp[i], "lines") + offset_varnames[1],
                          x = unit(-0.5, "lines"),
                          just = "right", gp = gp_varnames[[4]])
              } else {
                grid.text(var, y = unit(1, "npc") + unit(1 + tsp - labsp[i], "lines") + offset_varnames[1],
                          x = unit(1, "npc") + unit(0.5, "lines"),
                          just = "left", gp = gp_varnames[[2]])
              }
            } else {
              if (labbl_varnames[i]) {
                grid.text(var, y = unit(labsp[i], "lines") + -1 * offset_varnames[3],
                          x = unit(-0.5, "lines"), just = "right",
                          gp = gp_varnames[[4]])
              } else {
                grid.text(var, y = unit(labsp[i], "lines") + -1 * offset_varnames[3],
                          x = unit(1, "npc") + unit(0.5, "lines"),
                          just = "left", gp = gp_varnames[[2]])
              }
            }
          } else {
            if (tl_labels[i]) {
              if (labbl_varnames[i]) {
                grid.text(var, x = unit(lsp - 1 - labsp[i], "lines") + -1 * offset_varnames[4],
                          y = unit(-0.5, "lines"), just = "right", rot = 90,
                          gp = gp_varnames[[4]])
              } else {
                grid.text(var, x = unit(lsp - 1 - labsp[i], "lines") + -1 * offset_varnames[4],
                          y = unit(1, "npc") + unit(0.5, "lines"),
                          just = "left", rot = 90, gp = gp_varnames[[2]])
              }
            } else {
              if (labbl_varnames[i]) {
                grid.text(var, x = unit(1, "npc") + unit(labsp[i], "lines") + offset_varnames[2],
                          y = unit(-0.5, "lines"),
                          just = "right", rot = 90, gp = gp_varnames[[4]])
              } else {
                grid.text(var, x = unit(1, "npc") + unit(labsp[i], "lines") + offset_varnames[2],
                          y = unit(1, "npc") + unit(0.5, "lines"),
                          just = "left", rot = 90, gp = gp_varnames[[2]])
              }
            }
          }
        }
      }
    }

    ## draw labels
    split <- function(vind = 1, root = paste(prefix, "cell:", sep = ""),
                      left = TRUE, right = TRUE, top = TRUE, bottom = TRUE) {
      n <- d[[vind]]
      vl <- length(n)
      sp <- split_vertical[vind]
      labseq <- seq_along(n)
      if (!sp) labseq <- rev(labseq)

      for (labind in labseq) {
        mlab <- paste(root, dn[vind], "=", n[labind], sep = "")
        if (labels[vind] && (rep[vind] || !printed[[vind]][labind])) {
          lab <- if (!is.null(set_labels) && !is.null(set_labels[[dn[vind]]]))
            set_labels[[dn[vind]]][labind]
          else
            labs[[vind]][labind]
          if (labels_varnames[vind])
            lab <- if (!is.null(set_varnames) && !is.na(set_varnames[dn[vind]]))
              paste(set_varnames[dn[vind]], lab, sep = sep)
          else
              paste(dn[vind], lab, sep = sep)
          if (sp) {
            if (tl_labels[vind]) {
              if (top) {
                seekViewport(mlab)
                if (clip[vind])
                  pushViewport(viewport(height = unit(1, "npc") + 2 * offset_labels[1] +
                                        unit(2 * (2 + tsp - labsp[vind]), "lines"),
                                        clip = "on"))
                if (boxes[vind])
                  grid.rect(height = unit(0.8, "lines"),
                            y = unit(1, "npc") + offset_labels[1] +
                                unit(1 + tsp - labsp[vind] - (2 + as.numeric(offset_labels[1]) + tsp - labsp[vind]) * clip[vind], "lines"),
                            gp = gpar(fill = fill_boxes[[vind]][labind]))
                grid.text(lab,
                          y = unit(1, "npc") + offset_labels[1] +
                              unit(1 + tsp - labsp[vind] - (2 + as.numeric(offset_labels[1]) + tsp - labsp[vind]) * clip[vind], "lines"),
                          x = unit(0.15 * switch(pos_labels[1], left =, bottom = 1, center =, centre = 0, -1) * boxes[vind], "lines") +
                              unit(switch(pos_labels[1], left =, bottom = 0, center =, centre = 0.5, 1), "npc"),
                          rot = rot_labels[1], just = just_labels[1],
                          gp = gp_labels[[vind]])
                if (clip[vind]) popViewport()
                printed[[vind]][labind] <<- TRUE
              }
            } else {
              if (bottom) {
                seekViewport(mlab)
                if (clip[vind])
                  pushViewport(viewport(height = unit(1, "npc") + 2 * offset_labels[3] +
                                        unit(2 * (1 + abs(labsp[vind])), "lines"),
                                        clip = "on"))
###
                if (boxes[vind])
                  grid.rect(height = unit(0.8, "lines"),
                            y = -1 * offset_labels[3] + unit(labsp[vind] + (1 + as.numeric(offset_labels[3]) + abs(labsp[vind])) * clip[vind], "lines"),
                            gp = gpar(fill = fill_boxes[[vind]][labind]))

                grid.text(lab,
                          y = -1 * offset_labels[3] + unit(labsp[vind] + (1 + as.numeric(offset_labels[3]) + abs(labsp[vind])) * clip[vind], "lines"),
                          x = unit(0.15 * switch(pos_labels[3], left =, bottom = 1, center =, centre = 0, -1) * boxes[vind], "lines") +
                              unit(switch(pos_labels[3], left =, bottom = 0, center =, centre = 0.5, 1), "npc"),
                          rot = rot_labels[3], just = just_labels[3],
                          gp = gp_labels[[vind]])
                if (clip[vind]) popViewport()
                printed[[vind]][labind] <<- TRUE
              }
            }
          } else {
            if (tl_labels[vind]) {
              if (left) {
                seekViewport(mlab)
                if (clip[vind])
                  pushViewport(viewport(width = unit(1, "npc") + 2 * offset_labels[4] +
                                        unit(2 * (2 - lsp + labsp[vind]), "lines"),
                                        clip = "on"))
                if (boxes[vind])
                  grid.rect(width = unit(0.8, "lines"),
                            x = -1 * offset_labels[4] + unit(lsp - 1 - labsp[vind] + (2 - lsp + as.numeric(offset_labels[4]) + labsp[vind]) * clip[vind], "lines"),
                            gp = gpar(fill = fill_boxes[[vind]][labind]))

                grid.text(lab,
                          x = -1 * offset_labels[4] + unit(lsp - 1 - labsp[vind] + (2 - lsp + as.numeric(offset_labels[4]) + labsp[vind]) * clip[vind], "lines"),
                          y = unit(0.15 * switch(pos_labels[4], left =, bottom = 1, centre = 0, -1) * boxes[vind], "lines") +
                          unit(switch(pos_labels[4], left =, bottom = 0, center =, centre = 0.5, 1), "npc"),
                          rot = rot_labels[4], just = just_labels[4],
                          gp = gp_labels[[vind]])

                if (clip[vind]) popViewport()
                printed[[vind]][labind] <<- TRUE
              }
            } else {
              if (right) {
                seekViewport(mlab)
                if (clip[vind])
                  pushViewport(viewport(width = unit(1, "npc") + 2 * offset_labels[2] +
                                        unit(2 * (1 + abs(labsp[vind])), "lines"),
                                        clip = "on"))
                if (boxes[vind])
                  grid.rect(width = unit(0.8, "lines"),
                            x = offset_labels[2] + unit(1, "npc") +
                                unit(labsp[vind] - (1 + as.numeric(offset_labels[2]) + abs(labsp[vind])) * clip[vind], "lines"),
                            gp = gpar(fill = fill_boxes[[vind]][labind]))
                grid.text(lab,
                          x = offset_labels[2] + unit(1, "npc") + unit(0.1, "lines") +
                              unit(labsp[vind] - (1 + as.numeric(offset_labels[2]) + abs(labsp[vind])) * clip[vind], "lines"),
                          y = unit(0.15 * switch(pos_labels[2], left =, bottom = 1, center =, centre = 0, -1) * boxes[vind], "lines") +
                              unit(switch(pos_labels[2], left =, bottom = 0, center =, centre = 0.5, 1), "npc"),
                          rot = rot_labels[2], just = just_labels[2],
                          gp = gp_labels[[vind]])

                if (clip[vind]) popViewport()
                printed[[vind]][labind] <<- TRUE
              }
            }
          }
        }

        if (vind < ld) Recall(vind + 1, paste(mlab, ",", sep = ""),
                              if (sp) left && labind == 1 else left,
                              if (sp) right && labind == vl else right,
                              if (!sp) top && labind == 1 else top,
                              if (!sp) bottom && labind == vl else bottom)
      }
    }
    ## patch for alternating labels, part 1
    if (any(alternate_labels)) {
      ## save set_labels
      set_labels_hold <- set_labels

      ## create vanilla set_labels-object
      set_labels <- d

      ## copy old set_labels
      if (!is.null(set_labels_hold))
        set_labels[names(set_labels_hold)] <- set_labels_hold

      ## mask half of the labels
      for (i in which(alternate_labels))
        if (length(d[[i]]) > 1)
          set_labels[[i]][seq(2, length(d[[i]]), 2)] <- ""
    }

    split()

    ## patch for alternating labels, part 2
    if (any(alternate_labels)) {
      ## create again vanilla set_labels-object
      set_labels <- d

      ## copy again old set_labels
      if (!is.null(set_labels_hold))
        set_labels[names(set_labels_hold)] <- set_labels_hold

      ## clear all non-alternated labels
      labels[!alternate_labels] <- FALSE

      ## mask other half of alternated labels
      for (i in which(alternate_labels))
        set_labels[[i]][seq(1, length(d[[i]]), 2)] <- ""

      ## invert tl_labels and labsp
      tl_labels <- ! tl_labels
      labsp <- -labsp

      ## label again
      split()
    }
    seekViewport(paste(prefix, "base", sep = ""))
    upViewport(1)
  }
}
class(labeling_border) <- "grapcon_generator"

labeling_doubledecker <- function(lab_pos = c("bottom", "top"),
                                  dep_varname = TRUE,
                                  boxes = NULL,
                                  clip = NULL,
                                  labbl_varnames = FALSE,
                                  rot_labels = rep.int(0, 4),
                                  pos_labels = c("left", "center", "left", "center"),
                                  just_labels = c("left", "left", "left", "center"),
                                  varnames = NULL,
                                  gp_varnames = gpar(fontsize = 12, fontface = 2),
                                  offset_varnames = c(0, -0.6, 0, 0),
                                  tl_labels = NULL,
                                  ...) {
  lab_pos <- match.arg(lab_pos)

  if (inherits(gp_varnames, "gpar"))
      gp_varnames <- list(gp_varnames)
  gp_varnames <- pexpand(gp_varnames, 4, list(gpar(fontsize = 12, fontface = 2)), c("top", "right", "bottom", "left"))

  function(d, split_vertical, condvars, prefix = "") {
    if (is.table(d))
      d <- dimnames(d)
    ld <- length(d)
    dn <- names(d)

    ## expand dimension parameters
    boxes <- pexpand(boxes, ld, c(rep.int(TRUE, ld - 1), FALSE), dn)
    clip <- pexpand(clip, ld, c(rep.int(TRUE, ld - 1), FALSE), dn)
    varnames <- pexpand(varnames, ld, c(rep.int(TRUE, ld - 1), FALSE), dn)
    tl_labels <- pexpand(tl_labels, ld, c(rep.int(lab_pos == "top", ld - 1), FALSE), dn)
    if (!is.null(labbl_varnames))
        labbl_varnames <- pexpand(labbl_varnames, ld, FALSE, dn)

    ## expand side parameters
    rot_labels <- pexpand(rot_labels, 4, c(0, 0, 0, 0),
                          c("top", "right", "bottom", "left"))
    pos_labels <- pexpand(pos_labels, 4,
                          c("left", "center", "left", "center"),
                          c("top", "right", "bottom", "left"),
                          c("left", "center", "right"))
    just_labels <- pexpand(just_labels, 4,
                           c("left", "left", "left", "center"),
                           c("top", "right", "bottom", "left"),
                           c("left", "center", "right"))
    offset_varnames <- if (!is.unit(offset_varnames))
      unit(pexpand(offset_varnames, 4,
                   c(0, -0.6, 0, 0),
                   c("top","right","bottom","left")), "lines")
    else
      rep(offset_varnames, length.out = 4)

    labeling_border(boxes = boxes,
                    clip = clip,
                    labbl_varnames = labbl_varnames,
                    rot_labels = rot_labels,
                    pos_labels = pos_labels,
                    just_labels = just_labels,
                    varnames = varnames,
                    gp_varnames = gp_varnames,
                    offset_varnames = offset_varnames,
                    tl_labels = tl_labels,
                    ...
                    )(d, split_vertical, condvars, prefix)
    if (!(is.logical(dep_varname) && !dep_varname)) {
        if (is.null(dep_varname) || is.logical(dep_varname))
            dep_varname <- names(d)[length(d)]
        seekViewport(paste(prefix, "margin_right", sep = ""))
        grid.text(dep_varname,
                  x = unit(0.5, "lines"), y = unit(1, "npc"), just = c("left","top"),
                  gp = gp_varnames[[2]])
    }
  }
}
class(labeling_doubledecker) <- "grapcon_generator"

labeling_left <- function(rep = FALSE, pos_varnames = "left",
                        pos_labels = "left", just_labels = "left", ...)
  labeling_border(rep = rep, pos_varnames = pos_varnames,
              pos_labels = pos_labels, just_labels = just_labels, ...)
class(labeling_left) <- "grapcon_generator"

labeling_left2 <- function(tl_labels = TRUE, clip = TRUE, pos_varnames = "left",
                        pos_labels = "left", just_labels = "left", ...)
  labeling_border(tl_labels = tl_labels, clip = clip, pos_varnames = pos_varnames,
              pos_labels = pos_labels, just_labels = just_labels, ...)
class(labeling_left2) <- "grapcon_generator"

labeling_cboxed <- function(tl_labels = TRUE, boxes = TRUE, clip = TRUE, pos_labels = "center", ...)
  labeling_border(tl_labels = tl_labels, boxes = boxes, clip = clip, pos_labels = pos_labels, ...)
class(labeling_cboxed) <- "grapcon_generator"

labeling_lboxed <- function(tl_labels = FALSE, boxes = TRUE, clip = TRUE, pos_labels = "left", just_labels = "left", labbl_varnames = FALSE, ...)
  labeling_border(tl_labels = tl_labels, boxes = boxes, clip = clip, pos_labels = pos_labels, labbl_varnames = labbl_varnames, just_labels = just_labels, ...)
class(labeling_lboxed) <- "grapcon_generator"

labeling_values <-
function(value_type = c("observed", "expected", "residuals"),
         suppress = NULL, digits = 1, clip_cells = FALSE, ...)
{
   value_type <- match.arg(value_type)
   if (value_type == "residuals" && is.null(suppress))
       suppress <- 2
   if (is.null(suppress))
       suppress <- 0
   if (length(suppress) == 1)
       suppress <- c(-suppress, suppress)

   function(d, split_vertical, condvars, prefix) {
       lookup <- if (value_type == "observed") "x" else value_type
       if (!exists(lookup, envir = parent.frame()))
           stop(paste("Could not find", dQuote(value_type), "object."))

       values <- get(lookup, envir = parent.frame())
       values <- ifelse((values > suppress[2]) | (values < suppress[1]),
                        round(values, digits), NA)
       labeling_border(...)(d, split_vertical, condvars, prefix)
       labeling_cells(text = values, clip_cells = clip_cells, ...)(d, split_vertical, condvars, prefix)
   }
}
class(labeling_values) <- "grapcon_generator"

labeling_residuals <-
function(suppress = NULL, digits = 1, clip_cells = FALSE, ...)
    labeling_values(value_type = "residuals", suppress = suppress,
                    digits = digits, clip_cells = clip_cells, ...)
class(labeling_residuals) <- "grapcon_generator"
