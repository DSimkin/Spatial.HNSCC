
.limit <- function(v) {    
  m = 10
  Sign = sign(v)
  v = abs(v)
  v.lim = Sign * (floor(m * v) / m)
  return(v.lim)
}

.limits <- function(v, symmetric = TRUE) {
  stopifnot(length(v) == 2, all(sapply(v, class) == 'numeric'))
  
  if (!symmetric) {
    return(sapply(v, .limit))
  }
  
  limit = .limit(min(abs(v)))
  return(c(-1 * limit, 1 * limit))
}


.axis.spacer <- function(breaks, labels, limits, levels = NULL) {
  if (!is.null(labels) & !is.null(levels)) {
    breaks = levels %in% labels
    labels = levels[breaks]
  }
  if (is.null(breaks)) {
    breaks = seq(limits[[1]], limits[[2]], limits[[2]])
  }
  if (is.null(labels)) {
    labels = breaks
  }
  return(list(breaks = breaks, labels = labels))
}

.limit <- function(v) {    
  m = 10
  Sign = sign(v)
  v = abs(v)
  v.lim = Sign * (floor(m * v) / m)
  return(v.lim)
}

.limits <- function(v, symmetric = TRUE) {
  stopifnot(length(v) == 2, all(sapply(v, class) == 'numeric'))
  
  if (!symmetric) {
    return(sapply(v, .limit))
  }
  
  limit = .limit(min(abs(v)))
  return(c(-1 * limit, 1 * limit))
}


.axis.spacer <- function(breaks, labels, limits, levels = NULL) {
  if (!is.null(labels) & !is.null(levels)) {
    breaks = levels %in% labels
    labels = levels[breaks]
  }
  if (is.null(breaks)) {
    breaks = seq(limits[[1]], limits[[2]], limits[[2]])
  }
  if (is.null(labels)) {
    labels = breaks
  }
  return(list(breaks = breaks, labels = labels))
}



gmap <- function(dat,
                 x,
                 y,
                 fill = 1,
                 type = "continuous",
                 geom = c("tile", "raster"),
                 limits = c(-0.5, 0.5),
                 lim.find = F,
                 lim.sym = T,
                 midpoint = 0,
                 x.name = NULL,
                 y.name = NULL,
                 angle = NULL,
                 axis.rel = 1,
                 title.rel = 1.1,
                 title = NULL,
                 subtitle = NULL,
                 caption = NULL,
                 text.size = 12,
                 ratio = NULL,
                 tile.size = 0.1,
                 tile.col = "gray",
                 cols = NULL,
                 col = NULL,
                 na.value = 'white',
                 legend.position = 'right',
                 legend.height = 0.4,
                 legend.width = 0.6,
                 legend.rel = 0.8,
                 legend.colour = 'black',
                 ticks.linewidth = 0.5,
                 breaks = waiver(),
                 labels = waiver(),
                 x.breaks = waiver(),
                 y.breaks = waiver(),
                 x.labels = waiver(),
                 y.labels = waiver(),
                 num = F,
                 y.num = num, 
                 x.num = num,
                 legend.breaks = NULL,
                 legend.labels = NULL,
                 legend.title = NULL,
                 legend.title.rel = 0.8,
                 legend.title.hjust = 0.5,
                 expand = c(0,0),
                 minimal = FALSE,
                 magma_pal = FALSE) {
  
  x = rlang::enquo(x)
  y = rlang::enquo(y)
  fill = rlang::enquo(fill)
  xname = rlang::quo_name(x)
  yname = rlang::quo_name(y)
  
  if (class(breaks) != 'waiver') {
    x.breaks = breaks
    y.breaks = breaks
  }
  
  if (class(labels) != 'waiver') {
    x.labels = labels
    y.labels = labels
  }
  
  if (x.num | num | class(dat %>% pull(!!x)) == 'numeric') {
    dat = dat %>% dplyr::mutate(!!xname := as.numeric(!!x))
    x.scale.FUN = ggplot2::scale_x_continuous
  } else {
    x.scale.FUN = ggplot2::scale_x_discrete
  }
  if (y.num | num | class(dat %>% pull(!!y)) == 'numeric') {
    dat = dat %>% dplyr::mutate(!!yname := as.numeric(!!y))
    y.scale.FUN = ggplot2::scale_y_continuous
  } else {
    y.scale.FUN = ggplot2::scale_y_discrete
  }
  
  x.scale = quote(x.scale.FUN(expand = expand,
                              breaks = x.breaks,
                              labels = x.labels))
  y.scale = quote(y.scale.FUN(expand = expand,
                              breaks = y.breaks,
                              labels = y.labels))
  
  
  if (!is.null(col)) {
    cols = col
    legend.position = 'none'
  }
  
  if (is.null(cols) & type == "discrete") {
    cols = readRDS("/home/labs/tirosh/dorsi/Spatial_HNSCC/hotmap.rds")
  }
  
  if (is.null(cols) & limits[[1]] >= 0 & type == "continuous") {
    cols = RColorBrewer::brewer.pal(9, 'YlOrRd')
  } else if (is.null(cols) & type == "continuous") {
    cols = c("#053061",
             "#2166ac",
             "#4393c3",
             "#92c5de", 
             "#d1e5f0",
             "#f7f7f7",
             "#fddbc7",
             "#f4a582",
             "#d6604d",
             "#b2182b",
             "#67001f")
  }
  
  if (magma_pal) {
    library(RColorBrewer)
    library(viridis)
    cols <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
  }
  
  if (lim.find) {
    v = dat %>% dplyr::select(!!fill) %>% range
    limits = .limits(v = v, symmetric = lim.sym)
  }
  
  if (isTRUE(angle)) {
    angle = 45
  }
  
  if (is.numeric(angle)) {
    angle = ggpubr::rotate_x_text(angle = angle, vjust = 1)
  }
  
  legend = .axis.spacer(breaks = legend.breaks, labels = legend.labels, limits = limits)
  legend.breaks = legend$breaks
  legend.labels = legend$labels
  
  # geomtile = geom_tile(col = tile.col, size = tile.size)
  # if (any(sapply(list(tile.size, tile.col), is.null))) {
  #   geomtile = ggplot2::geom_tile()
  # }
  
  if(type == "continuous" & length(midpoint) == 1 && midpoint == 0) {
    grad_steps <- (max(limits) - min(limits)) / (length(cols) - 1)
    col_vals <- seq(from = min(limits), to = max(limits), by = grad_steps)
  } else if (type == "continuous" & length(midpoint) == 1 && midpoint != 0) {
    zero_col_val <- round(length(cols) / 2)
    low_grad_steps <- (midpoint - min(limits)) / (zero_col_val - 1)
    high_grad_steps <- (max(limits) - midpoint) / (zero_col_val - 1)
    col_vals <- c(seq(from = min(limits), to = midpoint, by = low_grad_steps), seq(from = midpoint, to = max(limits), by = high_grad_steps)[-1])
  } else if (type == "continuous" & length(midpoint) > 1) {
    zero_col_val <- round(length(cols) / 2)
    cols <- c(cols[1:(zero_col_val - 1)], rep(cols[zero_col_val], 2), cols[(zero_col_val + 1):length(cols)])
    low_grad_steps <- (min(midpoint) - min(limits)) / (zero_col_val - 1)
    high_grad_steps <- (max(limits) - max(midpoint)) / (zero_col_val - 1)
    col_vals <- c(seq(from = min(limits), to = min(midpoint), by = low_grad_steps), 
                  seq(from = max(midpoint), to = max(limits), by = high_grad_steps))
  }
  
  if(type == "discrete") {
    scale_fill <- ggplot2::scale_fill_gradient(2,low="#FFFFFF", high="#0c2a50ff",
                                               limits = limits,
                                               expand = expand,
                                               oob = scales::squish,
                                               breaks = breaks,
                                               labels = labels,
                                               name = legend.title,
                                               guide = guide_legend(frame.colour = 'black', label=TRUE))
  } else if(type == "continuous") {
    scale_fill <- ggplot2::scale_fill_gradientn(colors = cols,
                                                values = scales::rescale(col_vals),
                                                limits = limits,
                                                expand = expand,
                                                oob = scales::squish,
                                                breaks = legend.breaks,
                                                labels = legend.breaks,
                                                name = legend.title,
                                                na.value = na.value,
                                                guide = ggplot2::guide_colorbar(frame.colour='black',
                                                                                ticks.colour='black',
                                                                                title.position='top',
                                                                                title.hjust=legend.title.hjust,
                                                                                barwidth=legend.width,
                                                                                barheight=legend.height))
  }
  
  geom_selec <- switch(match.arg(geom),
                       tile = geom_tile(color = tile.col, size = tile.size),
                       raster = geom_raster())
  
  if (minimal) {
    G <- ggplot2::ggplot(dat, aes(x = !!x, y = !!y, fill = !!fill, group = 1)) +
      geom_selec +
      eval(scale_fill) +
      ggplot2::labs(x = x.name,
                    y = y.name,
                    title = title,
                    subtitle = subtitle,
                    caption = caption) 
  } else {
    G <- ggplot2::ggplot(dat, aes(x = !!x, y = !!y, fill = !!fill, group = 1)) +
      geom_selec +
      eval(scale_fill) +
      ggplot2::labs(x = x.name,
                    y = y.name,
                    title = title,
                    subtitle = subtitle,
                    caption = caption) +
      ggplot2::theme_bw(base_size = text.size) +
      angle +
      ggplot2::theme(aspect.ratio = ratio,
                     panel.grid = ggplot2::element_blank(),
                     title = ggplot2::element_text(size = ggplot2::rel(title.rel)),
                     axis.title = ggplot2::element_text(size = ggplot2::rel(axis.rel)),
                     legend.position = legend.position,
                     legend.text = ggplot2::element_text(size = ggplot2::rel(legend.rel),
                                                         colour = legend.colour,
                                                         hjust = legend.title.hjust),
                     legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.rel)),
                     legend.margin = margin(t = -0.5, unit='cm'),
                     legend.key.height = grid::unit(legend.height, "cm"),
                     legend.key.width = grid::unit(legend.width, "cm")) +
      eval(x.scale) +
      eval(y.scale)  
  }
  
  G
}





gmap2 <- function(dat,
                 x,
                 y,
                 fill = 1,
                 type = c("continuous", "discrete"),
                 geom = c("tile", "raster"),
                 limits = c(-0.5, 0.5),
                 lim.find = F,
                 lim.sym = T,
                 midpoint = 0,
                 x.name = NULL,
                 y.name = NULL,
                 angle = NULL,
                 axis.rel = 1,
                 title.rel = 1.1,
                 title = NULL,
                 subtitle = NULL,
                 caption = NULL,
                 text.size = 12,
                 ratio = NULL,
                 tile.size = 0.1,
                 tile.col = "gray",
                 cols = NULL,
                 col = NULL,
                 na.value = 'white',
                 legend.position = 'right',
                 legend.height = 0.4,
                 legend.width = 0.6,
                 legend.rel = 0.8,
                 legend.colour = 'black',
                 ticks.linewidth = 0.5,
                 breaks = waiver(),
                 labels = waiver(),
                 x.breaks = waiver(),
                 y.breaks = waiver(),
                 x.labels = waiver(),
                 y.labels = waiver(),
                 num = F,
                 y.num = num, 
                 x.num = num,
                 legend.breaks = NULL,
                 legend.labels = NULL,
                 legend.title = NULL,
                 legend.title.rel = 0.8,
                 expand = c(0,0),
                 minimal = FALSE,
                 magma_pal = FALSE) {
  
  x = rlang::enquo(x)
  y = rlang::enquo(y)
  fill = rlang::enquo(fill)
  xname = rlang::quo_name(x)
  yname = rlang::quo_name(y)
  
  if (class(breaks) != 'waiver') {
    x.breaks = breaks
    y.breaks = breaks
  }
  
  if (class(labels) != 'waiver') {
    x.labels = labels
    y.labels = labels
  }
  
  if (x.num | num | class(dat %>% pull(!!x)) == 'numeric') {
    dat = dat %>% dplyr::mutate(!!xname := as.numeric(!!x))
    x.scale.FUN = ggplot2::scale_x_continuous
  } else {
    x.scale.FUN = ggplot2::scale_x_discrete
  }
  if (y.num | num | class(dat %>% pull(!!y)) == 'numeric') {
    dat = dat %>% dplyr::mutate(!!yname := as.numeric(!!y))
    y.scale.FUN = ggplot2::scale_y_continuous
  } else {
    y.scale.FUN = ggplot2::scale_y_discrete
  }
  
  x.scale = quote(x.scale.FUN(expand = expand,
                              breaks = x.breaks,
                              labels = x.labels))
  y.scale = quote(y.scale.FUN(expand = expand,
                              breaks = y.breaks,
                              labels = y.labels))
  
  
  if (!is.null(col)) {
    cols = col
    legend.position = 'none'
  }
  
  if (is.null(cols) & type == "discrete") {
    cols = readRDS("/home/labs/tirosh/dorsi/Spatial_HNSCC/hotmap.rds")
  }
  
  if (is.null(cols) & limits[[1]] >= 0 & type == "continuous") {
    cols = RColorBrewer::brewer.pal(9, 'YlOrRd')
  } else if (is.null(cols) & type == "continuous") {
    cols = c("#053061",
             "#2166ac",
             "#4393c3",
             "#92c5de", 
             "#d1e5f0",
             "#f7f7f7",
             "#fddbc7",
             "#f4a582",
             "#d6604d",
             "#b2182b",
             "#67001f")
  }
  
  if (magma_pal) {
    library(RColorBrewer)
    library(viridis)
    cols <- c(grDevices::colorRampPalette(c(colorspace::lighten("antiquewhite", .75), "antiquewhite", rev(viridis::magma(30, begin = .25, end = .9))[1]))(15), rev(viridis::magma(15, begin = .25, end = .9)))
  }
  
  if (lim.find) {
    v = dat %>% dplyr::select(!!fill) %>% range
    limits = .limits(v = v, symmetric = lim.sym)
  }
  
  if (isTRUE(angle)) {
    angle = 45
  }
  
  if (is.numeric(angle)) {
    angle = ggpubr::rotate_x_text(angle = angle, vjust = 1)
  }
  
  legend = .axis.spacer(breaks = legend.breaks, labels = legend.labels, limits = limits)
  legend.breaks = legend$breaks
  legend.labels = legend$labels
  
  # geomtile = geom_tile(col = tile.col, size = tile.size)
  # if (any(sapply(list(tile.size, tile.col), is.null))) {
  #   geomtile = ggplot2::geom_tile()
  # }
  
  if(type == "continuous" & length(midpoint) == 1 && midpoint == 0) {
    grad_steps <- (max(limits) - min(limits)) / (length(cols) - 1)
    col_vals <- seq(from = min(limits), to = max(limits), by = grad_steps)
  } else if (type == "continuous" & length(midpoint) == 1 && midpoint != 0) {
    zero_col_val <- round(length(cols) / 2)
    low_grad_steps <- (midpoint - min(limits)) / (zero_col_val - 1)
    high_grad_steps <- (max(limits) - midpoint) / (zero_col_val - 1)
    col_vals <- c(seq(from = min(limits), to = midpoint, by = low_grad_steps), seq(from = midpoint, to = max(limits), by = high_grad_steps)[-1])
  } else if (type == "continuous" & length(midpoint) > 1) {
    zero_col_val <- round(length(cols) / 2)
    cols <- c(cols[1:(zero_col_val - 1)], rep(cols[zero_col_val], 2), cols[(zero_col_val + 1):length(cols)])
    low_grad_steps <- (min(midpoint) - min(limits)) / (zero_col_val - 1)
    high_grad_steps <- (max(limits) - max(midpoint)) / (zero_col_val - 1)
    col_vals <- c(seq(from = min(limits), to = min(midpoint), by = low_grad_steps), 
                  seq(from = max(midpoint), to = max(limits), by = high_grad_steps))
  }
  
  if(type == "discrete") {
    scale_fill <- ggplot2::scale_fill_gradient(2,low="#FFFFFF", high="#0c2a50ff",
                                               limits = limits,
                                               expand = expand,
                                               oob = scales::squish,
                                               breaks = breaks,
                                               labels = labels,
                                               name = legend.title,
                                               guide = guide_legend(frame.colour = 'black', label=TRUE))
  } else if(type == "continuous") {
    scale_fill <- ggplot2::scale_fill_gradientn(colors = cols,
                                                values = scales::rescale(col_vals),
                                                limits = limits,
                                                expand = expand,
                                                oob = scales::squish,
                                                breaks = legend.breaks,
                                                labels = legend.breaks,
                                                name = legend.title,
                                                na.value = na.value,
                                                guide = ggplot2::guide_colorbar(frame.colour='black',
                                                                                ticks.colour='black',
                                                                                title.position='top',
                                                                                title.hjust=0.5,
                                                                                barwidth=legend.width,
                                                                                barheight=legend.height))
  }
  
  geom_selec <- switch(match.arg(geom),
                       tile = geom_tile(color = tile.col, size = tile.size),
                       raster = geom_raster())
  
  if (minimal) {
    G <- ggplot2::ggplot(dat, aes(x = !!x, y = !!y, fill = !!fill, group = 1)) +
      geom_selec +
      eval(scale_fill) +
      ggplot2::labs(x = x.name,
                    y = y.name,
                    title = title,
                    subtitle = subtitle,
                    caption = caption) 
  } else {
    G <- ggplot2::ggplot(dat, aes(x = !!x, y = !!y, fill = !!fill, group = 1)) +
      geom_selec +
      eval(scale_fill) +
      ggplot2::labs(x = x.name,
                    y = y.name,
                    title = title,
                    subtitle = subtitle,
                    caption = caption) +
      ggplot2::theme_bw(base_size = text.size) +
      angle +
      ggplot2::theme(aspect.ratio = ratio,
                     panel.grid = ggplot2::element_blank(),
                     title = ggplot2::element_text(size = ggplot2::rel(title.rel)),
                     axis.title = ggplot2::element_text(size = ggplot2::rel(axis.rel)),
                     legend.position = legend.position,
                     legend.text = ggplot2::element_text(size = ggplot2::rel(legend.rel),
                                                         colour = legend.colour,
                                                         hjust = 0.5),
                     legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.rel)),
                     legend.margin = margin(t = -0.5, unit='cm'),
                     legend.key.height = grid::unit(legend.height, "cm"),
                     legend.key.width = grid::unit(legend.width, "cm")) +
      eval(x.scale) +
      eval(y.scale)  
  }
  
  G
}
