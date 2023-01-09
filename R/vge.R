

## empirical variogram function
vge <- function(xi, s_var, val,
                max_range, lag_buffer, overlap_factor,
                valid_n = 30,
                plot = TRUE){

  sample_var <- var(xi[, val])
  pair_val <- comb_2(xi[, val])

  pair_val[, dist_s := as.numeric(dist(xi[, s_var]))] # add variables[dist_s] on the exist table
  lag_d_s <- lag_dist(pair_val$dist_s, buffer = lag_buffer, m_factor = overlap_factor)
  pair_val[, cut_s :=
             cut(dist_s,
                 c(unique(unlist(lag_d_s)),
                   pair_val[, ceiling(max(dist_s)+lag_buffer)]),
                 right = F)]

  pair_val[, n := .(n = length(V1)),
           keyby = .(cut_s)]

  # invalid(= n group <= 30) pair filtering
  invalid_pair <- which(pair_val$n <= valid_n)

  # experimental variogram calculation
  vg_e = pair_val[-invalid_pair, .(
    cov = cov(V1, V2, use = "everything"),
    # conf_int_cor = cor.test(V1, V2)$conf.int,
    sd_x = sd(V1, na.rm = T),
    sd_y = sd(V2, na.rm = T)),
    keyby = .(cut_s)]
  vg_e[, gamma := sample_var - cov]
  vg_e[, sd_xy := sd_x * sd_y]
  vg_e[, conf_int_cov := sd_x * sd_y]
  vg_e[, lag_s := mean(
    as.numeric(
      gsub("\\[|\\)", "", unlist(str_split(cut_s, ",")))
    )
  ),
  by = .(cut_s)]

  # 2 row confidence interval to 1 row
  vg_e[, .(conf_cov_l = conf_int_cov[2], conf_cov_u = conf_int_cov[1]),
       by = .(cut_s, cov, sd_x, sd_y, sd_xy, lag_s)]

  # confirm empirical variogram
  vg_e_lim <- vg_e %>%
    filter(lag_s <= max_range * max(lag_s))

  # plotting option
  if(plot == TRUE){

    print(
      plot(vg_e_lim$lag_s, vg_e_lim$gamma, pch = 16)
    )


  }

  return(vg_e_lim)


}



vge_st <- function(xi, s_var, t_var, max_range, lag_buffer, overlap_factor, plot = TRUE){

  sample_var <- var(xi$val)
  pair_val <- comb_2(xi$val)

  pair_val[, dist_s := as.numeric(dist(xi[, s_var]))] # add variables[dist_s] on the exist table
  pair_val[, dist_t := as.numeric(dist(xi[, t_var]))] # add variables[dist_t] on the exist table

  lag_d_s <- lag_dist(pair_val$dist_s, buffer = lag_buffer, m_factor = overlap_factor)
  lag_d_t <- lag_dist(pair_val$dist_t, buffer = lag_buffer, m_factor = overlap_factor)


  #cut_s / cut_t: separation for s and t
  pair_val[, cut_s :=
             cut(dist_s,
                 c(unique(unlist(lag_d_s)),
                   pair_val[, ceiling(max(dist_s)+lag_buffer)]),
                 right = F)]
  pair_val[, cut_t :=
             cut(dist_t,
                 c(unique(unlist(lag_d_t)),
                   pair_val[, ceiling(max(dist_t)+lag_buffer)]),
                 right = F)]


  # invalid(= n group <= 30) pair filtering
  pair_val[, n := .(n = length(V1)),
           keyby = .(cut_s, cut_t)]
  invalid_pair <- which(pair_val$n <= 30)


  # experimental variogram calculation
  vg_e = pair_val[-invalid_pair, .(
    cov = cov(V1, V2, use = "everything"),
    # conf_int_cor = cor.test(V1, V2)$conf.int,
    sd_x = sd(V1, na.rm = T),
    sd_y = sd(V2, na.rm = T)),
    keyby = .(cut_s, cut_t)]
  vg_e[, gamma := sample_var - cov]
  vg_e[, sd_xy := sd_x * sd_y]
  vg_e[, conf_int_cov := sd_x * sd_y]
  vg_e[, lag_s := mean(as.numeric(gsub("\\[|\\)", "",unlist(str_split(cut_s, ","))))), by = .(cut_s)]
  vg_e[, lag_t := mean(as.numeric(gsub("\\[|\\)", "",unlist(str_split(cut_t, ","))))), by = .(cut_t)]


  # 2 row confidence interval to 1 row
  vg_e[, .(conf_cov_l = conf_int_cov[2], conf_cov_u = conf_int_cov[1]),
       by = .(cut_s, cut_t, cov, sd_x, sd_y, sd_xy, lag_s, lag_t)]


  # confirm empirical variogram
  vg_e_lim <- vg_e %>%
    filter(lag_s <= max_range * max(lag_s), lag_t <= max_range * max(lag_t))


  # plotting option
  if(plot == TRUE){

    print(
      lattice::wireframe(gamma ~ lag_s * lag_t, vg_e_lim,
                         outer = TRUE, shade = TRUE,
                         scale = list(arrows = FALSE, cex = 1.5),
                         zlab = list(label = "Semi-variance", cex = 1.5, rot = 90),
                         xlab = list(label = "Distance(h)", cex = 1.5, rot = 30),
                         ylab = list(label = "Distance(t)", cex = 1.5, rot = -40))
    )


  }


  return(vg_e_lim)

}



# anisotropic experimental variogram function for xyzt
vge_aniso <- function(xi, s_var, val, aniso = NULL,
                      max_range, lag_buffer, overlap_factor,
                      valid_n = 30,
                      plot = TRUE){


  aniso <- aniso
  aniso_dist <- paste0("dist_", aniso)
  aniso_cut <- paste0("cut_", aniso)
  aniso_lag <- paste0("lag_", aniso)
  aniso_conf <- paste0(c("conf_cov_l_", "conf_cov_u_"),
                       rep(aniso, each = 2))

  sample_var <- var(xi[, val])
  pair_val <- comb_2(xi[, val])

  # add variables[dist_x, y, z] on the exist table
  pair_val[, (aniso_dist) :=
             lapply(aniso, function(x) as.numeric(dist(xi[, x])))]

  # lag_d_s <- lag_dist(pair_val$dist_x, buffer = lag_buffer, m_factor = overlap_factor)

  lag_d_list <- lapply(aniso_dist,
                       function(x){

                         lag_dist(pair_val[[x]],
                                  buffer = lag_buffer,
                                  m_factor = overlap_factor)

                       })

  names(lag_d_list) <- aniso_dist

  pair_val[, (aniso_cut) :=
             lapply(aniso_dist,
                    function(x){
                      cut(pair_val[[x]],
                          c(unique(unlist(lag_d_list[[x]])),
                            ceiling(max(pair_val[[x]])) + lag_buffer),
                          right = F)
                    })
  ]


  pair_val[, n := .(n = length(V1)),
           keyby = aniso_cut]


  # invalid(= n group <= 30) pair filtering
  invalid_pair <- which(pair_val$n <= 10)

  # experimental variogram calculation
  vg_e = pair_val[-invalid_pair, .(
    cov = cov(V1, V2, use = "everything"),
    # conf_int_cor = cor.test(V1, V2)$conf.int,
    sd_x = sd(V1, na.rm = T),
    sd_y = sd(V2, na.rm = T)),
    keyby = aniso_cut]
  vg_e[, gamma := sample_var - cov]
  vg_e[, conf_int_gamma := sample_var - sd_x * sd_y]


  vg_e[, (aniso_lag) :=
         lapply(.SD, function(x){
           mean(
             as.numeric(
               gsub("\\[|\\)", "", unlist(str_split(x, ",")))
             )
           )
         }),
       .SDcols = aniso_cut,
       keyby = aniso_cut]


  # covariance confidence interval(Not Recommended)
  agg_cols <- c(aniso_cut, "cov", "sd_x", "sd_y", aniso_lag)
  vg_e[, (aniso_conf) :=
         lapply(.SD, function(x){
           c(x[1], x[2])
         }), .SDcols = "conf_int_gamma",
       keyby = agg_cols]


  # confirm empirical variogram
  vg_e_lim <- vg_e %>%
    filter(if_all(starts_with("lag_"), ~ . <= max_range * max(vg_e$lag_x)))


  # plotting option
  if(plot == TRUE){

    op <- par(no.readonly = T)
    par(mfrow = c(2, 2))

    for(pp in aniso_lag){

      plot(vg_e_lim[[pp]], vg_e_lim$gamma, pch = 16,
           main = paste0("Gamma For ", pp))

    }

    par(op)

  }


  return(vg_e_lim)


}




