

library(STK)
library(tidyverse)
library(lubridate)


# using sample data
data("no3_4d_sample", package = "STK")
data("target_coord_sample", package = "STK")


ex_vgm <- vge_st(
  xi = no3_4d_sample, s_var = 1:3, t_var = 4,
  max_range = 0.5, lag_buffer = 0.1, overlap_factor = 2,
  plot = TRUE
)


# step 3: fitting theoretical variogram using function `vgt()`

vg_fit <- vgt_st(
  ex_vgm = ex_vgm,
  vgm_list = list(vgm_s = "Sph", vgm_t = "Sph", vgm_jt = "Sph", vgm_st = "SumMetric"),
  par0 = c(100, 300, 0.2, 100, 300, 0.3, 1.0, 1.0, 0.2, 0.1),
  lower = c(0.01, 0.01, 0.05, 0.01, 0.01, 0.05, 0.01, 0.01, 0.05, 0.0001),
  upper = c(400, 600, 2, 300, 500, 2, 300, 300, 2, 300),
  plot = TRUE
)



stk_res <- st_krige(xi = no3_4d_sample, x0 = target_coord_sample,
                    s_var = 1:3, t_var = 4,
                    vgm_st = "SumMetric",
                    vgm_options = list(vgm_s = "Sph",
                                       vgm_t = "Sph",
                                       vgm_jt = "Sph",
                                       pars = vg_fit$par),
                    step_size = 10000, n_cores = 10)




