
## theoretical variogram fitting function
vgt <- function(ex_vgm, vgm_list, par0, lower, upper, plot = TRUE){

  match_funs <- lapply(vgm_list, match.fun)
  fit <- optim(par0, match_funs[['vgm_s']],
               d = ex_vgm$lag_s, gamma = ex_vgm$gamma,
               return = "error",
               gr=NULL, method = "L-BFGS-B", lower = lower, upper = upper)

  new_s <- seq(0, max(ex_vgm$lag_s), length.out = 100)

  fitted_gamma <- match_funs[['vgm_s']](d = new_s, gamma = NULL,
                                        pars = fit$par,
                                        return = "value")
  fitted_vgm <- data.frame(new_s, fitted_gamma)

  if(plot == TRUE){

    # plotting fitted variogram
    # compare with empirical variogram
    plot(ex_vgm$lag_s, ex_vgm$gamma, pch = 16)
    lines(fitted_vgm$new_s, fitted_vgm$fitted_gamma)

  }

  res <- list(par = fit$par, fitted_vgm = fitted_vgm)

  return(res)


}




## theoretical variogram fitting function for spatio-temporal kriging

vgt_st <- function(ex_vgm, vgm_list, par0, lower, upper, plot = TRUE){

  # Spatio-Temporal Variogram Fitting
  match_funs <- lapply(vgm_list, match.fun)

  fit <- optim(par0, match_funs[['vgm_st']],
               ds = ex_vgm$lag_s, dt = ex_vgm$lag_t, gamma = ex_vgm$gamma,
               vgm_s = match_funs[['vgm_s']],
               vgm_t = match_funs[['vgm_t']],
               vgm_jt = match_funs[['vgm_jt']],
               return = "error",
               gr=NULL, method = "L-BFGS-B", lower = lower, upper = upper)

  new_s <- seq(0, max(ex_vgm$lag_s), length.out = 100)
  new_t <- seq(0, max(ex_vgm$lag_t), length.out = 100)

  new_st <- expand.grid(lag_s = new_s, lag_t = new_t)


  fitted_gamma <- match_funs[['vgm_st']](ds = new_st$lag_s, dt = new_st$lag_t, gamma = NULL,
                                         vgm_s = "Sph", vgm_t = "Sph", vgm_jt = "Sph",
                                         pars = fit$par,
                                         return = "value")
  fitted_vgm <- data.frame(new_st, fitted_gamma)


  if(plot == TRUE){

    # plotting fitted variogram
    # compare with empirical variogram
    exp_vgm <- lattice::wireframe(gamma ~ lag_s * lag_t, ex_vgm,
                                  outer = TRUE, shade = TRUE,
                                  scale=list(arrows = FALSE, cex = 1.5),
                                  main = "Experimental Variogram",
                                  zlab = list(label = "Semi-variance", cex = 1.5, rot = 90),
                                  xlab = list(label = "Distance(h)", cex = 1.5, rot = 30),
                                  ylab = list(label = "Distance(t)", cex = 1.5, rot = -40), more = T)

    fit_vgm <- lattice::wireframe(fitted_gamma ~ lag_s * lag_t, fitted_vgm,
                                  outer = TRUE, shade = TRUE,
                                  scale=list(arrows = FALSE, cex = 1.5),
                                  main = "Fitted Variogram",
                                  zlim = range(ex_vgm$gamma),
                                  zlab = list(label = "Semi-variance", cex = 1.5, rot = 90),
                                  xlab = list(label = "Distance(h)", cex = 1.5, rot = 30),
                                  ylab = list(label = "Distance(t)", cex = 1.5, rot = -40))

    print(
      gridExtra::grid.arrange(exp_vgm, fit_vgm, ncol = 2)
    )

  }


  res <- list(par = fit$par, err = fit$value, fitted_vgm = fitted_vgm)

  return(res)

}
