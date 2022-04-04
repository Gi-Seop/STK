
# kriging core
st_krige <- function(xi, x0, s_var, t_var,
                     vgm_st = "SumMetric", vgm_options,
                     step_size = 10000, n_cores = 8,
                     progress = T, silence = F){


  obs_h <- crossdist(as.matrix(xi[, s_var]), as.matrix(xi[, s_var]))
  obs_t <- crossdist(as.matrix(xi[, t_var]), as.matrix(xi[, t_var]))

  if(vgm_st != "SumMetric"){

    message("Only SumMetric st variogram is provided yet.")

  }

  vgm_st <- match.fun(vgm_st)
  xs <- as.matrix(rbind(-xi[, -which(colnames(xi) == "val")], 0))
  xs_t <- t(-xs)
  mat_zero <- matrix(0, nrow = nrow(xs_t), ncol = ncol(xs))

  if(silence == F){
    message("\npreparing LHS matrix including matrix inverse...\n please wait...")
  }

  mat_L_inv <- vgm_st(ds = obs_h, dt = obs_t, gamma = NULL,
                      vgm_s = vgm_options$vgm_s,
                      vgm_t = vgm_options$vgm_t,
                      vgm_jt = vgm_options$vgm_jt,
                      pars = vgm_options$pars,
                      return = "value") %>%
    matrix(., nrow = dim(obs_h)[1], ncol = dim(obs_h)[2]) %T>%
    {assign(x = "ss", value = max(., na.rm = T), envir = .GlobalEnv)} %>%
    {ss - .} %>%
    `diag<-`(., ss) %>%
    cbind(., rep(-1, dim(obs_h)[1])) %>%
    rbind(c(rep(1, dim(obs_h)[2]), 0)) %>%
    cbind(., xs) %>%
    rbind(., cbind(xs_t, mat_zero)) %>%
    solve

  if(silence == F){
    message("LHS matrix is prepared..\n calculating for target coordinates")
  }

  # step_size <- step_size
  step_size <- ifelse(step_size > nrow(x0),
                      nrow(x0),
                      step_size)
  max_val = floor(nrow(x0)/step_size)

  if(progress == T){
    pb <- txtProgressBar(min = 0, max = max_val, style = 3)
  }

  res <- rbindlist(

    lapply(1:max_val,
           function(kk){

             if(kk != max_val){
               sample_list = (1+(kk - 1) * step_size):(kk * step_size)
             }else{
               sample_list = (1+(kk - 1) * step_size):nrow(x0)
             }


             target_h <- crossdist(as.matrix(x0[sample_list, s_var]),
                                   as.matrix(xi[, s_var]))
             target_t <- crossdist(as.matrix(x0[sample_list,t_var]),
                                   as.matrix(xi[,t_var]))


             mat_t <- vgm_st(ds = target_h, dt = target_t, gamma = NULL,
                             vgm_s = vgm_options$vgm_s,
                             vgm_t = vgm_options$vgm_t,
                             vgm_jt = vgm_options$vgm_jt,
                             pars = vgm_options$pars,
                             return = "value") %>%
               matrix(., nrow = dim(target_h)[1], ncol = dim(target_h)[2]) %>%
               {ss - .} %>%
               cbind(., 1) %>%
               cbind(., x0[sample_list,]) %>%
               t


             coefs <- eigenMapMatMult2(mat_L_inv , mat_t, n_cores = n_cores)
             gc(reset = T)


             lambda_id <- 1:nrow(xi)
             omega_id <- (nrow(xi)+1):dim(coefs)[1]
             lambda <- coefs[lambda_id,]
             omega <- coefs[omega_id,]


             est_var1 = diag(eigenMapMatMult2(t(coefs[lambda_id,]),
                                              mat_t[lambda_id,],
                                              n_cores = n_cores))
             est_var2 = diag(eigenMapMatMult2(t(coefs[omega_id,]),
                                              mat_t[omega_id,],
                                              n_cores = n_cores))

             if(progress == T){
               setTxtProgressBar(pb, kk)
             }

             return(

               data.table(est_val = as.numeric(xi[,which(colnames(xi) == "val")] %*% lambda),
                          est_var =  ss - est_var1 + est_var2 )

             )


           }
    )

  )

  return(res)

}
