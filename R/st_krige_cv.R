
# cross-validation

st_krige_cv <- function(xi, method = "nfold", n = 10,
                        s_var = 1:2, t_var = 3,
                        vgm_st = "SumMetric", vgm_options,
                        step_size = 10000, n_cores = 8){

  xi <- xi

  # shuffled sample index
  fold_id <- cut(sample(1:nrow(xi)), breaks = 10, labels = F)

  pb <- txtProgressBar(min = 0, max = n, style = 3)

  fold_res <- rbindlist(

    lapply(1:n,

           function(ff){

             target_coord <- xi[which(fold_id == ff), 1:4]
             temp_res1 <- st_krige(xi = xi[-which(fold_id == ff), ],
                                   x0 = target_coord,
                                   s_var = 1:2, t_var = 3,
                                   vgm_st = vgm_st,
                                   vgm_options = vgm_options,
                                   step_size = step_size, n_cores = n_cores,
                                   progress = F, silence = T)
             org_val <- xi[which(fold_id == ff), 5]
             temp_res2 <- data.table(temp_res1, org_val)

             gc(reset = T)

             setTxtProgressBar(pb, ff)

             return(temp_res2)

           }

    )

  )

  return(fold_res)

}
