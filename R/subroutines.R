

cut_to_lagd <- function(x1){
  lapply(strsplit(gsub("\\(|\\]", "", as.character(x1)), ","),
         function(x2) mean(as.numeric(x2))) %>% unlist
}



comb_2 <- function(x){

  ind <- combi2inds(x)
  return(setDT(list(x[ind$xid], x[ind$yid])))

}

comb_diff <- function(x){

  ind <- combi2inds(x)
  return(setDT(list(x[ind$xid] - x[ind$yid])))

}




# min-max normalize function
scale_minmax <- function(x) {

  (x-min(x))/(max(x)-min(x))

}

# min-max normalize function: fixed scale
fixed_scale_minmax <- function(x, fixed_min, fixed_max){

  (x - fixed_min)/(fixed_max - fixed_min)

}


# lag distance overlapping function
lag_dist <- function(d, buffer, m_factor){

  starts <- seq(0, max(d), by = buffer/m_factor)
  ends <- dplyr::lead(starts, m_factor)

  lag_interval <- data.frame(lag_st = starts, lag_ed = ends)[1:length(which(!is.na(ends))), ]
  return(lag_interval)

}




p_dec <- function(t_seq, y, p, m){

  SV <- as.numeric(y)
  Tj <- p/m
  wj <- 2*pi/Tj

  constant_term <- rep(1, length(t_seq))
  cos_term <- t(apply(as.matrix(t_seq), 1, function(y) cos(wj*y)))
  sin_term <- t(apply(as.matrix(t_seq), 1, function(y) sin(wj*y)))

  if(length(m) == 1){
    cos_term <- t(cos_term)
    sin_term <- t(sin_term)
  }

  colnames(cos_term) <- paste0("coef_", "cos", m)
  colnames(sin_term) <- paste0("coef_", "sin", m)

  PP <- cbind(constant_term, cos_term, sin_term)

  CC <- solve(t(PP)%*%PP)%*%(t(PP)%*%SV)

  eSV <- PP%*%CC

  res <- list(coefs = CC[,1], fitted = as.numeric(eSV),
              residuals = as.numeric(y - eSV))
  return(res)

}

