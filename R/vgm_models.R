


# variogram model for spacial ans temporal covariance
Exp <- function(d, gamma, pars, return = "error"){
  
  v_model <- ifelse(d <= pars[3],
                    pars[1] + pars[2]*(1- exp(-3*(d/pars[3]))),
                    pars[1] + pars[2])
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}



Gau <- function(d, gamma, pars, return = "error"){
  
  v_model <- ifelse(d <= pars[3],
                    pars[1] + pars[2]*(1- exp(-3*(d/pars[3])^2)),
                    pars[1] + pars[2])
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}



Sph <- function(d, gamma, pars, return = "value"){
  
  v_model <- ifelse(d <= pars[3],
                    pars[1] + pars[2]*(1.5*(d/pars[3]) -0.5*(d/pars[3])^3),
                    pars[1] + pars[2]) 
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}



Cos <- function(d, gamma, pars, return = "error"){
  
  v_model <- pars[1]*(1-cos(pi*d/pars[2]))
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}



Gslib <- function(d, gamma, pars, return = "error"){
  
  v_model <- pars[1]*
    (1-exp(-3*d/pars[2])*
       cos(pi*d/pars[3]))
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}



SumMetric <- function(ds, dt, gamma,
                      vgm_s = "Sph", vgm_t = "Sph", vgm_jt = "Sph",
                      pars,
                      return = "value"){
  
  par_s <- pars[1:3]
  par_t <- pars[4:6]
  par_st <- pars[7:9]
  STaniso <- pars[10]
  
  vgm_s <- match.fun(vgm_s)
  vgm_t <- match.fun(vgm_t)
  vgm_jt <- match.fun(vgm_jt)
  
  gamma_s <- vgm_s(d = ds, gamma = gamma, pars = par_s, return = "value")
  gamma_t <- vgm_t(d = dt, gamma = gamma, pars = par_t, return = "value")
  
  new_st <- sqrt((ds^2+(STaniso*dt)^2))
  gamma_st <- vgm_jt(d = new_st, gamma = gamma, pars = par_st, return = "value")
  v_model <- gamma_s + gamma_t + gamma_st
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}


SumMetric2 <- function(ds, dt, gamma,
                       vgm_s = "Sph", vgm_t = "Sph", vgm_jt = "Sph",
                       pars,
                       return = "value"){
  
  if(vgm_s == "Zhou"){
    
    par_s <- pars[1:5]
    par_t <- pars[6:8]
    par_st <- pars[9:11]
    STaniso <- pars[12]
    
  }else{
    
    par_s <- pars[1:3]
    par_t <- pars[4:6]
    par_st <- pars[7:9]
    STaniso <- pars[10]
    
  }
  
  vgm_s <- match.fun(vgm_s)
  vgm_t <- match.fun(vgm_t)
  vgm_jt <- match.fun(vgm_jt)
  
  gamma_s <- vgm_s(d = ds, gamma = gamma, pars = par_s, return = "value")
  gamma_t <- vgm_t(d = dt, gamma = gamma, pars = par_t, return = "value")
  
  new_st <- sqrt((ds^2+(STaniso*dt)^2))
  gamma_st <- vgm_jt(d = new_st, gamma = gamma, pars = par_st, return = "value")
  v_model <- gamma_s + gamma_t + gamma_st
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
}



Zhou <- function(d, gamma, pars, return = "value"){
  
  v_model <- ifelse(d[, 1] == 0 & d[, 2] == 0 & d[, 3] == 0,
                    pars[1] + pars[2],
                    pars[1] +
                      pars[2]*exp(-sqrt((d[, 1]/par[3])^2 +
                                          (d[, 2]/par[4])^2 +
                                          (d[, 3]/par[5])^2))) 
  
  if(return == "value"){
    return(v_model)
  }else if(return == "error"){
    return(sum((gamma - v_model)^2))
  }
  
  
}







