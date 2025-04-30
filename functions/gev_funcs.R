
pred_gev <- function(mod, newdata, conf_int = 0.95){
  q_conf <- (1-conf_int)/2
  q_conf <- qnorm(q_conf, lower.tail = FALSE)


  ## add fit and se.fit on the **link** scale
  pred_est <- predict(mod, newdata, se.fit = TRUE, type = "link")

  est_mgcv <- newdata %>%
    bind_cols(as.data.frame(pred_est[[1]])) %>%
    rename(est_location = V1) %>%
    rename(est_scale = V2) %>%
    rename(est_shape = V3) %>%
      bind_cols(as.data.frame(pred_est[[2]])) %>%
      rename(se_location = V1) %>%
      rename(se_scale = V2) %>%
      rename(se_shape = V3)

  ## create the interval and backtransform
      est_mgcv <- est_mgcv %>%
        mutate(est_location  = est_location,
            location_upr = est_location + (q_conf * se_location),
            location_lwr = est_location - (q_conf * se_location))

      est_mgcv <- est_mgcv %>%
            mutate(est_scale  = est_scale,
                scale_upr = est_scale + (q_conf * se_scale),
                scale_lwr = est_scale - (q_conf * se_scale)) %>%
            mutate(est_scale = exp(est_scale),
                  scale_upr = exp(scale_upr),
                  scale_lwr = exp(scale_lwr))

        est_mgcv <- est_mgcv %>%
              mutate(est_shape  = est_shape,
                shape_upr = est_shape + (q_conf * se_shape),
                shape_lwr = est_shape - (q_conf * se_shape))

    return(est_mgcv)
}




### Function to calculate mean of a GEV distribution
gev_mean <- function(loc, scale, shape){
  if(shape >= 1){
    mean_est <- inf
  } else if (shape == 0) {
    mean_est <- loc + scale * exp(1)
  } else if (shape >= 0 & shape < 1) {
    g_1 <- gamma(1-shape)
    mean_est <- loc + (scale *(g_1-1))/shape
  }
  return(mean_est )
}


### Function to calculate variance of a GEV distribution
gev_var <- function(loc, scale, shape){
  if(shape >= 0.5){
    var_est <- inf
  } else if (shape == 0) {
    var_est <- (scale^2)*((pi^2)/6)
  } else if (shape >= 0 & shape < 0.5) {
    g_1 <- gamma(1-shape)
    g_2 <- gamma(1-2*shape)
    var_est <- g_2 - (g_1^2)
    var_est <- ((scale_x^2)*var_est)/(shape_x^2)
  }
  return(var_est )
}




draw_from_gev <- function(mod, n, newdata){

  a <- nrow(newdata)

  Cg <- predict(mod, newdata, type = "lpmatrix")
  Vb <- vcov(mod)
  sims <- rmvn(n, mu = coef(mod), sig = Vb)

  int_cols <- which(startsWith(colnames(sims), "(Int") == TRUE)
  
  loc_cols <- which(startsWith(colnames(sims), "s(") == TRUE)
  loc_cols <- c(int_cols[[1]], loc_cols)
  
  scale_cols <- which(startsWith(colnames(sims), "s.1(") == TRUE)
  scale_cols <- c(int_cols[[2]], scale_cols)
  
  shape_cols <- which(startsWith(colnames(sims), "s.2(") == TRUE)
  shape_cols <- c(int_cols[[3]], shape_cols)

  fits_loc <- Cg[,loc_cols] %*% t(sims)[loc_cols,]
  fits_scale <- Cg[,scale_cols] %*% t(sims)[scale_cols,]
  fits_shape <- Cg[,shape_cols] %*% t(sims)[shape_cols,]

  for(j in seq(1,n)){
    pred_j <- newdata %>%
      mutate(draw = j) %>%
      mutate(loc_link = fits_loc[,j]) %>%
      mutate(scale_link = fits_scale[,j]) %>%
      mutate(shape_link = fits_shape[,j])
    if (j == 1){
      pred_df <- pred_j
    } else {
      pred_df <- pred_df %>% bind_rows(pred_j)
    }
  }

  pred_df <- pred_df %>%
    mutate(loc = loc_link) %>%
    mutate(scale = exp(scale_link)) %>%
    mutate(shape = shape_link)

  return(pred_df)
}



draw_from_gev_linear <- function(mod, n, newdata){

  a <- nrow(newdata)

  Cg <- predict(mod, newdata, type = "lpmatrix")
  Vb <- vcov(mod)
  sims <- rmvn(n, mu = coef(mod), sig = Vb)

  int_cols <- which(startsWith(colnames(sims), "(Int") == TRUE)
  
  loc_cols <- which(colnames(sims) == "year")
  loc_cols <- c(int_cols[[1]], loc_cols)
  
  scale_cols <- which(colnames(sims) == "year.1")
  scale_cols <- c(int_cols[[2]], scale_cols)
  
  shape_cols <- which(colnames(sims) == "year.2")
  shape_cols <- c(int_cols[[3]], shape_cols)

  fits_loc <- Cg[,loc_cols] %*% t(sims)[loc_cols,]
  fits_scale <- Cg[,scale_cols] %*% t(sims)[scale_cols,]
  fits_shape <- Cg[,shape_cols] %*% t(sims)[shape_cols,]

  for(j in seq(1,n)){
    pred_j <- newdata %>%
      mutate(draw = j) %>%
      mutate(loc_link = fits_loc[,j]) %>%
      mutate(scale_link = fits_scale[,j]) %>%
      mutate(shape_link = fits_shape[,j])
    if (j == 1){
      pred_df <- pred_j
    } else {
      pred_df <- pred_df %>% bind_rows(pred_j)
    }
  }

  pred_df <- pred_df %>%
    mutate(loc = loc_link) %>%
    mutate(scale = exp(scale_link)) %>%
    mutate(shape = shape_link)

  return(pred_df)
}



summarize_gev_draw <- function(model_draws, var = "quantile", level = 0.5, conf_int = 0.95){

  ### Check in case there are columns with these names
  model_draws <- model_draws %>%
    select_if(!names(.) %in% c('est', 'slope_est', 'delta_est', 'delta_per'))


  if(var == "quantile"){
    model_draws$est <- NA_real_
    for (j in seq(1,dim(model_draws)[1])){
      model_draws$est[[j]] <-  evd::qgev(level, loc = model_draws$loc[[j]], scale = model_draws$scale[[j]], shape  = model_draws$shape[[j]])
    }
  } else{
    model_draws <- model_draws %>%
      rename("est" = all_of(var))
  }


  model_draws <- model_draws %>%
    mutate(delta_est = est - lag(est)) %>%
    mutate(delta_per = year - lag(year)) %>%
    mutate(slope_est = delta_est / delta_per)

    q_conf <- (1-conf_int)/2

  summary_df <- model_draws %>%
    #drop_na(year, est, slope_est) %>%
    group_by(year) %>%
    summarize(est_median = quantile(est, 0.5, na.rm=TRUE),
      est_lower = quantile(est, q_conf, na.rm=TRUE),
      est_upper = quantile(est, 1-q_conf, na.rm=TRUE),
      slope_median = quantile(slope_est, 0.5, na.rm=TRUE),
      slope_lower = quantile(slope_est, q_conf, na.rm=TRUE),
      slope_upper = quantile(slope_est, 1-q_conf, na.rm=TRUE)) %>%
    ungroup()

    return(list(summary = summary_df, model_draws = model_draws))
}
