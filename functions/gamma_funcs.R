disp_to_shape <- function(disp){
1/ exp(-7 + log(1 + exp(disp)))
}




pred_gammals <- function(mod, newdata, conf_int = 0.95){
  q_conf <- (1-conf_int)/2
  q_conf <- qnorm(q_conf, lower.tail = FALSE)


  ## add fit and se.fit on the **link** scale
  pred_est <- predict(mod, newdata, se.fit = TRUE, type = "link")

  est_mgcv <- newdata %>%
    bind_cols(as.data.frame(pred_est[[1]])) %>%
    rename(est_mean = V1) %>%
    rename(est_dispersion = V2) %>%
      bind_cols(as.data.frame(pred_est[[2]])) %>%
      rename(se_mean = V1) %>%
      rename(se_dispersion = V2)

  ## create the interval and backtransform
      est_mgcv <- est_mgcv %>%
        mutate(est_mean  = est_mean,
            mean_upr = est_mean + (q_conf * se_mean),
            mean_lwr = est_mean - (q_conf * se_mean)) %>%
        mutate(est_mean = exp(est_mean),
              mean_upr = exp(mean_upr),
              mean_lwr = exp(mean_lwr))

      est_mgcv <- est_mgcv %>%
            mutate(est_dispersion  = est_dispersion,
                dispersion_upr = est_dispersion + (q_conf * se_dispersion),
                dispersion_lwr = est_dispersion - (q_conf * se_dispersion)) %>%
            mutate(est_shape = disp_to_shape(est_dispersion),
                  shape_upr = disp_to_shape(dispersion_upr),
                  shape_lwr = disp_to_shape(dispersion_lwr))

    return(est_mgcv)
}



draw_from_gammals <- function(mod, n, newdata){

  a <- nrow(newdata)

  Cg <- predict(mod, newdata, type = "lpmatrix")
  Vb <- vcov(mod)
  sims <- rmvn(n, mu = coef(mod), sig = Vb)

  int_cols <- which(startsWith(colnames(sims), "(Int") == TRUE)

  mean_cols <- int_cols[1]:(int_cols[2]-1)
  disp_cols <- int_cols[2]:length(colnames(sims))

  fits_mean <- Cg[,mean_cols] %*% t(sims)[mean_cols,]
  fits_disp <- Cg[,disp_cols] %*% t(sims)[disp_cols,]

  for(j in seq(1,n)){
    pred_j <- newdata %>%
      mutate(draw = j) %>%
      mutate(mean_link = fits_mean[,j]) %>%
      mutate(disp_link = fits_disp[,j])
    if (j == 1){
      pred_df <- pred_j
    } else {
      pred_df <- pred_df %>% bind_rows(pred_j)
    }
  }

  pred_df <- pred_df %>%
    mutate(mean = exp(mean_link)) %>%
    mutate(shape = disp_to_shape(disp_link)) %>%
  	mutate(scale = mean/shape)

  return(pred_df)
}


summarize_gammals_draw <- function(model_draws, var = "quantile", level = 0.5, conf_int = 0.95){

  ### Check in case there are columns with these names
  model_draws <- model_draws %>%
    select_if(!names(.) %in% c('est', 'slope_est', 'delta_est', 'delta_per'))


  if(var == "quantile"){
    model_draws$est <- NA_real_
    for (j in seq(1,dim(model_draws)[1])){
      model_draws$est[[j]] <-  qgamma(level, scale = model_draws$scale[[j]], shape  = model_draws$shape[[j]])
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
