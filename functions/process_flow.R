
###########################################################################
### Function to process flow
###########################################################################

process_flow <- function(flow){

  ### expand to make sure it contains all dates
  flow <- flow %>%
    drop_na(flow)
  dates_j <- data.frame(date = seq(as.Date(min(flow$date)), as.Date(max(flow$date)), by = "1 day"))
  flow <- dates_j %>%
    left_join(flow, by = "date")

  ### Calculate the daily difference in flow
  flow <- flow %>%
  	mutate(jdate = yday(date), month = month(date), year = year(date)) %>%
  	mutate(prev_flow = lag(flow, 1)) %>%
  	mutate(diff = flow - prev_flow) %>%
    drop_na(flow) %>%
    drop_na(diff)  %>%
    mutate(prev_diff = lag(diff,1)) %>%
    mutate(inc = diff > 0, prev_inc = prev_diff > 0)

    ### Extract the increases
    inc_df <- flow %>%
      filter(diff > 0) %>%
      select(date, flow, jdate, month, year, prev_flow, diff)

    ### Extract the decreases
    dec_df <- flow %>%
        filter(diff <= 0) %>%
        select(date, flow, jdate, month, year, prev_flow, diff) %>%
        mutate(b_coef = -log(flow/prev_flow))

    ### Extract the transition likelihood
    trans_df <- flow %>%
        select(date, flow, jdate, month, year, inc, prev_inc)

return(list(inc = inc_df, dec = dec_df, trans = trans_df))

}




process_flow_runoff <- function(flow_mm_day){

  ### expand to make sure it contains all dates
  flow_mm_day <- flow_mm_day %>%
    drop_na(flow_mm_day)
  dates_j <- data.frame(date = seq(as.Date(min(flow_mm_day$date)), as.Date(max(flow_mm_day$date)), by = "1 day"))
  flow_mm_day <- dates_j %>%
    left_join(flow_mm_day, by = "date")

  ### Calculate the daily difference in flow_mm_day
  flow_mm_day <- flow_mm_day %>%
  	mutate(jdate = yday(date), month = month(date), year = year(date)) %>%
  	mutate(prev_flow = lag(flow_mm_day, 1)) %>%
  	mutate(diff = flow_mm_day - prev_flow) %>%
    drop_na(flow_mm_day) %>%
    drop_na(diff)  %>%
    mutate(prev_diff = lag(diff,1)) %>%
    mutate(inc = diff > 0, prev_inc = prev_diff > 0)

    ### Extract the increases
    inc_df <- flow_mm_day %>%
      filter(diff > 0) %>%
      select(date, flow_mm_day, jdate, month, year, prev_flow, diff)

    ### Extract the decreases
    dec_df <- flow_mm_day %>%
        filter(diff <= 0) %>%
        select(date, flow_mm_day, jdate, month, year, prev_flow, diff) %>%
        mutate(b_coef = -log(flow_mm_day/prev_flow))

    ### Extract the transition likelihood
    trans_df <- flow_mm_day %>%
        select(date, flow_mm_day, jdate, month, year, inc, prev_inc)

return(list(inc = inc_df, dec = dec_df, trans = trans_df))

}
