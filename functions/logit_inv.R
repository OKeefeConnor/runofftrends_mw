logit_inv <- function(est){
  1.5 * binomial()$linkinv(est) - 1
}
