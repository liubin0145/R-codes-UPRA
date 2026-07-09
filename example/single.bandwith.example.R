

source("my.function.single.R")
source("argmax.R")
source("covariance.est.R")
source("mosum.single.R")
library(parallel)

para <- list(
  type="sign",
  paral="F",
  n = 600,
  p = 20,
  B = 50,
  bandwith=80,
  degree = 3,
  light_tail = T,
  cp_exit = TRUE,
  q0_list = c(0.3, 0.5, 0.7),
  outlier=F
)

dat <- Data_gen_multicpt(
  n = para$n,
  p = para$p,
  q0_list = para$q0_list,
  cp_exit = para$cp_exit,
  light_tail = para$light_tail,
  degree = para$degree,
  outlier=para$outlier
)

out <- mosum(
  dat$X,
  bandwith = para$bandwith,
  B = para$B,
  type = para$type,
  paral = para$paral
)







