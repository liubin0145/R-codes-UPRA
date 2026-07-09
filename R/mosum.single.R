
build_local_windows <- function(cp, n, L_min = 60, L_max = 150, rho = 0.05) {
  if (length(cp) == 0) return(list())
  
  cp <- sort((cp))
  cp_full <- c(1, cp, n)
  out <- vector("list", length(cp))
  
  for (i in seq_along(cp)) {
    center <- cp[i]
    pos <- i + 1   # 
    
    left_all  <- cp_full[1:(pos - 1)]
    right_all <- cp_full[(pos + 1):length(cp_full)]
    
    left_dist  <- center - left_all
    right_dist <- right_all - center
    
    ## 
    left_ok  <- which(left_dist  >= L_min)
    right_ok <- which(right_dist >= L_min)
    
   # left_anchor  <- if (length(left_ok)  > 0) left_all[max(left_ok)]   else cp_full[pos - 1]
  #  right_anchor <- if (length(right_ok) > 0) right_all[min(right_ok)] else cp_full[pos + 1]
    
    left_anchor  <- if (length(left_ok)  > 0) left_all[max(left_ok)]   else max(1, center - L_min)
    right_anchor <- if (length(right_ok) > 0) right_all[min(right_ok)] else min(n, center + L_min) 
    
    
    raw_minus <- center - left_anchor
    raw_plus  <- right_anchor - center
    
    ## 
    L_minus <- min(raw_minus, L_max)
    L_plus  <- min(raw_plus,  L_max)
    
    ## 
    s1 <- max(1, floor(center - (1 - rho) * L_minus))
    e1 <- min(center - 1, floor(center - rho * L_minus))
    
    s2 <- max(center + 1, ceiling(center + rho * L_plus))
    e2 <- min(n, ceiling(center + (1 - rho) * L_plus))
    
    ## 
    if (e1 < s1) {
      s1 <- max(1, center - 1)
      e1 <- max(1, center - 1)
    }
    if (e2 < s2) {
      s2 <- min(n, center + 1)
      e2 <- min(n, center + 1)
    }
    
    ## 
    b_eff <- max(2, min(e1 - s1 + 1, e2 - s2 + 1))
    
    ## 
    search_left  <- max(2, min(floor(raw_minus / 2), L_max))
    search_right <- max(2, min(floor(raw_plus  / 2), L_max))
    
    out[[i]] <- list(
      center = center,
      left_anchor = left_anchor,
      right_anchor = right_anchor,
      raw_minus = raw_minus,
      raw_plus = raw_plus,
      L_minus = L_minus,
      L_plus = L_plus,
      s1 = s1, e1 = e1,
      s2 = s2, e2 = e2,
      b_eff = b_eff,
      search_left = search_left,
      search_right = search_right
    )
  }
  
  out
}


choose_c_by_ratio <- function(theta.vec, b_eff, d, n, c.grid) {
  abs_theta <- abs(theta.vec)
  ord <- order(abs_theta, decreasing = TRUE)
  abs_sorted <- abs_theta[ord]
  
  ratio_df <- data.frame(
    c_value = c.grid,
    w_plus = NA_real_,
    k_selected = NA_integer_,
    ratio_value = NA_real_
  )
  
  for (ll in seq_along(c.grid)) {
    c_val <- c.grid[ll]
    w_plus <- c_val * sqrt(log(d * n) / (2 * b_eff))
    
    selected_idx <- which(abs_theta > w_plus)
    k_sel <- length(selected_idx)
    
    ratio_df$w_plus[ll] <- w_plus
    ratio_df$k_selected[ll] <- k_sel
    
    ## 
    if (k_sel >= 1 && (2 * k_sel) <= length(abs_sorted)) {
      num <- sqrt(sum(abs_sorted[1:k_sel]^2))
      den <- sqrt(sum(abs_sorted[(k_sel + 1):(2 * k_sel)]^2))
      
      if (den > 0) {
        ratio_df$ratio_value[ll] <- num / den
      } else {
        ratio_df$ratio_value[ll] <- Inf
      }
    } else {
      ratio_df$ratio_value[ll] <- NA_real_
    }
  }
  
  ## 
  if (all(is.na(ratio_df$ratio_value))) {
    best_idx <- which.min(abs(c.grid - 1))  # 
    selected_idx <- ord[1]
    return(list(
      c_best = ratio_df$c_value[best_idx],
      w_plus_best = ratio_df$w_plus[best_idx],
      k_best = 1,
      gamma_idx = selected_idx,
      ratio_best = NA_real_,
      ratio_df = ratio_df
    ))
  }
  
  best_idx <- which.max(ratio_df$ratio_value)
  c_best <- ratio_df$c_value[best_idx]
  w_plus_best <- ratio_df$w_plus[best_idx]
  gamma_idx <- which(abs_theta > w_plus_best)
  k_best <- length(gamma_idx)
  
  ## 
  if (k_best == 0) {
    gamma_idx <- ord[1]
    k_best <- 1
  }
  
  list(
    c_best = c_best,
    w_plus_best = w_plus_best,
    k_best = k_best,
    gamma_idx = gamma_idx,
    ratio_best = ratio_df$ratio_value[best_idx],
    ratio_df = ratio_df
  )
}

get_true_change_points <- function(n, q0_list) {
  q0_list <- sort(unique(q0_list))
  cp_inner <- floor(q0_list * n)
  cp_inner <- unique(cp_inner[cp_inner > 1 & cp_inner < n])
  c(1, cp_inner, n)
}


linear.mosum=function(X, bandwith=80, B=100, paral="F"){
  source("my.function.single.R")
  n=nrow(X)
  d=ncol(X)
  h=sqrt(log(d)/n)
  b=bandwith
  sigma=rep(1,d) #variance of the kernel 
  K=seq(b+1,n-b-1,1) #search domain
  inner.dimen=function(j,K){
    mosum.cusum.inner=numeric()
    mosum.cusum.inner[1]=1/(2*b)^(3/2)*sum((outer(X[(K[1]+1):min(n,(K[1]+b)),j],X[max(1,(K[1]-b)):K[1],j],"-")))
    for(k in 1:(length(K)-1)){
      mosum.cusum.inner[k+1]=mosum.cusum.inner[k]+1/(2*b)^(3/2)*(-sum((outer(X[(K[k]+1):(K[k]+1),j],X[max(1,(K[k]-b)):max(1,(K[k]-b)),j],"-")))
                                                                 -sum((outer(X[(K[k]+2):min(n,(K[k]+b)),j],X[max(1,(K[k]-b)):max(1,(K[k]-b)),j],"-")))
                                                                 -sum((outer(X[(K[k]+1):(K[k]+1),j],X[max(1,(K[k]-b+1)):(K[k]),j],"-")))
                                                               +sum((outer(X[min(n,(K[k]+b+1)):min(n,(K[k]+b+1)),j],X[max(1,(K[k]-b+1)):(K[k]),j],"-")))
                                                                 +sum((outer(X[(K[k]+2):min(n,(K[k]+b)),j],X[(K[k]+1):(K[k]+1),j],"-")))    
                                                                 +sum((outer(X[min(n,(K[k]+b+1)):min(n,(K[k]+b+1)),j],X[(K[k]+1):(K[k]+1),j],"-"))))
   
      
      
      
       }
    return(mosum.cusum.inner)
  }
  
mosum.cusum=matrix(0,nrow=length(K),ncol=d)

temp= lapply(1:d,inner.dimen,K)
for(i in 1:d){mosum.cusum[,i]=temp[[i]]}
mosum.cusum=abs(mosum.cusum)
mosum.cusum=as.matrix(mosum.cusum)
mosum.aggre=apply(mosum.cusum,1,max)
mosum.sta=max(mosum.aggre)

mosum.boot.quantile=linear.mosum.boot(X,bandwith=b,B=B,paral=paral)

#abline(h=mosum.boot.quantile)
reject=as.numeric(mosum.sta>mosum.boot.quantile)
G=b
#==== obtain initial estimator===========================
  T_values=c(rep(0,G),mosum.aggre,rep(0,G+1))# 
  n <- length(T_values)  # 
  significant_windows <- list()  # 
  eta=0.15
  c_alpha_W=mosum.boot.quantile
  
  for (v_j in G:(n - G)) {
    for (w_j in (v_j + ceiling(eta * G)):(n - G)) {
    
      if (all(T_values[v_j:w_j] >= c_alpha_W)) {
        if ((v_j > G && T_values[v_j - 1] < c_alpha_W) && 
            (w_j < n - G && T_values[w_j + 1] < c_alpha_W)) {
          significant_windows <- append(significant_windows, list(c(v_j, w_j)))
          break  
        }
      }
    }
  }
  

  change_points <- sapply(significant_windows, function(window) {
    v_j <- window[1]
    w_j <- window[2]
    v_j + which.max(T_values[v_j:w_j]) - 1  # 
  })
  
 
cpt.est.initial=c(1,sort(change_points),n)

cpt.est.actual.global <- cpt.est.initial[-c(1, length(cpt.est.initial))]
local.win.initial.global <- build_local_windows(
  cp = cpt.est.actual.global,
  n = n
)

#==== obtain the estimator for the coordiates
cpt.est.actual=cpt.est.initial[-c(1,length(cpt.est.initial))]
cpt.est.actual <- cpt.est.initial[-c(1, length(cpt.est.initial))]
local.win.initial <- local.win.initial.global

if(length(cpt.est.actual) >= 1){
  theta.had <- matrix(0, nrow = d, ncol = length(cpt.est.actual))
  Gamma.head <- vector("list", length(cpt.est.actual))
  
  c.selected <- rep(NA_real_, length(cpt.est.actual))
  w_plus.selected <- rep(NA_real_, length(cpt.est.actual))
  k.selected <- rep(NA_integer_, length(cpt.est.actual))
  ratio.best <- rep(NA_real_, length(cpt.est.actual))
  ratio.info <- vector("list", length(cpt.est.actual))
  
  for(i in 1:length(cpt.est.actual)){
    s1 <- local.win.initial[[i]]$s1
    e1 <- local.win.initial[[i]]$e1
    s2 <- local.win.initial[[i]]$s2
    e2 <- local.win.initial[[i]]$e2
    cat("s1,e1=",c(s1,e1),"\n")
    cat("s2,e2=",c(s2,e2),"\n")
    for(j in 1:d){
      theta.had[j, i] <- 1 / ((e1 - s1 + 1) * (e2 - s2 + 1)) *
        sum(outer(X[s2:e2, j], X[s1:e1, j], "-"))
    }
    
    
    
    abs_theta_i <- sort(abs(theta.had[, i]), decreasing = TRUE)
    base_i <- sqrt(log(d * n) / (2 * local.win.initial[[i]]$b_eff))
    
    K_min <- 2
    K_max <- floor(2*log(d))
    c_lower <- abs_theta_i[K_max + 1] / base_i   # 
    c_upper <- abs_theta_i[K_min] / base_i       # 
    choose_res <- choose_c_by_ratio(
      theta.vec = theta.had[, i],
      b_eff = local.win.initial[[i]]$b_eff,
      d = d,
      n = n,
      c.grid = seq(c_lower,c_upper,0.5)
    )
    Gamma.head[[i]] <- choose_res$gamma_idx
    c.selected[i] <- choose_res$c_best
    w_plus.selected[i] <- choose_res$w_plus_best
    k.selected[i] <- choose_res$k_best
    ratio.best[i] <- choose_res$ratio_best
    ratio.info[[i]] <- choose_res$ratio_df
  }} else {
    Gamma.head <- NULL
    theta.had <- NULL
    c.selected <- NULL
    w_plus.selected <- NULL
    k.selected <- NULL
    ratio.best <- NULL
    ratio.info <- NULL
  }

#==== local refinement =====     
cpt.est.refine <- numeric()
if(length(cpt.est.actual) >= 1){
  for(i in 1:length(cpt.est.actual)){
    left_search  <- local.win.initial[[i]]$search_left
    right_search <- local.win.initial[[i]]$search_right
    
    K.new <- seq(
      max(cpt.est.actual[i] - left_search + 1, 2),
      min(cpt.est.actual[i] + right_search - 1, n - 1),
      1
    )
    
    mosum.cusum.new <- matrix(0, nrow = length(K.new), ncol = d)
    temp <- lapply(1:d, inner.dimen, K.new)
    for(j in 1:d){
      mosum.cusum.new[, j] <- theta.had[j, i] * temp[[j]]
    }
    mosum.aggre.new <- apply(mosum.cusum.new[, Gamma.head[[i]], drop = FALSE], 1, sum)
    cpt.est.refine[i] <- K.new[which.max(mosum.aggre.new)]
  }
  cpt.est.refine <- c(1, sort(cpt.est.refine), n)
} else {
  cpt.est.refine <- c(1, n)
}


#===== obtain confidence interval ==============
source("argmax.R")
source("covariance.est.R")

cpt.est.refine.actual <- cpt.est.refine[-c(1, length(cpt.est.refine))]
local.win.refine <- build_local_windows(
  cp = cpt.est.refine.actual,
  n = n
)

confidence <- matrix(NA, nrow = length(cpt.est.refine.actual), ncol = 2)
plugin_argmax_quantiles <- vector("list", length(cpt.est.refine.actual))
argmax_probs <- quantile_probs_default()

for(i in 1:length(cpt.est.refine.actual)){
  s1 <- local.win.refine[[i]]$s1
  e1 <- local.win.refine[[i]]$e1
  s2 <- local.win.refine[[i]]$s2
  e2 <- local.win.refine[[i]]$e2
  cat("s1,e1=",c(s1,e1),"\n")
  cat("s2,e2=",c(s2,e2),"\n")
  X.before <- X[s1:e1, Gamma.head[[i]], drop = FALSE]
  X.after  <- X[s2:e2, Gamma.head[[i]], drop = FALSE]
  
  cov.est <- linear.covariance.est(
    X.before, X.after,
    theta = theta.had[Gamma.head[[i]], i]
  )
  
  plugin_argmax_quantiles[[i]] <- simulate_argmax_quantiles(
    theta  = theta.had[Gamma.head[[i]], i],
    Sigma1 = cov.est[[1]],
    Sigma2 = cov.est[[2]],
    Sigma3 = cov.est[[3]],
    Sigma4 = cov.est[[4]],
    probs = argmax_probs,
    nsim = 2000
  )
  
  idx_ci <- match(c(0.025, 0.975), argmax_probs)
  argmax.quantile <- plugin_argmax_quantiles[[i]][idx_ci]
  confidence[i, ] <- cpt.est.refine.actual[i] -
    1 * rev(argmax.quantile) * (norm(theta.had[Gamma.head[[i]], i], "2"))^(-2)
}


final.res=list(
mosum.sta=mosum.sta,
mosum.reject=reject,
cpt.est.initial=cpt.est.initial,
cpt.est.refine=cpt.est.refine,
cpt.conf=confidence
)

return(final.res)
}




sign.mosum=function(X, bandwith=80, B, paral="F"){
  source("my.function.single.R")
  n=nrow(X)
  d=ncol(X)
  h=sqrt(log(d)/n)
  b=bandwith#floor(2*h*n)   #bandwith
  sigma=rep(1,d) #variance of the kernel 
  K=seq(b+1,n-b-1,1) #search domain
  inner.dimen=function(j,K){
    mosum.cusum.inner=numeric()
    mosum.cusum.inner[1]=1/(2*b)^(3/2)*sum(sign(outer(X[(K[1]+1):min(n,(K[1]+b)),j],X[max(1,(K[1]-b)):K[1],j],"-")))
    for(k in 1:(length(K)-1)){
      mosum.cusum.inner[k+1]=mosum.cusum.inner[k]+1/(2*b)^(3/2)*(-sum(sign(outer(X[(K[k]+1):(K[k]+1),j],X[max(1,(K[k]-b)):max(1,(K[k]-b)),j],"-")))
                                                                 -sum(sign(outer(X[(K[k]+2):min(n,(K[k]+b)),j],X[max(1,(K[k]-b)):max(1,(K[k]-b)),j],"-")))
                                                                 -sum(sign(outer(X[(K[k]+1):(K[k]+1),j],X[max(1,(K[k]-b+1)):(K[k]),j],"-")))
                                                                 +sum(sign(outer(X[min(n,(K[k]+b+1)):min(n,(K[k]+b+1)),j],X[max(1,(K[k]-b+1)):(K[k]),j],"-")))
                                                                 +sum(sign(outer(X[(K[k]+2):min(n,(K[k]+b)),j],X[(K[k]+1):(K[k]+1),j],"-")))
                                                                 +sum(sign(outer(X[min(n,(K[k]+b+1)):min(n,(K[k]+b+1)),j],X[(K[k]+1):(K[k]+1),j],"-"))))
    }
    return(mosum.cusum.inner)
  }
  
  mosum.cusum=matrix(0,nrow=length(K),ncol=d)
  
  temp= lapply(1:d,inner.dimen,K)
  for(i in 1:d){mosum.cusum[,i]=temp[[i]]}
  mosum.cusum=abs(mosum.cusum)
  mosum.cusum=as.matrix(mosum.cusum)
  mosum.aggre=apply(mosum.cusum,1,max)
  mosum.sta=max(mosum.aggre)
  #plot(K,mosum.aggre)
  mosum.boot.quantile=sign.mosum.boot(X,bandwith=b,B=B,paral=paral)
  #abline(h=mosum.boot.quantile)
  reject=as.numeric(mosum.sta>mosum.boot.quantile)
  G=b
  #==== obtain initial estimator===========================
  T_values=c(rep(0,G),mosum.aggre,rep(0,G+1))# 
  n <- length(T_values)  # 
  significant_windows <- list()  # 
  eta=0.15
  c_alpha_W=mosum.boot.quantile
  for (v_j in G:(n - G)) {
    for (w_j in (v_j + ceiling(eta * G)):(n - G)) {
      if (all(T_values[v_j:w_j] >= c_alpha_W)) {
        if ((v_j > G && T_values[v_j - 1] < c_alpha_W) && 
            (w_j < n - G && T_values[w_j + 1] < c_alpha_W)) {
          significant_windows <- append(significant_windows, list(c(v_j, w_j)))
          break  
        }
      }
    }
  }
  
  change_points <- sapply(significant_windows, function(window) {
    v_j <- window[1]
    w_j <- window[2]
    v_j + which.max(T_values[v_j:w_j]) - 1  
  })
  
  
  cpt.est.initial=c(1,sort(change_points),n)
  
  
  cpt.est.actual.global <- cpt.est.initial[-c(1, length(cpt.est.initial))]
  local.win.initial.global <- build_local_windows(
    cp = cpt.est.actual.global,
    n = n
  )
  
  
  #==== obtain the estimator for the coordiates
  cpt.est.actual <- cpt.est.initial[-c(1, length(cpt.est.initial))]
  w_plus=sqrt(log(d*n)/(2*b))
  local.win.initial <- local.win.initial.global
  if(length(cpt.est.actual) >= 1){
    theta.had <- matrix(0, nrow = d, ncol = length(cpt.est.actual))
    Gamma.head <- vector("list", length(cpt.est.actual))
    
    c.selected <- rep(NA_real_, length(cpt.est.actual))
    w_plus.selected <- rep(NA_real_, length(cpt.est.actual))
    k.selected <- rep(NA_integer_, length(cpt.est.actual))
    ratio.best <- rep(NA_real_, length(cpt.est.actual))
    ratio.info <- vector("list", length(cpt.est.actual))
    
    for(i in 1:length(cpt.est.actual)){
      s1 <- local.win.initial[[i]]$s1
      e1 <- local.win.initial[[i]]$e1
      s2 <- local.win.initial[[i]]$s2
      e2 <- local.win.initial[[i]]$e2
      
      cat("s1,e1=",c(s1,e1),"\n")
      cat("s2,e2=",c(s2,e2),"\n")
      
      
      for(j in 1:d){
        theta.had[j, i] <- 1 / ((e1 - s1 + 1) * (e2 - s2 + 1)) *
          sum(sign(outer(X[s2:e2, j], X[s1:e1, j], "-")))
      }
      
      abs_theta_i <- sort(abs(theta.had[, i]), decreasing = TRUE)
      base_i <- sqrt(log(d * n) / (2 * local.win.initial[[i]]$b_eff))
      K_min <- 2
      K_max <- floor(2*log(d))
      c_lower <- abs_theta_i[K_max + 1] / base_i   
      c_upper <- abs_theta_i[K_min] / base_i       
      choose_res <- choose_c_by_ratio(
        theta.vec = theta.had[, i],
        b_eff = local.win.initial[[i]]$b_eff,
        d = d,
        n = n,
        c.grid = seq(c_lower,c_upper,0.5)
      )
      Gamma.head[[i]] <- choose_res$gamma_idx
      c.selected[i] <- choose_res$c_best
      w_plus.selected[i] <- choose_res$w_plus_best
      k.selected[i] <- choose_res$k_best
      ratio.best[i] <- choose_res$ratio_best
      ratio.info[[i]] <- choose_res$ratio_df
    }} else {
      Gamma.head <- NULL
      theta.had <- NULL
      c.selected <- NULL
      w_plus.selected <- NULL
      k.selected <- NULL
      ratio.best <- NULL
      ratio.info <- NULL
    }
  #==== local refinement =====    
  cpt.est.refine <- numeric()
  if(length(cpt.est.actual) >= 1){
    for(i in 1:length(cpt.est.actual)){
      left_search  <- local.win.initial[[i]]$search_left
      right_search <- local.win.initial[[i]]$search_right
      
      K.new <- seq(
        max(cpt.est.actual[i] - left_search + 1, 2),
        min(cpt.est.actual[i] + right_search - 1, n - 1),
        1
      )
      
      mosum.cusum.new <- matrix(0, nrow = length(K.new), ncol = d)
      temp <- lapply(1:d, inner.dimen, K.new)
      for(j in 1:d){
        mosum.cusum.new[, j] <- theta.had[j, i] * temp[[j]]
      }
      mosum.aggre.new <- apply(mosum.cusum.new[, Gamma.head[[i]], drop = FALSE], 1, sum)
      cpt.est.refine[i] <- K.new[which.max(mosum.aggre.new)]
    }
    cpt.est.refine <- c(1, sort(cpt.est.refine), n)
  } else {
    cpt.est.refine <- c(1, n)
  }
  
  #===== obtain confidence interval ==============
  source("argmax.R")
  source("covariance.est.R")
  cpt.est.refine.actual <- cpt.est.refine[-c(1, length(cpt.est.refine))]
  local.win.refine <- build_local_windows(
    cp = cpt.est.refine.actual,
    n = n
  )
  
  confidence <- matrix(NA, nrow = length(cpt.est.refine.actual), ncol = 2)
  plugin_argmax_quantiles <- vector("list", length(cpt.est.refine.actual))
  argmax_probs <- quantile_probs_default()
  
  for(i in 1:length(cpt.est.refine.actual)){
    s1 <- local.win.refine[[i]]$s1
    e1 <- local.win.refine[[i]]$e1
    s2 <- local.win.refine[[i]]$s2
    e2 <- local.win.refine[[i]]$e2
    cat("s1,e1=",c(s1,e1),"\n")
    cat("s2,e2=",c(s2,e2),"\n")
    X.before <- X[s1:e1, Gamma.head[[i]], drop = FALSE]
    X.after  <- X[s2:e2, Gamma.head[[i]], drop = FALSE]
   
    cov.est <- sign.covariance.est(
      X.before, X.after,
      theta = theta.had[Gamma.head[[i]], i]
    )
    
    plugin_argmax_quantiles[[i]] <- simulate_argmax_quantiles(
      theta  = theta.had[Gamma.head[[i]], i],
      Sigma1 = cov.est[[1]],
      Sigma2 = cov.est[[2]],
      Sigma3 = cov.est[[3]],
      Sigma4 = cov.est[[4]],
      probs = argmax_probs,
      nsim = 2000
    )
    
    idx_ci <- match(c(0.025, 0.975), argmax_probs)
    argmax.quantile <- plugin_argmax_quantiles[[i]][idx_ci]
    confidence[i, ] <- cpt.est.refine.actual[i] -
      1 * rev(argmax.quantile) * (norm(theta.had[Gamma.head[[i]], i], "2"))^(-2)
  }
  
 
  final.res=list(
    mosum.sta=mosum.sta,
    mosum.reject=reject,
    cpt.est.initial=cpt.est.initial,
    cpt.est.refine=cpt.est.refine,
    cpt.conf=confidence
  )
  return(final.res)
}

  
  
  
  
  
  


mosum=function(X, bandwith=80, B=100, type="linear", paral="F"){
  if(type=="linear"){
    res=linear.mosum(X, bandwith, B, paral=paral)
  }else if(type=="sign"){
    res=sign.mosum(X, bandwith, B, paral=paral)
  }
  return(res)
}




