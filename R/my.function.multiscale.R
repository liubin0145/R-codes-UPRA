# linear.mosum is the function using the linear kernels
# mosum(X) is the function using the sign kernels

#library(foreach)
#library(doParallel)
#library(MASS)
#source("my.function.R")

# cp is the estimated change point, e.g, c(180,300,420)
# n is the sample size 
# L_min is the minimum sample size condition
#L_max is the maximum sample size condition
# rho is the removing parameter
build_local_windows <- function(cp, n, L_min = 60, L_max = 150, rho = 0.05) {
  if (length(cp) == 0) return(list())
  
  cp <- sort((cp))
  cp_full <- c(1, cp, n)
  out <- vector("list", length(cp))
  
  for (i in seq_along(cp)) {
    center <- cp[i]
    pos <- i + 1   
    
    left_all  <- cp_full[1:(pos - 1)]
    right_all <- cp_full[(pos + 1):length(cp_full)]
    
    left_dist  <- center - left_all
    right_dist <- right_all - center
    
   
    left_ok  <- which(left_dist  >= L_min)
    right_ok <- which(right_dist >= L_min)
    
    # left_anchor  <- if (length(left_ok)  > 0) left_all[max(left_ok)]   else cp_full[pos - 1]
    #  right_anchor <- if (length(right_ok) > 0) right_all[min(right_ok)] else cp_full[pos + 1]
    
    left_anchor  <- if (length(left_ok)  > 0) left_all[max(left_ok)]   else max(1, center - L_min)
    right_anchor <- if (length(right_ok) > 0) right_all[min(right_ok)] else min(n, center + L_min) 
    
    
    raw_minus <- center - left_anchor
    raw_plus  <- right_anchor - center
    
    
    L_minus <- min(raw_minus, L_max)
    L_plus  <- min(raw_plus,  L_max)
    
   
    s1 <- max(1, floor(center - (1 - rho) * L_minus))
    e1 <- min(center - 1, floor(center - rho * L_minus))
    
    s2 <- max(center + 1, ceiling(center + rho * L_plus))
    e2 <- min(n, ceiling(center + (1 - rho) * L_plus))
    
   
    if (e1 < s1) {
      s1 <- max(1, center - 1)
      e1 <- max(1, center - 1)
    }
    if (e2 < s2) {
      s2 <- min(n, center + 1)
      e2 <- min(n, center + 1)
    }
    
    
    b_eff <- max(2, min(e1 - s1 + 1, e2 - s2 + 1))
    
   
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
  
 
  if (all(is.na(ratio_df$ratio_value))) {
    best_idx <- which.min(abs(c.grid - 1))  
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






interval_overlap <- function(v1, w1, v2, w2) {
  !(w1 < v2 || w2 < v1)
}

make_inner_dimen <- function(X, b, kernel_type = c("linear", "sign")) {
  kernel_type <- match.arg(kernel_type)
  n <- nrow(X)
  
  if (kernel_type == "linear") {
    inner.dimen <- function(j, K) {
      mosum.cusum.inner <- numeric()
      mosum.cusum.inner[1] <- 1 / (2 * b)^(3/2) *
        sum(outer(X[max(1, (K[1] - b)):K[1], j],
                  X[(K[1] + 1):min(n, (K[1] + b)), j], "-"))
      for (k in 1:(length(K) - 1)) {
        mosum.cusum.inner[k + 1] <- mosum.cusum.inner[k] + 1 / (2 * b)^(3/2) * (
          -sum(outer(X[max(1, (K[k] - b)):max(1, (K[k] - b)), j],
                     X[(K[k] + 1):(K[k] + 1), j], "-"))
          -sum(outer(X[max(1, (K[k] - b)):max(1, (K[k] - b)), j],
                     X[(K[k] + 2):min(n, (K[k] + b)), j], "-"))
          -sum(outer(X[max(1, (K[k] - b + 1)):(K[k]), j],
                     X[(K[k] + 1):(K[k] + 1), j], "-"))
          +sum(outer(X[max(1, (K[k] - b + 1)):(K[k]), j],
                     X[min(n, (K[k] + b + 1)):min(n, (K[k] + b + 1)), j], "-"))
          +sum(outer(X[(K[k] + 1):(K[k] + 1), j],
                     X[(K[k] + 2):min(n, (K[k] + b)), j], "-"))
          +sum(outer(X[(K[k] + 1):(K[k] + 1), j],
                     X[min(n, (K[k] + b + 1)):min(n, (K[k] + b + 1)), j], "-"))
        )
      }
      mosum.cusum.inner
    }
  } else {
    inner.dimen <- function(j, K) {
      mosum.cusum.inner <- numeric()
      mosum.cusum.inner[1] <- 1 / (2 * b)^(3/2) *
        sum(sign(outer(X[max(1, (K[1] - b)):K[1], j],
                       X[(K[1] + 1):min(n, (K[1] + b)), j], "-")))
      for (k in 1:(length(K) - 1)) {
        mosum.cusum.inner[k + 1] <- mosum.cusum.inner[k] + 1 / (2 * b)^(3/2) * (
          -sum(sign(outer(X[max(1, (K[k] - b)):max(1, (K[k] - b)), j],
                          X[(K[k] + 1):(K[k] + 1), j], "-")))
          -sum(sign(outer(X[max(1, (K[k] - b)):max(1, (K[k] - b)), j],
                          X[(K[k] + 2):min(n, (K[k] + b)), j], "-")))
          -sum(sign(outer(X[max(1, (K[k] - b + 1)):(K[k]), j],
                          X[(K[k] + 1):(K[k] + 1), j], "-")))
          +sum(sign(outer(X[max(1, (K[k] - b + 1)):(K[k]), j],
                          X[min(n, (K[k] + b + 1)):min(n, (K[k] + b + 1)), j], "-")))
          +sum(sign(outer(X[(K[k] + 1):(K[k] + 1), j],
                          X[(K[k] + 2):min(n, (K[k] + b)), j], "-")))
          +sum(sign(outer(X[(K[k] + 1):(K[k] + 1), j],
                          X[min(n, (K[k] + b + 1)):min(n, (K[k] + b + 1)), j], "-")))
        )
      }
      mosum.cusum.inner
    }
  }
  inner.dimen
}

extract_candidates_one_bandwidth <- function(X, bandwidth, B = 100, paral = "F",
                                             kernel_type = c("linear", "sign"),
                                             eta = 0.15) {
  kernel_type <- match.arg(kernel_type)
  n <- nrow(X)
  d <- ncol(X)
  b <- bandwidth
  
  if (2 * b >= n) {
    return(list(
      bandwidth = bandwidth,
      mosum.sta = NA_real_,
      boot.quantile = NA_real_,
      reject = 0,
      candidates = data.frame()
    ))
  }
  
  K <- seq(b + 1, n - b - 1, 1)
  inner.dimen <- make_inner_dimen(X, b, kernel_type)
  
  mosum.cusum <- matrix(0, nrow = length(K), ncol = d)
  temp <- lapply(1:d, inner.dimen, K)
  for (i in 1:d) mosum.cusum[, i] <- temp[[i]]
  
  mosum.cusum <- abs(mosum.cusum)
  mosum.aggre <- apply(mosum.cusum, 1, max)
  mosum.sta <- max(mosum.aggre)
  
  boot.quantile <- if (kernel_type == "linear") {
    linear.mosum.boot(X, bandwith = b, B = B, paral = paral)
  } else {
    sign.mosum.boot(X, bandwith = b, B = B, paral = paral)
  }
  
  reject <- as.numeric(mosum.sta > boot.quantile)
  
  T_values <- c(rep(0, b), mosum.aggre, rep(0, b + 1))
  TT <- length(T_values)
  significant_windows <- list()
  
  for (v_j in b:(TT - b)) {
    start_w <- v_j + ceiling(eta * b)
    if (start_w > (TT - b)) next
    for (w_j in start_w:(TT - b)) {
      if (all(T_values[v_j:w_j] >= boot.quantile)) {
        if ((v_j > b && T_values[v_j - 1] < boot.quantile) &&
            (w_j < TT - b && T_values[w_j + 1] < boot.quantile)) {
          significant_windows <- append(significant_windows, list(c(v_j, w_j)))
          break
        }
      }
    }
  }
  
  if (length(significant_windows) == 0) {
    cand <- data.frame(
      cp = numeric(),
      v = integer(),
      w = integer(),
      stat = numeric(),
      stat_std = numeric(),
      bandwidth = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    cand <- do.call(rbind, lapply(significant_windows, function(win) {
      v_j <- win[1]
      w_j <- win[2]
      local_vals <- T_values[v_j:w_j]
      cp_j <- v_j + which.max(local_vals) - 1
      stat_j <- max(local_vals)
      data.frame(
        cp = cp_j,
        v = v_j,
        w = w_j,
        stat = stat_j,
        stat_std = stat_j / max(boot.quantile, 1e-8),
        bandwidth = bandwidth,
        stringsAsFactors = FALSE
      )
    }))
  }
  
  list(
    bandwidth = bandwidth,
    mosum.sta = mosum.sta,
    boot.quantile = boot.quantile,
    reject = reject,
    candidates = cand
  )
}

dedup_one_scale <- function(cand, L0) {
  if (is.null(cand) || nrow(cand) <= 1) return(cand)
  cand <- cand[order(cand$cp), , drop = FALSE]
  grp <- integer(nrow(cand))
  grp[1] <- 1
  g <- 1
  for (i in 2:nrow(cand)) {
    if ((cand$cp[i] - cand$cp[i - 1]) < L0) {
      grp[i] <- g
    } else {
      g <- g + 1
      grp[i] <- g
    }
  }
  
  keep_idx <- unlist(tapply(seq_len(nrow(cand)), grp, function(idx) {
    idx[which.max(cand$stat[idx])]
  }))
  out <- cand[sort(keep_idx), , drop = FALSE]
  rownames(out) <- NULL
  out
}

build_anchors <- function(cand_list, L0) {
  if (length(cand_list) == 0) return(data.frame())
  
  anchors <- cand_list[[1]]
  if (is.null(anchors) || nrow(anchors) == 0) {
    anchors <- data.frame(
      cp = numeric(), v = integer(), w = integer(),
      stat = numeric(), stat_std = numeric(),
      bandwidth = numeric(), scale_id = integer(),
      stringsAsFactors = FALSE
    )
  }
  if (nrow(anchors) > 0) anchors$anchor_id <- seq_len(nrow(anchors))
  next_id <- nrow(anchors)
  
  if (length(cand_list) >= 2) {
    for (m in 2:length(cand_list)) {
      cand <- cand_list[[m]]
      if (is.null(cand) || nrow(cand) == 0) next
      
      for (i in seq_len(nrow(cand))) {
        overlap_or_close <- FALSE
        if (nrow(anchors) > 0) {
          overlap_or_close <- any(
            mapply(interval_overlap, cand$v[i], cand$w[i], anchors$v, anchors$w) |
              (abs(cand$cp[i] - anchors$cp) < L0)
          )
        }
        if (!overlap_or_close) {
          next_id <- next_id + 1
          new_row <- cand[i, , drop = FALSE]
          new_row$anchor_id <- next_id
          anchors <- rbind(anchors, new_row)
        }
      }
    }
  }
  
  anchors <- anchors[order(anchors$cp), , drop = FALSE]
  rownames(anchors) <- NULL
  anchors
}

assign_clusters_nearest <- function(all_candidates, anchors) {
  if (is.null(all_candidates) || nrow(all_candidates) == 0) return(all_candidates)
  if (is.null(anchors) || nrow(anchors) == 0) {
    all_candidates$cluster_id <- integer(nrow(all_candidates))
    return(all_candidates)
  }
  
  cluster_id <- integer(nrow(all_candidates))
  
  for (i in seq_len(nrow(all_candidates))) {
    d <- abs(all_candidates$cp[i] - anchors$cp)
    min_idx <- which(d == min(d))
    
    if (length(min_idx) > 1) {
      ov <- mapply(interval_overlap,
                   all_candidates$v[i], all_candidates$w[i],
                   anchors$v[min_idx], anchors$w[min_idx])
      if (any(ov)) {
        min_idx <- min_idx[which(ov)[1]]
      } else {
        min_idx <- min_idx[1]
      }
    }
    
    cluster_id[i] <- anchors$anchor_id[min_idx[1]]
  }
  
  all_candidates$cluster_id <- cluster_id
  all_candidates
}

merged_from_clusters <- function(clustered_candidates) {
  if (is.null(clustered_candidates) || nrow(clustered_candidates) == 0) return(data.frame())
  
  clist <- split(clustered_candidates, clustered_candidates$cluster_id)
  merged <- do.call(rbind, lapply(clist, function(df) {
    df <- df[order(-df$stat_std, -df$stat, df$bandwidth), , drop = FALSE]
    out <- df[1, , drop = FALSE]
    out$n_support <- nrow(df)
    out$n_scales <- length(unique(df$bandwidth))
    out
  }))
  
  merged <- merged[order(merged$cp), , drop = FALSE]
  rownames(merged) <- NULL
  merged
}

merge_multiscale_candidates <- function(cand_list, L0) {
  cand_list_dedup <- lapply(cand_list, dedup_one_scale, L0 = L0)
  
  anchors <- build_anchors(cand_list_dedup, L0 = L0)
  
  pooled <- do.call(rbind, Filter(function(x) !is.null(x) && nrow(x) > 0, cand_list_dedup))
  if (is.null(pooled) || nrow(pooled) == 0) {
    pooled <- data.frame()
    merged <- data.frame()
  } else {
    rownames(pooled) <- NULL
    pooled <- pooled[order(pooled$bandwidth, pooled$cp), , drop = FALSE]
    pooled <- assign_clusters_nearest(pooled, anchors)
    merged <- merged_from_clusters(pooled)
  }
  
  list(
    candidates_dedup = cand_list_dedup,
    anchors = anchors,
    pooled = pooled,
    merged = merged
  )
}



multiscale_stage1 <- function(X, bandwidths, B = 100, paral = "F",
                              kernel_type = c("linear", "sign"),
                              eta = 0.15, L0 = 20) {
  kernel_type <- match.arg(kernel_type)
  bandwidths <- sort(unique(bandwidths))
  
  scale_results <- lapply(seq_along(bandwidths), function(m) {
    out <- extract_candidates_one_bandwidth(
      X = X,
      bandwidth = bandwidths[m],
      B = B,
      paral = paral,
      kernel_type = kernel_type,
      eta = eta
    )
    if (nrow(out$candidates) > 0) {
      out$candidates$scale_id <- m
    } else {
      out$candidates$scale_id <- integer(0)
    }
    out
  })
  
  cand_list <- lapply(scale_results, function(x) x$candidates)
  merge_out <- merge_multiscale_candidates(cand_list, L0 = L0)
  
  scale_summary <- data.frame(
    bandwidth = sapply(scale_results, function(x) x$bandwidth),
    mosum_sta = sapply(scale_results, function(x) x$mosum.sta),
    boot_quantile = sapply(scale_results, function(x) x$boot.quantile),
    reject = sapply(scale_results, function(x) x$reject),
    n_raw_candidates = sapply(scale_results, function(x) nrow(x$candidates)),
    n_dedup_candidates = sapply(merge_out$candidates_dedup, function(x) if (is.null(x)) 0 else nrow(x)),
    stringsAsFactors = FALSE
  )
  
  overall_stat <- suppressWarnings(max(scale_summary$mosum_sta / pmax(scale_summary$boot_quantile, 1e-8), na.rm = TRUE))
  if (!is.finite(overall_stat)) overall_stat <- NA_real_
  overall_reject <- as.numeric(any(scale_summary$reject == 1, na.rm = TRUE))
  
  list(
    scale_results = scale_results,
    scale_summary = scale_summary,
    merge_out = merge_out,
    overall_stat = overall_stat,
    overall_reject = overall_reject
  )
}






#====================
Data_gen_multicpt <- function(n, p, q0_list=c(0.3,0.5,0.7), 
                              cp_exit = TRUE,
                              light_tail = TRUE,degree=3,outlier=F){
  ####num_q0: number of chang points
  ####Bin, I will let you design how to generate the simulated data
  library(MASS)
  X=matrix(0,nrow=n,ncol=p)
  mu1<- rep(0,p)
  mu2<- rep(0,p)
  mu3<- rep(0,p)
  mu4<- rep(0,p)
  if(light_tail){
    delta<- runif(5,12*sqrt(log(p)/n),12*sqrt(log(p)/n))*c(1,-1,1,-1,1)
    set1<- 1:5
    mu2[set1]=mu1[set1]+delta
    mu3=mu1
    mu4=mu2
  }else{
    delta<- runif(5,15*sqrt(log(p)/n),15*sqrt(log(p)/n))*c(1,-1,1,-1,1)
    set1<- 1:5
    mu2[set1]=mu1[set1]+delta
    mu3=mu1
    mu4=mu2
  }
  Q<- (1:n)/n
  #  cov.matrix=diag(rep(1,p))
  z<- mvrnorm(n,rep(0,p),block.matrix(p))
  # z<- mvrnorm(n,rep(0,p),cov.matrix)
  w<- rchisq(n,degree)
  w<- matrix(rep(w,p),ncol=p)
  constant=0
  if(light_tail){
    epsilon=z
    constant=dnorm(0)
  }else{
    epsilon=z/sqrt(w/degree)
    constant=dt(0,degree)
  }
  q0 <- q0_list[1]
  q1 <- q0_list[2]
  q2<- q0_list[3]
  if(cp_exit){
    for(i in 1:n){
      X[i,]<- mu1*as.numeric(Q[i]<=q0) +
        mu2*as.numeric(Q[i]<=q1) *as.numeric(Q[i]>q0) +
        mu3*as.numeric(Q[i]<=q2)*as.numeric(Q[i]>q1) + mu4*as.numeric(Q[i]>q2)+ epsilon[i,]
    }
  }else{
    for(i in 1:n){
      X[i,]<- mu1*as.numeric(Q[i]<=q0) +
        mu1*as.numeric(Q[i]<=q1) *as.numeric(Q[i]>q0) +
        mu1*as.numeric(Q[i]<=q2)*as.numeric(Q[i]>q1) + mu1*as.numeric(Q[i]>q2)+ epsilon[i,]
    }
  }
  if(outlier){
    outlier.percent=0.1
    outlier.loc=sample(1:n,floor(n*outlier.percent))
    X[outlier.loc,]=X[outlier.loc,]+matrix(sample(c(-5, 5), size = floor(n*outlier.percent)*p, replace = TRUE),ncol=p)
  }
  return(list(X=X,mu=cbind(mu1,mu2,mu3,mu4),mu.new=cbind(mu2-mu1,mu3-mu2,mu4-mu3)*1))
}


Data_gen_multicpt_var <- function(n, p, q0_list=c(0.25,0.5,0.75), 
                                  cp_exit = TRUE,
                                  light_tail = TRUE,degree=3){
  ####num_q0: number of chang points
  ####Bin, I will let you design how to generate the simulated dat
  library(MASS)
  X=matrix(0,nrow=n,ncol=p)
  set=1:5
  mu1<- rep(sqrt(0.3),p)
  mu2=mu1
  mu2[set]<- c(rep(sqrt(2),5))
  mu3<- mu1
  mu4<- mu2
  Q<- (1:n)/n
  z<- mvrnorm(n,rep(0,p),diag(rep(1,p)))
  #w<- rchisq(n,degree)
  w<- matrix(rt(n*p,degree),ncol=p)/sqrt(5/3)
  if(light_tail){
    X=z
  }else{
    X=w
  }
  q0 <- q0_list[1]
  q1 <- q0_list[2]
  q2<- q0_list[3]
  if(cp_exit){
    for(i in 1:n){
      X[i,]<- mu1*as.numeric(Q[i]<=q0)*X[i,] +
        mu2*as.numeric(Q[i]<=q1) *as.numeric(Q[i]>q0)*X[i,] +
        mu3*as.numeric(Q[i]<=q2)*as.numeric(Q[i]>q1)*X[i,] + mu4*as.numeric(Q[i]>q2)*X[i,]
    }
  }else{
    for(i in 1:n){
      X[i,]<- mu1*as.numeric(Q[i]<=q0)*X[i,] +
        mu1*as.numeric(Q[i]<=q1) *as.numeric(Q[i]>q0)*X[i,] +
        mu1*as.numeric(Q[i]<=q2)*as.numeric(Q[i]>q1)*X[i,] + mu1*as.numeric(Q[i]>q2)*X[i,]
    }
  }
  
  return(list(X=X))
  
}



# mosum.boot is the bootstrap function for the sign kernel
# linear.mosum.boot is the bootstrap function for the linear kernel
# X is the orignal data
# B is the number of bootstrap replications 


sign.mosum.boot=function(X,bandwith,B,paral="F"){
  library(parallel)
  boot.inner=function(bb){
    n=nrow(X)
    d=ncol(X)
    h=sqrt(log(d)/n)
    b=bandwith#floor(h*n)
    sigma=rep(1,d)
    K=seq(b+1,n-b-1,1)
    e=rnorm(n)
    mosum.cusum.boot=matrix(0,nrow=length(K),ncol=d)
    inner.dimen=function(j){
      mosum.cusum.inner=numeric()
    #  mosum.cusum.inner2=numeric()
    #  mosum.cusum.inner3=numeric() 
    #  theta.head=mean(sign(outer(X[(K[1]-b):K[1],j],X[(K[1]+1):(K[1]+b),j],"-")))
      mosum.cusum.inner[1]=1/(2*b)^(3/2)*sum(sign(outer(X[(K[1]-b):K[1],j],X[(K[1]+1):(K[1]+b),j],"-"))*outer(e[(K[1]-b):K[1]],e[(K[1]+1):(K[1]+b)],"+"))
    #  mosum.cusum.inner2[1]=1/(2*b)^(3/2)*sum(outer(e[(K[1]-b):K[1]],e[(K[1]+1):(K[1]+b)],"+"))
    #  mosum.cusum.inner3[1]=mosum.cusum.inner[1]-theta.head*mosum.cusum.inner2[1]
      for(k in 1:(length(K)-1)){
      #  theta.head=mean(sign(outer(X[(K[k]+1-b):(K[k]+1),j],X[(K[k]+1+1):(K[k]+1+b),j],"-")))
        mosum.cusum.inner[k+1]=mosum.cusum.inner[k]+1/(2*b)^(3/2)*(-sum(sign(outer(X[(K[k]-b):(K[k]-b),j],X[(K[k]+1):(K[k]+1),j],"-"))*outer(e[(K[k]-b):(K[k]-b)],e[(K[k]+1):(K[k]+1)],"+"))
                                                         -sum(sign(outer(X[(K[k]-b):(K[k]-b),j],X[(K[k]+2):(K[k]+b),j],"-"))*outer(e[(K[k]-b):(K[k]-b)],e[(K[k]+2):(K[k]+b)],"+"))
                                                         -sum(sign(outer(X[(K[k]-b+1):(K[k]),j],X[(K[k]+1):(K[k]+1),j],"-"))*outer(e[(K[k]-b+1):(K[k])],e[(K[k]+1):(K[k]+1)],"+"))
                                                         +sum(sign(outer(X[(K[k]-b+1):(K[k]),j],X[(K[k]+b+1):(K[k]+b+1),j],"-"))*outer(e[(K[k]-b+1):(K[k])],e[(K[k]+b+1):(K[k]+b+1)],"+"))
                                                         +sum(sign(outer(X[(K[k]+1):(K[k]+1),j],X[(K[k]+2):(K[k]+b),j],"-"))*outer(e[(K[k]+1):(K[k]+1)],e[(K[k]+2):(K[k]+b)],"+"))
                                                         +sum(sign(outer(X[(K[k]+1):(K[k]+1),j],X[(K[k]+b+1):(K[k]+b+1),j],"-"))*outer(e[(K[k]+1):(K[k]+1)],e[(K[k]+b+1):(K[k]+b+1)],"+")))
        
        
   #     mosum.cusum.inner2[k+1]= mosum.cusum.inner2[k]+1/(2*b)^(3/2)*(-sum(outer(e[(K[k]-b):(K[k]-b)],e[(K[k]+1):(K[k]+1)],"+"))
   #     -sum(outer(e[(K[k]-b):(K[k]-b)],e[(K[k]+2):(K[k]+b)],"+"))
   #     -sum(outer(e[(K[k]-b+1):(K[k])],e[(K[k]+1):(K[k]+1)],"+"))
   #     +sum(outer(e[(K[k]-b+1):(K[k])],e[(K[k]+b+1):(K[k]+b+1)],"+"))
   #     +sum(outer(e[(K[k]+1):(K[k]+1)],e[(K[k]+2):(K[k]+b)],"+"))
    #    +sum(outer(e[(K[k]+1):(K[k]+1)],e[(K[k]+b+1):(K[k]+b+1)],"+")))
   #     mosum.cusum.inner3[k+1]=mosum.cusum.inner[k+1]-theta.head*mosum.cusum.inner2[k+1]
      }
      return(mosum.cusum.inner)
    }
    
    temp=lapply(1:d,inner.dimen)
    for(i in 1:d){mosum.cusum.boot[,i]=temp[[i]]}
    mosum.cusum.boot=abs(mosum.cusum.boot)
    mosum.cusum.boot=as.matrix(mosum.cusum.boot)
    mosum.aggre.boot=apply(mosum.cusum.boot,1,max)
    return(mosum.sta.boot=max(mosum.aggre.boot))
  }
  
  #Boot=sapply(1:B,boot.inner)
 # mc.cores = min(B, parallel::detectCores() - 1)
 # print(mc.cores)
 # Boot <- unlist(parallel::mclapply(1:B, boot.inner, mc.cores = mc.cores))
  if(paral==T){
    mc.cores = min(B, parallel::detectCores() - 1)
    Boot <- unlist(parallel::mclapply(1:B, boot.inner, mc.cores = mc.cores))
  }else{Boot=sapply(1:B,boot.inner)}
  return(boot.quantile=quantile(Boot,0.95)) 
  
}



linear.mosum.boot=function(X,bandwith,B,paral="F"){
  library(parallel)
  boot.inner=function(bb){
    n=nrow(X)
    d=ncol(X)
    h=sqrt(log(d)/n)
    b=bandwith#floor(h*n)
    sigma=rep(1,d)
    K=seq(b+1,n-b-1,1)
    e=rnorm(n)
    mosum.cusum.boot=matrix(0,nrow=length(K),ncol=d)
    inner.dimen=function(j){
      mosum.cusum.inner=numeric()
      #  mosum.cusum.inner2=numeric()
      #  mosum.cusum.inner3=numeric() 
      #  theta.head=mean(sign(outer(X[(K[1]-b):K[1],j],X[(K[1]+1):(K[1]+b),j],"-")))
      mosum.cusum.inner[1]=1/(2*b)^(3/2)*sum((outer(X[(K[1]-b):K[1],j],X[(K[1]+1):(K[1]+b),j],"-"))*outer(e[(K[1]-b):K[1]],e[(K[1]+1):(K[1]+b)],"+"))
      #  mosum.cusum.inner2[1]=1/(2*b)^(3/2)*sum(outer(e[(K[1]-b):K[1]],e[(K[1]+1):(K[1]+b)],"+"))
      #  mosum.cusum.inner3[1]=mosum.cusum.inner[1]-theta.head*mosum.cusum.inner2[1]
      for(k in 1:(length(K)-1)){
        #  theta.head=mean(sign(outer(X[(K[k]+1-b):(K[k]+1),j],X[(K[k]+1+1):(K[k]+1+b),j],"-")))
        mosum.cusum.inner[k+1]=mosum.cusum.inner[k]+1/(2*b)^(3/2)*(-sum((outer(X[(K[k]-b):(K[k]-b),j],X[(K[k]+1):(K[k]+1),j],"-"))*outer(e[(K[k]-b):(K[k]-b)],e[(K[k]+1):(K[k]+1)],"+"))
                                                                   -sum((outer(X[(K[k]-b):(K[k]-b),j],X[(K[k]+2):(K[k]+b),j],"-"))*outer(e[(K[k]-b):(K[k]-b)],e[(K[k]+2):(K[k]+b)],"+"))
                                                                   -sum((outer(X[(K[k]-b+1):(K[k]),j],X[(K[k]+1):(K[k]+1),j],"-"))*outer(e[(K[k]-b+1):(K[k])],e[(K[k]+1):(K[k]+1)],"+"))
                                                                   +sum((outer(X[(K[k]-b+1):(K[k]),j],X[(K[k]+b+1):(K[k]+b+1),j],"-"))*outer(e[(K[k]-b+1):(K[k])],e[(K[k]+b+1):(K[k]+b+1)],"+"))
                                                                   +sum((outer(X[(K[k]+1):(K[k]+1),j],X[(K[k]+2):(K[k]+b),j],"-"))*outer(e[(K[k]+1):(K[k]+1)],e[(K[k]+2):(K[k]+b)],"+"))
                                                                   +sum((outer(X[(K[k]+1):(K[k]+1),j],X[(K[k]+b+1):(K[k]+b+1),j],"-"))*outer(e[(K[k]+1):(K[k]+1)],e[(K[k]+b+1):(K[k]+b+1)],"+")))
        
        
        #     mosum.cusum.inner2[k+1]= mosum.cusum.inner2[k]+1/(2*b)^(3/2)*(-sum(outer(e[(K[k]-b):(K[k]-b)],e[(K[k]+1):(K[k]+1)],"+"))
        #     -sum(outer(e[(K[k]-b):(K[k]-b)],e[(K[k]+2):(K[k]+b)],"+"))
        #     -sum(outer(e[(K[k]-b+1):(K[k])],e[(K[k]+1):(K[k]+1)],"+"))
        #     +sum(outer(e[(K[k]-b+1):(K[k])],e[(K[k]+b+1):(K[k]+b+1)],"+"))
        #     +sum(outer(e[(K[k]+1):(K[k]+1)],e[(K[k]+2):(K[k]+b)],"+"))
        #    +sum(outer(e[(K[k]+1):(K[k]+1)],e[(K[k]+b+1):(K[k]+b+1)],"+")))
        #     mosum.cusum.inner3[k+1]=mosum.cusum.inner[k+1]-theta.head*mosum.cusum.inner2[k+1]
      }
      return(mosum.cusum.inner)
    }
    
    temp=lapply(1:d,inner.dimen)
    for(i in 1:d){mosum.cusum.boot[,i]=temp[[i]]}
    mosum.cusum.boot=abs(mosum.cusum.boot)
    mosum.cusum.boot=as.matrix(mosum.cusum.boot)
    mosum.aggre.boot=apply(mosum.cusum.boot,1,max)
    return(mosum.sta.boot=max(mosum.aggre.boot))
  }
if(paral==T){
   mc.cores = min(B, parallel::detectCores() - 1)
   Boot <- unlist(parallel::mclapply(1:B, boot.inner, mc.cores = mc.cores))
}else{Boot=sapply(1:B,boot.inner)}
 # t3=Sys.time()
  return(boot.quantile=quantile(Boot,0.95)) 
  
}






rmse=function(x,tau0){
  sqrt(mean((x-tau0)^2))
}

s0.p=function(x,s0,p){
  d=length(x)
  y=sort(abs(x),decreasing = F)
  z= sum(y[(d-s0+1):d]^p)^{1/p}
  return(z)
  
}

block.matrix=function(dd){
  d=dd
  k=floor(d/5)
  NN=seq(1,d,5)
  A=matrix(0,nrow=d,ncol=d)
  for(i in 1:k){
    for(j in NN[i]:(NN[i]+4)){
      for(l in j:(NN[i]+4)){
        A[j,l]=0.6
      }
    }
  }
  A=A+t(A)
  diag(A)=runif(d,1,1)
  A
}

band.matrix=function(dd){
  d=dd
  A=matrix(0,nrow=d,ncol=d)
  for(i in 1:d){
    for(j in 1:i){
      A[i,j]=0.6^abs(i-j)
    }
  }
  A=A+t(A)-diag(rep(1,d))
  A
}  

equi.matrix=function(d){
  a=matrix(0.2,nrow=d,ncol=d)
  diag(a)=rep(1,d)
  a
}


hausdorff.distance=function(v1,v2,n){
  p1=length(v1)
  p2=length(v2)
  if ( p1*p2==0){ dis=n} else{
    distance.mat=matrix(0,nrow=p1,ncol=p2)
    for( i in 1: p1){
      for( j in 1: p2){
        distance.mat[i,j]=abs(v1[i]-v2[j])
        
      }
    }
    dis=max(max( apply(distance.mat, 1, min)),max( apply(distance.mat, 2, min)))/n
  }
  return(dis)
  
}

accuracy=function(x,tau0){
  length(which(x==tau0))/length(x)
}

rand.index=function(cpt1=c(1,300,500,700,1000),cpt2,n=1000){
  
  partition1<- list()
  partition2<- list()
  for(i in 1:(length(cpt1)-1)  ){
    partition1[[i]]<- seq((cpt1[i]+1),cpt1[i+1],1)   
  }
  
  for(i in 1:(length(cpt2)-1)){
    partition2[[i]]<- seq( (cpt2[i]+1),cpt2[i+1],1)   
  }
  
  a<- 0
  b<- 0
  c<- 0
  d<- 0
  
  
  for(i in 2:n){
    for(j in 1:(i-1)){
      temp1=0
      temp2=0
      for(k1 in 1:length(partition1)){
        if(as.numeric(is.element(i,partition1[[k1]]))==1&  as.numeric(is.element(j,partition1[[k1]]))==1){temp1=1}
      }
      
      for(k2 in 1:length(partition2)){
        if(as.numeric(is.element(i,partition2[[k2]]))==1&  as.numeric(is.element(j,partition2[[k2]]))==1){temp2=1}
      }
      
      if(temp1==1&temp2==1){a=a+1}
      if(temp1==0&temp2==0){b=b+1}
      if(temp1==1&temp2==0){c=c+1}
      if(temp1==0&temp2==1){d=d+1}
    }
  }
  rand.index<- (a+b)/choose(n,2)
  return(rand.index)
}



covariance.est=function(x,y,theta){
  n.x=nrow(x)
  n.y=nrow(y)
  
  
  n.x.half=floor(n.x/2)
  x.before=x[1:floor(n.x/2),] # to estimate the kernel function 
  x.after=x[(1+floor(n.x/2)):n.x,]
  
  h1.underline=function(z){
    n.z=nrow(z)
    d.z=ncol(z)
    
    
    res=matrix(NA,nrow=n.z,ncol=d.z)
    for(i in 1:n.z){
      temp=matrix(rep(z[i,],n.x.half),nrow=n.x.half,byrow = T)
      res[i,]=apply(sign(x.before-temp),2,mean)
    }
    return(res)
  }
  
  
  
  
  n.y.half=floor(n.y/2)
  y.before=y[1:floor(n.y/2),] # to estimate the kernel function 
  y.after=y[(1+floor(n.y/2)):n.y,]
  h1.overline=function(z){
    n.z=nrow(z)
    d.z=ncol(z)
    
    
    res=matrix(NA,nrow=n.z,ncol=d.z)
    for(i in 1:n.z){
      temp=matrix(rep(z[i,],n.y.half),nrow=n.y.half,byrow = T)
      res[i,]=apply(sign(y.before-temp),2,mean)
    }
    return(res)
  }
  
  
  h1.prime=function(z){
    n.z=nrow(z)
    d.z=ncol(z)
    
    res=matrix(NA,nrow=n.z,ncol=d.z)
    for(i in 1:n.z){
      temp=matrix(rep(z[i,],n.y),nrow=n.y,byrow = T)
      res[i,]=apply((sign(y-temp)-matrix(rep(theta,n.y),nrow=n.y,byrow = T)),2,mean)
    }
    return(res)
  }
  
  h2.prime=function(z){
    n.z=nrow(z)
    d.z=ncol(z)
    
    res=matrix(NA,nrow=n.z,ncol=d.z)
    for(i in 1:n.z){
      temp=matrix(rep(z[i,],n.x),nrow=n.x,byrow = T)
      res[i,]=apply((sign(temp-x)-matrix(rep(theta,n.x),nrow=n.x,byrow = T)),2,mean)
    }
    return(res)
  }
  
  cov.h1.prime=cov(h1.prime(x))  
  cov.h2.prime=cov(h2.prime(y))
  cov.h1.overline.h2.prime=cov(h1.overline(y.after)-h2.prime(y.after))
  cov.h1.underline.h1.prime=cov(h1.underline(x.after)+h1.prime(x.after))
  res=list(cov.h1.prime=cov.h1.prime,cov.h2.prime=cov.h2.prime,cov.h1.underline.h1.prime=cov.h1.underline.h1.prime,cov.h1.overline.h2.prime=cov.h1.overline.h2.prime)
  return(res)
}

top_abs_indices <- function(x, n = 3) {
  return(head(order(abs(x), decreasing = TRUE), n))
}












