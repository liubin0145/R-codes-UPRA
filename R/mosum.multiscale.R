
linear.mosum.multiscale=function(X,bandwidths=80,B=100,paral="F",L0=20,eta=0.15){
  source("my.function.multiscale.R")
  n=nrow(X)
  d=ncol(X)
  b=min(sort(unique(bandwidths)))

  stage1 <- multiscale_stage1(
    X = X,
    bandwidths = sort(unique(bandwidths)),
    B = B,
    paral = paral,
    kernel_type = "linear",
    eta = eta,
    L0 = L0
  )

  mosum.sta <- stage1$overall_stat
  reject <- stage1$overall_reject

  merged.initial <- stage1$merge_out$merged
  if (nrow(merged.initial) > 0) {
    cpt.est.initial.actual <- sort(unique(merged.initial$cp))
  } else {
    cpt.est.initial.actual <- numeric(0)
  }
  cpt.est.initial <- c(1, cpt.est.initial.actual, n)

  cpt.est.actual.global <- cpt.est.initial.actual
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
    cat("legnth left=",e1-s1,"\n")
    cat("s2,e2=",c(s2,e2),"\n")
    cat("legnth right=",e2-s2,"\n")
    for(j in 1:d){
      theta.had[j, i] <- 1 / ((e1 - s1 + 1) * (e2 - s2 + 1)) *
        sum(outer(X[s1:e1, j], X[s2:e2, j], "-"))
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
    b_refine_i <- max(2, min(local.win.initial[[i]]$search_left,
                             local.win.initial[[i]]$search_right))
    inner_local <- make_inner_dimen(X, b_refine_i, "linear")
    temp <- lapply(1:d, inner_local, K.new)
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

for(i in 1:length(cpt.est.refine.actual)){
  s1 <- local.win.refine[[i]]$s1
  e1 <- local.win.refine[[i]]$e1
  s2 <- local.win.refine[[i]]$s2
  e2 <- local.win.refine[[i]]$e2
  cat("s1,e1=",c(s1,e1),"\n")
  cat("legnth left=",e1-s1,"\n")
  cat("s2,e2=",c(s2,e2),"\n")
  cat("legnth right=",e2-s2,"\n")
  X.before <- X[s1:e1, Gamma.head[[i]], drop = FALSE]
  X.after  <- X[s2:e2, Gamma.head[[i]], drop = FALSE]
  
  cov.est <- linear.covariance.est(
    X.before, X.after,
    theta = theta.had[Gamma.head[[i]], i]
  )
  
  sim.argmax <- numeric()
  for(s in 1:2000){
    sim.argmax[s] <- argmax(
      theta  = theta.had[Gamma.head[[i]], i],
      Sigma1 = cov.est[[1]],
      Sigma2 = cov.est[[2]],
      Sigma3 = cov.est[[3]],
      Sigma4 = cov.est[[4]]
    )
  }
  
  argmax.quantile <- quantile(sim.argmax, c(0.025, 0.975))
  confidence[i, ] <- cpt.est.refine.actual[i] -
    1 * rev(argmax.quantile) * (norm(theta.had[Gamma.head[[i]], i], "2"))^(-2)
}




final.res=list(
mosum.sta=mosum.sta,
mosum.reject=reject,
cpt.est.initial=cpt.est.initial,
cpt.est.refine=cpt.est.refine,
cpt.conf=confidence,
anchors = stage1$merge_out$anchors,
pooled.candidates = stage1$merge_out$pooled,
merged.initial = merged.initial
)
return(final.res)
}



sign.mosum.multiscale=function(X,bandwidths=80,B=100,paral="F",L0=20,eta=0.15){
  source("my.function.multiscale.R")
  n=nrow(X)
  d=ncol(X)
  b=min(sort(unique(bandwidths)))
  
  stage1 <- multiscale_stage1(
    X = X,
    bandwidths = sort(unique(bandwidths)),
    B = B,
    paral = paral,
    kernel_type = "sign",
    eta = eta,
    L0 = L0
  )
  
  mosum.sta <- stage1$overall_stat
  reject <- stage1$overall_reject
  
  merged.initial <- stage1$merge_out$merged
  if (nrow(merged.initial) > 0) {
    cpt.est.initial.actual <- sort(unique(merged.initial$cp))
  } else {
    cpt.est.initial.actual <- numeric(0)
  }
  cpt.est.initial <- c(1, cpt.est.initial.actual, n)
  
  cpt.est.actual.global <- cpt.est.initial.actual
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
            sum(sign(outer(X[s1:e1, j], X[s2:e2, j], "-")))
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
        b_refine_i <- max(2, min(local.win.initial[[i]]$search_left,
                                 local.win.initial[[i]]$search_right))
        inner_local <- make_inner_dimen(X, b_refine_i, "sign")
        temp <- lapply(1:d, inner_local, K.new)
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
      
      sim.argmax <- numeric()
      for(s in 1:2000){
        sim.argmax[s] <- argmax(
          theta  = theta.had[Gamma.head[[i]], i],
          Sigma1 = cov.est[[1]],
          Sigma2 = cov.est[[2]],
          Sigma3 = cov.est[[3]],
          Sigma4 = cov.est[[4]]
        )
      }
      
      argmax.quantile <- quantile(sim.argmax, c(0.025, 0.975))
      confidence[i, ] <- cpt.est.refine.actual[i] -
        1 * rev(argmax.quantile) * (norm(theta.had[Gamma.head[[i]], i], "2"))^(-2)
    }
 
    final.res=list(
      mosum.sta=mosum.sta,
      mosum.reject=reject,
      cpt.est.initial=cpt.est.initial,
      cpt.est.refine=cpt.est.refine,
      cpt.conf=confidence,
      anchors = stage1$merge_out$anchors,
      pooled.candidates = stage1$merge_out$pooled,
      merged.initial = merged.initial
    )
  return(final.res)
}



mosum.multiscale=function(X,bandwidths=80,B=100,type="linear",paral="F",L0=20,eta=0.15){
  if(type=="linear"){
    res=linear.mosum.multiscale(X,bandwidths,B,paral=paral,L0=L0,eta=eta)
  }else if(type=="sign"){
    res=sign.mosum.multiscale(X,bandwidths,B,paral=paral,L0=L0,eta=eta)
  }
  return(res)
}




