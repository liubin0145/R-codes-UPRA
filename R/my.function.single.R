
#====================
Data_gen_multicpt <- function(n, p, q0_list=c(0.3,0.5,0.7), 
                              cp_exit = TRUE,
                              light_tail = TRUE,degree=3,outlier=F){
  source("my.function.single.R")
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
      mosum.cusum.inner[1]=1/(2*b)^(3/2)*sum(sign(outer(X[(K[1]+1):(K[1]+b),j],X[(K[1]-b):K[1],j],"-"))*outer(e[(K[1]+1):(K[1]+b)],e[(K[1]-b):K[1]],"+"))
    #  mosum.cusum.inner2[1]=1/(2*b)^(3/2)*sum(outer(e[(K[1]-b):K[1]],e[(K[1]+1):(K[1]+b)],"+"))
    #  mosum.cusum.inner3[1]=mosum.cusum.inner[1]-theta.head*mosum.cusum.inner2[1]
      for(k in 1:(length(K)-1)){
      #  theta.head=mean(sign(outer(X[(K[k]+1-b):(K[k]+1),j],X[(K[k]+1+1):(K[k]+1+b),j],"-")))
        mosum.cusum.inner[k+1]=mosum.cusum.inner[k]+1/(2*b)^(3/2)*(-sum(sign(outer(X[(K[k]+1):(K[k]+1),j],X[(K[k]-b):(K[k]-b),j],"-"))*outer(e[(K[k]+1):(K[k]+1)],e[(K[k]-b):(K[k]-b)],"+"))
                                                         -sum(sign(outer(X[(K[k]+2):(K[k]+b),j],X[(K[k]-b):(K[k]-b),j],"-"))*outer(e[(K[k]+2):(K[k]+b)],e[(K[k]-b):(K[k]-b)],"+"))
                                                         -sum(sign(outer(X[(K[k]+1):(K[k]+1),j],X[(K[k]-b+1):(K[k]),j],"-"))*outer(e[(K[k]+1):(K[k]+1)],e[(K[k]-b+1):(K[k])],"+"))
                                                         +sum(sign(outer(X[(K[k]+b+1):(K[k]+b+1),j],X[(K[k]-b+1):(K[k]),j],"-"))*outer(e[(K[k]+b+1):(K[k]+b+1)],e[(K[k]-b+1):(K[k])],"+"))
                                                         +sum(sign(outer(X[(K[k]+2):(K[k]+b),j],X[(K[k]+1):(K[k]+1),j],"-"))*outer(e[(K[k]+2):(K[k]+b)],e[(K[k]+1):(K[k]+1)],"+"))
                                                         +sum(sign(outer(X[(K[k]+b+1):(K[k]+b+1),j],X[(K[k]+1):(K[k]+1),j],"-"))*outer(e[(K[k]+b+1):(K[k]+b+1)],e[(K[k]+1):(K[k]+1)],"+")))
        
        
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
      mosum.cusum.inner[1]=1/(2*b)^(3/2)*sum((outer(X[(K[1]+1):(K[1]+b),j],X[(K[1]-b):K[1],j],"-"))*outer(e[(K[1]+1):(K[1]+b)],e[(K[1]-b):K[1]],"+"))
      #  mosum.cusum.inner2[1]=1/(2*b)^(3/2)*sum(outer(e[(K[1]-b):K[1]],e[(K[1]+1):(K[1]+b)],"+"))
      #  mosum.cusum.inner3[1]=mosum.cusum.inner[1]-theta.head*mosum.cusum.inner2[1]
      for(k in 1:(length(K)-1)){
        #  theta.head=mean(sign(outer(X[(K[k]+1-b):(K[k]+1),j],X[(K[k]+1+1):(K[k]+1+b),j],"-")))
        mosum.cusum.inner[k+1]=mosum.cusum.inner[k]+1/(2*b)^(3/2)*(-sum((outer(X[(K[k]+1):(K[k]+1),j],X[(K[k]-b):(K[k]-b),j],"-"))*outer(e[(K[k]+1):(K[k]+1)],e[(K[k]-b):(K[k]-b)],"+"))
                                                                   -sum((outer(X[(K[k]+2):(K[k]+b),j],X[(K[k]-b):(K[k]-b),j],"-"))*outer(e[(K[k]+2):(K[k]+b)],e[(K[k]-b):(K[k]-b)],"+"))
                                                                   -sum((outer(X[(K[k]+1):(K[k]+1),j],X[(K[k]-b+1):(K[k]),j],"-"))*outer(e[(K[k]+1):(K[k]+1)],e[(K[k]-b+1):(K[k])],"+"))
                                                                   +sum((outer(X[(K[k]+b+1):(K[k]+b+1),j],X[(K[k]-b+1):(K[k]),j],"-"))*outer(e[(K[k]+b+1):(K[k]+b+1)],e[(K[k]-b+1):(K[k])],"+"))
                                                                   +sum((outer(X[(K[k]+2):(K[k]+b),j],X[(K[k]+1):(K[k]+1),j],"-"))*outer(e[(K[k]+2):(K[k]+b)],e[(K[k]+1):(K[k]+1)],"+"))
                                                                   +sum((outer(X[(K[k]+b+1):(K[k]+b+1),j],X[(K[k]+1):(K[k]+1),j],"-"))*outer(e[(K[k]+b+1):(K[k]+b+1)],e[(K[k]+1):(K[k]+1)],"+")))
        
        
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
  
  #如何判断i与j落在相同的partition还是不同的partition呢？
  
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



