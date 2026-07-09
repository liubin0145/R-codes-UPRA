
library(MASS)
# x is the data before cpt
# y is the data after cpt
# theta is the estimated difference of the kernel
sign.covariance.est=function(x,y,theta){
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


linear.covariance.est=function(x,y,theta){
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
      res[i,]=apply((x.before-temp),2,mean)
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
      res[i,]=apply((y.before-temp),2,mean)
    }
    return(res)
  }
  
  
  h1.prime=function(z){
    n.z=nrow(z)
    d.z=ncol(z)
    
    res=matrix(NA,nrow=n.z,ncol=d.z)
    for(i in 1:n.z){
      temp=matrix(rep(z[i,],n.y),nrow=n.y,byrow = T)
      res[i,]=apply(((y-temp)-matrix(rep(theta,n.y),nrow=n.y,byrow = T)),2,mean)
    }
    return(res)
  }
  
  h2.prime=function(z){
    n.z=nrow(z)
    d.z=ncol(z)
    
    res=matrix(NA,nrow=n.z,ncol=d.z)
    for(i in 1:n.z){
      temp=matrix(rep(z[i,],n.x),nrow=n.x,byrow = T)
      res[i,]=apply(((temp-x)-matrix(rep(theta,n.x),nrow=n.x,byrow = T)),2,mean)
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