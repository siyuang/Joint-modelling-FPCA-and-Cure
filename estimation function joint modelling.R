# Estimation functions for joint modelling of longitudinal data and time to event outcome using FPCA and mixture cure model #

# load library #
library(orthogonalsplinebasis)
library(mvtnorm)
library(MASS)
library(boot)
library(Matrix)
library(lme4)
library(statmod)
library(nlme)
library(pROC)
library(survival)

# b-spline functions

bspline.basis<-function(t.min,t.max,n.knots,n.order){
  each.time<-t.max
  breaks<-seq(t.min,t.max,length.out = n.knots)
  # degree = order-1
  n.basis<-length(breaks)+n.order-2
  full.knots<- expand.knots(breaks, order=n.order) 
  ob<-OrthogonalSplineBasis(full.knots,order=n.order)
  D.temp<-OuterProdSecondDerivative(ob)
  return(list('ob'=ob,D.temp=D.temp,'n.basis'=n.basis))
  #return(nu.standard)
}

b.basis<-function(ob,t){
  evaluate(ob,t)
}



# Functions of each component

# Function of X(t)

x.t.fun<-function(t,gamma,v.obs,theta.x,phi.xi,b.basis,ob){
  basis.temp<-b.basis(ob,t=t)
  x.t<-as.vector(t(gamma)%*%t(v.obs))+as.vector(t(theta.x)%*%t(basis.temp))+t(phi.xi%*%t( basis.temp))
  return(t(x.t))
}

x.t.fun1<-function(t,gamma,v.obs,theta.x,phi.xi,basis.Y.vec){
  mean.vec<-as.vector(t(gamma)%*%t(v.obs))
  mean.temp<-matrix(rep(as.vector(t(theta.x)%*%t(basis.Y.vec)),nrow(phi.xi)),ncol=length(Y.vec),byrow = TRUE)
  rf.temp<-phi.xi%*%t(basis.Y.vec)
  x.t<-t((mean.vec+mean.temp)+rf.temp)
  return(x.t)
}

# hazard function of T

h.t.fun<-function(t,h0.num.hat,Y.vec,beta.hat,gamma,v.obs,v.obs.base,theta.x,phi.xi,b.basis,ob.x){
  beta1<-beta.hat[1]
  if(length(beta.hat)>1){
    beta.base<-beta.hat[2:length(beta.hat)]
  } else {beta.base<-0}
  x.t.temp<-x.t.fun(t,gamma,v.obs,theta.x,phi.xi,b.basis,ob.x)
  pos<-vapply(t,function(u){which(u==Y.vec)},FUN.VALUE = c(1))
  h.t.temp<-h0.num.hat[pos]*t(exp(beta1*x.t.temp+as.numeric(beta.base%*%v.obs.base)))
  return(t(h.t.temp))
}

# Survival function of T

S.t.fun<-function(t,h0.num.hat,Y.vec,beta,gamma,v.obs,v.obs.base,theta.x,phi.xi,b.basis,ob.x){
  if(t==0){
    S.t.out<-rep(1,q.sample)
  } else {
    pos<-which(t>=Y.vec)
    h0.temp<-h.t.fun(t=Y.vec[pos],h0.num.hat,Y.vec,beta,gamma,v.obs,v.obs.base,theta.x,phi.xi,b.basis,ob.x)
    S.t.temp<-exp(-rowSums(h0.temp))
    S.t.out<-S.t.temp
  }
  return(S.t.out)
}

# density function of T

f.t.fun<-function(t,h0.num.hat,Y.vec,beta,gamma,v.obs,v.obs.base,theta.x,phi.xi,b.basis,ob.x){
  #Y.vec.temp<-Y.vec[max(t)>=Y.vec]
  if(t==0){
    f.t.out<-rep(0,q.sample)
  } else {
    h0.temp<-h.t.fun(t,h0.num.hat,Y.vec,beta,gamma,v.obs,v.obs.base,theta.x,phi.xi,b.basis,ob.x)
    S.t.temp<-S.t.fun(t,h0.num.hat,Y.vec,beta,gamma,v.obs,v.obs.base,theta.x,phi.xi,b.basis,ob.x)
    f.t.out<-S.t.temp*h0.temp
  }
  return(as.vector(f.t.out))
}


# Common components of expectation components #
# Monte-carlo integration was used #
# phi.xi.sample is generated from pre-specified distribution

f_obs_delta1<-function(phi.xi.sample,
                       Y.obs,v.obs,v.obs.base,
                       h0.num.hat,Y.vec,beta.hat,
                       gamma.hat,theta.x.hat,sigma.hat,
                       b.basis,ob.x){
  n<-length(Y.obs)
  
  f_t<-NULL
  
  for(i in 1:n){
    f_t.temp<-f.t.fun(Y.obs[i],h0.num.hat,Y.vec,beta.hat,gamma.hat,v.obs[i,],v.obs.base[i,],theta.x.hat,phi.xi.sample[[i]],
                      b.basis,ob.x)
    f_t<-cbind(f_t,f_t.temp)
    
  } # 1*n #
  

  out<-f_t
  return(out) # 1*n #
}


f_obs0_int<-function(phi.xi.sample,
                     Y.obs,v.obs,v.obs.base,
                     h0.num.hat,Y.vec,beta.hat,
                     gamma.hat,theta.x.hat,sigma.hat,
                     b.basis,ob.x){
  n<-length(Y.obs)
  
  # current survival function #
  S_t<-NULL
  
  for(i in 1:n){
    S_t.temp<-S.t.fun(Y.obs[i],h0.num.hat,Y.vec,beta.hat,gamma.hat,v.obs[i,],v.obs.base[i,],theta.x.hat,phi.xi.sample[[i]],
                      b.basis,ob.x)
    S_t<-cbind(S_t,S_t.temp)
  } # 1*n #
  
  out1<-S_t
  return(out1) # 1*n #
}

######################################################################
# --- The following functions are used for maximization procedure ---#
#                                                                    #
#                                                                    #
######################################################################

# parameters of the incident rate #
# par is the coefficient alpha    #
exp_inc<-function(par,v.inc.1,v.inc.0,out_obs1.1,out_deno_obs1,out_k1_obs0.0,out_k0_obs0.0,f_obs0_k1_fin){
  c1.1<-as.vector(par%*%t(v.inc.1))-log(1+as.vector(exp(par%*%t(v.inc.1))))
  c1.0<-f_obs0_k1_fin*as.vector(par%*%t(v.inc.0))-log(1+as.vector(exp(par%*%t(v.inc.0))))
  c1.1[is.na(c1.1)]<--10000
  c1.0[is.na(c1.0)]<--10000
  ll<--(sum(c1.1)+sum(c1.0))
  return(ll)
}

exp_inc_deriv<-function(par,v.inc.1,v.inc.0,f_obs0_k1_fin){
  temp.1<-exp(par%*%t(v.inc.1))/(1+exp(par%*%t(v.inc.1)))
  temp.0<-exp(par%*%t(v.inc.0))/(1+exp(par%*%t(v.inc.0)))
  c1.1<-as.vector(1-temp.1)*v.inc.1
  c1.0<-as.vector(f_obs0_k1_fin-temp.0)*v.inc.0
  ll<-colSums(c1.1)+colSums(c1.0)
  return(ll)
}

# This is for standard error of alpha
alpha.se.fun<-function(par,v.inc.1,v.inc.0,out_obs1.1,out_deno_obs1,out_k1_obs0.0,
                       out_k0_obs0.0,f_obs0_k1_fin){
  res<-matrix(0,ncol(v.inc.1),ncol(v.inc.1))
  for(i in 1:nrow(v.inc.1)){
    c1.1<-as.numeric(1-inv.logit(par%*%t(v.inc.1[i,])))
    c1.1.square<-c1.1^2*(v.inc.1[i,]%*%t(v.inc.1[i,]))
    res<-res+c1.1.square
  }
  
  for(i in 1:nrow(v.inc.0)){
    c1.0.temp1<-f_obs0_k1_fin[i]^2
    c1.0.temp2<--2*f_obs0_k1_fin[i]*inv.logit(par%*%t(v.inc.0[i,]))
    c1.0.temp3<-inv.logit(par%*%t(v.inc.0[i,]))^2
    c1.0.square<-as.numeric(c1.0.temp1+c1.0.temp2+c1.0.temp3)*(v.inc.0[i,]%*%t(v.inc.0[i,]))
    res<-res+c1.0.square
  }
  
  return(res)
  
}


# components in longitudinal continuous predictors #
# update measurement error variance sigma^2 #
# sigma_hat_ker is for each subject #
sigma_hat_ker<-function(phi.xi.sample.1,phi.xi.sample.0,
                        Y.obs.1,x.obs.1,x.obs.t.1,vx.obs.1,n1,
                        Y.obs.0,x.obs.0,x.obs.t.0,vx.obs.0,n0,
                        gamma.hat,theta.x.hat,
                        out_obs1.1,out_deno_obs1,
                        out_k1_obs0.0,out_k0_obs0.0,f_obs0_k1_fin,f_obs0_k0_fin,
                        out_deno_obs0_k1,out_deno_obs0_k0,
                        x.i.standard.1,x.i.standard.0){
  
  out<-NULL
  for(i in 1:n1){
    mean.temp<-as.vector(gamma.hat%*%t(vx.obs.1[i,]))+t(theta.x.hat)%*%t(x.i.standard.1[[i]]) # 1*Mi
    fpc.temp<-phi.xi.sample.1[[i]]%*%t(x.i.standard.1[[i]]) # q.sample*Mi
    out.temp<-as.vector(x.obs.1[[i]]-mean.temp)-t(fpc.temp) # Mi*q.sample
    out.temp1<-colSums(t(out.temp*out.temp)*(out_obs1.1[i,]/out_deno_obs1[i]))                     # q.sample
    out.temp2<-out.temp1      # q.sample
    out<-c(out,sum(out.temp2))
  }
  
  for(i in 1:n0){
    mean.temp<-as.vector(gamma.hat%*%t(vx.obs.0[i,]))+t(theta.x.hat)%*%t(x.i.standard.0[[i]]) # 1*Mi
    fpc.temp<-phi.xi.sample.0[[i]]%*%t(x.i.standard.0[[i]]) # q.sample*Mi
    out.temp<-as.vector(x.obs.0[[i]]-mean.temp)-t(fpc.temp) # Mi*q.sample
    out.temp1<-colSums(t(out.temp*out.temp)*
                         (out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] +
                            out_k0_obs0.0*f_obs0_k0_fin[i]) )                            # q.sample
    out.temp2<-out.temp1
    out<-c(out,sum(out.temp2))
  }
  res<-sum(out)/(length(unlist(x.obs.1))+length(unlist(x.obs.0)))
  return(res)
}



# survival function components # 
first.deriv.fun<-function(beta.hat,Y.vec,v.obs.base.1,v.obs.base.0,
                          phi.xi.sample.1,phi.xi.sample.0,
                          Y.obs.1,x.obs.1,x.obs.t.1,
                          vx.obs.1,ind.obs.1,n1,
                          Y.obs.0,x.obs.0,x.obs.t.0,
                          vx.obs.0,ind.obs.0,n0,
                          gamma.hat,theta.x.hat,sigma.hat,
                          out_obs1.1,out_deno_obs1,
                          out_k1_obs0.0,out_k0_obs0.0,
                          f_obs0_k1_fin,f_obs0_k0_fin,
                          out_deno_obs0_k1,out_deno_obs0_k0,
                          x.i.standard.1,x.i.standard.0,
                          basis.Y.vec.x){
  beta.x<-beta.hat[1]
  if(length(beta.hat)>1){
    beta.base<-beta.hat[2:length(beta.hat)]
  } else {beta.base<-0}
  
  Y.obs.mat<-matrix(rep(Y.obs,length(Y.vec)),ncol=length(Y.vec))
  Y.obs.mat.1<-Y.obs.mat[ind.obs==1,]
  Y.obs.mat.0<-Y.obs.mat[ind.obs==0,]
  Y.obs.ind<-t(apply(Y.obs.mat,1,function(vec){vec>=Y.vec}))
  Y.obs.ind.1<-Y.obs.ind[ind.obs==1,]
  Y.obs.ind.0<-Y.obs.ind[ind.obs==0,]
  risk.t<-colSums(Y.obs.ind)
  ind.obs.mat<-matrix(rep(ind.obs,length(Y.vec)),ncol=length(Y.vec))
  ind.obs.mat.1<-ind.obs.mat[ind.obs==1,]
  ind.obs.mat.0<-ind.obs.mat[ind.obs==0,]
  event.t.mat<-as.matrix(I(t(Y.obs.mat)==Y.vec))*t(ind.obs.mat)
  event.t.mat.1<-t(event.t.mat)[ind.obs==1,]
  event.t.mat.0<-t(event.t.mat)[ind.obs==0,]
  event.t.mat.mod<-rbind(event.t.mat.1,event.t.mat.0)
  event.t<-rowSums(event.t.mat)
  event.t.1<-rowSums(t(event.t.mat.1))
  event.t.0<-rowSums(t(event.t.mat.0))
  E.x.s1<-rep(0,length(Y.vec))
  E.x.s1.1.1<-rep(0,length(Y.vec))
  E.x.s1.1.base<-rep(0,length(Y.vec))
  E.x.s1.2.1<-rep(0,length(Y.vec))
  E.x.s1.2.base<-rep(0,length(Y.vec))
  E.x.s1.3.1<-rep(0,length(Y.vec))
  E.x.s1.3.base<-rep(0,length(Y.vec))
  
  
  
  for(i in 1:n1){
    #delta_obs==1
    
    x.t.temp<-x.t.fun1(Y.vec,gamma.hat,vx.obs.1[i,],theta.x.hat,phi.xi.sample.1[[i]],basis.Y.vec.x) # 1*nt
    v.obs.base.mat.temp<-matrix(rep(v.obs.base.1[i,],length(Y.vec)),nrow=length(Y.vec),ncol=length(v.obs.base.1[i,]),byrow = TRUE)
    # denominators
    eta.temp<-beta.x*x.t.temp+as.numeric(t(beta.base)%*%(v.obs.base.1[i,]))
    temp<-t(exp(eta.temp)*Y.obs.ind.1[i,])*(out_obs1.1[i,]/out_deno_obs1[i]) # n1*Y.vec
    E.x.s1<-E.x.s1+colSums(temp)
    
    
    # for first derivative # 
    # numerator
    temp.1.1<-(t(temp)*x.t.temp) # n1*n.vec
    #temp.1.base<-(t(temp)*v.obs.base.1[i,]) 
    temp.1.base<-vapply(v.obs.base.1[i,],function(u){rowSums(t(temp)*u)},FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s1.1.1<-E.x.s1.1.1+rowSums(temp.1.1)
    E.x.s1.1.base<-E.x.s1.1.base+temp.1.base
    
    # first component
    temp.2.1<-t(event.t.mat.1[i,]*x.t.temp)*out_obs1.1[i,]/out_deno_obs1[i] # n1*vec
    #temp.2.base<-t(event.t.mat.1[i,]*v.obs.base.mat.temp)*out_obs1.1[i,]/out_deno_obs1[i] # n1*vec
    temp.2.base<-vapply(v.obs.base.1[i,],
                        function(u){colSums(matrix(rep(t(event.t.mat.1[i,]*u),length(out_obs1.1[i,])),ncol=length(Y.vec),nrow=length(out_obs1.1[i,]),byrow=TRUE)*out_obs1.1[i,]/out_deno_obs1[i])},
                        FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s1.2.1<-E.x.s1.2.1+colSums(temp.2.1)
    E.x.s1.2.base<-E.x.s1.2.base+ temp.2.base
    
    # for second derivative
    # first numerator
    temp.3.1<-(t(temp)*x.t.temp*x.t.temp) # n1*n.vec
    #temp.3.base<-(t(temp)* v.obs.base.mat.temp* v.obs.base.mat.temp) # n1*n.vec 
    temp.3.base<-vapply(v.obs.base.1[i,],
                        function(u){rowSums(t(temp)*u*u)},
                        FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s1.3.1<-E.x.s1.3.1+rowSums(temp.3.1)
    E.x.s1.3.base<-E.x.s1.3.base+temp.3.base
    
  }
  #delta_obs=0
  
  E.x.s0<-rep(0,length(Y.vec))
  E.x.s0.1.1<-rep(0,length(Y.vec))
  E.x.s0.1.base<-rep(0,length(Y.vec))
  E.x.s0.2.1<-rep(0,length(Y.vec))
  E.x.s0.2.base<-rep(0,length(Y.vec))
  E.x.s0.3.1<-rep(0,length(Y.vec))
  E.x.s0.3.base<-rep(0,length(Y.vec))
  
  for(i in 1:n0){ 
    x.t.temp0<-x.t.fun1(Y.vec,gamma.hat,vx.obs.0[i,],theta.x.hat,phi.xi.sample.0[[i]],basis.Y.vec.x) # 1*nt
    v.obs.base.mat.temp0<-matrix(rep(v.obs.base.0[i,],length(Y.vec)),nrow=length(Y.vec),ncol=length(v.obs.base.0[i,]),byrow = TRUE)
    
    
    
    # denominators
    eta.temp0<-beta.x*x.t.temp0+as.numeric(t(beta.base)%*%(v.obs.base.0[i,]))
    temp0<-t(exp(eta.temp0)*Y.obs.ind.0[i,])*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] # n1*Y.vec
    E.x.s0<-E.x.s0+colSums(temp0)
    
    # first derivative
    # numerators
    temp0.1.1<-(t(temp0)*x.t.temp0) # n1*n.vec
    temp0.1.base<-vapply(v.obs.base.0[i,],function(u){rowSums(t(temp0)*u)},FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s0.1.1<-E.x.s0.1.1+rowSums(temp0.1.1)
    E.x.s0.1.base<-E.x.s0.1.base+temp0.1.base
    
    # first component
    temp0.2.1<-t(event.t.mat.0[i,]*x.t.temp0*Y.obs.ind.0[i,])*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] # n1*vec
    #temp0.2.base<-t(event.t.mat.0[i,]*s.t.temp0*Y.obs.ind.0[i,])*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] # n1*vec
    temp0.2.base<-vapply(v.obs.base.0[i,],
                         function(u){colSums(matrix(rep(t(event.t.mat.0[i,]*u),length(out_k1_obs0.0[i,])),ncol=length(Y.vec),nrow=length(out_k1_obs0.0[i,]),byrow=TRUE)*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i])},
                         FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s0.2.1<-E.x.s0.2.1+colSums(temp0.2.1)
    E.x.s0.2.base<-E.x.s0.2.base+temp0.2.base
    
    # for second derivative
    temp0.3.1<-(t(temp0)*x.t.temp0*x.t.temp0) # n1*n.vec
    
    #temp0.3.base<-(t(temp0)*s.t.temp0*x.t.temp0) # n1*n.vec 
    temp0.3.base<-vapply(v.obs.base.0[i,],
                         function(u){rowSums(t(temp0)*u*u)},
                         FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s0.3.1<-E.x.s0.3.1+rowSums(temp0.3.1)
    E.x.s0.3.base<-E.x.s0.3.base+temp0.3.base
    
  }
  
  deno.temp<-E.x.s1+E.x.s0
  h0.num.new<-event.t/deno.temp
  #h0.new<-stepfun(Y.vec,c(0,h0.num.new))
  
  # beta of survival function #
  # first derivative
  first.dev.1.x<-sum(E.x.s1.2.1,E.x.s0.2.1)
  first.dev.1.base<-colSums(E.x.s1.2.base+E.x.s0.2.base)
  first.dev.1<-c(first.dev.1.x,first.dev.1.base)
  
  first.dev.num.x<-E.x.s1.1.1+E.x.s0.1.1
  first.dev.num.base<-E.x.s1.1.base+E.x.s0.1.base
  first.dev.num<-event.t*cbind(first.dev.num.x,first.dev.num.base)
  
  first.dev.dno<-deno.temp
  
  first.dev<-first.dev.1-colSums(first.dev.num/first.dev.dno)
  
  
  
  return(first.dev)
}

second.deriv.fun<-function(beta.hat,Y.vec,v.obs.base.1,v.obs.base.0,
                           phi.xi.sample.1,phi.xi.sample.0,
                           Y.obs.1,x.obs.1,x.obs.t.1,
                           vx.obs.1,ind.obs.1,n1,
                           Y.obs.0,x.obs.0,x.obs.t.0,
                           vx.obs.0,ind.obs.0,n0,
                           gamma.hat,theta.x.hat,sigma.hat,
                           out_obs1.1,out_deno_obs1,
                           out_k1_obs0.0,out_k0_obs0.0,
                           f_obs0_k1_fin,f_obs0_k0_fin,
                           out_deno_obs0_k1,out_deno_obs0_k0,
                           x.i.standard.1,x.i.standard.0,
                           basis.Y.vec.x){
  beta.x<-beta.hat[1]
  if(length(beta.hat)>1){
    beta.base<-beta.hat[2:length(beta.hat)]
  } else {beta.base<-0}
  
  
  E.x.s1.3.1.base<-rep(0,length(Y.vec))
  Y.obs.mat<-matrix(rep(Y.obs,length(Y.vec)),ncol=length(Y.vec))
  Y.obs.mat.1<-Y.obs.mat[ind.obs==1,]
  Y.obs.mat.0<-Y.obs.mat[ind.obs==0,]
  Y.obs.ind<-t(apply(Y.obs.mat,1,function(vec){vec>=Y.vec}))
  Y.obs.ind.1<-Y.obs.ind[ind.obs==1,]
  Y.obs.ind.0<-Y.obs.ind[ind.obs==0,]
  risk.t<-colSums(Y.obs.ind)
  ind.obs.mat<-matrix(rep(ind.obs,length(Y.vec)),ncol=length(Y.vec))
  ind.obs.mat.1<-ind.obs.mat[ind.obs==1,]
  ind.obs.mat.0<-ind.obs.mat[ind.obs==0,]
  event.t.mat<-as.matrix(I(t(Y.obs.mat)==Y.vec))*t(ind.obs.mat)
  event.t.mat.1<-t(event.t.mat)[ind.obs==1,]
  event.t.mat.0<-t(event.t.mat)[ind.obs==0,]
  event.t.mat.mod<-rbind(event.t.mat.1,event.t.mat.0)
  event.t<-rowSums(event.t.mat)
  event.t.1<-rowSums(t(event.t.mat.1))
  event.t.0<-rowSums(t(event.t.mat.0))
  E.x.s1<-rep(0,length(Y.vec))
  E.x.s1.1.1<-rep(0,length(Y.vec))
  E.x.s1.1.base<-rep(0,length(Y.vec))
  E.x.s1.2.1<-rep(0,length(Y.vec))
  E.x.s1.2.base<-rep(0,length(Y.vec))
  E.x.s1.3.1<-rep(0,length(Y.vec))
  E.x.s1.3.base<-rep(0,length(Y.vec))
  
  
  for(i in 1:n1){
    #delta_obs==1
    
    x.t.temp<-x.t.fun1(Y.vec,gamma.hat,vx.obs.1[i,],theta.x.hat,phi.xi.sample.1[[i]],basis.Y.vec.x) # 1*nt
    v.obs.base.mat.temp<-matrix(rep(v.obs.base.1[i,],length(Y.vec)),nrow=length(Y.vec),ncol=length(v.obs.base.1[i,]),byrow = TRUE)
    # denominators
    eta.temp<-beta.x*x.t.temp+as.numeric(t(beta.base)%*%(v.obs.base.1[i,]))
    temp<-t(exp(eta.temp)*Y.obs.ind.1[i,])*(out_obs1.1[i,]/out_deno_obs1[i]) # n1*Y.vec
    E.x.s1<-E.x.s1+colSums(temp)
    
    
    # for first derivative # 
    # numerator
    temp.1.1<-(t(temp)*x.t.temp) # n1*n.vec
    #temp.1.base<-(t(temp)*v.obs.base.1[i,]) 
    temp.1.base<-vapply(v.obs.base.1[i,],function(u){rowSums(t(temp)*u)},FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s1.1.1<-E.x.s1.1.1+rowSums(temp.1.1)
    E.x.s1.1.base<-E.x.s1.1.base+temp.1.base
    
    # first component
    temp.2.1<-t(event.t.mat.1[i,]*x.t.temp)*out_obs1.1[i,]/out_deno_obs1[i] # n1*vec
    #temp.2.base<-t(event.t.mat.1[i,]*v.obs.base.mat.temp)*out_obs1.1[i,]/out_deno_obs1[i] # n1*vec
    temp.2.base<-vapply(v.obs.base.1[i,],
                        function(u){colSums(matrix(rep(t(event.t.mat.1[i,]*u),length(out_obs1.1[i,])),ncol=length(Y.vec),nrow=length(out_obs1.1[i,]),byrow=TRUE)*out_obs1.1[i,]/out_deno_obs1[i])},
                        FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s1.2.1<-E.x.s1.2.1+colSums(temp.2.1)
    E.x.s1.2.base<-E.x.s1.2.base+ temp.2.base
    
    # for second derivative
    # first numerator
    temp.3.1<-rowSums(t(temp)*x.t.temp*x.t.temp) # n1*n.vec
    #temp.3.base<-(t(temp)* v.obs.base.mat.temp* v.obs.base.mat.temp) # n1*n.vec 
    temp.3.base<-vapply(rowSums(t(temp)),function(u){u*(v.obs.base.1[i,])%*%t(v.obs.base.1[i,])},
                        FUN.VALUE = matrix(0,length(v.obs.base.1[i,]),length(v.obs.base.1[i,])))
    
    temp.3.1.base<-vapply(rowSums(t(temp)*x.t.temp),function(u){u*v.obs.base.1[i,]},
                          FUN.VALUE = rep(0,length(v.obs.base.1[i,])))
    
    E.x.s1.3.1<-E.x.s1.3.1+temp.3.1
    E.x.s1.3.base<-E.x.s1.3.base+temp.3.base
    E.x.s1.3.1.base<-E.x.s1.3.1.base+temp.3.1.base
    
  }
  #delta_obs=0
  
  E.x.s0<-rep(0,length(Y.vec))
  E.x.s0.1.1<-rep(0,length(Y.vec))
  E.x.s0.1.base<-rep(0,length(Y.vec))
  E.x.s0.2.1<-rep(0,length(Y.vec))
  E.x.s0.2.base<-rep(0,length(Y.vec))
  E.x.s0.3.1<-rep(0,length(Y.vec))
  E.x.s0.3.base<-rep(0,length(Y.vec))
  E.x.s0.3.1.base<-rep(0,length(Y.vec))
  
  for(i in 1:n0){ 
    x.t.temp0<-x.t.fun1(Y.vec,gamma.hat,vx.obs.0[i,],theta.x.hat,phi.xi.sample.0[[i]],basis.Y.vec.x) # 1*nt
    v.obs.base.mat.temp0<-matrix(rep(v.obs.base.0[i,],length(Y.vec)),nrow=length(Y.vec),ncol=length(v.obs.base.0[i,]),byrow = TRUE)
    
    
    
    # denominators
    eta.temp0<-beta.x*x.t.temp0+as.numeric(t(beta.base)%*%(v.obs.base.0[i,]))
    temp0<-t(exp(eta.temp0)*Y.obs.ind.0[i,])*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] # n1*Y.vec
    E.x.s0<-E.x.s0+colSums(temp0)
    
    # first derivative
    # numerators
    temp0.1.1<-(t(temp0)*x.t.temp0) # n1*n.vec
    temp0.1.base<-vapply(v.obs.base.0[i,],function(u){rowSums(t(temp0)*u)},FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s0.1.1<-E.x.s0.1.1+rowSums(temp0.1.1)
    E.x.s0.1.base<-E.x.s0.1.base+temp0.1.base
    
    # first component
    temp0.2.1<-t(event.t.mat.0[i,]*x.t.temp0*Y.obs.ind.0[i,])*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] # n1*vec
    #temp0.2.base<-t(event.t.mat.0[i,]*s.t.temp0*Y.obs.ind.0[i,])*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] # n1*vec
    temp0.2.base<-vapply(v.obs.base.0[i,],
                         function(u){colSums(matrix(rep(t(event.t.mat.0[i,]*u),length(out_k1_obs0.0[i,])),ncol=length(Y.vec),nrow=length(out_k1_obs0.0[i,]),byrow=TRUE)*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i])},
                         FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s0.2.1<-E.x.s0.2.1+colSums(temp0.2.1)
    E.x.s0.2.base<-E.x.s0.2.base+temp0.2.base
    
    # for second derivative
    
    # first numerator
    temp0.3.1<-rowSums(t(temp)*x.t.temp0*x.t.temp0) # n1*n.vec
    #temp.3.base<-(t(temp)* v.obs.base.mat.temp* v.obs.base.mat.temp) # n1*n.vec 
    temp0.3.base<-vapply(rowSums(t(temp0)),function(u){u*(v.obs.base.0[i,])%*%t(v.obs.base.0[i,])},
                         FUN.VALUE = matrix(0,length(v.obs.base.0[i,]),length(v.obs.base.0[i,])))
    
    temp0.3.1.base<-vapply(rowSums(t(temp0)*x.t.temp0),function(u){u*v.obs.base.0[i,]},
                           FUN.VALUE = rep(0,length(v.obs.base.0[i,])))
    
    E.x.s0.3.1<-E.x.s0.3.1+temp0.3.1
    E.x.s0.3.base<-E.x.s0.3.base+temp0.3.base
    E.x.s0.3.1.base<-E.x.s0.3.1.base+temp0.3.1.base
    
  }
  
  deno.temp<-E.x.s1+E.x.s0
  h0.num.new<-event.t/deno.temp
  #h0.new<-stepfun(Y.vec,c(0,h0.num.new))
  
  # beta of survival function #
  # first derivative
  first.dev.1.x<-sum(E.x.s1.2.1,E.x.s0.2.1)
  first.dev.1.base<-colSums(E.x.s1.2.base+E.x.s0.2.base)
  first.dev.1<-c(first.dev.1.x,first.dev.1.base)
  
  first.dev.num.x<-E.x.s1.1.1+E.x.s0.1.1
  first.dev.num.base<-E.x.s1.1.base+E.x.s0.1.base
  first.dev.num<-event.t*cbind(first.dev.num.x,first.dev.num.base)
  
  first.dev.dno<-deno.temp
  
  first.dev<-first.dev.1-colSums(first.dev.num/first.dev.dno)
  
  #second derivative
  # first part
  second.dev.num.1.x<-E.x.s1.3.1+E.x.s0.3.1   # 1*n.vec
  second.dev.num.1.base<-E.x.s1.3.base+E.x.s0.3.base   # n.base*n.base*n.vec
  second.dev.num.1.x.base<-E.x.s1.3.1.base+E.x.s0.3.1.base   # n.base*n.vec
  
  second.dev.dno.1<-first.dev.dno
  second.dev.1.x<-event.t*second.dev.num.1.x/second.dev.dno.1
  second.dev.1.base<-array(dim=dim(second.dev.num.1.base))
  for(i in 1:length(Y.vec)){
    second.dev.1.base[,,i]<-event.t[i]*second.dev.num.1.base[,,i]/second.dev.dno.1[i]
  }
  second.dev.1.x.base<-event.t*t(second.dev.num.1.x.base)/second.dev.dno.1
  
  second.first<-matrix(0,dim(second.dev.num.1.base)[1]+1,dim(second.dev.num.1.base)[1]+1)
  second.first[1,1]<--sum(second.dev.1.x)
  second.first[2:nrow(second.first),2:nrow(second.first)]<--rowSums(second.dev.1.base,dims=2)
  second.first[1,2:nrow(second.first)]<-second.first[2:nrow(second.first),1]<--colSums(second.dev.1.x.base)
  
  # second part
  second.dev.num.2.x<- first.dev.num.x* first.dev.num.x
  second.dev.num.2.base<-array(dim=dim(second.dev.num.1.base))
  for(i in 1:length(Y.vec)){
    second.dev.num.2.base[,,i]<-first.dev.num.base[i,]%*%t(first.dev.num.base[i,])
  }
  second.dev.num.2.x.base<-first.dev.num.x*first.dev.num.base
  
  second.dev.dno.2<-first.dev.dno^2
  second.dev.2.x<-event.t*second.dev.num.2.x/second.dev.dno.2
  second.dev.2.base<-array(dim=dim(second.dev.num.2.base))
  for(i in 1:length(Y.vec)){
    second.dev.2.base[,,i]<-event.t[i]*second.dev.num.2.base[,,i]/second.dev.dno.2[i]
  }
  second.dev.2.x.base<-event.t*second.dev.num.2.x.base/second.dev.dno.2
  
  second.second<-matrix(0,dim(second.dev.num.2.base)[1]+1,dim(second.dev.num.2.base)[1]+1)
  second.second[1,1]<-sum(second.dev.2.x)
  second.second[2:nrow(second.second),2:nrow(second.second)]<-rowSums(second.dev.2.base,dims=2)
  second.second[1,2:nrow(second.second)]<-second.second[2:nrow(second.second),1]<-colSums(second.dev.2.x.base)
  
  second.dev<-second.first+second.second
  
  return(second.dev)
}

# baseline hazard function 
h0.new.fun<-function(par,Y.vec,v.obs.base.1,v.obs.base.0,
                     phi.xi.sample.1,phi.xi.sample.0,
                     Y.obs.1,x.obs.1,x.obs.t.1,vx.obs.1,ind.obs.1,n1,
                     Y.obs.0,x.obs.0,x.obs.t.0,vx.obs.0,ind.obs.0,n0,
                     gamma.hat,theta.x.hat,sigma.hat,
                     out_obs1.1,out_deno_obs1,
                     out_k1_obs0.0,out_k0_obs0.0,f_obs0_k1_fin,f_obs0_k0_fin,
                     out_deno_obs0_k1,out_deno_obs0_k0,
                     x.i.standard.1,x.i.standard.0,
                     basis.Y.vec.x){
  beta.x<-par[1]
  if(length(par)>1){
    beta.base<-par[2:length(par)]
  } else {beta.base<-0}
  
  
  
  Y.obs.mat<-matrix(rep(Y.obs,length(Y.vec)),ncol=length(Y.vec))
  Y.obs.mat.1<-Y.obs.mat[ind.obs==1,]
  Y.obs.mat.0<-Y.obs.mat[ind.obs==0,]
  Y.obs.ind<-t(apply(Y.obs.mat,1,function(vec){vec>=Y.vec}))
  Y.obs.ind.1<-Y.obs.ind[ind.obs==1,]
  Y.obs.ind.0<-Y.obs.ind[ind.obs==0,]
  risk.t<-colSums(Y.obs.ind)
  ind.obs.mat<-matrix(rep(ind.obs,length(Y.vec)),ncol=length(Y.vec))
  ind.obs.mat.1<-ind.obs.mat[ind.obs==1,]
  ind.obs.mat.0<-ind.obs.mat[ind.obs==0,]
  event.t.mat<-as.matrix(I(t(Y.obs.mat)==Y.vec))*t(ind.obs.mat)
  event.t.mat.1<-t(event.t.mat)[ind.obs==1,]
  event.t.mat.0<-t(event.t.mat)[ind.obs==0,]
  event.t.mat.mod<-rbind(event.t.mat.1,event.t.mat.0)
  event.t<-rowSums(event.t.mat)
  event.t.1<-rowSums(t(event.t.mat.1))
  event.t.0<-rowSums(t(event.t.mat.0))
  E.x.s1<-rep(0,length(Y.vec))
  E.x.s1.1.1<-rep(0,length(Y.vec))
  E.x.s1.1.base<-rep(0,length(Y.vec))
  E.x.s1.2.1<-rep(0,length(Y.vec))
  E.x.s1.2.base<-rep(0,length(Y.vec))
  E.x.s1.3.1<-rep(0,length(Y.vec))
  E.x.s1.3.base<-rep(0,length(Y.vec))
  
  
  
  for(i in 1:n1){
    #delta_obs==1
    
    x.t.temp<-x.t.fun1(Y.vec,gamma.hat,vx.obs.1[i,],theta.x.hat,phi.xi.sample.1[[i]],basis.Y.vec.x) # 1*nt
    v.obs.base.mat.temp<-matrix(rep(v.obs.base.1[i,],length(Y.vec)),nrow=length(Y.vec),ncol=length(v.obs.base.1[i,]),byrow = TRUE)
    # denominators
    eta.temp<-beta.x*x.t.temp+as.numeric(t(beta.base)%*%(v.obs.base.1[i,]))
    temp<-t(exp(eta.temp)*Y.obs.ind.1[i,])*(out_obs1.1[i,]/out_deno_obs1[i]) # n1*Y.vec
    E.x.s1<-E.x.s1+colSums(temp)
    
    
    # for first derivative # 
    # numerator
    temp.1.1<-(t(temp)*x.t.temp) # n1*n.vec
    #temp.1.base<-(t(temp)*v.obs.base.1[i,]) 
    temp.1.base<-vapply(v.obs.base.1[i,],function(u){rowSums(t(temp)*u)},FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s1.1.1<-E.x.s1.1.1+rowSums(temp.1.1)
    E.x.s1.1.base<-E.x.s1.1.base+temp.1.base
    
    # first component
    temp.2.1<-t(event.t.mat.1[i,]*x.t.temp)*out_obs1.1[i,]/out_deno_obs1[i] # n1*vec
    #temp.2.base<-t(event.t.mat.1[i,]*v.obs.base.mat.temp)*out_obs1.1[i,]/out_deno_obs1[i] # n1*vec
    temp.2.base<-vapply(v.obs.base.1[i,],
                        function(u){colSums(matrix(rep(t(event.t.mat.1[i,]*u),length(out_obs1.1[i,])),ncol=length(Y.vec),nrow=length(out_obs1.1[i,]),byrow=TRUE)*out_obs1.1[i,]/out_deno_obs1[i])},
                        FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s1.2.1<-E.x.s1.2.1+colSums(temp.2.1)
    E.x.s1.2.base<-E.x.s1.2.base+ temp.2.base
    
    # for second derivative
    # first numerator
    temp.3.1<-(t(temp)*x.t.temp*x.t.temp) # n1*n.vec
    #temp.3.base<-(t(temp)* v.obs.base.mat.temp* v.obs.base.mat.temp) # n1*n.vec 
    temp.3.base<-vapply(v.obs.base.1[i,],
                        function(u){rowSums(t(temp)*u*u)},
                        FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s1.3.1<-E.x.s1.3.1+rowSums(temp.3.1)
    E.x.s1.3.base<-E.x.s1.3.base+temp.3.base
    
  }
  #delta_obs=0
  
  E.x.s0<-rep(0,length(Y.vec))
  E.x.s0.1.1<-rep(0,length(Y.vec))
  E.x.s0.1.base<-rep(0,length(Y.vec))
  E.x.s0.2.1<-rep(0,length(Y.vec))
  E.x.s0.2.base<-rep(0,length(Y.vec))
  E.x.s0.3.1<-rep(0,length(Y.vec))
  E.x.s0.3.base<-rep(0,length(Y.vec))
  
  for(i in 1:n0){ 
    x.t.temp0<-x.t.fun1(Y.vec,gamma.hat,vx.obs.0[i,],theta.x.hat,phi.xi.sample.0[[i]],basis.Y.vec.x) # 1*nt
    v.obs.base.mat.temp0<-matrix(rep(v.obs.base.0[i,],length(Y.vec)),nrow=length(Y.vec),ncol=length(v.obs.base.0[i,]),byrow = TRUE)
    
    
    
    # denominators
    eta.temp0<-beta.x*x.t.temp0+as.numeric(t(beta.base)%*%(v.obs.base.0[i,]))
    temp0<-t(exp(eta.temp0)*Y.obs.ind.0[i,])*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] # n1*Y.vec
    E.x.s0<-E.x.s0+colSums(temp0)
    
    # first derivative
    # numerators
    temp0.1.1<-(t(temp0)*x.t.temp0) # n1*n.vec
    temp0.1.base<-vapply(v.obs.base.0[i,],function(u){rowSums(t(temp0)*u)},FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s0.1.1<-E.x.s0.1.1+rowSums(temp0.1.1)
    E.x.s0.1.base<-E.x.s0.1.base+temp0.1.base
    
    # first component
    temp0.2.1<-t(event.t.mat.0[i,]*x.t.temp0*Y.obs.ind.0[i,])*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] # n1*vec
    #temp0.2.base<-t(event.t.mat.0[i,]*s.t.temp0*Y.obs.ind.0[i,])*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i] # n1*vec
    temp0.2.base<-vapply(v.obs.base.0[i,],
                         function(u){colSums(matrix(rep(t(event.t.mat.0[i,]*u),length(out_k1_obs0.0[i,])),ncol=length(Y.vec),nrow=length(out_k1_obs0.0[i,]),byrow=TRUE)*out_k1_obs0.0[i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i])},
                         FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s0.2.1<-E.x.s0.2.1+colSums(temp0.2.1)
    E.x.s0.2.base<-E.x.s0.2.base+temp0.2.base
    
    # for second derivative
    temp0.3.1<-(t(temp0)*x.t.temp0*x.t.temp0) # n1*n.vec
    
    #temp0.3.base<-(t(temp0)*s.t.temp0*x.t.temp0) # n1*n.vec 
    temp0.3.base<-vapply(v.obs.base.0[i,],
                         function(u){rowSums(t(temp0)*u*u)},
                         FUN.VALUE = rep(0,length=length(Y.vec)))
    
    E.x.s0.3.1<-E.x.s0.3.1+rowSums(temp0.3.1)
    E.x.s0.3.base<-E.x.s0.3.base+temp0.3.base
    
  }
  
  deno.temp<-E.x.s1+E.x.s0
  h0.num.new<-event.t/deno.temp
  
  
  return(h0.num.new)
}


# conditional distribution of phi^T*xi

phi.xi.cond.dist<-function(x.obs,gamma.hat,vx.obs,theta.x.hat,phi.hat,rho.hat,x.i.standard,sigma.hat,Theta.hat){
  n<-length(x.obs)
  mean.all<-NULL
  var.all<-list()
  for(i in 1:n){
    #cov.x<-x.i.standard[[i]]%*%Theta.hat%*%t(x.i.standard[[i]])+diag(sigma.hat,length(x.obs[[i]]),length(x.obs[[i]]))
    cov.x<-x.i.standard[[i]]%*%Theta.hat%*%t(x.i.standard[[i]])+sigma.hat[[i]]
    cov.x.xi<-x.i.standard[[i]]%*%Theta.hat
    e.x<-t(theta.x.hat)%*%t(x.i.standard[[i]])+as.vector(gamma.hat%*%vx.obs[i,])
    mean.temp<-t(cov.x.xi)%*%solve(cov.x)%*%t(x.obs[[i]]-e.x)
    var.temp<-Theta.hat-t(cov.x.xi)%*%solve(cov.x)%*%(cov.x.xi)
    mean.all<-cbind(mean.all,mean.temp)
    var.all[[i]]<-var.temp
    
  }
  return(list("phi.xi.cond.mean"=mean.all,"phi.xi.cond.var"=var.all))
}

# generate phi^T*xi 
phi.xi.gen<-function(n,q.sample,phi.xi.cond.mean,phi.xi.cond.var){
  phi.xi.sample<-list()
  for(i in 1:n){
    phi.xi.sample[[i]]<-rmvnorm(q.sample,mean=phi.xi.cond.mean[,i],sigma=phi.xi.cond.var[[i]])
  }
  return(phi.xi.sample)
}


# initial estimation of h0
h0.new.dep.fun<-function(Y.vec,Y.obs,v.obs.base,gamma.hat,vx.obs,beta.hat,xi.hat,theta.x.hat,phi.hat,
                         b.basis,ob.x){
  beta1<-beta.hat[1]
  if(length(beta.hat)>1){
    beta.base<-beta.hat[2:length(beta.hat)]
  } else {beta.base<-0}
  
  
  Y.obs.mat<-matrix(rep(Y.obs,length(Y.vec)),ncol=length(Y.vec))
  Y.obs.ind<-t(apply(Y.obs.mat,1,function(vec){vec>=Y.vec}))
  risk.t<-colSums(Y.obs.ind)
  ind.obs.mat<-matrix(rep(ind.obs,length(Y.vec)),ncol=length(Y.vec))
  event.t.mat<-as.matrix(I(t(Y.obs.mat)==Y.vec))*t(ind.obs.mat)
  event.t<-rowSums(event.t.mat)
  phi.xi.hat<-t(phi.hat%*%xi.hat)
  eta<-NULL
  for(t in Y.vec){
    eta.temp<-exp( beta1*x.t.fun(t,gamma.hat,vx.obs,theta.x.hat,phi.xi.hat,b.basis,ob.x)+t(as.vector(beta.base)%*%t(v.obs.base)))
    eta<-cbind(eta,eta.temp)
  }
  
  eta.t<-eta*Y.obs.ind
  deno.temp<-colSums(eta.t)
  h0.num.new<-event.t/deno.temp
  h0.new<-stepfun(Y.vec,c(0,h0.num.new))
  return(list("h0.new"=h0.new,"h0.num"=h0.num.new))
}



#zstar approximation for binary observations

multiple.update<-function(n,s.hat,z.i.standard,Psi.hat){
  var.zstar.hat<-list()
  R.hat<-list()
  S.hat<-list()
  for(i in 1:n){
    var.zstar.hat[[i]]<-diag(1/deriv.inv.logit.fun(as.vector(s.hat[[i]])))
    R.hat[[i]]<-var.zstar.hat[[i]]+z.i.standard[[i]]%*%Psi.hat%*%t(z.i.standard[[i]])
    S.hat[[i]]<-Psi.hat-Psi.hat%*%t(z.i.standard[[i]])%*%solve(R.hat[[i]])%*%z.i.standard[[i]]%*%Psi.hat
  }
  return(list("var.zstar.hat"=var.zstar.hat,"R.hat"=R.hat,"S.hat"=S.hat))
}


s.hat.fun.1<-function(n,Psi.hat,theta.s.new,tau.new,v.obs,z.i.standard,var.zstar.hat,z.star.hat,zeta.eta.hat){
  s.new<-list()
  for(i in 1:n){
    #zeta.eta.new[i,]<-Psi.hat%*%t(z.i.standard[[i]])%*%solve(R.new[[i]])%*%
    #  (t(z.star.hat[[i]])-z.i.standard[[i]]%*%theta.s.hat-as.numeric(tau.hat%*%v.obs[i,]))
    s.new[[i]]<-z.i.standard[[i]]%*%theta.s.new+as.numeric(tau.new%*%v.obs[i,])+
      t(zeta.eta.hat[i,]%*%t(z.i.standard[[i]]))
  }
  return(s.new)
}


zstar.update.new<-function(n1,n0,z.obs.1,z.obs.0,v.obs.1,v.obs.0,theta.s.hat,tau.hat,
                           z.i.standard.1,z.i.standard.0,Psi.hat,s.hat.1,s.hat.0,
                           var.zstar.hat.1,R.hat.1,S.hat.1,
                           var.zstar.hat.0,R.hat.0,S.hat.0,
                           n.basis){
  
  temp1.hat<-matrix(0,n.basis,n.basis)
  
  for(i in 1:n1){
    temp1.hat<-temp1.hat+t(z.i.standard.1[[i]])%*%solve(R.hat.1[[i]])%*%z.i.standard.1[[i]]
  }
  for(i in 1:n0){
    temp1.hat<-temp1.hat+t(z.i.standard.0[[i]])%*%solve(R.hat.0[[i]])%*%z.i.standard.0[[i]]
  }
  #A.new<-solve(temp1.hat+penalty)%*%temp1.hat%*%solve(temp1.hat+penalty)
  var.zeta.eta.hat.1<-list()
  z.star.new.1<-list()
  for(i in 1:n1){
    var.zeta.eta.hat.1[[i]]<-S.hat.1[[i]]+Psi.hat%*%t(z.i.standard.1[[i]])%*%solve(R.hat.1[[i]])%*%z.i.standard.1[[i]]%*%
      solve(temp1.hat)%*%t(z.i.standard.1[[i]])%*%solve(R.hat.1[[i]])%*%z.i.standard.1[[i]]%*%Psi.hat
    tmp<-diag((z.i.standard.1[[i]])%*%var.zeta.eta.hat.1[[i]]%*%t(z.i.standard.1[[i]]))
    z.star.new.1[[i]]<-matrix((1/deriv.inv.logit.fun(s.hat.1[[i]]))*(z.obs.1[[i]]-inv.logit(s.hat.1[[i]]))+
                                s.hat.1[[i]]-0.5*(1/deriv.inv.logit.fun(s.hat.1[[i]]))*second.deriv.inv.logit.fun(s.hat.1[[i]])*tmp,nrow=1)
  }
  var.zeta.eta.hat.0<-list()
  z.star.new.0<-list()
  for(i in 1:n0){
    var.zeta.eta.hat.0[[i]]<-S.hat.0[[i]]+Psi.hat%*%t(z.i.standard.0[[i]])%*%solve(R.hat.0[[i]])%*%z.i.standard.0[[i]]%*%
      solve(temp1.hat)%*%t(z.i.standard.0[[i]])%*%solve(R.hat.0[[i]])%*%z.i.standard.0[[i]]%*%Psi.hat
    tmp<-diag((z.i.standard.0[[i]])%*%var.zeta.eta.hat.0[[i]]%*%t(z.i.standard.0[[i]]))
    z.star.new.0[[i]]<-matrix((1/deriv.inv.logit.fun(s.hat.0[[i]]))*(z.obs.0[[i]]-inv.logit(s.hat.0[[i]]))+
                                s.hat.0[[i]]-0.5*(1/deriv.inv.logit.fun(s.hat.0[[i]]))*second.deriv.inv.logit.fun(s.hat.0[[i]])*tmp,nrow=1)
  }
  return(list("z.star.new.1"=z.star.new.1,"z.star.new.0"=z.star.new.0))
}


# compute the observed likelihood and AIC

obs.like.x.fun<-function(x.obs,x.obs.t,vx.obs,gamma.hat,theta.x.hat,phi.xi.hat,sigma.hat,
                         n.knots.sele.x,n.order,t.min,t.max,ob.x){
  
  s.t.temp<-as.vector(x.t.fun(x.obs.t,gamma.hat,vx.obs,theta.x.hat,phi.xi.hat,b.basis,ob.x))
  z.t.temp<-inv.logit(s.t.temp)
  f.z.temp<-sum(x.obs*log(z.t.temp)+(1-x.obs)*log(1-z.t.temp))
  #f.x.temp<-log(dmvnorm(x.obs,mean=s.t.temp,sigma=diag(sqrt(sigma.hat),length(x.obs))))
  return(f.z.temp)
}

obs.like.fun<-function(Y.obs.1,Y.obs.0,v.inc.1,v.inc.0,
                       alpha.hat,h0.num.hat,v.obs.base.1,v.obs.base.0,beta.hat,
                       x.obs.1,x.obs.0,x.obs.t.1,x.obs.t.0,vx.obs.1,vx.obs.0,
                       gamma.hat,theta.x.hat,phi.xi.hat.1,phi.xi.hat.0,sigma.hat,
                       n.knots.sele.x,n.order,t.min,t.max){
  n<-length(c(Y.obs.1,Y.obs.0))
  Y.vec<-sort(unique(Y.obs.1))
  Y.vec<-c(0,Y.vec)
  
  each.time<-t.max
  breaks.x<-seq(t.min,t.max,length.out = n.knots.sele.x)
  n.basis.x<-length(breaks.x)+n.order-2
  full.knots.x<- expand.knots(breaks.x, order=n.order) 
  ob.x<-OrthogonalSplineBasis(full.knots.x,order=n.order)
  
  f.x<-0
  ll.temp<-0
  for(i in 1:length(Y.obs.1)){
    # distribution of x observations
    
    f.x.temp<-obs.like.x.fun(x.obs.1[[i]],x.obs.t.1[[i]],vx.obs.1[i,],gamma.hat,theta.x.hat,
                             phi.xi.hat.1[,i],sigma.hat,
                             n.knots.sele.x,n.order,t.min,t.max,ob.x)
    f.x<-f.x+f.x.temp
    
    # x at Ti
    
    x.t.temp<-x.t.fun(Y.obs.1[i],gamma.hat,vx.obs.1[i,],theta.x.hat,phi.xi.hat.1[,i],b.basis,ob.x)
    # fail rate
    pi.temp<-as.numeric(inv.logit(t(v.inc.1[i,])%*%(alpha.hat)))
    
    # survival function
    
    temp<-f.t.fun(Y.obs.1[i],h0.num.hat,Y.vec,beta.hat,gamma.hat,vx.obs.1[i,],v.obs.base.1[i,],theta.x.hat,
                  phi.xi.hat.1[,i],b.basis,ob.x)
    
    ll.temp<-ll.temp+log(temp)
    
  }
  
  for(i in 1:length(Y.obs.0)){
    # distribution of x observations
    
    f.x.temp<-obs.like.x.fun(x.obs.0[[i]],x.obs.t.0[[i]],vx.obs.0[i,],gamma.hat,theta.x.hat,
                             phi.xi.hat.0[,i],sigma.hat,
                             n.knots.sele.x,n.order,t.min,t.max,ob.x)
    f.x<-f.x+f.x.temp
    
    # x at Ti
    
    x.t.temp<-x.t.fun(Y.obs.0[i],gamma.hat,vx.obs.0[i,],theta.x.hat,
                      phi.xi.hat.0[,i],b.basis,ob.x)
    # fail rate
    pi.temp<-as.numeric(inv.logit(t(v.inc.0[i,])%*%(alpha.hat)))
    
    # survival function
    
    temp<-S.t.fun(Y.obs.0[i],h0.num.hat,Y.vec,beta.hat,gamma.hat,vx.obs.0[i,],v.obs.base.0[i,],theta.x.hat,
                  phi.xi.hat.0[,i],b.basis,ob.x)
    
    ll.temp<-ll.temp+log(temp)
    
  }
  
  res<-f.x+ll.temp
  return(res)
}

AIC.fun<-function(Y.obs.1,Y.obs.0,v.inc.1,v.inc.0,alpha.hat,h0.hat,v.obs.base.1,v.obs.base.0,beta.hat,
                  x.obs.1,x.obs.0,x.obs.t.1,x.obs.t.0,vx.obs.1,vx.obs.0,
                  gamma.hat,theta.x.hat,phi.xi.hat.1,phi.xi.hat.0,sigma.hat,
                  n.knots.sele.x,n.order,t.min,t.max,n.phi){
  ll<-obs.like.fun(Y.obs.1,Y.obs.0,v.inc.1,v.inc.0,alpha.hat,h0.hat,v.obs.base.1,v.obs.base.0,beta.hat,
                   x.obs.1,x.obs.0,x.obs.t.1,x.obs.t.0,vx.obs.1,vx.obs.0,
                   gamma.hat,theta.x.hat,phi.xi.hat.1,phi.xi.hat.0,sigma.hat,
                   n.knots.sele.x,n.order,t.min,t.max)
  each.time<-t.max
  breaks.x<-seq(t.min,t.max,length.out = n.knots.sele.x)
  n.basis.x<-length(breaks.x)+n.order-2
  
  AIC.res<--2*ll+2*(n.basis.x*(n.phi+1)+length(alpha.hat)+length(gamma.hat)+1)
  return(AIC.res)
}


