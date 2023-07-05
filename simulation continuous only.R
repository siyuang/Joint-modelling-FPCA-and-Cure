# FPCA and Mixture cure model simulation #
# This program is used for the numerical simulation studies of the FPCA and Mixture cure model simulation #
# Loading functions #
rm(list=ls())

source("estimation function joint modelling.R")
source("estimation function continuous.R")
##############################################
# Data generation                            #
# Observations:                              #
# Observation time: Y=min(T,C)               #
# Longitudinal observation: X                #
# Scalar observation: V                      #
##############################################


# Define frames for the final results
# incidence
alpha.fin<-NULL
# components of the x observation
gamma.fin<-NULL
theta.x.fin<-list()
phi.fin<-list()
sigma.fin<-NULL
rho.fin<-list()

# components of the z observation
tau.fin<-NULL
theta.s.fin<-NULL
eta.fin<-list()
lambda.fin<-NULL

# survival function
beta.fin<-NULL
h0.fin<-list()

# dimension
n.phi.fin<-NULL
n.knots.fin<-NULL

args<-1
pos<-1

#########################
# repetition starts here#
#                       #
#########################

for(ite.fin in args*100+c(1:10)){
  set.seed(ite.fin)
  n<-600
  t.min<-0
  t.max<-10
  # First, generate the incident rate #
  v1<-rbinom(n,1,0.5)
  v2<-rbinom(n,1,0.7)
  v3<-rnorm(n,0,1)
  v.obs<-cbind(v1,v3)
  v.inc<-cbind(1,v.obs)
  v.obs.base<-cbind(v2,v3)
  vx.obs<-cbind(v1)
  alpha<-c(0.5,-0.2,1)
  k.prob<-exp(alpha%*%t(v.inc))/(1+exp(alpha%*%t(v.inc)))
  k.true<-rbinom(n,1,as.vector(k.prob))
  table(k.true)
  
  # for the event case ,we need to generate the functional processed and the survival time #
  # suppose domain is 0:10 #
  gamma<-0.1
  mu.t.x<-function(t){
    rep(0,length(t))
  }
  mu.x.fun<-function(t,gamma,v.obs){
    mu.t.x(t)+as.vector(gamma%*%t(v.obs))
  }
  phi1.fun<-function(t){
    #-cos(pi*t/10)/sqrt(5)
    sin(2*pi*t/10)/sqrt(5)
  }
  phi2.fun<-function(t){
    #sin(pi*t/10)/sqrt(5)
    -cos(2*pi*t/10)/sqrt(5)
  }
  lambda1.true<-5
  lambda2.true<-2
  xi1.true<-rnorm(n,0,sqrt(lambda1.true))
  xi2.true<-rnorm(n,0,sqrt(lambda2.true))
  
  x.true.fun<-function(t,mu.x.fun,phi1.fun,phi2.fun,xi1.true,xi2.true,gamma,v.obs){
    true<-mu.x.fun(t,gamma,v.obs)+xi1.true*phi1.fun(t)+xi2.true*phi2.fun(t)
    return(true)
  }
  
  # event time #
  # PH model is assumed #
  
  # baseline hazard #
  h0.fun<-function(t){
    t^2/100
  }
  
  # 
  beta1.true<-1
  beta.base.true<-c(1,-0.5)
  
  # survival function
  h.fun<-function(t,h0.fun,x.true.fun,beta1.true,beta.base.true,mu.x.fun,phi1.fun,phi2.fun,xi1.true,xi2.true,gamma,v.obs,v.obs.base){
    eta.temp<-beta1.true*x.true.fun(t,mu.x.fun,phi1.fun,phi2.fun,xi1.true,xi2.true,gamma,v.obs)+as.numeric(beta.base.true%*%v.obs.base)
    h0.fun(t)*exp(eta.temp)
  }
  Fbar.fun<-function(t,h0.fun,x.true.fun,beta1.true,beta.base.true,mu.x.fun,phi1.fun,phi2.fun,xi1.true,xi2.true,gamma,v.obs,v.obs.base){
    exp(-stats::integrate(h.fun,0,t,h0.fun=h0.fun,x.true.fun=x.true.fun,beta1.true=beta1.true,beta.base.true=beta.base.true,
                          mu.x.fun=mu.x.fun,phi1.fun=phi1.fun,phi2.fun=phi2.fun,
                          xi1.true=xi1.true,xi2.true=xi2.true,gamma=gamma,v.obs=v.obs,v.obs.base=v.obs.base)$value)
  }
  Fbar.fun.gen<-function(t,h0.fun,x.true.fun,beta1.true,beta.base.true,mu.x.fun,phi1.fun,phi2.fun,xi1.true,xi2.true,gamma,v.obs,v.obs.base,
                         Fbar.true){
    Fbar.fun(t,h0.fun,x.true.fun,beta1.true,beta.base.true,mu.x.fun,phi1.fun,phi2.fun,xi1.true,xi2.true,gamma,v.obs,v.obs.base)-Fbar.true
  }

  # true event time
  Fbar.true<-runif(n,0,1)
  t.true<-NULL
  t.true.ind<-NULL
  for(i in 1:n){
    t.true.temp<-try(uniroot(Fbar.fun.gen,interval=c(0,100),
                             h0.fun,x.true.fun,beta1.true,beta.base.true,mu.x.fun,phi1.fun,phi2.fun,xi1.true[i],xi2.true[i],
                             gamma,vx.obs[i,],v.obs.base[i,],Fbar.true[i])$root,TRUE)
    if(!inherits(t.true.temp, "try-error")){
      if(t.true.temp>0){
        t.true<-c(t.true,t.true.temp)
        t.true.ind<-c(t.true.ind,i)
      }
    }
  }
  length(t.true)
  summary(t.true)
  
  
  # consider cured #
  t.true[k.true==0]<-10
  summary(t.true)
  # censoring generate #
  cen.shape<-8
  cen.scale<-12
  cen.scale*gamma(1+1/cen.shape)
  cen.true<-rweibull(n,cen.shape,cen.scale)
  summary(cen.true)
  
  # temporary observation time #
  ind.obs.true<-as.numeric(I(t.true<=cen.true))
  table(ind.obs.true)
  t.obs.true<-t.true
  t.obs.true[ind.obs.true==0]<-cen.true[ind.obs.true==0]
  t.obs<-t.obs.true
  summary(t.obs)
  # censoring indicator, consider cure proportion #
  ind.obs<-ind.obs.true*k.true
  table(ind.obs)
  
  # truncate at 10 #
  ind.obs[t.obs>10]<-0
  t.obs[t.obs>10]<-10
  
  table(t.obs==10,ind.obs==1)
  summary(t.obs)
  
  # generate longitudinal continuous process observation #
  
  x.obs<-list()
  x.obs.t<-list()
  for(i in 1:n){
    n.sample<-sample(c(6:10),1,prob=rep(1/5,5))
    t.temp<-sort(runif(n.sample,t.min,t.obs[i]))
    #t.obs.temp<-t.temp[t.temp<=t.obs[i]] 
    t.obs.temp<-t.temp
    x.obs.t[[i]]<-t.obs.temp
    s.obs.temp<-x.true.fun(t.obs.temp,mu.x.fun,phi1.fun,phi2.fun,xi1.true[i],xi2.true[i],gamma,vx.obs[i,])
    #prob.obs.temp<-exp(s.obs.temp)/(1+exp(s.obs.temp))
    x.obs[[i]]<-s.obs.temp+rnorm(n.sample,0,sqrt(0.1))
  }
  
  # observation set #
  Y.obs<-t.obs
  ind.obs<-ind.obs
  x.obs<-x.obs
  x.obs.t<-x.obs.t
  vx.obs<-vx.obs
  
  
  ##############################################
  # Estimation:                                #
  # --- Observations:                      --- #
  # Observation time: Y=min(T,C)               #
  # Indicator delta=I(T<=C)                    #
  # Longitudinal observation: X                #
  # Scalar observation: V                      #
  # --- Paramters:                         --- #
  # FPCA of X(t)                               #
  # coefficient                                #
  #                                            #
  ##############################################
  
  # EM algorithm were applied
  # E-step: calculate expected log-complete function
  # settings
  
  
  Y.obs.1<-Y.obs[ind.obs==1]
  ind.obs.1<-ind.obs[ind.obs==1]
  v.obs.1<-as.matrix(v.obs[ind.obs==1,])
  v.inc.1<-v.inc[ind.obs==1,]
  x.obs.1<-x.obs[ind.obs==1]
  x.obs.t.1<-x.obs.t[ind.obs==1]
  v.obs.base.1<-as.matrix(v.obs.base[ind.obs==1,])
  vx.obs.1<-as.matrix(vx.obs[ind.obs==1,])

  
  Y.obs.0<-Y.obs[ind.obs==0]
  ind.obs.0<-ind.obs[ind.obs==0]
  v.obs.0<-as.matrix(v.obs[ind.obs==0,])
  v.inc.0<-v.inc[ind.obs==0,]
  x.obs.0<-x.obs[ind.obs==0]
  x.obs.t.0<-x.obs.t[ind.obs==0]
  v.obs.base.0<-as.matrix(v.obs.base[ind.obs==0,])
  vx.obs.0<-as.matrix(vx.obs[ind.obs==0,])
  
  
  n1<-length(ind.obs.1)
  n0<-length(ind.obs.0)
  n.phi.can<-c(1:2)
  n.knots.can<-seq(8:10)
  par.can<-expand.grid(n.phi.can,n.knots.can)
  
  
  alpha.fin.temp<-NULL
  
  gamma.fin.temp<-NULL
  theta.x.fin.temp<-list()
  eta.fin.temp<-list()
  lambda.fin.temp<-list()
  
  beta.fin.temp<-NULL
  h0.fin.temp<-list()
  
  n.phi.fin.temp<-NULL
  n.knots.fin.temp<-NULL
  AIC.res<-NULL
  
  
  
  for(cand in 1:nrow(par.can)){
    #------- Given initial estimates -------#
    # pre-specified parameters: B-spline bases, FPCs
    n.phi<-par.can[cand,1]
    n.knots.sele.x<-par.can[cand,2]
    n.order<-4     
    P<-0.9
    
    t.min<-0
    t.max<-10
    # alpha in the incidence
    alpha.hat<-rep(0,ncol(v.inc))
    
    # FPCA  

    tau1.sele.x<-2
    tau2.sele.x<-2
    
    
    
    # components of the z observation
    
    each.time<-t.max
    #breaks<-    1+  (0:(n.knots-1))*((each.time-1)/(n.knots-1))
    breaks.x<-seq(t.min,t.max,length.out = n.knots.sele.x)
    
    # degree = order-1
    n.basis.x<-length(breaks.x)+n.order-2
    
    ##### Threshold of selecting principal component
    
    # generate b-spline basis
    full.knots.x<- expand.knots(breaks.x, order=n.order) 
    ob.x<-OrthogonalSplineBasis(full.knots.x,order=n.order)
    X.i.standard.x<-list()
    for(i in 1:n){
      X.i.standard.x[[i]]<-evaluate(ob.x,x.obs.t[[i]]) 
    }
    x.i.standard.1<-X.i.standard.x[ind.obs==1]
    x.i.standard.0<-X.i.standard.x[ind.obs==0]
    # D.temp is penalty matrix
    D.temp.x<-OuterProdSecondDerivative(ob.x)
    
    
    estimate.x.res<-estimation.cont(x.obs,X.i.standard.x,tau1.sele.x,tau2.sele.x,D.temp.x,n.basis.x,vx.obs)
    X1.lambda.fin<-estimate.x.res$lambda
    X1.theta.fin<-estimate.x.res$theta
    X1.beta.fin<-estimate.x.res$beta
    X1.psi.fin<-estimate.x.res$psi.fin
    X1.sigma.fin<-estimate.x.res$var.err
    
    
    t.grid.standard<-seq(t.min,t.max,0.005)
    nu.standard.x<-evaluate(ob.x,t.grid.standard)
    pc1.X1.temp<-nu.standard.x%*%X1.theta.fin[,1]
    if(pc1.X1.temp[which.max(abs(pc1.X1.temp[1:(length(pc1.X1.temp)/2)]))]>0){
      X1.theta.fin[,1]<--X1.theta.fin[,1]
    }
    if(ncol(X1.lambda.fin)==2){ 
      pc2.X1.temp<-nu.standard.x%*%X1.theta.fin[,2]
      if(pc2.X1.temp[(length(pc2.X1.temp)/4):(length(pc2.X1.temp)*3/4)][which.max(abs(pc2.X1.temp[(length(pc2.X1.temp)/4):(length(pc2.X1.temp)*3/4)]))]<0){
        X1.theta.fin[,2]<--X1.theta.fin[,2]
      }
    }
    xi.hat<-xi.estimation.cont(log.likelihood.cont,n,x.obs,X1.beta.fin,X1.theta.fin,
                                X.i.standard.x,X1.lambda.fin,vx.obs)
    
    
    gamma.hat<-X1.beta.fin[(n.basis.x+1):length(X1.beta.fin)]
    theta.x.hat<-X1.beta.fin[1:n.basis.x]
    eta.hat<-X1.theta.fin
    lambda.hat<-X1.lambda.fin
    Theta.hat<-X1.psi.fin
    sigma.hat<-as.numeric(X1.sigma.fin)

    # survival function
    beta.hat<-c(0,rep(0,ncol(v.obs.base)))
    Y.vec<-sort(unique(Y.obs[ind.obs==1]))
    Y.vec<-c(0,Y.vec)
    h0.num.hat<-h0.new.dep.fun(Y.vec,Y.obs,v.obs.base,gamma.hat,vx.obs,beta.hat,
                               t(xi.hat),theta.x.hat,eta.hat,
                               b.basis,ob.x)$h0.num
    # generate samples for monte carlo integration
    q.sample<-50
    
    
    phi.xi.cond.dist.hat.1<-phi.xi.cond.dist(x.obs.1,gamma.hat,vx.obs.1,theta.x.hat,
                                               eta.hat,lambda.hat,x.i.standard.1,sigma.hat,Theta.hat)
    phi.xi.cond.mean.hat.1<-phi.xi.cond.dist.hat.1$phi.xi.cond.mean
    phi.xi.cond.var.hat.1<-phi.xi.cond.dist.hat.1$phi.xi.cond.var
    
    phi.xi.sample.1<-phi.xi.gen(n1,q.sample,phi.xi.cond.mean.hat.1,phi.xi.cond.var.hat.1)
    
    phi.xi.cond.dist.hat.0<-phi.xi.cond.dist(x.obs.0,gamma.hat,vx.obs.0,theta.x.hat,
                                               eta.hat,lambda.hat,x.i.standard.0,sigma.hat,Theta.hat)
    phi.xi.cond.mean.hat.0<-phi.xi.cond.dist.hat.0$phi.xi.cond.mean
    phi.xi.cond.var.hat.0<-phi.xi.cond.dist.hat.0$phi.xi.cond.var
    
    phi.xi.sample.0<-phi.xi.gen(n0,q.sample,phi.xi.cond.mean.hat.0,phi.xi.cond.var.hat.0)
    
    lambda.hat<-lambda.hat[1:n.phi,1:n.phi]
    eta.hat<-eta.hat[,1:n.phi]
    ###############################
    #--- iteration starts here ---#
    ###############################
    iteration<-1
    like.hat<-99999
    #while(iteration<=50|(iteration<=500&diff3>=abs(like.hat*0.01))){
    while(iteration<=50){
      #par.hat<-c(alpha.hat,gamma.hat,theta.x.hat,unlist(phi.hat),sigma.hat,diag(rho.hat),
      #           gamma.hat,unlist(theta.x.hat),unlist(eta.hat),diag(lambda.hat),beta.hat)
      #par.hat<-c(alpha.hat,gamma.hat,theta.x.hat,unlist(Theta.hat),sigma.hat,
      #           gamma.hat,unlist(theta.x.hat),unlist(Theta.hat),beta.hat)
      par.hat1<-c(alpha.hat,gamma.hat,beta.hat)
      par.hat2<-c(theta.x.hat,lambda.hat,eta.hat)
      
      q.sample<-50+floor(iteration/20)*5
      #q.sample<-10
      
      
      # common components received from the current estimates
      
      out_obs1<-f_obs_delta1(phi.xi.sample.1,
                             Y.obs.1,vx.obs.1,v.obs.base.1,
                             h0.num.hat,Y.vec,beta.hat,
                             gamma.hat,theta.x.hat,x.obs.1,
                             b.basis,ob.x)
      out_obs1.1<-t(out_obs1)/q.sample
      
      out_obs0<-f_obs0_int(phi.xi.sample.0,
                           Y.obs.0,vx.obs.0,v.obs.base.0,
                           h0.num.hat,Y.vec,beta.hat,
                           gamma.hat,theta.x.hat,x.obs.0,
                           b.basis,ob.x)
      
      out_k1_obs0<-t(out_obs0)
      #out_k1_obs0[Y.obs.0==10]<-0
      out_k0_obs0<-1
      
      out_k1_obs0.0<-out_k1_obs0/q.sample
      out_k0_obs0.0<-out_k0_obs0/q.sample
      
      out_deno_obs1<-rowSums(out_obs1.1)
      #out_deno_obs1[which(out_k1_obs0==0)]<-1
      out_deno_obs0_k1<-rowSums(out_k1_obs0.0)
      out_deno_obs0_k1[which(out_deno_obs0_k1==0)]<-1
      out_deno_obs0_k0<-1
      
      # delta=0
      # k1=1
      p_k1_nu<-alpha.hat%*%t(v.inc.0)
      p_k1<-as.vector(exp(p_k1_nu)/(1+exp(p_k1_nu))) # 1*n #
      f_obs0_k1<-rowSums(p_k1*out_k1_obs0.0)
      
      #delta=0
      # k1=0
      p_k0<-1-p_k1 # 1*n #
      f_obs0_k0<-p_k0
      
      f_obs0_k1_fin<-f_obs0_k1/(f_obs0_k1+f_obs0_k0)
      f_obs0_k1_fin[is.na(f_obs0_k1_fin)]<-0
      
      f_obs0_k0_fin<-f_obs0_k0/(f_obs0_k1+f_obs0_k0)
      f_obs0_k0_fin[is.na(f_obs0_k0_fin)]<-1
      head(f_obs0_k0_fin)
      summary(f_obs0_k0_fin)
      summary(out_deno_obs1)
      
      # ---- component 1: incident rate ---- #
      #alpha.new<-optimize(exp_inc,interval=c(-100,100),vx.obs.1,vx.obs.0,out_obs1.1,out_deno_obs1,out_k1_obs0.0,out_k0_obs0.0,f_obs0_k1_fin)$minimum
      par.ini<-alpha.hat
      alpha.new<-stats::optim(par.ini,exp_inc,v.inc.1=v.inc.1,v.inc.0=v.inc.0,out_obs1.1=out_obs1.1,
                              out_deno_obs1=out_deno_obs1,out_k1_obs0.0=out_k1_obs0.0,
                              out_k0_obs0.0=out_k0_obs0.0,f_obs0_k1_fin=f_obs0_k1_fin)$par
      
      
      # ---- component 2: mean functions of X and Z ---- #
      
      #  measurement error  # 
      sigma.new<-sigma_hat_ker(phi.xi.sample.1,phi.xi.sample.0,
                               Y.obs.1,x.obs.1,x.obs.t.1,vx.obs.1,n1,
                               Y.obs.0,x.obs.0,x.obs.t.0,vx.obs.0,n0,
                               gamma.hat,theta.x.hat,
                               out_obs1.1,out_deno_obs1,
                               out_k1_obs0.0,out_k0_obs0.0,f_obs0_k1_fin,f_obs0_k0_fin,
                               out_deno_obs0_k1,out_deno_obs0_k0,
                               x.i.standard.1,x.i.standard.0)
      # mean function of S(t)
     
      # mean function coefficients theta and gamma #
      left<-matrix(0,nrow=n.basis.x+ncol(vx.obs),ncol=n.basis.x+ncol(vx.obs))
      right<-matrix(0,nrow=n.basis.x+ncol(vx.obs),ncol=1)
      for(i in 1:n1){
        Ax.temp<-cbind(x.i.standard.1[[i]],vx.obs.1[i,])
        W.temp<-x.i.standard.1[[i]]%*%Theta.hat%*%t(x.i.standard.1[[i]])+sigma.new*diag(1,length(x.obs.1[[i]]))
        temp1<-t(Ax.temp)%*%solve(W.temp)%*%Ax.temp
        temp2<-t(Ax.temp)%*%solve(W.temp)%*%(x.obs.1[[i]])
        left<-left+temp1
        right<-right+temp2
      }
      for(i in 1:n0){
        Ax.temp<-cbind(x.i.standard.0[[i]],vx.obs.0[i,])
        W.temp<-x.i.standard.0[[i]]%*%Theta.hat%*%t(x.i.standard.0[[i]])+sigma.new*diag(1,length(x.obs.0[[i]]))
        temp1<-t(Ax.temp)%*%solve(W.temp)%*%Ax.temp
        temp2<-t(Ax.temp)%*%solve(W.temp)%*%(x.obs.0[[i]])
        left<-left+temp1
        right<-right+temp2
      }
      theta_tau_new<-solve(left)%*%right
      theta.x.new<-theta_tau_new[1:n.basis.x]
      gamma.new<-theta_tau_new[(n.basis.x+1):length(theta_tau_new)]
      
      
      
      
      
      
      phi.xi.cond.dist.hat.1<-phi.xi.cond.dist(x.obs.1,gamma.new,vx.obs.1,theta.x.new,
                                                 eta.hat,lambda.hat,x.i.standard.1,sigma.hat,Theta.hat)
      phi.xi.cond.mean.hat.1<-phi.xi.cond.dist.hat.1$phi.xi.cond.mean
      phi.xi.cond.var.hat.1<-phi.xi.cond.dist.hat.1$phi.xi.cond.var
      
      phi.xi.sample.1<-phi.xi.gen(n1,q.sample,phi.xi.cond.mean.hat.1,phi.xi.cond.var.hat.1)
      
      phi.xi.cond.dist.hat.0<-phi.xi.cond.dist(x.obs.0,gamma.new,vx.obs.0,theta.x.new,
                                                 eta.hat,lambda.hat,x.i.standard.0,sigma.hat,Theta.hat)
      phi.xi.cond.mean.hat.0<-phi.xi.cond.dist.hat.0$phi.xi.cond.mean
      phi.xi.cond.var.hat.0<-phi.xi.cond.dist.hat.0$phi.xi.cond.var
      
      phi.xi.sample.0<-phi.xi.gen(n0,q.sample,phi.xi.cond.mean.hat.0,phi.xi.cond.var.hat.0)
      

      out_obs1<-f_obs_delta1(phi.xi.sample.1,
                             Y.obs.1,vx.obs.1,v.obs.base.1,
                             h0.num.hat,Y.vec,beta.hat,
                             gamma.new,theta.x.new,sigma.hat,
                             b.basis,ob.x)
      out_obs1.1<-t(out_obs1)/q.sample
      
      out_obs0<-f_obs0_int(phi.xi.sample.0,
                           Y.obs.0,vx.obs.0,v.obs.base.0,
                           h0.num.hat,Y.vec,beta.hat,
                           gamma.new,theta.x.new,sigma.hat,
                           b.basis,ob.x)
      
      out_k1_obs0<-t(out_obs0)
      out_k0_obs0<-1
      out_k1_obs0.0<-out_k1_obs0/q.sample
      out_k0_obs0.0<-out_k0_obs0/q.sample
      
      out_deno_obs1<-rowSums(out_obs1.1)
      
      out_deno_obs0_k1<-rowSums(out_k1_obs0.0)
      out_deno_obs0_k1[which(out_deno_obs0_k1==0)]<-1
      out_deno_obs0_k0<-1
      
      # delta=0
      # k1=1
      p_k1_nu<-alpha.hat%*%t(v.inc.0)
      p_k1<-as.vector(exp(p_k1_nu)/(1+exp(p_k1_nu))) # 1*n #
      f_obs0_k1<-rowSums(p_k1*out_k1_obs0.0)
      
      #delta=0
      # k1=0
      p_k0<-1-p_k1 # 1*n #
      f_obs0_k0<-p_k0
      
      f_obs0_k1_fin<-f_obs0_k1/(f_obs0_k1+f_obs0_k0)
      f_obs0_k1_fin[is.na(f_obs0_k1_fin)]<-0
      
      f_obs0_k0_fin<-f_obs0_k0/(f_obs0_k1+f_obs0_k0)
      f_obs0_k0_fin[is.na(f_obs0_k0_fin)]<-1
      
      # --- components of continuous observation X --- #
      # update FPC components 
      
      # the coefficient of FPC

      
      
      Theta.temp<-matrix(0,n.basis.x,n.basis.x)
      gamma.xi.sample.1<-list()
      gamma.xi.sample.0<-list()
      for(i in 1:n1){
        gamma.xi.sample.1[[i]]<-t(phi.xi.sample.1[[i]])
        E.xi.temp1<-(gamma.xi.sample.1[[i]])%*%(t(gamma.xi.sample.1[[i]])*out_obs1.1[i,]/out_deno_obs1[i])
        Theta.temp<-Theta.temp+E.xi.temp1
        
      }
      for(i in 1:n0){
        gamma.xi.sample.0[[i]]<-t(phi.xi.sample.0[[i]])
        E.xi.temp1<-(gamma.xi.sample.0[[i]])%*%(t(gamma.xi.sample.0[[i]])*
                                                      (out_k1_obs0.0 [i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i]+
                                                         out_k0_obs0.0*f_obs0_k0_fin[i]))
        Theta.temp<-Theta.temp+E.xi.temp1
      }
      
      Theta.eigen.temp<-eigen(Theta.temp/n)
      lambda.new<-diag(Theta.eigen.temp$values[1:n.phi],ncol = n.phi,nrow=n.phi)
      eta.new<-as.matrix(Theta.eigen.temp$vectors[,1:n.phi])
      
      t.grid.standard<-seq(t.min,t.max,0.005)
      nu.standard.x<-evaluate(ob.x,t.grid.standard)
      
      
      pc1.X1.temp<-nu.standard.x%*%eta.new[,1]
      if(pc1.X1.temp[which.max(abs(pc1.X1.temp[1:(length(pc1.X1.temp)/2)]))]>0){
        X1.theta.fin[,1]<--X1.theta.fin[,1]
      }
      if(ncol(X1.lambda.fin)==2){ 
        pc2.X1.temp<-nu.standard.x%*%X1.theta.fin[,2]
        if(pc2.X1.temp[(length(pc2.X1.temp)/4):(length(pc2.X1.temp)*3/4)][which.max(abs(pc2.X1.temp[(length(pc2.X1.temp)/4):(length(pc2.X1.temp)*3/4)]))]<0){
          X1.theta.fin[,2]<--X1.theta.fin[,2]
        }
      }
      
      
      
      Theta.new<-eta.new%*%lambda.new%*%t(eta.new)

 
      phi.xi.cond.dist.hat.1<-phi.xi.cond.dist(x.obs.1,gamma.new,vx.obs.1,theta.x.new,
                                                 eta.new,lambda.new,x.i.standard.1,sigma.hat,Theta.new)
      phi.xi.cond.mean.hat.1<-phi.xi.cond.dist.hat.1$phi.xi.cond.mean
      phi.xi.cond.var.hat.1<-phi.xi.cond.dist.hat.1$phi.xi.cond.var
      
      phi.xi.sample.1<-phi.xi.gen(n1,q.sample,phi.xi.cond.mean.hat.1,phi.xi.cond.var.hat.1)
      
      phi.xi.cond.dist.hat.0<-phi.xi.cond.dist(x.obs.0,gamma.new,vx.obs.0,theta.x.new,
                                                 eta.new,lambda.new,x.i.standard.0,sigma.hat,Theta.new)
      phi.xi.cond.mean.hat.0<-phi.xi.cond.dist.hat.0$phi.xi.cond.mean
      phi.xi.cond.var.hat.0<-phi.xi.cond.dist.hat.0$phi.xi.cond.var
      
      phi.xi.sample.0<-phi.xi.gen(n0,q.sample,phi.xi.cond.mean.hat.0,phi.xi.cond.var.hat.0)
      
      xi.eta.hat.1<-t(phi.xi.cond.mean.hat.1)
      xi.eta.hat.0<-t(phi.xi.cond.mean.hat.0)
      
      out_obs1<-f_obs_delta1(phi.xi.sample.1,
                             Y.obs.1,vx.obs.1,v.obs.base.1,
                             h0.num.hat,Y.vec,beta.hat,
                             gamma.new,theta.x.new,sigma.hat,
                             b.basis,ob.x)
      out_obs1.1<-t(out_obs1)/q.sample
      
      out_obs0<-f_obs0_int(phi.xi.sample.0,
                           Y.obs.0,vx.obs.0,v.obs.base.0,
                           h0.num.hat,Y.vec,beta.hat,
                           gamma.new,theta.x.new,sigma.hat,
                           b.basis,ob.x)
      
      out_k1_obs0<-t(out_obs0)
      out_k0_obs0<-1
      
      out_k1_obs0.0<-out_k1_obs0/q.sample
      out_k0_obs0.0<-out_k0_obs0/q.sample
      
      out_deno_obs1<-rowSums(out_obs1.1)
      
      out_deno_obs0_k1<-rowSums(out_k1_obs0.0)
      out_deno_obs0_k1[which(out_deno_obs0_k1==0)]<-1
      out_deno_obs0_k0<-1
      
      # delta=0
      # k1=1
      p_k1_nu<-alpha.hat%*%t(v.inc.0)
      p_k1<-as.vector(exp(p_k1_nu)/(1+exp(p_k1_nu))) # 1*n #
      f_obs0_k1<-rowSums(p_k1*out_k1_obs0.0)
      
      #delta=0
      # k1=0
      p_k0<-1-p_k1 # 1*n #
      f_obs0_k0<-p_k0
      
      f_obs0_k1_fin<-f_obs0_k1/(f_obs0_k1+f_obs0_k0)
      f_obs0_k1_fin[is.na(f_obs0_k1_fin)]<-0
      
      f_obs0_k0_fin<-f_obs0_k0/(f_obs0_k1+f_obs0_k0)
      f_obs0_k0_fin[is.na(f_obs0_k0_fin)]<-1
      # survival function components #
      
      basis.Y.vec.z<-b.basis(ob.x,t=Y.vec)
      
      beta.new<-rootSolve::multiroot(first.deriv.fun,start=beta.hat,
                                     Y.vec=Y.vec,v.obs.base.1=v.obs.base.1,v.obs.base.0=v.obs.base.0,
                                     phi.xi.sample.1=phi.xi.sample.1,phi.xi.sample.0=phi.xi.sample.0,
                                     Y.obs.1=Y.obs.1,x.obs.1=x.obs.1,x.obs.t.1=x.obs.t.1,
                                     vx.obs.1=vx.obs.1,ind.obs.1=ind.obs.1,n1=n1,
                                     Y.obs.0=Y.obs.0,x.obs.0=x.obs.0,x.obs.t.0=x.obs.t.0,
                                     vx.obs.0=vx.obs.0,ind.obs.0=ind.obs.0,n0=n0,
                                     gamma.hat=gamma.new,theta.x.hat=theta.x.new,sigma.hat=NULL,
                                     out_obs1.1=out_obs1.1,out_deno_obs1=out_deno_obs1,
                                     out_k1_obs0.0=out_k1_obs0.0,out_k0_obs0.0=out_k0_obs0.0,
                                     f_obs0_k1_fin= f_obs0_k1_fin,f_obs0_k0_fin=f_obs0_k0_fin,
                                     out_deno_obs0_k1=out_deno_obs0_k1,out_deno_obs0_k0=out_deno_obs0_k0,
                                     x.i.standard.1=x.i.standard.1,x.i.standard.0=x.i.standard.0,
                                     basis.Y.vec.x=basis.Y.vec.z)$root
      
      
      
      h0.num.new<-h0.new.fun(beta.new,
                             Y.vec=Y.vec,v.obs.base.1=v.obs.base.1,v.obs.base.0=v.obs.base.0,
                             phi.xi.sample.1=phi.xi.sample.1,phi.xi.sample.0=phi.xi.sample.0,
                             Y.obs.1=Y.obs.1,x.obs.1=x.obs.1,x.obs.t.1=x.obs.t.1,
                             vx.obs.1=vx.obs.1,ind.obs.1=ind.obs.1,n1=n1,
                             Y.obs.0=Y.obs.0,x.obs.0=x.obs.0,x.obs.t.0=x.obs.t.0,
                             vx.obs.0=vx.obs.0,ind.obs.0=ind.obs.0,n0=n0,
                             gamma.hat=gamma.new,theta.x.hat=theta.x.new,sigma.hat=NULL,
                             out_obs1.1=out_obs1.1,out_deno_obs1=out_deno_obs1,
                             out_k1_obs0.0=out_k1_obs0.0,out_k0_obs0.0=out_k0_obs0.0,
                             f_obs0_k1_fin= f_obs0_k1_fin,f_obs0_k0_fin=f_obs0_k0_fin,
                             out_deno_obs0_k1=out_deno_obs0_k1,out_deno_obs0_k0=out_deno_obs0_k0,
                             x.i.standard.1=x.i.standard.1,x.i.standard.0=x.i.standard.0,
                             basis.Y.vec.x=basis.Y.vec.z)
      
      phi.xi.hat.1<-phi.xi.cond.mean.hat.1
      phi.xi.hat.0<-phi.xi.cond.mean.hat.0
      
      vx.obs.1<-vx.obs.1
      vx.obs.0<-vx.obs.0
      
      like.new<-obs.like.fun(Y.obs.1,Y.obs.0,v.inc.1,v.inc.0,alpha,h0.num.hat,v.obs.base.1,v.obs.base.0,beta.hat,
                             x.obs.1,x.obs.0,x.obs.t.1,x.obs.t.0,vx.obs.1,vx.obs.0,
                             gamma.hat,theta.x.hat,phi.xi.hat.1,phi.xi.hat.0,sigma.hat,
                             n.knots.sele.x,n.order,t.min,t.max)
      
      # compare difference between new and hat outcomes ##
      
      par.new1<-c(alpha.new,gamma.new,beta.new)
      par.new2<-c(theta.x.new,lambda.new,eta.new)
      diff1<-max(abs(par.hat1-par.new1))
      diff2<-max(abs(par.hat2-par.new2))
      diff3<-max(abs(like.new-like.hat))
      print(paste("Diff:",c(diff1,diff2,diff3)))
      print(paste("iteration:",iteration))
      print(paste("Max Diff1:",which.max(abs(par.hat1-par.new1))))
      print(paste("Max Diff2:",which.max(abs(par.hat2-par.new2))))
      print(paste("beta:",beta.new))
      print(paste("lambda:",lambda.new))
      if(diff1<=0.0001 & diff2<=0.001) break
      if(iteration>=50&diff3<abs(like.hat)*0.01) break
      # --- set all new to hat --- # 
      alpha.hat<-alpha.new
      
      # components of the z observation
      gamma.hat<-gamma.new
      theta.x.hat<-theta.x.new
      eta.hat<-eta.new
      lambda.hat<-lambda.new
      
      like.hat<-like.new
      Theta.hat<-Theta.new
      
      
      # survival function
      
      
      beta.hat<-as.vector(beta.new)
      h0.num.hat<-h0.num.new
      
      iteration<-iteration+1
      
      
    }
  # compute the AIC for parameter selection
  phi.xi.hat.1<-phi.xi.cond.mean.hat.1
  phi.xi.hat.0<-phi.xi.cond.mean.hat.0
  

  AIC.temp<-AIC.fun(Y.obs.1,Y.obs.0,v.inc.1,v.inc.0,alpha,h0.num.hat,v.obs.base.1,v.obs.base.0,beta.hat,
                    x.obs.1,x.obs.0,x.obs.t.1,x.obs.t.0,vx.obs.1,vx.obs.0,
                    gamma.hat,theta.x.hat,phi.xi.hat.1,phi.xi.hat.0,sigma.hat,
                    n.knots.sele.x,n.order,t.min,t.max,n.phi)
  AIC.res<-c(AIC.res,AIC.temp)
  
  
  # --- incidence --- #
  alpha.fin.temp<-cbind(alpha.fin.temp,alpha.new)
  
  
  # components of the z observation
  gamma.fin.temp<-c(gamma.fin.temp,gamma.new)
  theta.x.fin.temp[[cand]]<-theta.x.new
  eta.fin.temp[[cand]]<-eta.new
  lambda.fin.temp[[cand]]<-diag(lambda.new)
  
  
  # survival function
  beta.fin.temp<-cbind(beta.fin.temp,as.vector(beta.new))
  h0.fin.temp[[cand]]<-cbind(Y.vec,h0.num.new)
  
  # dimension
  n.phi.fin.temp<-c(n.phi.fin.temp,n.phi)
  n.knots.fin.temp<-c(n.knots.fin.temp,n.knots.sele.x)
  
  
  
  }
  
  
  # results from iterations
  
  min.AIC<-which.min(AIC.res)
  
  # results from iterations
  
  # --- incidence --- #
  alpha.fin<-cbind(alpha.fin,alpha.fin.temp[,min.AIC])
  
  
  # components of the z observation
  gamma.fin<-c(gamma.fin,gamma.fin.temp[min.AIC])
  theta.x.fin[[pos]]<-theta.x.fin.temp[[min.AIC]]
  eta.fin[[pos]]<-eta.fin.temp[[min.AIC]]
  lambda.fin[[pos]]<-lambda.fin.temp[[min.AIC]]
  
  
  # survival function
  beta.fin<-cbind(beta.fin,beta.fin.temp[,min.AIC])
  h0.fin[[pos]]<-h0.fin.temp[[min.AIC]]
  
  n.phi.fin<-c(n.phi.fin,n.phi.fin.temp[min.AIC])
  n.knots.fin<-c(n.knots.fin,n.knots.fin.temp[min.AIC])
  
  pos<-pos+1
}



