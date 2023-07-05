# FPCA and Mixture cure model simulation #
# This program is used for the numerical simulation studies of the FPCA and Mixture cure model simulation #
# Loading functions #
rm(list=ls())

source("estimation function joint modelling.R")
source("estimation function binary.R")

##############################################
# Data generation                            #
# Observations:                              #
# Observation time: Y=min(T,C)               #
# Longitudinal observation:Z                 #
# Scalar observation: V                      #
##############################################


# Define frames for the final results
# incidence
alpha.fin<-NULL
# components of the x observation
gamma.fin<-NULL
theta.x.fin<-NULL
phi.fin<-list()
sigma.fin<-NULL
rho.fin<-NULL

# components of the z observation
tau.fin<-NULL
theta.s.fin<-list()
eta.fin<-list()
lambda.fin<-list()

# survival function
beta.fin<-NULL
h0.fin<-list()

n.knots.fin<-NULL
n.phi.fin<-NULL


pos<-1
args<-1
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
  vz.obs<-cbind(v1)
  alpha<-c(0.5,-0.2,1)
  k.prob<-exp(alpha%*%t(v.inc))/(1+exp(alpha%*%t(v.inc)))
  k.true<-rbinom(n,1,k.prob)
  table(k.true)
  
  # for the event case ,we need to generate the functional processed and the survival time #
  # suppose domain is 0:10 #
  gamma<-0.1
  mu.t.s<-function(t){
    rep(0,length(t))
  }
  mu.s.fun<-function(t,gamma,v.obs){
    mu.t.s(t)+as.vector(gamma%*%t(v.obs))
  }
  phi1.fun<-function(t){
    #-cos(pi*t/10)/sqrt(5)
    sin(2*pi*t/10)/sqrt(5)
  }
  phi2.fun<-function(t){
    #sin(pi*t/10)/sqrt(5)
    cos(2*pi*t/10)/sqrt(5)
  }
  lambda1.true<-5
  lambda2.true<-2
  zeta1.true<-rnorm(n,0,sqrt(lambda1.true))
  zeta2.true<-rnorm(n,0,sqrt(lambda2.true))
  
  s.true.fun<-function(t,mu.s.fun,phi1.fun,phi2.fun,zeta1.true,zeta2.true,gamma,v.obs){
    true<-mu.s.fun(t,gamma,v.obs)+zeta1.true*phi1.fun(t)+zeta2.true*phi2.fun(t)
    return(true)
  }
  
  # generate event time #
  # PH model is assumed #
  
  # baseline hazard #
  h0.fun<-function(t){
    t^2/100
  }
  
  beta1.true<-1
  beta.base.true<-c(1,-0.5)
 
  # survival function
  h.fun<-function(t,h0.fun,x.true.fun,beta1.true,beta.base.true,mu.s.fun,phi1.fun,phi2.fun,zeta1.true,zeta2.true,gamma,v.obs,v.obs.base){
    eta.temp<-beta1.true*x.true.fun(t,mu.s.fun,phi1.fun,phi2.fun,zeta1.true,zeta2.true,gamma,v.obs)+as.numeric(beta.base.true%*%v.obs.base)
    h0.fun(t)*exp(eta.temp)
  }
  Fbar.fun<-function(t,h0.fun,x.true.fun,beta1.true,beta.base.true,mu.s.fun,phi1.fun,phi2.fun,zeta1.true,zeta2.true,gamma,v.obs,v.obs.base){
    exp(-stats::integrate(h.fun,0,t,h0.fun=h0.fun,x.true.fun=x.true.fun,beta1.true=beta1.true,beta.base.true=beta.base.true,
                          mu.s.fun=mu.s.fun,phi1.fun=phi1.fun,phi2.fun=phi2.fun,
                          zeta1.true=zeta1.true,zeta2.true=zeta2.true,gamma=gamma,v.obs=v.obs,v.obs.base=v.obs.base)$value)
  }
  Fbar.fun.gen<-function(t,h0.fun,x.true.fun,beta1.true,beta.base.true,mu.s.fun,phi1.fun,phi2.fun,zeta1.true,zeta2.true,gamma,v.obs,v.obs.base,
                         Fbar.true){
    Fbar.fun(t,h0.fun,x.true.fun,beta1.true,beta.base.true,mu.s.fun,phi1.fun,phi2.fun,zeta1.true,zeta2.true,gamma,v.obs,v.obs.base)-Fbar.true
  }

  # true event time
  Fbar.true<-runif(n,0,1)
  t.true<-NULL
  t.true.ind<-NULL
  for(i in 1:n){
    t.true.temp<-try(uniroot(Fbar.fun.gen,interval=c(0,100),
                             h0.fun,s.true.fun,beta1.true,beta.base.true,mu.s.fun,phi1.fun,phi2.fun,zeta1.true[i],zeta2.true[i],
                             gamma,vz.obs[i,],v.obs.base[i,],Fbar.true[i])$root,TRUE)
    if(!inherits(t.true.temp, "try-error")){
      if(t.true.temp>0){
        t.true<-c(t.true,t.true.temp)
        t.true.ind<-c(t.true.ind,i)
      }
    }
  }
  length(t.true)
  summary(t.true)
  
  
  # adjust the observation time due to cured #
  t.true[k.true==0]<-10
  summary(t.true)
  # generate censoring time #
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
  # censoring indicator,  consider cure proportion #
  ind.obs<-ind.obs.true*k.true
  table(ind.obs)
  
  # truncate at 10 #
  ind.obs[t.obs>10]<-0
  t.obs[t.obs>10]<-10

  
  summary(t.obs)
  
  # generate longitudinal binary process observation #
  
  z.obs<-list()
  z.obs.t<-list()
  for(i in 1:n){
    n.sample<-sample(c(6:10),1,prob=rep(1/5,5))
    t.temp<-sort(runif(n.sample,t.min,t.obs[i]))
    #t.obs.temp<-t.temp[t.temp<=t.obs[i]] 
    t.obs.temp<-t.temp
    z.obs.t[[i]]<-t.obs.temp
    s.obs.temp<-s.true.fun(t.obs.temp,mu.s.fun,phi1.fun,phi2.fun,zeta1.true[i],zeta2.true[i],gamma,vz.obs[i,])
    prob.obs.temp<-exp(s.obs.temp)/(1+exp(s.obs.temp))
    z.obs[[i]]<-rbinom(length(t.obs.temp),1, prob.obs.temp)
  }
  
  # observation set #
  Y.obs<-t.obs
  ind.obs<-ind.obs
  z.obs<-z.obs
  z.obs.t<-z.obs.t
  vz.obs<-vz.obs
  
  
  ##############################################
  # Estimation:                                #
  # --- Observations:                      --- #
  # Observation time: Y=min(T,C)               #
  # Indicator delta=I(T<=C)                    #
  # Longitudinal observation: Z                #
  # Scalar observation: V                      #
  # --- Paramters:                         --- #
  # FPCA of S(t)                               #
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
  z.obs.1<-z.obs[ind.obs==1]
  z.obs.t.1<-z.obs.t[ind.obs==1]
  v.obs.base.1<-as.matrix(v.obs.base[ind.obs==1,])
  
  
  Y.obs.0<-Y.obs[ind.obs==0]
  ind.obs.0<-ind.obs[ind.obs==0]
  v.obs.0<-as.matrix(v.obs[ind.obs==0,])
  v.inc.0<-v.inc[ind.obs==0,]
  z.obs.0<-z.obs[ind.obs==0]
  z.obs.t.0<-z.obs.t[ind.obs==0]
  v.obs.base.0<-as.matrix(v.obs.base[ind.obs==0,])
  
  
  n1<-length(ind.obs.1)
  n0<-length(ind.obs.0)
  n.phi.can<-c(1:2)
  n.knots.can<-rep(8:9)
  par.can<-expand.grid(n.phi.can,n.knots.can)
  
  
  alpha.fin.temp<-NULL

  tau.fin.temp<-NULL
  theta.s.fin.temp<-list()
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
    n.knots.sele.z<-par.can[cand,2]
    n.order<-4     
    P<-0.9
    
    t.min<-0
    t.max<-10
    # alpha in the incidence
    alpha.hat<-rep(0,ncol(v.inc))
    
    # FPCA  
    vz.obs.1<-as.matrix(vz.obs[ind.obs==1,])
    vz.obs.0<-as.matrix(vz.obs[ind.obs==0,])
    
    tau1.sele.z<-2
    tau2.sele.z<-2
 
    
    
    # components of the z observation
    
    each.time<-t.max
    #breaks<-    1+  (0:(n.knots-1))*((each.time-1)/(n.knots-1))
    breaks.z<-seq(t.min,t.max,length.out = n.knots.sele.z)
    
    # degree = order-1
    n.basis.z<-length(breaks.z)+n.order-2
    
    ##### Threshold of selecting principal component
    
    # generate b-spline basis
    full.knots.z<- expand.knots(breaks.z, order=n.order) 
    ob.z<-OrthogonalSplineBasis(full.knots.z,order=n.order)
    X.i.standard.z<-list()
    for(i in 1:n){
      X.i.standard.z[[i]]<-evaluate(ob.z,z.obs.t[[i]]) 
    }
    z.i.standard.1<-X.i.standard.z[ind.obs==1]
    z.i.standard.0<-X.i.standard.z[ind.obs==0]
    # D.temp is penalty matrix
    D.temp.z<-OuterProdSecondDerivative(ob.z)
    
    
    estimate.z.res<-estimation.bin(z.obs,X.i.standard.z,tau1.sele.z,tau2.sele.z,D.temp.z,n.basis.z,vz.obs)
    Z1.lambda.fin<-estimate.z.res$lambda
    Z1.theta.fin<-estimate.z.res$theta
    Z1.beta.fin<-estimate.z.res$beta
    Z1.psi.fin<-estimate.z.res$psi.fin
    
    t.grid.standard<-seq(t.min,t.max,0.005)
    nu.standard.z<-evaluate(ob.z,t.grid.standard)
    pc1.Z1.temp<-nu.standard.z%*%Z1.theta.fin[,1]
    if(pc1.Z1.temp[which.max(abs(pc1.Z1.temp[1:(length(pc1.Z1.temp)/2)]))]>0){
      Z1.theta.fin[,1]<--Z1.theta.fin[,1]
    }
    if(ncol(Z1.lambda.fin)==2){ 
      pc2.Z1.temp<-nu.standard.z%*%Z1.theta.fin[,2]
      if(pc2.Z1.temp[(length(pc2.Z1.temp)/4):(length(pc2.Z1.temp)*3/4)][which.max(abs(pc2.Z1.temp[(length(pc2.Z1.temp)/4):(length(pc2.Z1.temp)*3/4)]))]<0){
        Z1.theta.fin[,2]<--Z1.theta.fin[,2]
      }
    }
    zeta.hat<-xi.estimation.bin(log.likelihood.bin,n,z.obs,Z1.beta.fin,Z1.theta.fin,
                                X.i.standard.z,Z1.lambda.fin,vz.obs)
    
    
    tau.hat<-Z1.beta.fin[(n.basis.z+1):length(Z1.beta.fin)]
    theta.s.hat<-Z1.beta.fin[1:n.basis.z]
    eta.hat<-Z1.theta.fin
    lambda.hat<-Z1.lambda.fin
    Psi.hat<-Z1.psi.fin
    
    
    
    
    # update zstar
    z.star.hat<-estimate.z.res$W.star.fin
    
    s.hat<-estimate.z.res$eta.fin
    Psi.hat<-estimate.z.res$psi.fin

    
    z.star.hat.0<-z.star.hat[ind.obs==0]

    s.hat.0<-s.hat[ind.obs==0]
    
    z.star.hat.1<-z.star.hat[ind.obs==1]

    s.hat.1<-s.hat[ind.obs==1]
    
    var.zstar.hat.1<-list()
    for(i in 1:n1){
      var.zstar.hat.1[[i]]<-diag(as.vector(1/deriv.inv.logit.fun(s.hat.1[[i]])))
    }  
    var.zstar.hat.0<-list()
    for(i in 1:n0){
      var.zstar.hat.0[[i]]<-diag(as.vector(1/deriv.inv.logit.fun(s.hat.0[[i]])))
    }
    zeta.eta.hat.1<-t(estimate.z.res$gamma.fin[,ind.obs==1])
    zeta.eta.hat.0<-t(estimate.z.res$gamma.fin[,ind.obs==0])
    
    
    Psi.hat<-eta.hat%*%lambda.hat%*%t(eta.hat)
    
    update.ini.1<-multiple.update(n1,s.hat.1,z.i.standard.1,Psi.hat)
    R.hat.1<-update.ini.1$R.hat
    S.hat.1<-update.ini.1$S.hat
    
    update.ini.0<-multiple.update(n0,s.hat.0,z.i.standard.0,Psi.hat)
    R.hat.0<-update.ini.0$R.hat
    S.hat.0<-update.ini.0$S.hat
    # survival function
    beta.hat<-c(0,rep(0,ncol(v.obs.base)))
    Y.vec<-sort(unique(Y.obs[ind.obs==1]))
    Y.vec<-c(0,Y.vec)
    h0.num.hat<-h0.new.dep.fun(Y.vec,Y.obs,v.obs.base,tau.hat,vz.obs,beta.hat,
                               t(zeta.hat),theta.s.hat,eta.hat,
                               b.basis,ob.z)$h0.num
    # generate samples for monte carlo integration
    q.sample<-50
    
    
    eta.zeta.cond.dist.hat.1<-phi.xi.cond.dist(z.star.hat.1,tau.hat,vz.obs.1,theta.s.hat,
                                               eta.hat,lambda.hat,z.i.standard.1,var.zstar.hat.1,Psi.hat)
    eta.zeta.cond.mean.hat.1<-eta.zeta.cond.dist.hat.1$phi.xi.cond.mean
    eta.zeta.cond.var.hat.1<-eta.zeta.cond.dist.hat.1$phi.xi.cond.var
    
    eta.zeta.sample.1<-phi.xi.gen(n1,q.sample,eta.zeta.cond.mean.hat.1,eta.zeta.cond.var.hat.1)
    
    eta.zeta.cond.dist.hat.0<-phi.xi.cond.dist(z.star.hat.0,tau.hat,vz.obs.0,theta.s.hat,
                                               eta.hat,lambda.hat,z.i.standard.0,var.zstar.hat.0,Psi.hat)
    eta.zeta.cond.mean.hat.0<-eta.zeta.cond.dist.hat.0$phi.xi.cond.mean
    eta.zeta.cond.var.hat.0<-eta.zeta.cond.dist.hat.0$phi.xi.cond.var
    
    eta.zeta.sample.0<-phi.xi.gen(n0,q.sample,eta.zeta.cond.mean.hat.0,eta.zeta.cond.var.hat.0)
    
    ###############################
    #--- iteration starts here ---#
    ###############################
    iteration<-1
    like.hat<-99999
    while(iteration<=200){
      par.hat1<-c(alpha.hat,tau.hat,beta.hat)
      par.hat2<-c(theta.s.hat,lambda.hat,eta.hat)
      
      q.sample<-50+floor(iteration/20)*5

      
      
      # common components received from the current estimates
      
      out_obs1<-f_obs_delta1(eta.zeta.sample.1,
                             Y.obs.1,vz.obs.1,v.obs.base.1,
                             h0.num.hat,Y.vec,beta.hat,
                             tau.hat,theta.s.hat,z.star.hat.1,
                             b.basis,ob.z)
      out_obs1.1<-t(out_obs1)/q.sample
      
      out_obs0<-f_obs0_int(eta.zeta.sample.0,
                           Y.obs.0,vz.obs.0,v.obs.base.0,
                           h0.num.hat,Y.vec,beta.hat,
                           tau.hat,theta.s.hat,z.star.hat.0,
                           b.basis,ob.z)
      
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
      head(f_obs0_k0_fin)
      summary(f_obs0_k0_fin)
      summary(out_deno_obs1)
      
      # ---- component 1: incident rate ---- #
     
      par.ini<-alpha.hat
      alpha.new<-stats::optim(par.ini,exp_inc,v.inc.1=v.inc.1,v.inc.0=v.inc.0,out_obs1.1=out_obs1.1,
                              out_deno_obs1=out_deno_obs1,out_k1_obs0.0=out_k1_obs0.0,
                              out_k0_obs0.0=out_k0_obs0.0,f_obs0_k1_fin=f_obs0_k1_fin)$par
      
      # ---- component 2: mean functions of X and Z ---- #
      
      
      # mean function of S(t)
      
      sigma.zstar.new.1<-list()
      for(i in 1:n1){
        sigma.zstar.new.1[[i]]<-var.zstar.hat.1[[i]]
      }
      sigma.zstar.new.0<-list()
      for(i in 1:n0){
        sigma.zstar.new.0[[i]]<-var.zstar.hat.0[[i]]
      }
      
      # mean function coefficients theta and gamma #
      left<-matrix(0,nrow=n.basis.z+ncol(vz.obs),ncol=n.basis.z+ncol(vz.obs))
      right<-matrix(0,nrow=n.basis.z+ncol(vz.obs),ncol=1)
      for(i in 1:n1){

        Ax.temp<-cbind(z.i.standard.1[[i]],vz.obs.1[i,])
        W.temp<-z.i.standard.1[[i]]%*%Psi.hat%*%t(z.i.standard.1[[i]])+sigma.zstar.new.1[[i]]
        temp1<-t(Ax.temp)%*%solve(W.temp)%*%Ax.temp
        temp2<-t(Ax.temp)%*%solve(W.temp)%*%t(z.star.hat.1[[i]])
        left<-left+temp1
        right<-right+temp2
      }
      for(i in 1:n0){

        Ax.temp<-cbind(z.i.standard.0[[i]],vz.obs.0[i,])
        W.temp<-z.i.standard.0[[i]]%*%Psi.hat%*%t(z.i.standard.0[[i]])+sigma.zstar.new.0[[i]]
        temp1<-t(Ax.temp)%*%solve(W.temp)%*%Ax.temp
        temp2<-t(Ax.temp)%*%solve(W.temp)%*%t(z.star.hat.0[[i]])
        left<-left+temp1
        right<-right+temp2
      }
      theta_tau_new<-solve(left)%*%right
      theta.s.new<-theta_tau_new[1:n.basis.z]
      tau.new<-theta_tau_new[(n.basis.z+1):length(theta_tau_new)]
      
      
      ##############################################################
      # update zstar and conditional distributions 
      ###############################################################
      
      s.hat.1<-s.hat.fun.1(n1,Psi.hat,theta.s.new,tau.new,vz.obs.1,z.i.standard.1,var.zstar.hat.1,
                           z.star.hat.1,zeta.eta.hat.1)
      s.hat.0<-s.hat.fun.1(n0,Psi.hat,theta.s.new,tau.new,vz.obs.0,z.i.standard.0,var.zstar.hat.0,
                           z.star.hat.0,zeta.eta.hat.0)
      zeta.eta.new<-matrix(0,n,n.basis.z)
      
      zstar.update<-zstar.update.new(n1,n0,z.obs.1,z.obs.0,vz.obs.1,vz.obs.0,theta.s.new,tau.new,
                                     z.i.standard.1,z.i.standard.0,Psi.hat,s.hat.1,s.hat.0,
                                     var.zstar.hat.1,R.hat.1,S.hat.1,
                                     var.zstar.hat.0,R.hat.0,S.hat.0,n.basis.z)
      
      z.star.new.1<-zstar.update$z.star.new.1
      z.star.new.0<-zstar.update$z.star.new.0
      
      update.hat.1<-multiple.update(n1,s.hat.1,z.i.standard.1,Psi.hat)
      var.zstar.new.1<-update.hat.1$var.zstar.hat
      R.hat.1<-update.hat.1$R.hat
      S.hat.1<-update.hat.1$S.hat
      
      update.hat.0<-multiple.update(n0,s.hat.0,z.i.standard.0,Psi.hat)
      var.zstar.new.0<-update.hat.0$var.zstar.hat
      R.hat.0<-update.hat.0$R.hat
      S.hat.0<-update.hat.0$S.hat
      
      zeta.eta.new.1<-matrix(0,n1,n.basis.z)
      for(i in 1:n1){
        zeta.eta.new.1[i,]<-Psi.hat%*%t(z.i.standard.1[[i]])%*%solve(R.hat.1[[i]])%*%
          (t(z.star.new.1[[i]])-z.i.standard.1[[i]]%*%theta.s.new-as.numeric(tau.new%*%vz.obs.1[i,]))
      }
      zeta.eta.new.0<-matrix(0,n0,n.basis.z)
      for(i in 1:n0){
        zeta.eta.new.0[i,]<-Psi.hat%*%t(z.i.standard.0[[i]])%*%solve(R.hat.0[[i]])%*%
          (t(z.star.new.0[[i]])-z.i.standard.0[[i]]%*%theta.s.new-as.numeric(tau.new%*%vz.obs.0[i,]))
      }
      zeta.eta.hat.1<-zeta.eta.new.1
      zeta.eta.hat.0<-zeta.eta.new.0
      
      
      z.star.hat.1<-z.star.new.1
      z.star.hat.0<-z.star.new.0
      var.zstar.hat.1<-var.zstar.new.1
      var.zstar.hat.0<-var.zstar.new.0
      
      
      
      eta.zeta.cond.dist.hat.1<-phi.xi.cond.dist(z.star.hat.1,tau.new,vz.obs.1,theta.s.new,
                                                 eta.hat,lambda.hat,z.i.standard.1,var.zstar.hat.1,Psi.hat)
      eta.zeta.cond.mean.hat.1<-eta.zeta.cond.dist.hat.1$phi.xi.cond.mean
      eta.zeta.cond.var.hat.1<-eta.zeta.cond.dist.hat.1$phi.xi.cond.var
      
      eta.zeta.sample.1<-phi.xi.gen(n1,q.sample,eta.zeta.cond.mean.hat.1,eta.zeta.cond.var.hat.1)
      
      eta.zeta.cond.dist.hat.0<-phi.xi.cond.dist(z.star.hat.0,tau.new,vz.obs.0,theta.s.new,
                                                 eta.hat,lambda.hat,z.i.standard.0,var.zstar.hat.0,Psi.hat)
      eta.zeta.cond.mean.hat.0<-eta.zeta.cond.dist.hat.0$phi.xi.cond.mean
      eta.zeta.cond.var.hat.0<-eta.zeta.cond.dist.hat.0$phi.xi.cond.var
      
      eta.zeta.sample.0<-phi.xi.gen(n0,q.sample,eta.zeta.cond.mean.hat.0,eta.zeta.cond.var.hat.0)
      
      out_obs1<-f_obs_delta1(eta.zeta.sample.1,
                             Y.obs.1,vz.obs.1,v.obs.base.1,
                             h0.num.hat,Y.vec,beta.hat,
                             tau.new,theta.s.new,var.zstar.hat.1,
                             b.basis,ob.z)
      out_obs1.1<-t(out_obs1)/q.sample
      
      out_obs0<-f_obs0_int(eta.zeta.sample.0,
                           Y.obs.0,vz.obs.0,v.obs.base.0,
                           h0.num.hat,Y.vec,beta.hat,
                           tau.new,theta.s.new,var.zstar.hat.0,
                           b.basis,ob.z)
      
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
      
      
      # --- components of binary observation Z --- #
      # we consider to use quasi likelihood method 
      # under the current estimate, we estimate the approximated observation
      # use the approximated z.star.new and s.hat 
      # variance is 1/H'(s.hat)

      
      Psi.temp<-matrix(0,n.basis.z,n.basis.z)
      gamma.zeta.sample.1<-list()
      gamma.zeta.sample.0<-list()
      for(i in 1:n1){
        gamma.zeta.sample.1[[i]]<-t(eta.zeta.sample.1[[i]])
        E.zeta.temp1<-(gamma.zeta.sample.1[[i]])%*%(t(gamma.zeta.sample.1[[i]])*out_obs1.1[i,]/out_deno_obs1[i])
        Psi.temp<-Psi.temp+E.zeta.temp1
        
      }
      for(i in 1:n0){
        gamma.zeta.sample.0[[i]]<-t(eta.zeta.sample.0[[i]])
        E.zeta.temp1<-(gamma.zeta.sample.0[[i]])%*%(t(gamma.zeta.sample.0[[i]])*
                                                      (out_k1_obs0.0 [i,]*f_obs0_k1_fin[i]/out_deno_obs0_k1[i]+
                                                         out_k0_obs0.0*f_obs0_k0_fin[i]))
        Psi.temp<-Psi.temp+E.zeta.temp1
      }
      
      Psi.eigen.temp<-eigen(Psi.temp/n)
      lambda.new<-diag(Psi.eigen.temp$values[1:n.phi],ncol = n.phi,nrow=n.phi)
      eta.new<-as.matrix(Psi.eigen.temp$vectors[,1:n.phi])
      
      t.grid.standard<-seq(t.min,t.max,0.005)
      nu.standard.z<-evaluate(ob.z,t.grid.standard)
      
      
      pc1.Z1.temp<-nu.standard.z%*%eta.new[,1]
      if(pc1.Z1.temp[which.max(abs(pc1.Z1.temp[1:(length(pc1.Z1.temp)/2)]))]>0){
        Z1.theta.fin[,1]<--Z1.theta.fin[,1]
      }
      if(ncol(Z1.lambda.fin)==2){ 
        pc2.Z1.temp<-nu.standard.z%*%Z1.theta.fin[,2]
        if(pc2.Z1.temp[(length(pc2.Z1.temp)/4):(length(pc2.Z1.temp)*3/4)][which.max(abs(pc2.Z1.temp[(length(pc2.Z1.temp)/4):(length(pc2.Z1.temp)*3/4)]))]<0){
          Z1.theta.fin[,2]<--Z1.theta.fin[,2]
        }
      }
      
      
      
      Psi.new<-eta.new%*%lambda.new%*%t(eta.new)


      
      #############################################
      # update zstar and others another time      #
      #############################################
      #############################################
      # update zstar and conditional distributions#
      #############################################
      s.hat.1<-s.hat.fun.1(n1,Psi.new,theta.s.new,tau.new,vz.obs.1,z.i.standard.1,var.zstar.hat.1,
                           z.star.hat.1,zeta.eta.hat.1)
      s.hat.0<-s.hat.fun.1(n0,Psi.new,theta.s.new,tau.new,vz.obs.0,z.i.standard.0,var.zstar.hat.0,
                           z.star.hat.0,zeta.eta.hat.0)
      zstar.update<-zstar.update.new(n1,n0,z.obs.1,z.obs.0,vz.obs.1,vz.obs.0,theta.s.new,tau.new,
                                     z.i.standard.1,z.i.standard.0,Psi.new,s.hat.1,s.hat.0,
                                     var.zstar.hat.1,R.hat.1,S.hat.1,
                                     var.zstar.hat.0,R.hat.0,S.hat.0,n.basis.z)
      
      z.star.new.1<-zstar.update$z.star.new.1
      z.star.new.0<-zstar.update$z.star.new.0
      
      update.hat.1<-multiple.update(n1,s.hat.1,z.i.standard.1,Psi.new)
      var.zstar.new.1<-update.hat.1$var.zstar.hat
      R.hat.1<-update.hat.1$R.hat
      S.hat.1<-update.hat.1$S.hat
      
      update.hat.0<-multiple.update(n0,s.hat.0,z.i.standard.0,Psi.new)
      var.zstar.new.0<-update.hat.0$var.zstar.hat
      R.hat.0<-update.hat.0$R.hat
      S.hat.0<-update.hat.0$S.hat
      
      
      zeta.eta.new.1<-matrix(0,n1,n.basis.z)
      for(i in 1:n1){
        zeta.eta.new.1[i,]<-Psi.new%*%t(z.i.standard.1[[i]])%*%solve(R.hat.1[[i]])%*%
          (t(z.star.new.1[[i]])-z.i.standard.1[[i]]%*%theta.s.new-as.numeric(tau.new%*%vz.obs.1[i,]))
      }
      zeta.eta.new.0<-matrix(0,n0,n.basis.z)
      for(i in 1:n0){
        zeta.eta.new.0[i,]<-Psi.new%*%t(z.i.standard.0[[i]])%*%solve(R.hat.0[[i]])%*%
          (t(z.star.new.0[[i]])-z.i.standard.0[[i]]%*%theta.s.new-as.numeric(tau.new%*%vz.obs.0[i,]))
      }
      zeta.eta.hat.1<-zeta.eta.new.1
      zeta.eta.hat.0<-zeta.eta.new.0
      
      z.star.hat.1<-z.star.new.1
      z.star.hat.0<-z.star.new.0
      var.zstar.hat.1<-var.zstar.new.1
      var.zstar.hat.0<-var.zstar.new.0
      
      
      eta.zeta.cond.dist.hat.1<-phi.xi.cond.dist(z.star.hat.1,tau.new,vz.obs.1,theta.s.new,
                                                 eta.new,lambda.new,z.i.standard.1,var.zstar.hat.1,Psi.new)
      eta.zeta.cond.mean.hat.1<-eta.zeta.cond.dist.hat.1$phi.xi.cond.mean
      eta.zeta.cond.var.hat.1<-eta.zeta.cond.dist.hat.1$phi.xi.cond.var
      
      eta.zeta.sample.1<-phi.xi.gen(n1,q.sample,eta.zeta.cond.mean.hat.1,eta.zeta.cond.var.hat.1)
      
      eta.zeta.cond.dist.hat.0<-phi.xi.cond.dist(z.star.hat.0,tau.new,vz.obs.0,theta.s.new,
                                                 eta.new,lambda.new,z.i.standard.0,var.zstar.hat.0,Psi.new)
      eta.zeta.cond.mean.hat.0<-eta.zeta.cond.dist.hat.0$phi.xi.cond.mean
      eta.zeta.cond.var.hat.0<-eta.zeta.cond.dist.hat.0$phi.xi.cond.var
      
      eta.zeta.sample.0<-phi.xi.gen(n0,q.sample,eta.zeta.cond.mean.hat.0,eta.zeta.cond.var.hat.0)
      
      zeta.eta.hat.1<-t(eta.zeta.cond.mean.hat.1)
      zeta.eta.hat.0<-t(eta.zeta.cond.mean.hat.0)
      
      out_obs1<-f_obs_delta1(eta.zeta.sample.1,
                             Y.obs.1,vz.obs.1,v.obs.base.1,
                             h0.num.hat,Y.vec,beta.hat,
                             tau.new,theta.s.new,var.zstar.hat.1,
                             b.basis,ob.z)
      out_obs1.1<-t(out_obs1)/q.sample
      
      out_obs0<-f_obs0_int(eta.zeta.sample.0,
                           Y.obs.0,vz.obs.0,v.obs.base.0,
                           h0.num.hat,Y.vec,beta.hat,
                           tau.new,theta.s.new,var.zstar.hat.0,
                           b.basis,ob.z)
      
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
      
      basis.Y.vec.z<-b.basis(ob.z,t=Y.vec)
      
      beta.new<-rootSolve::multiroot(first.deriv.fun,start=beta.hat,
                                     Y.vec=Y.vec,v.obs.base.1=v.obs.base.1,v.obs.base.0=v.obs.base.0,
                                     phi.xi.sample.1=eta.zeta.sample.1,phi.xi.sample.0=eta.zeta.sample.0,
                                     Y.obs.1=Y.obs.1,x.obs.1=z.star.hat.1,x.obs.t.1=z.obs.t.1,
                                     vx.obs.1=vz.obs.1,ind.obs.1=ind.obs.1,n1=n1,
                                     Y.obs.0=Y.obs.0,x.obs.0=z.star.hat.0,x.obs.t.0=z.obs.t.0,
                                     vx.obs.0=vz.obs.0,ind.obs.0=ind.obs.0,n0=n0,
                                     gamma.hat=tau.new,theta.x.hat=theta.s.new,sigma.hat=NULL,
                                     out_obs1.1=out_obs1.1,out_deno_obs1=out_deno_obs1,
                                     out_k1_obs0.0=out_k1_obs0.0,out_k0_obs0.0=out_k0_obs0.0,
                                     f_obs0_k1_fin= f_obs0_k1_fin,f_obs0_k0_fin=f_obs0_k0_fin,
                                     out_deno_obs0_k1=out_deno_obs0_k1,out_deno_obs0_k0=out_deno_obs0_k0,
                                     x.i.standard.1=z.i.standard.1,x.i.standard.0=z.i.standard.0,
                                     basis.Y.vec.x=basis.Y.vec.z)$root
      
      
      
      h0.num.new<-h0.new.fun(beta.new,
                             Y.vec=Y.vec,v.obs.base.1=v.obs.base.1,v.obs.base.0=v.obs.base.0,
                             phi.xi.sample.1=eta.zeta.sample.1,phi.xi.sample.0=eta.zeta.sample.0,
                             Y.obs.1=Y.obs.1,x.obs.1=z.star.hat.1,x.obs.t.1=z.obs.t.1,
                             vx.obs.1=vz.obs.1,ind.obs.1=ind.obs.1,n1=n1,
                             Y.obs.0=Y.obs.0,x.obs.0=z.star.hat.0,x.obs.t.0=z.obs.t.0,
                             vx.obs.0=vz.obs.0,ind.obs.0=ind.obs.0,n0=n0,
                             gamma.hat=tau.new,theta.x.hat=theta.s.new,sigma.hat=NULL,
                             out_obs1.1=out_obs1.1,out_deno_obs1=out_deno_obs1,
                             out_k1_obs0.0=out_k1_obs0.0,out_k0_obs0.0=out_k0_obs0.0,
                             f_obs0_k1_fin= f_obs0_k1_fin,f_obs0_k0_fin=f_obs0_k0_fin,
                             out_deno_obs0_k1=out_deno_obs0_k1,out_deno_obs0_k0=out_deno_obs0_k0,
                             x.i.standard.1=z.i.standard.1,x.i.standard.0=z.i.standard.0,
                             basis.Y.vec.x=basis.Y.vec.z)
      
      eta.zeta.hat.1<-eta.zeta.cond.mean.hat.1
      eta.zeta.hat.0<-eta.zeta.cond.mean.hat.0
      
      vz.obs.1<-vz.obs.1
      vz.obs.0<-vz.obs.0
      
      like.new<-obs.like.fun(Y.obs.1,Y.obs.0,v.inc.1,v.inc.0,alpha,h0.num.hat,v.obs.base.1,v.obs.base.0,beta.hat,
                             z.obs.1,z.obs.0,z.obs.t.1,z.obs.t.0,vz.obs.1,vz.obs.0,
                             tau.hat,theta.s.hat,eta.zeta.hat.1,eta.zeta.hat.0,sigma.hat,
                             n.knots.sele.z,n.order,t.min,t.max)
      
      # compare difference between new and hat outcomes ##
      
      par.new1<-c(alpha.new,tau.new,beta.new)
      par.new2<-c(theta.s.new,lambda.new,eta.new)
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
      tau.hat<-tau.new
      theta.s.hat<-theta.s.new
      eta.hat<-eta.new
      lambda.hat<-lambda.new
      
      like.hat<-like.new
      Psi.hat<-Psi.new
      
      
      # survival function
      
      
      beta.hat<-as.vector(beta.new)
      h0.num.hat<-h0.num.new
      
      iteration<-iteration+1
      
      
    }
    
    # compute the AIC for parameter selection
    eta.zeta.hat.1<-eta.zeta.cond.mean.hat.1
    eta.zeta.hat.0<-eta.zeta.cond.mean.hat.0
    
    vz.obs.1<-vz.obs.1
    vz.obs.0<-vz.obs.0
    
    AIC.temp<-AIC.fun(Y.obs.1,Y.obs.0,v.inc.1,v.inc.0,alpha,h0.num.hat,v.obs.base.1,v.obs.base.0,beta.hat,
                      z.obs.1,z.obs.0,z.obs.t.1,z.obs.t.0,vz.obs.1,vz.obs.0,
                      tau.hat,theta.s.hat,eta.zeta.hat.1,eta.zeta.hat.0,sigma.hat,
                      n.knots.sele.z,n.order,t.min,t.max,n.phi)
    
    AIC.res<-c(AIC.res,AIC.temp)
    
    
    # --- incidence --- #
    alpha.fin.temp<-cbind(alpha.fin.temp,alpha.new)
    
    
    # components of the z observation
    tau.fin.temp<-c(tau.fin.temp,tau.new)
    theta.s.fin.temp[[cand]]<-theta.s.new
    eta.fin.temp[[cand]]<-eta.new
    lambda.fin.temp[[cand]]<-diag(lambda.new)
    
    
    # survival function
    beta.fin.temp<-cbind(beta.fin.temp,as.vector(beta.new))
    h0.fin.temp[[cand]]<-cbind(Y.vec,h0.num.new)
    
    # dimension
    n.phi.fin.temp<-c(n.phi.fin.temp,n.phi)
    n.knots.fin.temp<-c(n.knots.fin.temp,n.knots.sele.z)
    
    
  }
  
  # select the results with the minimum AIC
  min.AIC<-which.min(AIC.res)
  
  # results from iterations
  
  # --- incidence --- #
  alpha.fin<-cbind(alpha.fin,alpha.fin.temp[,min.AIC])
  
  
  # components of the z observation
  tau.fin<-c(tau.fin,tau.fin.temp[min.AIC])
  theta.s.fin[[pos]]<-theta.s.fin.temp[[min.AIC]]
  eta.fin[[pos]]<-eta.fin.temp[[min.AIC]]
  lambda.fin[[pos]]<-lambda.fin.temp[[min.AIC]]
  
  
  # survival function
  beta.fin<-cbind(beta.fin,beta.fin.temp[,min.AIC])
  h0.fin[[pos]]<-h0.fin.temp[[min.AIC]]
  
  n.phi.fin<-c(n.phi.fin,n.phi.fin.temp[min.AIC])
  n.knots.fin<-c(n.knots.fin,n.knots.fin.temp[min.AIC])
  
  pos<-pos+1
}

