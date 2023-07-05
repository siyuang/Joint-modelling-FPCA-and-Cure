# --- This program is used for marginal estimation of the parameters of the FPCA with continuous observations --- #



library(orthogonalsplinebasis)
library(mvtnorm)
library(MASS)
library(boot)
library(Matrix)
library(lme4)
library(statmod)
library(nlme)
library(pROC)

# --- Initial update of the parameters --- #


multi.update.ini.cont<-function(n,Y.obs,beta.current=beta.ini,X.i.standard,gamma.current=gamma.ini,
                                psi.current=psi.ini,var.err.ini=var.err.ini,penalty1,penalty2,n.basis,v.obs){
  eta.current<-list()
  sigma.current<-list()
  R.current<-list()
  S.current<-list()
  temp1.current<-matrix(0,n.basis,n.basis)
  for(i in 1:n){
    eta.current[[i]]<-(beta.current)%*%t(cbind(X.i.standard[[i]],v.obs[i,]))+t(gamma.current[,i])%*%t(X.i.standard[[i]])
    sigma.current[[i]]<-diag(var.err.ini,length(Y.obs[[i]]))
    R.current[[i]]<-sigma.current[[i]]+X.i.standard[[i]]%*%psi.current%*%t(X.i.standard[[i]])
    S.current[[i]]<-psi.current-psi.current%*%t(X.i.standard[[i]])%*%solve(R.current[[i]])%*%X.i.standard[[i]]%*%psi.current
    temp1.current<-temp1.current+t(X.i.standard[[i]])%*%solve(R.current[[i]])%*%X.i.standard[[i]]
  }
  A.current<-solve(temp1.current+penalty1)%*%temp1.current%*%solve(temp1.current+penalty1)
  var.gamma.current<-list()
  Y.obs<-Y.obs
  for(i in 1:n){
    var.gamma.current[[i]]<-S.current[[i]]+psi.current%*%t(X.i.standard[[i]])%*%solve(R.current[[i]])%*%X.i.standard[[i]]%*%
      A.current%*%t(X.i.standard[[i]])%*%solve(R.current[[i]])%*%X.i.standard[[i]]%*%psi.current
    
  }
  return(list("eta.current"=eta.current,"Sigma.current"=sigma.current,"R.current"=R.current,
              "S.current"=S.current,"temp1.current"=temp1.current,"var.gamma.current"=var.gamma.current,
              "Y.obs"=Y.obs,"A.current"=A.current))
}


# --- Update of the parameters within iterations --- #


multi.update.cont<-function(n,Y.obs,beta.current,X.i.standard,gamma.current,
                            psi.current,var.err.current,A.current,penalty1,penalty2,n.basis,v.obs){
  eta.current<-list()
  sigma.current<-list()
  R.current<-list()
  S.current<-list()
  temp1.current<-matrix(0,n.basis,n.basis)
  for(i in 1:n){
    eta.current[[i]]<-(beta.current)%*%t(cbind(X.i.standard[[i]],v.obs[i,]))+t(gamma.current[,i])%*%t(X.i.standard[[i]])
    sigma.current[[i]]<-diag(as.numeric(var.err.current),length(Y.obs[[i]]))
    R.current[[i]]<-sigma.current[[i]]+X.i.standard[[i]]%*%psi.current%*%t(X.i.standard[[i]])
    S.current[[i]]<-psi.current-psi.current%*%t(X.i.standard[[i]])%*%solve(R.current[[i]])%*%X.i.standard[[i]]%*%psi.current
    temp1.current<-temp1.current+t(X.i.standard[[i]])%*%solve(R.current[[i]])%*%X.i.standard[[i]]
  }
  A.new<-solve(temp1.current+penalty1)%*%temp1.current%*%solve(temp1.current+penalty1)
  var.gamma.current<-list()
  Y.obs<-Y.obs
  for(i in 1:n){
    var.gamma.current[[i]]<-S.current[[i]]+psi.current%*%t(X.i.standard[[i]])%*%solve(R.current[[i]])%*%X.i.standard[[i]]%*%
      A.current%*%t(X.i.standard[[i]])%*%solve(R.current[[i]])%*%X.i.standard[[i]]%*%psi.current
    
  }
  return(list("eta.current"=eta.current,"Sigma.current"=sigma.current,"R.current"=R.current,
              "S.current"=S.current,"temp1.current"=temp1.current,"var.gamma.current"=var.gamma.current,
              "Y.obs"=Y.obs,"A.new"=A.new))
}

# --- Update the variance of the measurement error --- #

var.err.update.cont<-function(n,Y.obs,beta.current,psi.current,X.i.standard,var.err.current,v.obs){
  var.err.new<-0
  N<-0
  for(i in 1:n){
    R.current.temp<-diag(length(Y.obs[[i]]))+X.i.standard[[i]]%*%psi.current%*%t(X.i.standard[[i]])/as.numeric(var.err.current)
    var.err.temp<-(Y.obs[[i]]-(beta.current)%*%t(cbind(X.i.standard[[i]],v.obs[i,])))%*%solve(R.current.temp)%*%t(Y.obs[[i]]-(beta.current)%*%t(cbind(X.i.standard[[i]],v.obs[i,])))
    var.err.new<- var.err.new+var.err.temp
    N<-N+length(Y.obs[[i]])
  }
  var.err.new<-var.err.new/N
  return( var.err.new)
}

# Func

estimation.cont<-function(Y.obs,X.i.standard,tau1,tau2,D.temp,n.basis,v.obs){
  n.temp<-length(Y.obs)
  penalty1<-tau1*D.temp
  penalty2<-tau2*D.temp
  # --- initial values --- #
  beta.ini<-rep(0,n.basis+ncol(v.obs))
  beta.ini.matrix<-matrix(rep(beta.ini,n.temp),n.temp,byrow=T)
  psi.ini<-diag(1,n.basis)
  # using psi to calculate initial J, theta,lambda,gamma,eta
  eigen.decom.ini<-eigen(psi.ini)
  eigen.val.ini<-eigen.decom.ini$values
  eigen.val.cumprob.ini<-cumsum(eigen.val.ini)/sum(eigen.val.ini)
  J.ini<-min(which(eigen.val.cumprob.ini>=P))
  theta.ini<-eigen.decom.ini$vectors[1:J.ini,]
  lambda.ini<-diag(eigen.val.ini[1:J.ini],J.ini,J.ini)
  gamma.ini<-matrix(0,n.basis,n.temp)
  var.err.ini<-1
  
  multi.ini<-multi.update.ini.cont(n.temp,Y.obs,beta.ini,X.i.standard,gamma.ini,psi.ini,var.err.ini,penalty1,penalty2,n.basis,v.obs)
  eta.ini<-multi.ini$eta.current
  sigma.ini<-multi.ini$Sigma.current
  R.ini<-multi.ini$R.current
  S.ini<-multi.ini$S.current
  temp1<-multi.ini$temp1.current
  A<-multi.ini$A.current
  
  # --- set all initial to current --- #
  
  beta.current<-beta.ini
  lambda.current<-lambda.ini
  psi.current<-psi.ini
  gamma.current<-gamma.ini
  eta.current<-eta.ini
  sigma.current<-sigma.ini
  R.current<-R.ini
  S.current<-S.ini
  var.err.current<-var.err.ini
  temp1.current<-temp1
  A.current<-A
  
  # --- beta.new --- #
  # --- update eta.current1,Y.obs1,S.current1,R.current1,sigma.current1
  left<-matrix(0,nrow=n.basis+ncol(v.obs),ncol=n.basis+ncol(v.obs))
  temp2<-matrix(0,n.basis+ncol(v.obs),1)
  for(i in 1:n.temp){
    left<-left+t(cbind(X.i.standard[[i]],v.obs[i,]))%*%solve(R.current[[i]])%*%cbind(X.i.standard[[i]],v.obs[i,])
    temp2<-temp2+t(cbind(X.i.standard[[i]],v.obs[i,]))%*%solve(R.current[[i]])%*%(Y.obs[[i]])
  }
  beta.new<-as.vector(solve(left)%*%temp2)
  
  multiple.update1<-multi.update.ini.cont(n.temp,Y.obs,beta.new,X.i.standard,gamma.current,psi.current,var.err.current,penalty1,penalty2,n.basis,v.obs)
  eta.current1<-multiple.update1$eta.current
  Sigma.current1<-multiple.update1$sigma.current
  R.current1<-multiple.update1$R.current
  S.current1<-multiple.update1$S.current
  temp.current1<-multiple.update1$temp1.current
  var.gamma.current1<-multiple.update1$var.gamma.current
  
  
  # --- variance of measurement error term ---  #
  var.err.new<-var.err.update.cont(n.temp,Y.obs,beta.new,psi.current,X.i.standard,var.err.current,v.obs)
  
  # --- gamma.new --- #
  
  gamma.new<-matrix(0,n.basis,n.temp)
  for(i in 1:n.temp){
    gamma.new[,i]<-psi.current%*%t(X.i.standard[[i]])%*%solve(R.current1[[i]])%*%((Y.obs[[i]])-cbind(X.i.standard[[i]],v.obs[i,])%*%beta.new)
  }
  
  # --- psi.current2 --- #
  psi.current2.temp<-matrix(0,n.basis,n.basis)
  for(i in 1:n.temp){
    psi.current2.temp<-psi.current2.temp+gamma.new[,i]%*%t(gamma.new[,i])+S.current1[[i]]
  }
  psi.current2<-psi.current2.temp/n.temp
  
  # --- J.new, theta.new, lambda.new --- #
  # --- psi.new --- #
  # --- update eta.new,Y.obs.new,S.new,R.new,sigma.new --- #
  eigen.decom.new<-eigen(psi.current2,symmetric = TRUE)
  eigen.val.new<-eigen.decom.new$values
  eigen.val.cumprob.new<-cumsum(eigen.val.new)/sum(eigen.val.new)
  J.new<-min(which(eigen.val.cumprob.new>=P))
  #J.new<-n.phi
  if(J.new==1){
    theta.new<-matrix(eigen.decom.new$vectors[,1])
  } else {theta.new<-eigen.decom.new$vectors[,1:J.new]}
  lambda.new<-diag(eigen.val.new[1:J.new],J.new,J.new)
  
  
  
  theta.new.star.temp<-matrix(0,n.basis,n.basis)
  for(i in 1:n.temp){
    if(length(Y.obs[[i]])>1){
      tmp<-diag(as.numeric(var.err.current),length(Y.obs[[i]]))
    } else{
      tmp<-as.numeric(var.err.current)
    }
    theta.new.star.temp<-theta.new.star.temp+t(X.i.standard[[i]])%*%solve(tmp)%*%X.i.standard[[i]]
  }
  theta.new.star<-matrix(0,J.new,n.basis)
  for(i in 1:J.new){
    theta.new.star[i,]<-solve(theta.new.star.temp*lambda.new[i,i]+penalty2)%*%
      (theta.new.star.temp*lambda.new[i,i])%*%theta.new[,i]
    
  }
  psi.new<-t(theta.new.star)%*%lambda.new%*%theta.new.star
  
  multi.update.new<-multi.update.cont(n.temp,Y.obs,beta.new,X.i.standard,gamma.new,psi.new,var.err.new,A.current,penalty1,penalty2,n.basis,v.obs)
  eta.new<-multi.update.new$eta.current
  sigma.new<-multi.update.new$Sigma.current
  R.new<-multi.update.new$R.current
  S.new<-multi.update.new$S.current
  temp1.new<-multi.update.new$temp1.current
  var.gamma.new<-multi.update.new$var.gamma.current
  
  A.new<-multi.update.new$A.new
  
  # --- set all new to current --- #
  
  beta.current<-beta.new
  lambda.current<-lambda.new
  theta.current<-theta.new
  psi.current<-psi.new
  gamma.current<-gamma.new
  eta.current<-eta.new
  sigma.current<-sigma.new
  R.current<-R.new
  S.current<-S.new
  
  temp1.current<-temp1.new
  A.current<-A.new
  J.current<-J.new
  var.err.current<-var.err.new
  
  # --- iteration starts here --------------------------------------------------------------------#
  
  j<-1
  threshold<-100
  while (((threshold>0.001)&(j<1000))) {
    
    
    
    # --- beta.new --- #
    # --- update eta.current1,Y.obs1,S.current1,R.current1,sigma.current1
    left<-matrix(0,nrow=n.basis+ncol(v.obs),ncol=n.basis+ncol(v.obs))
    temp2<-matrix(0,n.basis+ncol(v.obs),1)
    for(i in 1:n.temp){
      left<-left+t(cbind(X.i.standard[[i]],v.obs[i,]))%*%solve(R.current[[i]])%*%cbind(X.i.standard[[i]],v.obs[i,])
      temp2<-temp2+t(cbind(X.i.standard[[i]],v.obs[i,]))%*%solve(R.current[[i]])%*%(Y.obs[[i]])
    }
    beta.new<-as.vector(solve(left)%*%temp2)
    
    
    multiple.update1<-multi.update.cont(n.temp,Y.obs,beta.new,X.i.standard,gamma.current,psi.current,var.err.current,A.current,penalty1,penalty2,n.basis,v.obs)
    eta.current1<-multiple.update1$eta.current
    Sigma.current1<-multiple.update1$sigma.current
    R.current1<-multiple.update1$R.current
    S.current1<-multiple.update1$S.current
    temp1.current1<-multiple.update1$temp1.current
    var.gamma.current1<-multiple.update1$var.gamma.current
    
    # --- variance of measurement error term ---  #
    var.err.new<-var.err.update.cont(n.temp,Y.obs,beta.new,psi.current,X.i.standard,var.err.current,v.obs)
    
    
    # --- gamma.new --- #
    gamma.new<-matrix(0,n.basis,n.temp)
    for(i in 1:n.temp){
      gamma.new[,i]<-psi.current%*%t(X.i.standard[[i]])%*%solve(R.current1[[i]])%*%((Y.obs[[i]])-cbind(X.i.standard[[i]],v.obs[i,])%*%beta.new)
    }
    
    # --- psi.current2 --- #
    psi.current2.temp<-matrix(0,n.basis,n.basis)
    for(i in 1:n.temp){
      psi.current2.temp<-psi.current2.temp+gamma.new[,i]%*%t(gamma.new[,i])+S.current1[[i]]
    }
    psi.current2<-psi.current2.temp/n.temp
    
    # --- J.new, theta.new, lambda.new --- #
    # --- psi.new --- #
    # --- update eta.new,Y.obs.new,S.new,R.new,sigma.new --- #
    eigen.decom.new<-eigen(psi.current2,symmetric = TRUE)
    eigen.val.new<-eigen.decom.new$values
    eigen.val.cumprob.new<-cumsum(eigen.val.new)/sum(eigen.val.new)
    J.new<-min(which(eigen.val.cumprob.new>=P))
    #J.new<-n.phi
    if(J.new==1){
      theta.new<-matrix(eigen.decom.new$vectors[,1])
    } else {theta.new<-eigen.decom.new$vectors[,1:J.new]}
    lambda.new<-diag(eigen.val.new[1:J.new],J.new,J.new)
    
    
    
    theta.new.star.temp<-matrix(0,n.basis,n.basis)
    for(i in 1:n.temp){
      if(length(Y.obs[[i]])>1){
        tmp<-diag(as.numeric(var.err.current),length(Y.obs[[i]]))
      }else{
        tmp<-as.numeric(var.err.current)
      }
      theta.new.star.temp<-theta.new.star.temp+t(X.i.standard[[i]])%*%solve(tmp)%*%X.i.standard[[i]]
    }
    theta.new.star<-matrix(0,J.new,n.basis)
    for(i in 1:J.new){
      theta.new.star[i,]<-solve(theta.new.star.temp*lambda.new[i,i]+penalty2)%*%
        (theta.new.star.temp*lambda.new[i,i])%*%theta.new[,i]
      
    }
    psi.new<-t(theta.new.star)%*%lambda.new%*%theta.new.star
    
    multi.update.new<-multi.update.cont(n.temp,Y.obs,beta.new,X.i.standard,gamma.new,psi.new,var.err.new,A.current,penalty1,penalty2,n.basis,v.obs)
    eta.new<-multi.update.new$eta.current
    sigma.new<-multi.update.new$Sigma.current
    R.new<-multi.update.new$R.current
    S.new<-multi.update.new$S.current
    temp1.new<-multi.update.new$temp1.current
    var.gamma.new<-multi.update.new$var.gamma.current
    
    A.new<-multi.update.new$A.new
    
    # --- check threshold --- #
    res.temp.current<-c(beta.current,unlist(theta.current),unlist(lambda.current))
    res.temp.new<-c(beta.new,unlist(theta.new),unlist(lambda.new))
    if( length(as.numeric(res.temp.current) ) != length(as.numeric(res.temp.new) ))  threshold<-99
    else{
      threshold<- max(abs(as.numeric(res.temp.current)-as.numeric(res.temp.new))  )
    }
    
    # --- set all new to current --- #
    
    beta.current<-beta.new
    lambda.current<-lambda.new
    theta.current<-theta.new
    psi.current<-psi.new
    gamma.current<-gamma.new
    eta.current<-eta.new
    sigma.current<-sigma.new
    R.current<-R.new
    S.current<-S.new
    
    temp1.current<-temp1.new
    A.current<-A.new
    J.current<-J.new
    var.err.current<-var.err.new
    
    j<-j+1
    print(j)
    print(threshold)
  }
  return(list("beta.fin"=beta.new,"lambda.fin"=lambda.new,"theta.fin"=theta.new,"J.fin"=J.new,
              "gamma.fin"=gamma.new,"eta.fin"=eta.new,"var.err"=var.err.new,"psi.fin"=psi.new))
}


# marginal likelihood for the observations

marginal.likelihood.inner.cont<-function(beta.fin,X.i.standard,J.fin,Y.obs,theta.fin,lambda.fin,var.err,xi.fin,v.obs){
  n<-length(Y.obs)
  likelihood.full<-NULL
  if(J.fin==1){for(i in 1:n)
  {
    eta.temp<-(beta.fin)%*%t(cbind(X.i.standard[[i]],v.obs[i,]))+xi.fin[i,]%*%t(theta.fin)%*%t(X.i.standard[[i]])
    temp<-prod(as.numeric(1/sqrt((2*pi*var.err)))*exp(-(Y.obs[[i]]-eta.temp)^2/(2*as.numeric(var.err))))
    likelihood.full<-c(likelihood.full,temp)
  }
  } 
  if(J.fin==2){for(i in 1:n)
  {
    eta.temp<-(beta.fin)%*%t(cbind(X.i.standard[[i]],v.obs[i,]))+xi.fin[i,]%*%t(theta.fin)%*%t(X.i.standard[[i]])
    temp<-prod(as.numeric(1/sqrt((2*pi*var.err)))*exp(-(Y.obs[[i]]-eta.temp)^2/(2*as.numeric(var.err))))
    likelihood.full<-c(likelihood.full,temp)
  }
    
  }
  if(J.fin==3){for(i in 1:n)
  {
    eta.temp<-(beta.fin)%*%t(cbind(X.i.standard[[i]],v.obs[i,]))+
      xi.fin[i,]%*%t(theta.fin)%*%t(X.i.standard[[i]])
    temp<-prod(as.numeric(1/sqrt((2*pi*var.err)))*exp(-(Y.obs[[i]]-eta.temp)^2/(2*as.numeric(var.err))))
    likelihood.full<-c(likelihood.full,temp)
  }
  }
  return(likelihood.full)
}




# --- Use MLE to estimate the FPC scores --- #

log.likelihood.cont<-function(par,Y.obs,beta.fin,theta.fin,X.i.standard,lambda.fin,v.obs){
  mu.fin<-cbind(X.i.standard,v.obs)%*%beta.fin
  if(nrow(lambda.fin)==1){
    phi1.fin<-X.i.standard%*%theta.fin[,1]
    gamma.fin1<-par[1]*phi1.fin
    gamma.fin<-gamma.fin1
    lse<-sum((Y.obs-gamma.fin-mu.fin)^2)
    #logl1<-log(prob.fin)*(Y.obs)+log(1-prob.fin)*(1-Y.obs)
    #logl<-sum(logl1)+log(dmvnorm(x=par,mean=rep(0,length(par)),sigma = sqrt(lambda.fin)))
    #logl<-sum(logl1)+log(dnorm(par[1],0,sqrt(lambda.fin[1,1])))
  } else if(nrow(lambda.fin)==2){
    phi1.fin<-X.i.standard%*%theta.fin[,1]
    phi2.fin<-X.i.standard%*%theta.fin[,2]
    gamma.fin1<-par[1]*phi1.fin
    gamma.fin2<-par[2]*phi2.fin
    gamma.fin<-gamma.fin1+gamma.fin2
    lse<-sum((Y.obs-gamma.fin-mu.fin)^2)
    #prob.fin<-inv.logit(mu.fin+gamma.fin)
    #logl1<-log(prob.fin)*(Y.obs)+log(1-prob.fin)*(1-Y.obs)
    #logl<-sum(logl1)+log(dmvnorm(x=par,mean=rep(0,length(par)),sigma = sqrt(lambda.fin)))
    #logl<-sum(logl1)+log(dnorm(par[1],0,sqrt(lambda.fin[1,1])))+log(dnorm(par[2],0,sqrt(lambda.fin[2,2])))
  } else if (nrow(lambda.fin)==3){
    phi1.fin<-X.i.standard%*%theta.fin[,1]
    phi2.fin<-X.i.standard%*%theta.fin[,2]
    phi3.fin<-X.i.standard%*%theta.fin[,3]
    gamma.fin1<-par[1]*phi1.fin
    gamma.fin2<-par[2]*phi2.fin
    gamma.fin3<-par[3]*phi3.fin
    gamma.fin<-gamma.fin1+gamma.fin2+gamma.fin3
    lse<-sum((Y.obs-gamma.fin-mu.fin)^2)
    #prob.fin<-inv.logit(mu.fin+gamma.fin)
    #logl1<-log(prob.fin)*(Y.obs)+log(1-prob.fin)*(1-Y.obs)
    #logl<-sum(logl1)+log(dmvnorm(x=par,mean=rep(0,length(par)),sigma = sqrt(lambda.fin)))
    #logl<-sum(logl1)+log(dnorm(par[1],0,sqrt(lambda.fin[1,1])))+log(dnorm(par[2],0,sqrt(lambda.fin[2,2])))+
    #log(dnorm(par[3],0,sqrt(lambda.fin[3,3])))
  }
  lse
}


xi.estimation.cont<-function(log.likelihood,n,Y.obs,beta.fin,theta.fin,X.i.standard,lambda.fin,v.obs){
  xi.fin<-NULL
  for(i in 1:n){
    if(nrow(lambda.fin)==1){
      test.optim<-optimize(log.likelihood.cont,interval=c(-100,100),Y.obs=Y.obs[[i]],beta.fin = beta.fin,theta.fin=theta.fin,
                           X.i.standard = X.i.standard[[i]],lambda.fin = lambda.fin,v.obs=v.obs[i,]  )$minimum
      
    } else if (nrow(lambda.fin)==2){
      par<-c(0.5,0.5)
      test.optim<-stats::optim(par,log.likelihood.cont,Y.obs=Y.obs[[i]],beta.fin = beta.fin,theta.fin=theta.fin,
                               X.i.standard = X.i.standard[[i]],lambda.fin = lambda.fin ,v.obs=v.obs[i,] )$par
    } else if (nrow(lambda.fin)==3){
      par<-c(0.5,0.5,0.5)
      test.optim<-stats::optim(par,log.likelihood.cont,Y.obs=Y.obs[[i]],beta.fin = beta.fin,theta.fin=theta.fin,
                               X.i.standard = X.i.standard[[i]],lambda.fin = lambda.fin,v.obs=v.obs[i,]  )$par
    }
    xi.fin<-rbind(xi.fin,test.optim)
  }
  return(xi.fin)
}

# Use AIC to select the number of FPC and B-spline basis


marginal.likelihood.cont<-function(beta.fin,X.i.standard,J.fin,Y.obs,theta.fin,lambda.fin,var.err,xi.fin,v.obs){
  if(J.fin==1){
    temp.vec<-marginal.likelihood.inner.cont(beta.fin,X.i.standard,J.fin,Y.obs,theta.fin,lambda.fin,var.err,xi.fin,v.obs)
    res<--sum(log(temp.vec))}
  if(J.fin==2){
    temp.vec<-
      marginal.likelihood.inner.cont(beta.fin,X.i.standard,J.fin,Y.obs,theta.fin,lambda.fin,var.err,xi.fin,v.obs)
    res<--sum(log(temp.vec))}
  if(J.fin==3){
    temp.vec<-
      marginal.likelihood.inner.cont(beta.fin,X.i.standard,J.fin,Y.obs,theta.fin,lambda.fin,var.err,xi.fin,v.obs)
    res<--sum(log(temp.vec))}
  return(res)
}



AIC.estimation.cont<-function(n,Y.cont,X.i.standard,tau1=2,tau2=2,D.temp,n.basis,ob,v.obs){
  estimation.temp<-estimation.cont(Y.obs=Y.cont,X.i.standard=X.i.standard,
                                   tau1=tau1,tau2=tau2,D.temp=D.temp,n.basis=n.basis,v.obs=v.obs)
  beta.fin<-estimation.temp$beta.fin
  J.fin<-estimation.temp$J.fin
  theta.fin<-estimation.temp$theta.fin
  lambda.fin<-estimation.temp$lambda.fin
  var.err<-estimation.temp$var.err
  
  
  xi.fin<-xi.estimation.cont(log.likelihood.cont,n,Y.cont,beta.fin,theta.fin,
                             X.i.standard,lambda.fin,v.obs)
  
  likelihood.temp<-marginal.likelihood.cont(beta.fin,X.i.standard,J.fin,Y.cont,theta.fin,lambda.fin,var.err,xi.fin,v.obs)
  AIC.temp<-2*likelihood.temp+2*((J.fin+1)*n.basis+J.fin+1)
  return(AIC.temp)
}


estimation.tau.sele.cont<-function(Y.cont,t.grid,n,n.knots.sele,n.order,t.min,t.max,tau1.cand,tau2.cand,v.obs){
  ##### knot
  each.time<-t.max
  #breaks<-    1+  (0:(n.knots-1))*((each.time-1)/(n.knots-1))
  breaks<-seq(t.min,t.max,length.out = n.knots.sele)
  
  # degree = order-1
  n.basis<-length(breaks)+n.order-2
  
  ##### Threshold of selecting principal component
  
  # generate b-spline basis
  full.knots<- expand.knots(breaks, order=n.order) 
  ob<-OrthogonalSplineBasis(full.knots,order=n.order)
  X.i.standard<-list()
  for(i in 1:n){
    X.i.standard[[i]]<-evaluate(ob,t.grid[[i]]) 
  }
  # D.temp is penalty matrix
  D.temp<-OuterProdSecondDerivative(ob)
  
  tau.cand<-expand.grid(tau1.cand,tau2.cand)
  AIC.res<-NULL
  for(cand in 1:nrow(tau.cand)){
    res<-AIC.estimation.cont(n,Y.cont,X.i.standard,tau.cand[cand,1],tau.cand[cand,2],D.temp,n.basis,ob,v.obs)
    AIC.res<-c(AIC.res,res)
  }
  tau1.sele<-tau.cand[which.min(AIC.res),1]
  tau2.sele<-tau.cand[which.min(AIC.res),2]
  return(list("tau1"=tau1.sele,"tau2"=tau2.sele,"AIC.res"=AIC.res))
}

estimation.basis.sele.cont<-function(Y.cont,t.grid,n,n.knots,n.order,t.min,t.max,
                                     tau1,tau2,v.obs){
  
  AIC.res<-NULL
  for(cand in 1:length(n.knots)){
    ##### knot
    each.time<-t.max
    #breaks<-    1+  (0:(n.knots-1))*((each.time-1)/(n.knots-1))
    breaks<-seq(t.min,t.max,length.out = n.knots[cand])
    # degree = order-1
    n.basis<-length(breaks)+n.order-2
    ##### Threshold of selecting principal component
    # generate b-spline basis
    full.knots<- expand.knots(breaks, order=n.order) 
    ob<-OrthogonalSplineBasis(full.knots,order=n.order)
    X.i.standard<-list()
    for(i in 1:n){
      X.i.standard[[i]]<-evaluate(ob,t.grid[[i]]) 
    }
    # D.temp is penalty matrix
    D.temp<-OuterProdSecondDerivative(ob)
    
    res<-AIC.estimation.cont(n,Y.cont,X.i.standard,tau1,tau2,D.temp,n.basis,ob,v.obs)
    AIC.res<-c(AIC.res,res)
  }
  n.knots.sele<-n.knots[which.min(AIC.res)]
  return(list("n.knots.sele"=n.knots.sele,"AIC.res"=AIC.res))
}


