rm(list=ls())

#load libraries
library(fda)
library(ks)
library(mvtnorm)

#############################
###### Load functions #######
#############################
source("FSBS.R")
source("Auxiliar_Functions_FSBS.R")
source("Distances.R")
source("Competitors_Functions.R")




############################################### FSBS Scenario 4 #########################################


Tb=200
r=1
C=3
#the responses lies in dimension 3
n=10
#the observations lies in dimension 1
d=2
#we define the n_ts
nts<-matrix(NaN,nrow=Tb, ncol=1)
for(i in 1:Tb){
  nts[i]=n
}
N_t=sum(nts)
N_t



#we create our functions f*s
f1_star=function(x){
  return(0)
}
f2_star=function(x){
  aux=3*x[1]*x[2]
  return(aux)
}

f3_star=function(x){
  aux=0
  return(aux)
}


repetitions=100
sum_scale_hd=0
sum_hd=0
sum_kd=0
K_d=matrix(nrow=0,ncol=1)
scale_hd=matrix(nrow=0,ncol=1)
hd=matrix(nrow=0,ncol=1)
results_FSBS=matrix(nrow=0,ncol=2)
count_iter=0

for(rep in 1:repetitions){
  count_iter=count_iter+1
  #we create the xs
  X_s<-matrix(NaN,nrow=N_t, ncol=d)
  set.seed(rep)
  for(i in 1:N_t){
    for(j in 1:d)
    {
      X_s[i,j]=runif(1,0,1)
    }
  }
  #X_s
  
  
  #we create the \xi_s--------------------------------------------
  val_in_0_1=0.5
  number_basis=50
  set.seed(rep)
  b_t_js_0=matrix( rnorm(number_basis*(Tb+1),mean=0,sd=1), number_basis, Tb+1) 
  xi_t_matrix=matrix(NaN,n,Tb)
  aux=0
  for (t in 1:Tb) {
    aux=nts[t]+aux
    for (i in 1:n) {
      #xti is X_s[aux-nts[t]+i,]
      X_s[aux-nts[t]+i,]
      #creating xi_0 at xti to be aux_2
      b_t_0=b_t_js_0[,1]
      aux_2=0
      for (j in 1:number_basis) {
        #creating h_j at x_ti to be aux_1
        aux_1=1
        aux_1=(sqrt(2)/2*pi)*sin(j*X_s[aux-nts[t]+i,1])*(aux_1)*(1/45)*(sqrt(2)/2*pi)*cos((number_basis-j)*X_s[aux-nts[t]+i,2])
        aux_2=aux_2+b_t_0[j]*(aux_1)
      }
      #creating xi_t at xti to be aux_3
      aux_3=val_in_0_1^t*(aux_2)
      for(l in 1:t){
        #we need xi_t_aux(x,i) to be aux_4
        b_t_0_js=b_t_js_0[,l]
        aux_4=0
        for (j in 1:number_basis) {
          aux_1=1
          aux_1=(sqrt(2)/2*pi)*sin(j*X_s[aux-nts[t]+i,1])*(aux_1)*(1/45)*(sqrt(2)/2*pi)*cos((number_basis-j)*X_s[aux-nts[t]+i,2])
          aux_4=aux_4+b_t_0_js[j]*(aux_1)
        }
        #aux=aux+val_in_0_1^(t-i)*b_t_js[i]*h_j(x,i)
        aux_3=aux_3+val_in_0_1^(t-l)*(aux_4)
      }
      xi_t_matrix[i,t]=aux_3
    }
  }
  #we create the delta_ts--------------------------------------
  delta_s<-matrix(NaN,nrow=N_t, ncol=1)
  m=max(nts)
  #m
  #ruido gaussiano
  sigma=matrix(NaN,nrow=m,ncol=m)
  for(i in 1:m){
    for(j in 1:m){
      if(i==j){
        sigma[i,j]=0.5
      }
      else{
        sigma[i,j]=0
      }
    }
  }
  mean=matrix(NaN,nrow=m,ncol=1)
  for(i in 1:m){
    mean[i]=0
  }
  set.seed(rep)
  normal_noise <- rmvnorm(Tb, mean=mean, sigma=sigma)
  #normal_noise
  normal_noise=t(normal_noise)
  #normal_noise
  #matrix A
  A=diag(0.3,m,m)
  
  #delta_s_0=rmvnorm(1, mean=mean, sigma=diag(1,m,m))
  set.seed(rep)
  delta_s_0=rnorm(m,0,1)
  Delta=A%*%delta_s_0+normal_noise[,1]
  Delta=as.matrix(Delta)
  #Delta
  for(t in 2:Tb)
  {
    vaux=A%*%Delta[,t-1]+normal_noise[,t]
    Delta=cbind(Delta,vaux)
  }
  #Delta
  
  
  
  #we create the y_i_ts----------------------------------------
  y_t_is=matrix(NaN,nrow=N_t,ncol=1)
  aux=0
  for(t in 1:Tb)
  {
    aux=nts[t]+aux
    
    if (t<=100)
    {
      for(i in 1:nts[t])
      {
        y_t_is[aux-nts[t]+i]=f1_star(X_s[aux-nts[t]+i,])+xi_t_matrix[i,t]+Delta[i,t]
      }
      
    }
    else if(100<t&&t<=150)
      for(i in 1:nts[t])
      {
        y_t_is[aux-nts[t]+i]=f2_star(X_s[aux-nts[t]+i,])+xi_t_matrix[i,t]+Delta[i,t]
      }
    
    
    else if(t>150)
    {
      for(i in 1:nts[t])
      {
        y_t_is[aux-nts[t]+i]=f3_star(X_s[aux-nts[t]+i,])+xi_t_matrix[i,t]+Delta[i,t]
      }
    }
  }
  y_t_is
  
  
  #we create the data matrix containing x_t_is, y_t_is and labels
  labels_Ts=matrix(NaN,nrow=N_t,ncol=1)
  aux1=0
  for (i in 1:Tb) {
    aux1=nts[i]+aux1
    for (j in 1:nts[i]) {
      labels_Ts[aux1-nts[i]+j]=i
    }
  }
  #labels_Ts
  dat=matrix(NaN,nrow=N_t)
  dat=cbind(labels_Ts,X_s)
  dat=cbind(dat,y_t_is)
  dat
  
  
  
  
  #using the data---------------------------------------------------
  #we create the seeded intervals
  s.inter<-seeded.intervals(Tb,C)
  s.inter
  
  #s.inter_train<-seeded.intervals(Tb/2,C)
  #s.inter_train
  
  
  #the choice of bandwidth h_bar
  H_bar_stimated=Hpi.diag(X_s)
  H_bar_stimated
  h_bar_aux=H_bar_stimated
  h_bar=sqrt(mean(diag(h_bar_aux)))
  
  
  #we define tau and h
  #tau_1=6.416
  #h=1
  #tau=(7)*log(Tb)*sqrt(1/(n*(h^2))+1)
  #tau
  #we call our algorithm
  #start_time <- Sys.time()
  #S_1=seedBS(dat, s.inter, h, h_bar, tau, d,n )
  #end_time <- Sys.time()
  #end_time - start_time
  #S_1
  
  
  
  
  ######################################################
  ######################### CV #########################
  ######################################################
  #creating testing and training sets.
  dat_train=dat[(n+1):(2*n),]
  dat_test=dat[1:n,]
  for(i in 1:(Tb/2-1)){
    dat_test=rbind(dat_test,dat[(2*i*n+1):((2*i+1)*n),])
    dat_train=rbind(dat_train,dat[((2*i+1)*n+1):((2*(i+1))*n),])
  }
  #dat_test
  #dat_train
  dat_train[,1]=ceiling(dat_train[,1]/2)
  #dat_train
  
  #we create the seeded intervals to work with training data
  s.inter_train<-seeded.intervals(Tb/2,C)
  #s.inter_train
  
  
  #we create the estimator functions 
  g_hat_i=function(eta1,eta2,h,x,dat_train,h_bar,phat)
  {
    res=0
    for (i in (eta1+1):eta2) {
      res=res+statistic(x,dat_train,i,h,h_bar,d,phat)
    }
    res=res/(eta2-eta1)
    return(res)
  }
  
  #we create the possible values for h and tau
  h_int=seq(0.99,1.01,0.01)
  l_h_int=length(h_int)
  tau_int=c(.445,.45,.455)
  l_tau_int=length(tau_int)
  
  #we compute errors of estimation
  
  errors=matrix(NaN,l_h_int,l_tau_int)
  #start_time <- Sys.time()
  for (ind1 in 1:l_h_int) {
    for (ind2 in 1:l_tau_int) {
      S=seedBS(dat_train, s.inter_train, h_int[ind1], h_bar, tau_int[ind2], d,n )
      if(is.null(S)){
        errors[ind1,ind2]=10000000000
        next
      }
      else{
        S=as.vector(S)
        S=append(S,0)
        S=append(S,Tb/2)
        S=sort(S)
        error=0
        for(j in 1:(length(S)-1)){
          for(t in (S[j]+1):S[j+1]){
            for(i in 1:n){
              indices=which(dat_test[,1]==(2*t-1))
              X_t=dat_test[indices,2:(d+1)]
              Y_t=dat_test[indices,(d+2):(d+2)]
              phat=p_hat(X_t[i],dat_train,h_bar,d)
              error=error+(g_hat_i(S[j],S[j+1],h_int[ind1],X_t[i],dat_train,h_bar,phat)-Y_t[i])^2
            }
          }
        }
        errors[ind1,ind2]=error
      }  
    }
  }
  min_error=min(errors)
  h_tau=which(errors==min_error,arr.ind = TRUE)
  h_min=min(h_tau[,1])
  h_min=h_int[h_min]
  tau_min=max(h_tau[,2])
  tau_min=tau_int[tau_min]
  #using the data---------------------------------------------------
  #we create the seeded intervals
  s.inter<-seeded.intervals(Tb,C)
  S_1=seedBS(dat, s.inter, h_min, h_bar, tau_min, d,n )
  aux_matrix=matrix(NaN,length(S_1),2)
  aux_matrix[,1]=count_iter
  aux_matrix[,2]=sort(S_1)
  results_FSBS=rbind(results_FSBS,aux_matrix)
  vec2=c(100,150)
  dis=Hausdorff(S_1,vec2) 
  hd=rbind(hd,dis)
  scale_hd=rbind(scale_hd,dis/200)
  sum_hd=sum_hd+dis
  sum_scale_hd=sum_scale_hd+(dis/200)
  
  dist=K_distance(S_1,vec2)
  K_d=rbind(K_d,dist)
  sum_kd=sum_kd+dist
}
