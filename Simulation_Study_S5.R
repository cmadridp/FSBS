rm(list=ls())

#Load libraries
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

#############################
# Simulation setting
#reptime <- 1000 #for the monte carlo
d <- 3 # number of PC used in the change point algorithm
no_grids <- 50 # number of equal-spaced grid on [0,1]
num <- 20 # number of basis for smoothing
N <-200
burnin=30
#############################
# Change size setting
# Combination 1
ts_type <- c('BM') # Brownian Bridge
# Combination 2
# ts_type_sets <- c('BM') # Brownian Motion

change_type <- c('Sint')
sim_type<- c('TwoCP_Dep_GKernel')

#############################
# Critical values (for monte carlo)
# Zhang and Shao 2011 EJS (This method will fail when there is non-monotonic change)
crit90_sz_ejs <- c(29.6, 56.5, 81.5, 114.7, 150,183.8,223.5, 267.1, 308.5, 360)
crit95_sz_ejs <- c(40.1, 73.7, 103.6, 141.5, 182.7, 218.8, 267.3, 317.9, 360.7, 420.5)
crit99_sz_ejs <- c(68.6, 117.7, 160.0, 209.7, 265.8, 318.3, 368.0, 432.5, 483.6, 567.2)
crit_sz_ejs <- c(crit90_sz_ejs[d], crit95_sz_ejs[d], crit99_sz_ejs[d])

# BGHK is BerkesGabrysHorvathKokoszka(2009) JRSSB. (This method does not allow temporal dependence.)
# HK is HormannKokoszka(2010)AoS. Same critical value (This allows temporal dependence.)
crit90_bghk <- c(0.345,0.6068,0.8426,1.0653,1.279713,1.4852,1.690773,1.897365,2.096615)
crit95_bghk <- c(0.4605,0.7488,1.0014,1.2397,1.469,1.6847,1.8956,2.1242,2.323)
crit99_bghk <- c(0.740,1.0721,1.3521,1.6267,1.866702,2.12595,2.342252,2.589244,2.80978)
crit_bghk <- c(crit90_bghk[d], crit95_bghk[d], crit99_bghk[d])





############################################### FSBS Scenario 5 #########################################

#Competitor results matrices
#EJS
results_EJS_90=matrix(nrow=0,ncol=2)
results_EJS_95=matrix(nrow=0,ncol=2)
results_EJS_99=matrix(nrow=0,ncol=2)


sum_scale_hd_EJS_90=0
sum_hd_EJS_90=0
sum_kd_EJS_90=0
K_d_EJS_90=matrix(nrow=0,ncol=1)
scale_hd_EJS_90=matrix(nrow=0,ncol=1)
hd_EJS_90=matrix(nrow=0,ncol=1)

sum_scale_hd_EJS_95=0
sum_hd_EJS_95=0
sum_kd_EJS_95=0
K_d_EJS_95=matrix(nrow=0,ncol=1)
scale_hd_EJS_95=matrix(nrow=0,ncol=1)
hd_EJS_95=matrix(nrow=0,ncol=1)

sum_scale_hd_EJS_99=0
sum_hd_EJS_99=0
sum_kd_EJS_99=0
K_d_EJS_99=matrix(nrow=0,ncol=1)
scale_hd_EJS_99=matrix(nrow=0,ncol=1)
hd_EJS_99=matrix(nrow=0,ncol=1)

#BGHK
results_BGHK_90=matrix(nrow=0,ncol=2)
results_BGHK_95=matrix(nrow=0,ncol=2)
results_BGHK_99=matrix(nrow=0,ncol=2)


sum_scale_hd_BGHK_90=0
sum_hd_BGHK_90=0
sum_kd_BGHK_90=0
K_d_BGHK_90=matrix(nrow=0,ncol=1)
scale_hd_BGHK_90=matrix(nrow=0,ncol=1)
hd_BGHK_90=matrix(nrow=0,ncol=1)

sum_scale_hd_BGHK_95=0
sum_hd_BGHK_95=0
sum_kd_BGHK_95=0
K_d_BGHK_95=matrix(nrow=0,ncol=1)
scale_hd_BGHK_95=matrix(nrow=0,ncol=1)
hd_BGHK_95=matrix(nrow=0,ncol=1)

sum_scale_hd_BGHK_99=0
sum_hd_BGHK_99=0
sum_kd_BGHK_99=0
K_d_BGHK_99=matrix(nrow=0,ncol=1)
scale_hd_BGHK_99=matrix(nrow=0,ncol=1)
hd_BGHK_99=matrix(nrow=0,ncol=1)

#HK
results_HK_90=matrix(nrow=0,ncol=2)
results_HK_95=matrix(nrow=0,ncol=2)
results_HK_99=matrix(nrow=0,ncol=2)

sum_scale_hd_HK_90=0
sum_hd_HK_90=0
sum_kd_HK_90=0
K_d_HK_90=matrix(nrow=0,ncol=1)
scale_hd_HK_90=matrix(nrow=0,ncol=1)
hd_HK_90=matrix(nrow=0,ncol=1)

sum_scale_hd_HK_95=0
sum_hd_HK_95=0
sum_kd_HK_95=0
K_d_HK_95=matrix(nrow=0,ncol=1)
scale_hd_HK_95=matrix(nrow=0,ncol=1)
hd_HK_95=matrix(nrow=0,ncol=1)

sum_scale_hd_HK_99=0
sum_hd_HK_99=0
sum_kd_HK_99=0
K_d_HK_99=matrix(nrow=0,ncol=1)
scale_hd_HK_99=matrix(nrow=0,ncol=1)
hd_HK_99=matrix(nrow=0,ncol=1)


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
  # Monotonic change
  change_size <- cbind(0, sin((0:no_grids)/no_grids), 2*sin((0:no_grids)/no_grids), 3*sin((0:no_grids)/no_grids))
  dep_type <- 'ARH1_GKernel'
  ts_dep <- T
  cp <- round(N*c(0,1/3,2/3,1))
  no_seg <- length(cp)-1
  change <- c()
  for(seg_index in 1:no_seg){
    seg_tmp <- matrix(rep(change_size[,seg_index], cp[seg_index+1]-cp[seg_index]), no_grids+1, cp[seg_index+1]-cp[seg_index])
    change <- cbind(change, seg_tmp)
  }
  change
  Psi <- matrix(0,no_grids+1,no_grids+1)
  for(i in 0:no_grids){
    for(j in 0:no_grids){
      Psi[i+1,j+1] <- psi1(i/no_grids,j/no_grids)
    }
  }
  set.seed(rep)
  xr <- matrix(rnorm(no_grids*(N+burnin),0,1),no_grids,(N+burnin))
  xr <- rbind(0,xr) # BM starts at 0
  BM <- apply(xr,2,cumsum)/sqrt(no_grids)
  sim.data <- BM
  data <- matrix(0,no_grids+1,(N+burnin+1))
  for(i in 1:(N+burnin)){
    X <- (1/no_grids)*Psi%*%data[,i]+sim.data[,i]
    data[,i+1] <- X
  }
  AR <- data[,(burnin+2):(N+burnin+1)]
  fts=AR+change
  
  time3 <- proc.time()
  EJS_result <- SN_BS_functional(data=fts, d=d, num_basis=num, rangeval=c(0,1), no_grids=no_grids,
                                 start=1, end=N, critical_value=crit_sz_ejs)
  time4 <- proc.time()
  EJS_result <- post_cp_extract(EJS_result)
  EJS_result <- list(EJS_result90=EJS_result$cp_result90, EJS_result95=EJS_result$cp_result95,
                     EJS_result99=EJS_result$cp_result99, EJS_time=as.numeric(time4-time3)[1])
  EJS_result
  
  # BGHK
  time5 <- proc.time()
  BGHK_result <- BGHK_BS_functional(data=fts, d=d, num_basis=num, rangeval=c(0,1), no_grids=no_grids,
                                    start=1, end=N, critical_value=crit_bghk)
  time6 <- proc.time()
  BGHK_result <- post_cp_extract(BGHK_result)
  BGHK_result <- list(BGHK_result90=BGHK_result$cp_result90, BGHK_result95=BGHK_result$cp_result95,
                      BGHK_result99=BGHK_result$cp_result99, BGHK_time=as.numeric(time6-time5)[1])
  BGHK_result
  # HK
  time7 <- proc.time()
  HK_result <- HK_BS_functional(data=fts, d=d, num_basis=num, rangeval=c(0,1), no_grids=no_grids,
                                start=1, end=N, critical_value=crit_bghk)
  time8 <- proc.time()
  HK_result <- post_cp_extract(HK_result)
  HK_result <- list(HK_result90=HK_result$cp_result90, HK_result95=HK_result$cp_result95,
                    HK_result99=HK_result$cp_result99, HK_time=as.numeric(time8-time7)[1])
  HK_result 
  
  
  
  ####################################### FSBS ####################################### 
  values <- ((0:no_grids)/no_grids)
  X_s=matrix(NaN,(no_grids+1)*N,1)
  aux=0
  for (i in 1:N) {
    aux=no_grids+1+aux
    for (j in 1:(no_grids+1)) {
      X_s[aux-no_grids-1+j]=values[j]
    }
  }
  X_s
  
  y_t_is=matrix(NaN,(no_grids+1)*N,1)
  aux=0
  for (i in 1:N) {
    aux=no_grids+1+aux
    for (j in 1:(no_grids+1)) {
      y_t_is[aux-no_grids-1+j]=fts[j,i]
    }
  }
  y_t_is
  
  Tb=N
  r=1
  C=3
  #The responses lies in dimension 3
  n=no_grids+1
  #The observations lies in dimension 1
  d_FSBS=1
  #we define the n_ts
  nts<-matrix(NaN,nrow=Tb, ncol=1)
  for(i in 1:Tb){
    nts[i]=n
  }
  N_t=sum(nts)
  N_t
  #we create the data matrix containing x_t_is, y_t_is and labels
  labels_Ts=matrix(NaN,nrow=N_t,ncol=1)
  aux1=0
  for (i in 1:Tb) {
    aux1=nts[i]+aux1
    for (j in 1:nts[i]) {
      labels_Ts[aux1-nts[i]+j]=i
    }
  }
  labels_Ts
  dat=matrix(NaN,nrow=N_t)
  dat=cbind(labels_Ts,X_s)
  dat=cbind(dat,y_t_is)
  #dat
  
  #The choice of bandwidth h_bar
  H_bar_stimated=hpi(X_s)
  #H_bar_stimated
  h_bar=H_bar_stimated
  
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
  dat_train[,1]=ceiling(dat_train[,1]/2)

  
  #We create the seeded intervals to work with training data
  s.inter_train<-seeded.intervals(Tb/2,C)
  
  
  #We create the estimator functions 
  g_hat_i=function(eta1,eta2,h,x,dat_train,h_bar,phat)
  {
    res=0
    for (i in (eta1+1):eta2) {
      res=res+statistic(x,dat_train,i,h,h_bar,d_FSBS,phat)
    }
    res=res/(eta2-eta1)
    return(res)
  }
  
  #We create the possible values for h and tau
  h_int=c(0.00349,0.0035,0.00351)
  l_h_int=length(h_int)
  tau_int=c(0.0069,0.007,0.0071)
  l_tau_int=length(tau_int)
  
  #We compute errors of estimation
  errors=matrix(NaN,l_h_int,l_tau_int)
  for (ind1 in 1:l_h_int) {
    for (ind2 in 1:l_tau_int) {
      S=seedBS(dat_train, s.inter_train, h_int[ind1], h_bar, tau_int[ind2], d_FSBS,n )
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
              X_t=dat_test[indices,2:(d_FSBS+1)]
              Y_t=dat_test[indices,(d_FSBS+2):(d_FSBS+2)]
              phat=p_hat(X_t[i],dat_train,h_bar,d_FSBS)
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

  #Using the data---------------------------------------------------
  #We create the seeded intervals
  s.inter<-seeded.intervals(Tb,C)
  S_1=seedBS(dat, s.inter, h_min, h_bar, tau_min, d_FSBS,n )

  aux_matrix=matrix(NaN,length(S_1),2)
  aux_matrix[,1]=count_iter
  aux_matrix[,2]=sort(S_1)
  results_FSBS=rbind(results_FSBS,aux_matrix)
  vec2=c(68,134)
  dis=Hausdorff(S_1,vec2) 
  hd=rbind(hd,dis)
  scale_hd=rbind(scale_hd,dis/200)
  sum_hd=sum_hd+dis
  sum_scale_hd=sum_scale_hd+(dis/200)
  
  dist=K_distance(S_1,vec2)
  K_d=rbind(K_d,dist)
  sum_kd=sum_kd+dist
  
  
  #Competitors errors
  #Results_EJS_90
  aux_matrix_EJS_90=matrix(NaN,length(EJS_result[[1]]),2)
  aux_matrix_EJS_90[,1]=count_iter
  aux_matrix_EJS_90[,2]=sort(EJS_result[[1]])
  results_EJS_90=rbind(results_EJS_90,aux_matrix_EJS_90)
  
  vec2=c(68,134)
  dis_EJS_90=Hausdorff(EJS_result[[1]],vec2) 
  hd_EJS_90=rbind(hd_EJS_90,dis_EJS_90)
  scale_hd_EJS_90=rbind(scale_hd_EJS_90,dis_EJS_90/200)
  sum_hd_EJS_90=sum_hd_EJS_90+dis_EJS_90
  sum_scale_hd_EJS_90=sum_scale_hd_EJS_90+(dis_EJS_90/200)
  
  dist_EJS_90=K_distance(EJS_result[[1]],vec2)
  K_d_EJS_90=rbind(K_d_EJS_90,dist_EJS_90)
  sum_kd_EJS_90=sum_kd_EJS_90+dist_EJS_90
  
  #Results_EJS_95
  aux_matrix_EJS_95=matrix(NaN,length(EJS_result[[2]]),2)
  aux_matrix_EJS_95[,1]=count_iter
  aux_matrix_EJS_95[,2]=sort(EJS_result[[2]])
  results_EJS_95=rbind(results_EJS_95,aux_matrix_EJS_95)
  vec2=c(68,134)
  dis_EJS_95=Hausdorff(EJS_result[[2]],vec2) 
  hd_EJS_95=rbind(hd_EJS_95,dis_EJS_95)
  scale_hd_EJS_95=rbind(scale_hd_EJS_95,dis_EJS_95/200)
  sum_hd_EJS_95=sum_hd_EJS_95+dis_EJS_95
  sum_scale_hd_EJS_95=sum_scale_hd_EJS_95+(dis_EJS_95/200)
  
  dist_EJS_95=K_distance(EJS_result[[2]],vec2)
  K_d_EJS_95=rbind(K_d_EJS_95,dist_EJS_95)
  sum_kd_EJS_95=sum_kd_EJS_95+dist_EJS_95
  
  #Results_EJS_99
  aux_matrix_EJS_99=matrix(NaN,length(EJS_result[[3]]),2)
  aux_matrix_EJS_99[,1]=count_iter
  aux_matrix_EJS_99[,2]=sort(EJS_result[[3]])
  results_EJS_99=rbind(results_EJS_99,aux_matrix_EJS_99)
  
  vec2=c(68,134)
  dis_EJS_99=Hausdorff(EJS_result[[3]],vec2) 
  hd_EJS_99=rbind(hd_EJS_99,dis_EJS_99)
  scale_hd_EJS_99=rbind(scale_hd_EJS_99,dis_EJS_99/200)
  sum_hd_EJS_99=sum_hd_EJS_99+dis_EJS_99
  sum_scale_hd_EJS_99=sum_scale_hd_EJS_99+(dis_EJS_99/200)
  
  dist_EJS_99=K_distance(EJS_result[[3]],vec2)
  K_d_EJS_99=rbind(K_d_EJS_99,dist_EJS_99)
  sum_kd_EJS_99=sum_kd_EJS_99+dist_EJS_99
  
  #Results_BGHK_90
  aux_matrix_BGHK_90=matrix(NaN,length(BGHK_result[[1]]),2)
  aux_matrix_BGHK_90[,1]=count_iter
  aux_matrix_BGHK_90[,2]=sort(BGHK_result[[1]])
  results_BGHK_90=rbind(results_BGHK_90,aux_matrix_BGHK_90)
  
  vec2=c(68,134)
  dis_BGHK_90=Hausdorff(BGHK_result[[1]],vec2) 
  hd_BGHK_90=rbind(hd_BGHK_90,dis_BGHK_90)
  scale_hd_BGHK_90=rbind(scale_hd_BGHK_90,dis_BGHK_90/200)
  sum_hd_BGHK_90=sum_hd_BGHK_90+dis_BGHK_90
  sum_scale_hd_BGHK_90=sum_scale_hd_BGHK_90+(dis_BGHK_90/200)
  
  dist_BGHK_90=K_distance(BGHK_result[[1]],vec2)
  K_d_BGHK_90=rbind(K_d_BGHK_90,dist_BGHK_90)
  sum_kd_BGHK_90=sum_kd_BGHK_90+dist_BGHK_90
  
  #Results_BGHK_95
  aux_matrix_BGHK_95=matrix(NaN,length(BGHK_result[[2]]),2)
  aux_matrix_BGHK_95[,1]=count_iter
  aux_matrix_BGHK_95[,2]=sort(BGHK_result[[2]])
  results_BGHK_95=rbind(results_BGHK_95,aux_matrix_BGHK_95)
  
  vec2=c(68,134)
  dis_BGHK_95=Hausdorff(BGHK_result[[2]],vec2) 
  hd_BGHK_95=rbind(hd_BGHK_95,dis_BGHK_95)
  scale_hd_BGHK_95=rbind(scale_hd_BGHK_95,dis_BGHK_95/200)
  sum_hd_BGHK_95=sum_hd_BGHK_95+dis_BGHK_95
  sum_scale_hd_BGHK_95=sum_scale_hd_BGHK_95+(dis_BGHK_95/200)
  
  dist_BGHK_95=K_distance(BGHK_result[[2]],vec2)
  K_d_BGHK_95=rbind(K_d_BGHK_95,dist_BGHK_95)
  sum_kd_BGHK_95=sum_kd_BGHK_95+dist_BGHK_95
  
  #Results_BGHK_99
  aux_matrix_BGHK_99=matrix(NaN,length(BGHK_result[[3]]),2)
  aux_matrix_BGHK_99[,1]=count_iter
  aux_matrix_BGHK_99[,2]=sort(BGHK_result[[3]])
  results_BGHK_99=rbind(results_BGHK_99,aux_matrix_BGHK_99)
  
  vec2=c(68,134)
  dis_BGHK_99=Hausdorff(BGHK_result[[3]],vec2) 
  hd_BGHK_99=rbind(hd_BGHK_99,dis_BGHK_99)
  scale_hd_BGHK_99=rbind(scale_hd_BGHK_99,dis_BGHK_99/200)
  sum_hd_BGHK_99=sum_hd_BGHK_99+dis_BGHK_99
  sum_scale_hd_BGHK_99=sum_scale_hd_BGHK_99+(dis_BGHK_99/200)
  
  dist_BGHK_99=K_distance(BGHK_result[[3]],vec2)
  K_d_BGHK_99=rbind(K_d_BGHK_99,dist_BGHK_99)
  sum_kd_BGHK_99=sum_kd_BGHK_99+dist_BGHK_99
  
  #Results_HK_90
  aux_matrix_HK_90=matrix(NaN,length(HK_result[[1]]),2)
  aux_matrix_HK_90[,1]=count_iter
  aux_matrix_HK_90[,2]=sort(HK_result[[1]])
  results_HK_90=rbind(results_HK_90,aux_matrix_HK_90)
  
  
  vec2=c(68,134)
  dis_HK_90=Hausdorff(HK_result[[1]],vec2) 
  hd_HK_90=rbind(hd_HK_90,dis_HK_90)
  scale_hd_HK_90=rbind(scale_hd_HK_90,dis_HK_90/200)
  sum_hd_HK_90=sum_hd_HK_90+dis_HK_90
  sum_scale_hd_HK_90=sum_scale_hd_HK_90+(dis_HK_90/200)
  
  dist_HK_90=K_distance(HK_result[[1]],vec2)
  K_d_HK_90=rbind(K_d_HK_90,dist_HK_90)
  sum_kd_HK_90=sum_kd_HK_90+dist_HK_90
  
  #Results_HK_95
  aux_matrix_HK_95=matrix(NaN,length(HK_result[[2]]),2)
  aux_matrix_HK_95[,1]=count_iter
  aux_matrix_HK_95[,2]=sort(HK_result[[2]])
  results_HK_95=rbind(results_HK_95,aux_matrix_HK_95)
  
  vec2=c(68,134)
  dis_HK_95=Hausdorff(HK_result[[2]],vec2) 
  hd_HK_95=rbind(hd_HK_95,dis_HK_95)
  scale_hd_HK_95=rbind(scale_hd_HK_95,dis_HK_95/200)
  sum_hd_HK_95=sum_hd_HK_95+dis_HK_95
  sum_scale_hd_HK_95=sum_scale_hd_HK_95+(dis_HK_95/200)
  
  dist_HK_95=K_distance(HK_result[[2]],vec2)
  K_d_HK_95=rbind(K_d_HK_95,dist_HK_95)
  sum_kd_HK_95=sum_kd_HK_95+dist_HK_95
  
  #Results_HK_99
  aux_matrix_HK_99=matrix(NaN,length(HK_result[[3]]),2)
  aux_matrix_HK_99[,1]=count_iter
  aux_matrix_HK_99[,2]=sort(HK_result[[3]])
  results_HK_99=rbind(results_HK_99,aux_matrix_HK_99)
  
  vec2=c(68,134)
  dis_HK_99=Hausdorff(HK_result[[3]],vec2) 
  hd_HK_99=rbind(hd_HK_99,dis_HK_99)
  scale_hd_HK_99=rbind(scale_hd_HK_99,dis_HK_99/200)
  sum_hd_HK_99=sum_hd_HK_99+dis_HK_99
  sum_scale_hd_HK_99=sum_scale_hd_HK_99+(dis_HK_99/200)
  
  dist_HK_99=K_distance(HK_result[[3]],vec2)
  K_d_HK_99=rbind(K_d_HK_99,dist_HK_99)
  sum_kd_HK_99=sum_kd_HK_99+dist_HK_99
}




