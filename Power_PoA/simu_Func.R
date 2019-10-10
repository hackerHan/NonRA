#Function for calculating the discrepancy meansure of series system
discrep_series<-function(x)
{
  Ps<-x["Ps"]
  P1<-x["P1"]
  P2<-x["P2"]
  P3<-x["P3"]
  P4<-x["P4"]
  P5<-x["P5"]
  disc<-Ps/(P1*P2*P3*P4*P5)
  return(disc)
}

qt1=function(x)
{
  quantile(x,prob=0.05)
}
qt2=function(x)
{
  quantile(x,prob=0.95)
}


#nsim represents simulation times;
#M represents number of MC samples
#D1 represents true discrepancy;
#N represents sample size of system;
#Nc represents sample size of component of which the default is 100
#Dv1 represents the given range of discrepancy to calculate the PoA
simu_pw<-function(x,Pc=0.921,M=1000,nsim=10000,Nc=100,Dv1=c(0.95,1.05),D1,N)
{
  ns<-x[1] #Given sample size of system data
  d1<-x[2] #Given True discrepancy for simulating data
  #Subfunction
  simu<-function(D1,N)
  {
    D1=d1
    N=ns
    MC_comp<-data.frame(X=rbinom(5,size=Nc,prob=Pc),N=rep(Nc,5)) #Generate MC sample for component data
    MC_sys<-data.frame(X=rbinom(1,N,Pc^5*D1),N=N) #Generate MC sample for system data
    xc<-MC_comp$X #Pass number of componenent data
    xs<-MC_sys$X #Pass number of system data
    #Generate fractional random weights for system data
    rwt_sys<-rdirichlet(nsim,rep(1,N))*N 
    #Generate fractional random weights independently for component data
    rwt_c1<-rdirichlet(nsim,rep(1,Nc))*Nc
    rwt_c2<-rdirichlet(nsim,rep(1,Nc))*Nc
    rwt_c3<-rdirichlet(nsim,rep(1,Nc))*Nc
    rwt_c4<-rdirichlet(nsim,rep(1,Nc))*Nc
    rwt_c5<-rdirichlet(nsim,rep(1,Nc))*Nc
    
    #Generate data for simulation
    sys<-matrix(sample(c(rep(1,xs),rep(0,(N-xs))),N),N,1)
    X1<-matrix(sample(c(rep(1,xc[1]),rep(0,(Nc-xc[1]))),Nc),Nc,1)
    X2<-matrix(sample(c(rep(1,xc[2]),rep(0,(Nc-xc[2]))),Nc),Nc,1)
    X3<-matrix(sample(c(rep(1,xc[3]),rep(0,(Nc-xc[3]))),Nc),Nc,1)
    X4<-matrix(sample(c(rep(1,xc[4]),rep(0,(Nc-xc[4]))),Nc),Nc,1)
    X5<-matrix(sample(c(rep(1,xc[5]),rep(0,(Nc-xc[5]))),Nc),Nc,1)
    
    #Generate simulated reliability for system and components using fractional random weights
    simu_rb<-data.frame(index=1:nsim,Ps=rwt_sys%*%sys/length(sys),
                        P1=rwt_c1%*%X1/length(X1),
                        P2=rwt_c2%*%X2/length(X2),
                        P3=rwt_c3%*%X3/length(X3),
                        P4=rwt_c4%*%X4/length(X4),
                        P5=rwt_c5%*%X5/length(X5))
    simu_disc<-apply(simu_rb,1,discrep_series)
    pl<-sum(simu_disc<1)/nsim
    pu<-sum(simu_disc>1)/nsim
    P_2s<-2*min(pl,pu) #two-sided P value
    d1_est<-quantile(simu_disc,0.5) #Estimated discrepancy (Median value)
    d1_hwci<-(quantile(simu_disc,0.975)-quantile(simu_disc,0.025))/2 #Half width of confidence interval at α=0.05
    #Probability of agreement
    poa1=mean(simu_disc<=Dv1[2]&simu_disc>=Dv1[1])
    return(c(P_2s,d1_est,d1_hwci,poa1))
  }
  vec<-matrix(data=NA,nrow=M,ncol=4)
  colnames(vec)<-c("two-sided P value","Estimated discrepancy","Half width of CI","Probability of agreement")
  #Calculate P value & Estimated discrepancy & Half width of confidence interval at α=0.05 for 1000 MC samples with specific combination of sample size and true discrepancy
  #Calculate PoA
  temp<-apply(vec,1,simu)
  vec[,"two-sided P value"]<-temp[1,]
  vec[,"Estimated discrepancy"]<-temp[2,]
  vec[,"Half width of CI"]<-temp[3,]
  vec[,"Probability of agreement"]<-temp[4,]
  save(vec,file=paste("simu/s",ns,"d",d1,".rdata",sep=""))
  # Power at α=0.05
  pr0.05<-mean(temp[1,]<=0.05)
  # Power at a=0.1
  pr0.1<-mean(temp[1,]<=0.1)
  # Half width of credible interval
  cv.d1r<-mean(temp[3,])
  # Relative half Width of credible interval
  cv.d1<-cv.d1r/d1
  # Average PoA
  poa1_ave<-mean(temp[4,])
  # Lower PoA
  poal<-qt1(temp[4,])
  # Upper PoA
  poau<-qt2(temp[4,])
  return(c(pr0.05,pr0.1,cv.d1,cv.d1r,poa1_ave,poal,poau))
}

