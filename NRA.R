#First use Discrepancy measure 1(posted by Anderson Cook)

#Original data table for series system
real_series<-data.frame(data_type=c("system","c1","c2","c3","c4","c5"),X=c(17,111,96,76,98,79),N=c(25,120,100,80,100,80))
real_series$P=real$X/real$N
options(max.print=1e+06)

#Generate real data for simulation in the next section
set.seed(1)
sys<-sample(c(rep(1,17),rep(0,8)),25)
X1<-sample(c(rep(1,111),rep(0,9)),120)
X2<-sample(c(rep(1,96),rep(0,4)),100)
X3<-sample(c(rep(1,76),rep(0,4)),80)
X4<-sample(c(rep(1,98),rep(0,2)),100)
X5<-sample(c(rep(1,79),rep(0,1)),80)

#Define regular bootstrap resampling function(same with integer weight bootstrap resampling)
simu<-function(x)
{
  sim_data<-sample(x,length(x),replace=TRUE)
  targ<-length(sim_data[sim_data==1])
  return(targ)
}
#Define interval calculating function
quantf=function(x)
{
  quantile(x,prob=c(0.025,0.5,0.975))
}

#set bootstrap sample size 200000,generate 200000 regular resampling simulation data
set.seed(1)
sys_sim<-replicate(200000,simu(sys))
X1_sim<-replicate(200000,simu(X1))
X2_sim<-replicate(200000,simu(X2))
X3_sim<-replicate(200000,simu(X3))
X4_sim<-replicate(200000,simu(X4))
X5_sim<-replicate(200000,simu(X5))

#===========================================================================#
#First explore the scenario in series system

#generate simulation table suited for regular bootstrap
simu_table_series<-data.frame(index=1:200000,system=sys_sim,component1=X1_sim,component2=X2_sim,component3=X3_sim,
                              component4=X4_sim,component5=X5_sim,Ps=sys_sim/length(sys),P1=X1_sim/length(X1),P2=X2_sim/length(X2),
                              P3=X3_sim/length(X3),P4=X4_sim/length(X4),P5=X5_sim/length(X5))


#Define the discrepancy meansure of series system
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

#Generate discrepancy values for original data and simulated data
real_series_discrep<-real_series$P[1]/prod(real_series$P[2:6])
simu_series_discrep<-apply(simu_table_series,1,discrep_series)

#Generate histogram for discrepancy for series system

#Frequency histogram
par(mfrow=c(1,1))
hist(simu_series_discrep,freq=TRUE,breaks=seq(min(simu_series_discrep),max(simu_series_discrep),length.out=20),xlim=range(simu_series_discrep),main="Discrepancy measure for series system",xlab="Discrepancy measure")
abline(v=real_series_discrep,col="red",lty=2)
abline(v=quantf(simu_series_discrep),col="blue",lty=3)

#Density histogram
hist(simu_series_discrep,freq=FALSE,breaks=seq(min(simu_series_discrep),max(simu_series_discrep),length.out=20),xlim=range(simu_series_discrep),main="Discrepancy measure for series system",xlab="Discrepancy measure")
lines(density(simu_series_discrep),col="blue",lty=1)
abline(v=real_series_discrep,col="red",lty=2)
abline(v=quantf(simu_series_discrep),col="blue",lty=3)

#Kernal density plot for system discrepancy
plot(density(simu_series_discrep),type="l",main="density of discrepancy")
abline(v=real_series_discrep,col="red",lty=2)
text(x=real_series_discrep,y=0,labels=paste("real data discrep=",round(real_series_discrep,3)))
abline(v=quantf(simu_series_discrep),col="blue",lty=3)
text(x=median(simu_series_discrep),y=3.0,labels=paste("median=",round(median(simu_series_discrep),3)))
abline(v=1,col="red",lty=4)

#Calculate the cumulative probability of discrepancy in real data for series system
cum_disc_real_series<-diff(ecdf(simu_series_discrep)(c(0,real_series_discrep)))
diff(ecdf(simu_series_discrep)(c(0,1)))
#P_2sided=(1-0.91441)*2=0.17118>0.05 
#Accept the  Null hypothesis


#Generate histogram and density plot for system reliability and components
par(mfrow=c(2,3))
hist(simu_table_series$Ps,freq=TRUE,main="System reliability simulations",xlab="system reliability")
hist(simu_table_series$P1,freq=TRUE,main="Component1 reliability simulations",xlab="Component1 reliablity")
hist(simu_table_series$P2,freq=TRUE,main="Component2 reliability simulations",xlab="Component2 reliablity")
hist(simu_table_series$P3,freq=TRUE,main="Component3 reliability simulations",xlab="Component3 reliablity")
hist(simu_table_series$P4,freq=TRUE,main="Component4 reliability simulations",xlab="Component4 reliablity")
hist(simu_table_series$P5,freq=TRUE,main="Component5 reliability simulations",xlab="Component5 reliablity")

#Generate empirical distribution for discrepancy
par(mfrow=c(1,1))
plot(ecdf(simu_series_discrep),verticals=TRUE,do.p=FALSE,main="Discrepancy measure empirical distribution for series system",xlab="simulation discrepancy for series system")
abline(v=real_series_discrep,lty=4,col="blue")

#Generate empirical distribution for system reliability and components
par(mfrow=c(3,2))
plot(ecdf(simu_table_series$Ps),verticals=TRUE,do.p=TRUE,main="system reliability empirical distribution",xlab="system reliability")
plot(ecdf(simu_table_series$P1),verticals=TRUE,do.p=TRUE,main="Component1 reliability empirical distribution",xlab="component1 reliability")
plot(ecdf(simu_table_series$P2),verticals=TRUE,do.p=TRUE,main="Component2 reliability empirical distribution",xlab="component2 reliability")
plot(ecdf(simu_table_series$P3),verticals=TRUE,do.p=TRUE,main="Component3 reliability empirical distribution",xlab="component3 reliability")
plot(ecdf(simu_table_series$P4),verticals=TRUE,do.p=TRUE,main="Component4 reliability empirical distribution",xlab="component4 reliability")
plot(ecdf(simu_table_series$P5),verticals=TRUE,do.p=TRUE,main="Component5 reliability empirical distribution",xlab="component5 reliability")



#=====================================================================#
#Next explore the scenario in parallel system
#Original data table for series system
real_para<-data.frame(data_type=c("system","c1","c2","c3"),X=c(33,136,54,31),N=c(40,200,80,60))
real_para$P<-real_para$X/real_para$N

#Generate real data for simulation in the next section
set.seed(1)
sys<-sample(c(rep(1,33),rep(0,7)),40)
X1<-sample(c(rep(1,136),rep(0,64)),200)
X2<-sample(c(rep(1,54),rep(0,26)),80)
X3<-sample(c(rep(1,31),rep(0,29)),60)

#set bootstrap sample size 10000,generate 10000 regular resampling simulation data
set.seed(1)
sys_sim<-replicate(200000,simu(sys))
X1_sim<-replicate(200000,simu(X1))
X2_sim<-replicate(200000,simu(X2))
X3_sim<-replicate(200000,simu(X3))

#generate simulation table suited for regular bootstrap
simu_table_para<-data.frame(index=1:200000,system=sys_sim,component1=X1_sim,component2=X2_sim,component3=X3_sim,
                            Ps=sys_sim/length(sys),P1=X1_sim/length(X1),P2=X2_sim/length(X2),
                            P3=X3_sim/length(X3))


#Define the discrepancy measure in parallel system
discrep_para<-function(x)
{
  Ps<-x["Ps"]
  P1<-x["P1"]
  P2<-x["P2"]
  P3<-x["P3"]
  disc<-Ps/(1-(1-P1)*(1-P2)*(1-P3))
  return(disc)
}

#Generate discrepancy values for original data and simulated data
real_para_discrep<-real_para$P[1]/(1-(1-real_para$P[2])*(1-real_para$P[3])*(1-real_para$P[4]))
simu_para_discrep<-apply(simu_table_para,1,discrep_para)

#Generate histogram for discrepancy for parallel system
#Frequency histogram
par(mfrow=c(1,1))
hist(simu_para_discrep,freq=TRUE,breaks=seq(min(simu_para_discrep),max(simu_para_discrep),length.out=20),xlim=range(simu_para_discrep),main="Discrepancy measure for parallel system",xlab="Discrepancy measure")
abline(v=real_para_discrep,col="red",lty=2)
abline(v=quantf(simu_para_discrep),col="blue",lty=3)
#Density histogram
hist(simu_para_discrep,freq=FALSE,breaks=seq(min(simu_para_discrep),max(simu_para_discrep),length.out=20),xlim=range(simu_para_discrep),main="Discrepancy measure for parallel system",xlab="Discrepancy measure")
lines(density(simu_para_discrep),col="blue",lty=1)
abline(v=real_para_discrep,col="red",lty=2)
abline(v=quantf(simu_para_discrep),col="blue",lty=3)
#Kernal density plot for system discrepancy
plot(density(simu_para_discrep),type="l",main="density of discrepancy")
abline(v=real_para_discrep,col="red",lty=2)
text(x=real_para_discrep,y=0,labels=paste("real data discrep=",round(real_para_discrep,3)))
abline(v=quantf(simu_para_discrep),col="blue",lty=3)
text(x=median(simu_para_discrep),y=3.0,labels=paste("median=",round(median(simu_para_discrep),3)))
abline(v=1,col="red",lty=4)

#Generate empirical distribution for discrepancy
par(mfrow=c(1,1))
plot(ecdf(simu_para_discrep),verticals=TRUE,do.p=FALSE,main="Discrepancy measure empirical distribution for parallel system",xlab="simulation discrepancy for parallel system")
abline(v=real_para_discrep,lty=4,col="blue")

#Calculate the cumulative probability of discrepancy in real data
cum_disc_real_para<-diff(ecdf(simu_para_discrep)(c(0,real_para_discrep)))
diff(ecdf(simu_para_discrep)(c(0,1)))
#P_2sided=(1-0.986475)*2=0.02705
#Reject the Null hypothesis


#==========================================================================#
#Consider the complex system I
#Delta=Ps/{PA*[PB+PC-PB*PC]*[PD+PE-PD*PE]*PF}

#Original data table for complex system I
real_complex_I<-data.frame(data_type=c("system","A","B","C","D","E","F"),N=c(90,60,130,110,130,90,60),X=c(53,56,115,92,113,76,54))
real_complex_I$P=real_complex_I$X/real_complex_I$N

#Generate real data for simulation in the next section
set.seed(1)
sys<-sample(c(rep(1,53),rep(0,37)),90)
X1<-sample(c(rep(1,56),rep(0,4)),60)
X2<-sample(c(rep(1,115),rep(0,15)),130)
X3<-sample(c(rep(1,92),rep(0,18)),110)
X4<-sample(c(rep(1,113),rep(0,17)),130)
X5<-sample(c(rep(1,76),rep(0,14)),90)
X6<-sample(c(rep(1,54),rep(0,6)),60)

#set bootstrap sample size 10000,generate 10000 regular resampling simulation data
set.seed(1)
sys_sim<-replicate(200000,simu(sys))
X1_sim<-replicate(200000,simu(X1))
X2_sim<-replicate(200000,simu(X2))
X3_sim<-replicate(200000,simu(X3))
X4_sim<-replicate(200000,simu(X4))
X5_sim<-replicate(200000,simu(X5))
X6_sim<-replicate(200000,simu(X6))

#generate simulation table suited for regular bootstrap
simu_table_I<-data.frame(index=1:200000,system=sys_sim,component1=X1_sim,
                         component2=X2_sim,component3=X3_sim,
                         component4=X4_sim,component5=X5_sim,
                         Ps=sys_sim/length(sys),P1=X1_sim/length(X1),
                         P2=X2_sim/length(X2),P3=X3_sim/length(X3),
                         P4=X4_sim/length(X4),P5=X5_sim/length(X5),
                         P6=X6_sim/length(X6))

#Define discrepancy measure in complex system I
discrep_complex_I<-function(x)
{
  Ps<-x["Ps"]
  P1<-x["P1"]
  P2<-x["P2"]
  P3<-x["P3"]
  P4<-x["P4"]
  P5<-x["P5"]
  P6<-x["P6"]
  disc<-Ps/(P1*(P2+P3-P2*P3)*(P4+P5-P4*P5)*P6)
  return(disc)
}

#Generate discrepancy values for original data and simulated data
attach(real_complex_I)
real_I_discrep<-P[1]/(P[2]*(P[3]+P[4]-P[3]*P[4])*(P[5]+P[6]-P[5]*P[6])*P[7])
detach(real_complex_I)
simu_I_discrep<-apply(simu_table_I,1,discrep_complex_I)


#Generate histogram for discrepancy for complex system 1

#Frequency histogram
par(mfrow=c(1,1))
hist(simu_I_discrep,freq=TRUE,breaks=seq(min(simu_I_discrep),max(simu_I_discrep),length.out=20),xlim=range(simu_I_discrep),main="Discrepancy measure for complex system I",xlab="Discrepancy measure")
abline(v=real_I_discrep,col="red",lty=2)
abline(v=quantf(simu_I_discrep),col="blue",lty=3)

#Density histogram
hist(simu_I_discrep,freq=FALSE,breaks=seq(min(simu_I_discrep),max(simu_I_discrep),length.out=20),xlim=range(simu_I_discrep),main="Discrepancy measure for complex system I",xlab="Discrepancy measure")
lines(density(simu_I_discrep),col="blue",lty=1)
abline(v=real_I_discrep,col="red",lty=2)
abline(v=quantf(simu_I_discrep),col="blue",lty=3)

#Kernal density plot for system discrepancy
plot(density(simu_I_discrep),type="l",main="density of discrepancy")
abline(v=real_I_discrep,col="red",lty=2)
text(x=real_I_discrep,y=0,labels=paste("real data discrep=",round(real_I_discrep,3)))
abline(v=quantf(simu_I_discrep),col="blue",lty=3)
text(x=median(simu_I_discrep),y=3.0,labels=paste("median=",round(median(simu_I_discrep),3)))
abline(v=1,col='red',lty=4)

#Calculate the cumulative probability of discrepancy in real data for series system
cum_disc_real_I<-diff(ecdf(simu_I_discrep)(c(0,real_I_discrep)))
diff(ecdf(simu_I_discrep)(c(0,1)))
#P_2sided=(1-0.99932)*2=0.00136
#Reject the Null hypothesis



#===============================================#
#Consider the complex system II-a when Bridge C is disconnected
#Delta=Ps/{1-(1-PA*PB)*(1-PD*PE)}

#Original data table for complex system II-a
real_complex_II_a<-data.frame(data_type=c("system","A","B","C","D","E"),N=c(70,100,100,160,120,160),X=c(59,70,84,113,92,123))
real_complex_II_a$P=real_complex_II_a$X/real_complex_II_a$N

#Generate real data for simulation in the next section
set.seed(1)
sys<-sample(c(rep(1,59),rep(0,11)),70)
X1<-sample(c(rep(1,70),rep(0,30)),100)
X2<-sample(c(rep(1,84),rep(0,16)),100)
X3<-sample(c(rep(1,113),rep(0,47)),160)
X4<-sample(c(rep(1,92),rep(0,28)),120)
X5<-sample(c(rep(1,123),rep(0,37)),160)

#set bootstrap sample size 10000,generate 10000 regular resampling simulation data
set.seed(1)
sys_sim<-replicate(200000,simu(sys))
X1_sim<-replicate(200000,simu(X1))
X2_sim<-replicate(200000,simu(X2))
X3_sim<-replicate(200000,simu(X3))
X4_sim<-replicate(200000,simu(X4))
X5_sim<-replicate(200000,simu(X5))



#generate simulation table suited for regular bootstrap
simu_table_II_a<-data.frame(index=1:200000,system=sys_sim,component1=X1_sim,
                         component2=X2_sim,component3=X3_sim,
                         component4=X4_sim,component5=X5_sim,
                         Ps=sys_sim/length(sys),P1=X1_sim/length(X1),
                         P2=X2_sim/length(X2),P3=X3_sim/length(X3),
                         P4=X4_sim/length(X4),P5=X5_sim/length(X5))

#Define discrepancy measure in complex system II_a
discrep_complex_II_a<-function(x)
{
  Ps<-x["Ps"]
  P1<-x["P1"]
  P2<-x["P2"]
  P3<-x["P3"]
  P4<-x["P4"]
  P5<-x["P5"]
  disc<-Ps/(1-(1-P1*P2)*(1-P4*P5))
  return(disc)
}

#Generate discrepancy values for original data and simulated data
attach(real_complex_II_a)
real_II_a_discrep<-P[1]/(1-(1-P[2]*P[3])*(1-P[5]*P[6]))
detach(real_complex_II_a)
simu_II_a_discrep<-apply(simu_table_II_a,1,discrep_complex_II_a)


#Generate histogram for discrepancy for complex system II_a

#Frequency histogram
par(mfrow=c(1,1))
hist(simu_II_a_discrep,freq=TRUE,breaks=seq(min(simu_II_a_discrep),max(simu_II_a_discrep),length.out=20),xlim=range(simu_II_a_discrep),main="Discrepancy measure for complex system II-a",xlab="Discrepancy measure")
abline(v=real_II_a_discrep,col="red",lty=2)
abline(v=quantf(simu_II_a_discrep),col="blue",lty=3)

#Density histogram
hist(simu_II_a_discrep,freq=FALSE,breaks=seq(min(simu_II_a_discrep),max(simu_II_a_discrep),length.out=20),xlim=range(simu_II_a_discrep),main="Discrepancy measure for complex system II-a",xlab="Discrepancy measure")
lines(density(simu_II_a_discrep),col="blue",lty=1)
abline(v=real_II_a_discrep,col="red",lty=2)
abline(v=quantf(simu_II_a_discrep),col="blue",lty=3)

#Kernal density plot for system discrepancy
plot(density(simu_II_a_discrep),type="l",main="density of discrepancy")
abline(v=real_II_a_discrep,col="red",lty=2)
text(x=real_II_a_discrep,y=0,labels=paste("real data discrep=",round(real_II_a_discrep,3)))
abline(v=quantf(simu_II_a_discrep),col="blue",lty=3)
text(x=median(simu_II_a_discrep),y=3.0,labels=paste("median=",round(median(simu_II_a_discrep),3)))
abline(v=1,col="red",lty=4)

#Calculate the cumulative probability of discrepancy in real data for series system
cum_disc_real_II_a<-diff(ecdf(simu_II_a_discrep)(c(0,real_II_a_discrep)))
diff(ecdf(simu_II_a_discrep)(c(0,1)))
#P_2sided=2*0.396415=0.79283
#Accept the Null hypothesis

#====================================================#
#Consider the complex system II-b when Bridge C is connected
#Delta=Ps/{1-(1-PA*PB)*(1-PA*PE)*(1-PB*PD)*(1-PD*PE)}

#Original data table for complex system II-b
real_complex_II_b<-data.frame(data_type=c("system","A","B","C","D","E"),N=c(150,100,100,160,120,160),X=c(140,70,84,113,92,123))
real_complex_II_b$P=real_complex_II_b$X/real_complex_II_b$N

#Generate real data for simulation in the next section
set.seed(1)
sys<-sample(c(rep(1,140),rep(0,10)),150)
X1<-sample(c(rep(1,70),rep(0,30)),100)
X2<-sample(c(rep(1,84),rep(0,16)),100)
X3<-sample(c(rep(1,113),rep(0,47)),160)
X4<-sample(c(rep(1,92),rep(0,28)),120)
X5<-sample(c(rep(1,123),rep(0,37)),160)

#set bootstrap sample size 10000,generate 10000 regular resampling simulation data
set.seed(1)
sys_sim<-replicate(200000,simu(sys))
X1_sim<-replicate(200000,simu(X1))
X2_sim<-replicate(200000,simu(X2))
X3_sim<-replicate(200000,simu(X3))
X4_sim<-replicate(200000,simu(X4))
X5_sim<-replicate(200000,simu(X5))



#generate simulation table suited for regular bootstrap
simu_table_II_b<-data.frame(index=1:200000,system=sys_sim,component1=X1_sim,
                            component2=X2_sim,component3=X3_sim,
                            component4=X4_sim,component5=X5_sim,
                            Ps=sys_sim/length(sys),P1=X1_sim/length(X1),
                            P2=X2_sim/length(X2),P3=X3_sim/length(X3),
                            P4=X4_sim/length(X4),P5=X5_sim/length(X5))

#Define discrepancy measure in complex system II_b
discrep_complex_II_b<-function(x)
{
  Ps<-x["Ps"]
  P1<-x["P1"]
  P2<-x["P2"]
  P3<-x["P3"]
  P4<-x["P4"]
  P5<-x["P5"]
  disc<-Ps/(1-(1-P1*P2)*(1-P1*P5)*(1-P2*P4)*(1-P4*P5))
  return(disc)
}

#Generate discrepancy values for original data and simulated data
attach(real_complex_II_b)
real_II_b_discrep<-P[1]/(1-(1-P[2]*P[3])*(1-P[2]*P[6])*(1-P[3]*P[5])*(1-P[5]*P[6]))
detach(real_complex_II_b)
simu_II_b_discrep<-apply(simu_table_II_b,1,discrep_complex_II_b)

#Generate histogram for discrepancy for complex system II_a

#Frequency histogram
par(mfrow=c(1,1))
hist(simu_II_b_discrep,freq=TRUE,breaks=seq(min(simu_II_b_discrep),max(simu_II_b_discrep),length.out=20),xlim=range(simu_II_b_discrep),main="Discrepancy measure for complex system II-b",xlab="Discrepancy measure")
abline(v=real_II_b_discrep,col="red",lty=2)
abline(v=quantf(simu_II_b_discrep),col="blue",lty=3)

#Density histogram
hist(simu_II_b_discrep,freq=FALSE,breaks=seq(min(simu_II_b_discrep),max(simu_II_b_discrep),length.out=20),xlim=range(simu_II_a_discrep),main="Discrepancy measure for complex system II-b",xlab="Discrepancy measure")
lines(density(simu_II_b_discrep),col="blue",lty=1)
abline(v=real_II_b_discrep,col="red",lty=2)
abline(v=quantf(simu_II_b_discrep),col="blue",lty=3)

#Kernal density plot for system discrepancy
plot(density(simu_II_b_discrep),type="l",main="density of discrepancy")
abline(v=real_II_b_discrep,col="red",lty=2)
text(x=real_II_b_discrep,y=0,labels=paste("real data discrep=",round(real_II_a_discrep,3)))
abline(v=quantf(simu_II_b_discrep),col="blue",lty=3)
text(x=median(simu_II_b_discrep),y=3.0,labels=paste("median=",round(median(simu_II_b_discrep),3)))
abline(v=1,col="red",lty=4)

#Calculate the cumulative probability of discrepancy in real data for series system
cum_disc_real_II_b<-diff(ecdf(simu_II_b_discrep)(c(0,real_II_b_discrep)))
diff(ecdf(simu_II_b_discrep)(c(0,1)))
#P_2sided=2*(1-0.96404)=0.07192>0.05
#Accept the Null hypothesis if α=0.05
#Reject the Null hypothesis if α=0.1



#===========================================#
#Fractional random weights bootstrap
#Generate 10000 fractional random weights for system and components respectively from uniform dirichlet distributions
#sum of each group's weights are (sample size N)*1
library(MCMCpack)
set.seed(1)
wt_series_sys<-rdirichlet(200000,rep(1,25))*25
wt_series_c1<-rdirichlet(200000,rep(1,120))*120
wt_series_c2<-rdirichlet(200000,rep(1,100))*100
wt_series_c3<-rdirichlet(200000,rep(1,80))*80
wt_series_c4<-rdirichlet(200000,rep(1,100))*100
wt_series_c5<-rdirichlet(200000,rep(1,80))*80

dim(wt_series_sys)
dim(wt_series_c1)
dim(wt_series_c2)
dim(wt_series_c3)
dim(wt_series_c4)
dim(wt_series_c5)


real_series<-data.frame(data_type=c("system","c1","c2","c3","c4","c5"),X=c(17,111,96,76,98,79),N=c(25,120,100,80,100,80))
real_series$P=real$X/real$N

#Generate real data for simulation in the next section
set.seed(1)
sys<-matrix(sample(c(rep(1,17),rep(0,8)),25),25,1)
X1<-matrix(sample(c(rep(1,111),rep(0,9)),120),120,1)
X2<-matrix(sample(c(rep(1,96),rep(0,4)),100),100,1)
X3<-matrix(sample(c(rep(1,76),rep(0,4)),80),80,1)
X4<-matrix(sample(c(rep(1,98),rep(0,2)),100),100,1)
X5<-matrix(sample(c(rep(1,79),rep(0,1)),80),80,1)

#Generate simulated reliability for system and components by fractional random weights
simu_wt_series<-data.frame(index=1:200000,Ps=wt_series_sys%*%sys/length(sys),
                           P1=wt_series_c1%*%X1/length(X1),
                           P2=wt_series_c2%*%X2/length(X2),
                           P3=wt_series_c3%*%X3/length(X3),
                           P4=wt_series_c4%*%X4/length(X4),
                           P5=wt_series_c5%*%X5/length(X5))

#Generate discrepancy values for simulated data using fractional weight bootstrap
#We know for original data,the discrepancy is real_series_discrep
#real_series_discrep<-real_series$P[1]/prod(real_series$P[2:6])
disc_series_wt<-apply(simu_wt_series,1,discrep_series)

#Frequency histogram
par(mfrow=c(1,1))
hist(disc_series_wt,freq=TRUE,breaks=seq(min(disc_series_wt),max(disc_series_wt),length.out=20),xlim=range(disc_series_wt),main="Discrepancy measure for series system using fractional weight bootstrap",xlab="Discrepancy measure")
abline(v=real_series_discrep,col="red",lty=2)
abline(v=quantf(disc_series_wt),col="blue",lty=3)

#Density histogram
hist(disc_series_wt,freq=FALSE,breaks=seq(min(disc_series_wt),max(disc_series_wt),length.out=20),xlim=range(disc_series_wt),main="Discrepancy measure for series system using fractional weight bootstrap",xlab="Discrepancy measure")
lines(density(disc_series_wt),col="blue",lty=1)
abline(v=real_series_discrep,col="red",lty=2)
abline(v=quantf(disc_series_wt),col="blue",lty=3)

#Kernal density plot for system discrepancy
plot(density(disc_series_wt),type="l",main="density of discrepancy")
abline(v=real_series_discrep,col="red",lty=2)
text(x=real_series_discrep,y=0,labels=paste("real data discrep=",round(real_series_discrep,3)))
abline(v=quantf(disc_series_wt),col="blue",lty=3)
text(x=median(disc_series_wt),y=3.0,labels=paste("median=",round(median(disc_series_wt),3)))
abline(v=1,col="red",lty=4)

#Calculate the cumulative probability of discrepancy in real data for series system
cum_disc_real_series_wt<-diff(ecdf(disc_series_wt)(c(0,real_series_discrep)))
diff(ecdf(disc_series_wt)(c(0,1)))
#P_2sided=2*(1-0.92278)=0.15444
#Accept the Null hypothesis

#==================================================================================#
#Implement power test for series system
#Regular resampling
#Sample sizes:30、50、100
#Assuming system reliability is 0.7,component reliability 0.931

real_series_30<-data.frame(data_type=c("system","c1","c2","c3","c4","c5"),X=c(21,rep(28,5)),N=rep(30,6))
real_series_50<-data.frame(data_type=c("system","c1","c2","c3","c4","c5"),X=c(35,rep(47,5)),N=rep(50,6))
real_serie_100<-data.frame(data_type=c("system","c1","c2","c3","c4","c5"),X=c(70,rep(93,5)),N=rep(100,6))

#Generate real data for simulation in the next section
set.seed(1)
sys_30<-sample(c(rep(1,21),rep(0,9)),30)
X_30<-sample(c(rep(1,28),rep(0,2)),30)

sys_50<-sample(c(rep(1,35),rep(0,15)),50)
X_50<-sample(c(rep(1,47),rep(0,3)),50)

sys_100<-sample(c(rep(1,70),rep(0,30)),100)
X_100<-sample(c(rep(1,93),rep(0,7)),100)

#set bootstrap sample size 200000,generate 200000 regular resampling simulation data
set.seed(1)
sys_30_sim<-replicate(200000,simu(sys_30))
X_30_sim<-replicate(200000,simu(X_30))

sys_50_sim<-replicate(200000,simu(sys_50))
X_50_sim<-replicate(200000,simu(X_50))

sys_100_sim<-replicate(200000,simu(sys_100))
X_100_sim<-replicate(200000,simu(X_100))


simu_table_series_30<-data.frame(index=1:200000,system=sys_30_sim,component=X_30_sim,
                                 Ps=sys_30_sim/length(sys_30),Pc=X_30_sim/length(X_30))
simu_table_series_50<-data.frame(index=1:200000,system=sys_50_sim,component=X_50_sim,
                                 Ps=sys_50_sim/length(sys_50),Pc=X_50_sim/length(X_50))
simu_table_series_100<-data.frame(index=1:200000,system=sys_100_sim,component=X_100_sim,
                                 Ps=sys_100_sim/length(sys_100),Pc=X_100_sim/length(X_100))
#Generate discrepany values of different sample sizes for comparing the power later
simu_series_discrep_30<-simu_table_series_30$Ps/((simu_table_series_30$Pc)^5)
simu_series_discrep_50<-simu_table_series_50$Ps/((simu_table_series_50$Pc)^5)
simu_series_discrep_100<-simu_table_series_100$Ps/((simu_table_series_100$Pc)^5)

plot(density(simu_series_discrep_100),type="l",main="density of discrepancy")
text(x=real_series_discrep,y=0,labels=paste("real data discrep=",round(real_series_discrep,3)))
abline(v=quantf(simu_series_discrep_100),col="blue",lty=3)

#Power comparation
c1_30<-quantile(simu_series_discrep_30,0.025)
c2_30<-quantile(simu_series_discrep_30,0.975)
cm_30<-quantile(simu_series_discrep_30,0.5)
step_size<-seq(-2,2,0.01)
power_30<-as.numeric()

#Generate power values when sample size=30
for(i in 1:length(step_size))
{
  power_30[i]<-1-diff(ecdf(simu_series_discrep_30+step_size[i])(c(c1_30,c2_30)))
}

#Find 3 quantiles when sample size=50
c1_50<-quantile(simu_series_discrep_50,0.025)
c2_50<-quantile(simu_series_discrep_50,0.975)
cm_50<-quantile(simu_series_discrep_50,0.5)
step_size<-seq(-2,2,0.01)
power_50<-as.numeric()

#Generate power values for sample size=30
for(i in 1:length(step_size))
{
  power_50[i]<-1-diff(ecdf(simu_series_discrep_30+step_size[i])(c(c1_50,c2_50)))
}
c1_100<-quantile(simu_series_discrep_100,0.025)
c2_100<-quantile(simu_series_discrep_100,0.975)
cm_100<-quantile(simu_series_discrep_100,0.5)
step_size<-seq(-2,2,0.01)
power_100<-as.numeric()

#Generate power values for sample size=100
for(i in 1:length(step_size))
{
  power_100[i]<-1-diff(ecdf(simu_series_discrep_100+step_size[i])(c(c1_100,c2_100)))
}


pw_discrep<-c(cm_30+step_size,cm_50+step_size,cm_100+step_size)
pw_2_side<-data.frame(discrepancy=pw_discrep,power=c(power_30,power_50,power_100),sample_size=as.factor(c(rep(30,401),rep(50,401),rep(100,401))))
#Generate power curve for 3 different sample sizes
ggplot(pw_2_side)+
  geom_line(mapping=aes(x=discrepancy,y=power,group=sample_size,color=sample_size))






#===============================================================#
#Use random weight bootstrap to implement power test
#3 sample sizes:30 50 100
library(MCMCpack)
set.seed(1)
pw_wt_series_sys_30<-rdirichlet(200000,rep(1,30))*30
pw_wt_series_sys_50<-rdirichlet(200000,rep(1,50))*50
pw_wt_series_sys_100<-rdirichlet(200000,rep(1,100))*100
set.seed(2)
pw_wt_series_comp_30<-rdirichlet(200000,rep(1,30))*30
pw_wt_series_comp_50<-rdirichlet(200000,rep(1,50))*50
pw_wt_series_comp_100<-rdirichlet(200000,rep(1,100))*100

#Use previous original data  
#set.seed(1)
#sys_30<-sample(c(rep(1,21),rep(0,9)),30)
#X_30<-sample(c(rep(1,28),rep(0,2)),30)

#sys_50<-sample(c(rep(1,35),rep(0,15)),50)
#X_50<-sample(c(rep(1,47),rep(0,3)),50)

#sys_100<-sample(c(rep(1,70),rep(0,30)),100)
#X_100<-sample(c(rep(1,93),rep(0,7)),100)

pw_wt_series_30<-data.frame(index=1:200000,Ps=pw_wt_series_sys_30%*%sys_30/length(sys_30),
                           Pc=pw_wt_series_comp_30%*%X_30/length(X_30))
pw_wt_series_50<-data.frame(index=1:200000,Ps=pw_wt_series_sys_50%*%sys_50/length(sys_50),
                            Pc=pw_wt_series_comp_50%*%X_50/length(X_50))
pw_wt_series_100<-data.frame(index=1:200000,Ps=pw_wt_series_sys_100%*%sys_100/length(sys_100),
                            Pc=pw_wt_series_comp_100%*%X_100/length(X_100))

pw_series_discrep_30<-pw_wt_series_30$Ps/((pw_wt_series_30$Pc)^5)
pw_series_discrep_50<-pw_wt_series_50$Ps/((pw_wt_series_50$Pc)^5)
pw_series_discrep_100<-pw_wt_series_100$Ps/((pw_wt_series_100$Pc)^5)

plot(density(pw_series_discrep_50),type="l",main="density of discrepancy")


#Define α=0.05
pw_c1_30_1<-quantile(pw_series_discrep_30,0.025)
pw_c2_30_1<-quantile(pw_series_discrep_30,0.975)
pw_cm_30_1<-quantile(pw_series_discrep_30,0.5)
pw_c1_50_1<-quantile(pw_series_discrep_50,0.025)
pw_c2_50_1<-quantile(pw_series_discrep_50,0.975)
pw_cm_50_1<-quantile(pw_series_discrep_50,0.5)
pw_c1_100_1<-quantile(pw_series_discrep_100,0.025)
pw_c2_100_1<-quantile(pw_series_discrep_100,0.975)
pw_cm_100_1<-quantile(pw_series_discrep_100,0.5)
step_size<-seq(-2,2,0.01)
pw_30_1<-as.numeric()
pw_50_1<-as.numeric()
pw_100_1<-as.numeric()

#Generate power values when sample size=30
for(i in 1:length(step_size))
{
  pw_30_1[i]<-1-diff(ecdf(pw_series_discrep_30+step_size[i])(c(pw_c1_30,pw_c2_30)))
  pw_50_1[i]<-1-diff(ecdf(pw_series_discrep_50+step_size[i])(c(pw_c1_50,pw_c2_50)))
  pw_100_1[i]<-1-diff(ecdf(pw_series_discrep_100+step_size[i])(c(pw_c1_100,pw_c2_100)))
}


#Generate power plot
pw_wt_discrep_1<-c(pw_cm_30+step_size,pw_cm_50+step_size,pw_cm_100+step_size)
pw_wt_2_side_1<-data.frame(discrepancy=pw_wt_discrep,power=c(pw_30,pw_50,pw_100),sample_size=as.factor(c(rep(30,401),rep(50,401),rep(100,401))))
#Generate power curve for 3 different sample sizes
ggplot(pw_wt_2_side_1)+
  geom_line(mapping=aes(x=discrepancy,y=power,group=sample_size,color=sample_size))+
  labs(title="power curve at α=0.05")



#Define α=0.1
pw_c1_30_2<-quantile(pw_series_discrep_30,0.05)
pw_c2_30_2<-quantile(pw_series_discrep_30,0.95)
pw_cm_30_2<-quantile(pw_series_discrep_30,0.5)
pw_c1_50_2<-quantile(pw_series_discrep_50,0.05)
pw_c2_50_2<-quantile(pw_series_discrep_50,0.95)
pw_cm_50_2<-quantile(pw_series_discrep_50,0.5)
pw_c1_100_2<-quantile(pw_series_discrep_100,0.05)
pw_c2_100_2<-quantile(pw_series_discrep_100,0.95)
pw_cm_100_2<-quantile(pw_series_discrep_100,0.5)
step_size<-seq(-2,2,0.01)
pw_30_2<-as.numeric()
pw_50_2<-as.numeric()
pw_100_2<-as.numeric()

#Generate power values when sample size=30
for(i in 1:length(step_size))
{
  pw_30_2[i]<-1-diff(ecdf(pw_series_discrep_30+step_size[i])(c(pw_c1_30,pw_c2_30)))
  pw_50_2[i]<-1-diff(ecdf(pw_series_discrep_50+step_size[i])(c(pw_c1_50,pw_c2_50)))
  pw_100_2[i]<-1-diff(ecdf(pw_series_discrep_100+step_size[i])(c(pw_c1_100,pw_c2_100)))
}


#Generate power plot
pw_wt_discrep_2<-c(pw_cm_30+step_size,pw_cm_50+step_size,pw_cm_100+step_size)
pw_wt_2_side_2<-data.frame(discrepancy=pw_wt_discrep,power=c(pw_30,pw_50,pw_100),sample_size=as.factor(c(rep(30,401),rep(50,401),rep(100,401))))
#Generate power curve for 3 different sample sizes
ggplot(pw_wt_2_side_2)+
  geom_line(mapping=aes(x=discrepancy,y=power,group=sample_size,color=sample_size))+
  labs(title="power curve at α=0.1")

  
  





  




