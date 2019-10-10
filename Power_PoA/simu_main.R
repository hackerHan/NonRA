library(MCMCpack)
library(plotly)
library(ggplot2)
Ns<-seq(20,100,by=10) #System sample size
D1<-seq(0.5,1.5,by=0.1) #True discrepancy
Pc<-0.921

grid<-expand.grid(Ns,D1)
colnames(grid)<-c("sample size","True discrepancy")
grid$pr0.05<-0 #power at alpha 0.05
grid$pr0.1<-0 #power at alpha 0.1
grid$cv.d1<-0 #Relative half Width of credible interval
grid$cv.d1r<-0 #Half Width of Credible interval
grid$Ave.PoA<-0 #Average PoA
grid$PoAl<-0 #Lower PoA
grid$PoAu<-0 #Upper PoA

time1<-Sys.time()
MC_results<-apply(grid,1,simu_pw)
time2<-Sys.time()
time2-time1

grid$pr0.05<-MC_results[1,] 
grid$pr0.1<-MC_results[2,]
grid$cv.d1<-MC_results[3,]
grid$cv.d1r<-MC_results[4,]
grid$Ave.PoA<-MC_results[5,]
grid$PoAl<-MC_results[6,]
grid$PoAu<-MC_results[7,]
save(grid,file=paste("power",".rdata",sep=""))

# contour plot
attach(grid)
# contour plot of power
plot_ly(type="contour",x=`True discrepancy`,y=`sample size`,
                z=~pr0.05,colors=c("red","green"),
        contours=list(coloring="heatmap",showlabels=TRUE),
        hoverinfo="all")%>%
  layout(title="Power contour at alpha 0.05",xaxis=list(title="Discrepancy"),
         yaxis=list(title="System sample size"))
plot_ly(type="contour",x=`True discrepancy`,y=`sample size`,
        z=~pr0.1,colors=c("red","green"),
        contours=list(coloring="heatmap",showlabels=TRUE),
        hoverinfo="all",name="Power contour at alpha 0.1")%>%
  layout(title="Power contour at alpha 0.1",xaxis=list(title="Discrepancy"),
         yaxis=list(title="System sample size"))

#sliced profile plot
grid_D0.6<-grid[which(grid$`True discrepancy`==0.6),]
ggplot(data=grid_D0.6)+
  geom_point(mapping=aes(x=`sample size`,y=pr0.05))+
  geom_smooth(mapping=aes(y=pr0.05,x=`sample size`))+
  geom_hline(mapping=aes(yintercept=0.9),color="red")+
  ggtitle("Power at Delta=0.6,α=0.05")
ggplot(data=grid_D0.6)+
  geom_point(mapping=aes(x=`sample size`,y=pr0.1))+
  geom_smooth(mapping=aes(y=pr0.1,x=`sample size`))+
  geom_hline(mapping=aes(yintercept=0.9),color="red")+
  ggtitle("Power at Delta=0.6,α=0.1")

#CV
plot_ly(type="contour",x=`True discrepancy`,y=`sample size`,
        z=~cv.d1,colors=c("red","green"),
        contours=list(coloring="heatmap",showlabels=TRUE),
        hoverinfo="all")%>%
  layout(title="Relative Half width of Credible Interval",
         xaxis=list(title="Discrepancy",yaxis=list(title="System sample size")))
plot_ly(type="contour",x=`True discrepancy`,y=`sample size`,
        z=~cv.d1r,colors=c("red","green"),
        contours=list(coloring="heatmap",showlabels=TRUE),
        hoverinfo="all")%>%
  layout(title="Half width of Credible Interval",xaxis=list(title="Discrepancy"),
         yaxis=list(title="System sample size"))
#Average PoA
plot_ly(type="contour",x=`True discrepancy`,y=`sample size`,
        z=~Ave.PoA,colors=c("red","green"),
        contours=list(coloring="heatmap",showlabels=TRUE),
        hoverinfo="all")%>%
  layout(title="Bayesian probability of agreement on average",
         xaxis=list(title="Discrepancy",yaxis=list(title="System sample size")))
#lower and upper bound of PoA
plot_ly(type="contour",x=`True discrepancy`,y=`sample size`,
        z=~PoAl,colors=c("red","green"),
        contours=list(coloring="heatmap",showlabels=TRUE),
        hoverinfo="all")%>%
  layout(title="lower bound of Bayesian probability of agreement",
         xaxis=list(title="Discrepancy",yaxis=list(title="System sample size")))
plot_ly(type="contour",x=`True discrepancy`,y=`sample size`,
        z=~PoAu,colors=c("red","green"),
        contours=list(coloring="heatmap",showlabels=TRUE),
        hoverinfo="all")%>%
  layout(title="upper bound of Bayesian probability of agreement",
         xaxis=list(title="Discrepancy",yaxis=list(title="System sample size")))

detach(grid)


