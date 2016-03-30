library(MARSS)
library(reshape2)
#####Get Secchi TS Data#######
secchi <- read.delim("median_annual_secchi_22yrs.txt")
secchi = secchi[,c("lagoslakeid","Year","median_secchi")]

#####Get Cluster Data#########
cluster_pred <- read.csv("cluster_pred.csv", stringsAsFactors=FALSE)

#Transpose Data to wide formate
names(secchi) = c("id","variable","value")
testm = melt(secchi,id.vars = c("id","variable"),measure.vars = "value")
testw = dcast(secchi,id~variable)
testw = merge(cluster_pred[,c(1,2)],testw,by.x="Lagoslakeid",by.y="id")

#Merge in Clusters for plotting each cluster
clusters = sort(unique(testw$cluster))

#Standardize the data
for(i in 1:nrow(testw)){
  datavec = scale(as.numeric(testw[i,3:27]),center = TRUE,scale = TRUE)
  testw[i,3:27] = datavec
}

#Plot individual lake trends with common trend over the top
yrange = range(testw[,c(3:27)],na.rm=TRUE)

dev.new(width=173/25.4,height=173/25.4)
par(mfrow=c(4,2),oma=c(.5,1.5,.2,.2),mar=c(1,1,0,0) )
for(i in 1:length(clusters)){
  #Code to run DFA model for each cluster
   temp = data.matrix(testw[which(testw$cluster==clusters[i]),3:27])
#   model.list = list(m=1,R="diagonal and unequal")
#   cntr.list = list(maxit=10000, safe=TRUE)
#   dfa1 = MARSS(temp,model=model.list,z.score=FALSE,form="dfa",control=cntr.list)
  #Plotting all the data
  for (n in 1:nrow(temp)) {
    if(i == 1 | i == 3 | i==5) {
      if (n==1) {
        plot(x=c(1987:2011),y = temp[n,],ylim=yrange,type="l",col="grey",ylab="Annual Precip",xlab="",xaxt="n",yaxt="n")
        title(main=paste0("Cluster ",clusters[i]),line=-1)
        axis(side=1,labels = FALSE,tck=-0.02)
        axis(side=2,labels = FALSE,tck=-0.02)
        axis(side=2,line = -.75,lwd=0)
        box(lwd=1)
        if (i==3) mtext(side=2,"Standardized Secchi Depth",line=1.2,adj=4.5)
      }
    } else if (i == 7 & n == 1) {
        plot(x=c(1987:2011),y = temp[n,],ylim=yrange,type="l",col="grey",ylab="Annual Precip",xlab="",xaxt="n",yaxt="n")
        title(main=paste0("Cluster ",clusters[i]),line=-1)
        axis(side=1,labels = FALSE,tck=-0.02)
        axis(side=2,labels = FALSE,tck=-0.02)
        axis(side=1,line = -.75,lwd=0)
        axis(side=2,line = -.75,lwd=0)
        box(lwd=1)
    } else if (i == 8 & n == 1)  {
      plot(x=c(1987:2011),y = temp[n,],ylim=yrange,type="l",col="grey",ylab="Annual Precip",xlab="",xaxt="n",yaxt="n")
      title(main=paste0("Cluster ",clusters[i]),line=-1)
      axis(side=1,labels = FALSE,tck=-0.02)
      axis(side=2,labels = FALSE,tck=-0.02)
      axis(side=1,line = -.75,lwd=0)
      box(lwd=1)  
      } else if (n==1) {
      plot(x=c(1987:2011),y = temp[n,],ylim=yrange,type="l",col="grey",ylab="Annual Precip",xlab="",xaxt="n",yaxt="n")
      title(main=paste0("Cluster ",clusters[i]),line=-1)
      axis(side=1,labels = FALSE,tck=-0.02)
      axis(side=2,labels = FALSE,tck=-0.02)
      }
    if(n > 1) lines(x=c(1987:2011),y = temp[n,],col="grey")
  }
  #Code to plot the common DFA trend line
#   lines(x=c(1987:2011),y=as.numeric(dfa1$states),col="black",lwd=2)
#   lines(x=c(1987:2011),y=(as.numeric(dfa1$states)+2*as.numeric(dfa1$states.se)),col="black",lwd=2,lty=2)
#   lines(x=c(1987:2011),y=(as.numeric(dfa1$states)-2*as.numeric(dfa1$states.se)),col="black",lwd=2,lty=2)
}
