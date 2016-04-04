#########################
######Packages
#########################

library(maps)
library(maptools)
options(device="quartz")
library(R2WinBUGS)

#########################
######Functions
#########################

northarrow <- function(loc,size,bearing=0,cols,cex=1,...) {
  # checking arguments
  if(missing(loc))  stop("loc is missing")
  if(missing(size))  stop("size is missing")
  # default colors are white and black
  if(missing(cols)) cols <- rep(c("white","black"),8)
  # calculating coordinates of polygons
  radii <- rep(size/c(1,4,2,4),4)
  x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
  y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
  # drawing polygons
  for (i in 1:15) {
    x1 <- c(x[i],x[i+1],loc[1])
    y1 <- c(y[i],y[i+1],loc[2])
    polygon(x1,y1,col=cols[i])
  }
  # drawing the last polygon
  polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
  # drawing letters
  b <- c("E","N","W","S")
  for (i in 0:3) text((.3+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],(.3+par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
                      cex=cex,col="black")
}

#########################
######Get Data
#########################

# cluster_pred <- read.csv("~/Dropbox (Lottig Family)/CSI_LIMNO/CSI_LIMNO_Manuscripts-presentations/CSI_Secchi trends paper/cluster_pred.csv", stringsAsFactors=FALSE)
# lagos_lakes_10541 <- read.delim("~/Dropbox (Lottig Family)/CSI_LIMNO/CSI-LIMNO_DATA/LAGOSData/Version1.054.1/lagos_lakes_10541.txt", stringsAsFactors=FALSE)
# data = merge(cluster_pred,lagos_lakes_10541[,c(1,3,4)],by.x="Lagoslakeid",by.y="lagoslakeid")
# 
# ######Subset all data to just WISCONSIN
# WI_data = data
# WI_data = merge(WI_data,lagos_lakes_10541[,c("lagoslakeid","state_name")],by.x="Lagoslakeid",by.y="lagoslakeid")
# WI_data = WI_data[which(WI_data$state_name=="Wisconsin"),]
# 
# #######Look a SI Coefs to figure out what lakes to include
# # SecchiOut <- read.csv("~/Dropbox (Lottig Family)/CSI_LIMNO/CSI_LIMNO_Manuscripts-presentations/CSI_Secchi trends paper/Data/SecchiOut.csv",stringsAsFactors = FALSE)
# # SecchiOut = SecchiOut[,c("Lagoslakeid","cluster","SI_coef")]
# # SecchiOut = SecchiOut[which(SecchiOut$Lagoslakeid %in% WI_data$Lagoslakeid),]
# # SecchiOut = SecchiOut[order(SecchiOut$cluster,SecchiOut$SI_coef,decreasing = TRUE),]
# # best.lakes = data.frame(cluster=c(8,7,6,5,4,3,2,1), LagoslakeID=c(802,858,5215,5033,1667,1041,5501,4733))
# 
# ##########Get Secchi Raw Data
# secchi <- read.delim("~/Dropbox/CSI_LIMNO/CSI_LIMNO_Manuscripts-presentations/CSI_Secchi trends paper/Matlab code/data/median_annual_secchi_22yrs.txt")

# map.1 = WI_data[which(WI_data$Lagoslakeid==802),]
# points(map.1$nhd_long,map.1$nhd_lat)

# Load Noah's workspace
load('Figure_Data_TY.RData')
head(secchi)
secchi <- secchi[,c(1,2,3)]
dat <- secchi[secchi$lagoslakeid==802 |secchi$lagoslakeid==4733 | secchi$lagoslakeid==5501 | secchi$lagoslakeid==858, ]
unique(dat$lagoslakeid)
head(dat)
dim(dat)
####### Fill in missing years for all lakes with NA
dat2 <- expand.grid(lagoslakeid = unique(dat$lagoslakeid), Year = min(dat$Year):max(dat$Year))
head(dat2)
dim(dat2)

# dat3 has missing data for years with no data for each lake
dat3 <- merge(dat,dat2,all=TRUE)
dim(dat3)

dat <- dat3

#################################################################
########## BUGS CODE ############################################
#################################################################

# Define the model in the BUGS language and write a text file
sink("model2.txt")
cat("
    model {
    # Area-specific Model 2
    for (i in 1:N) {
      for (t in 1:T) {
        y2[i,t] ~ dnorm(temp2[i,t], tau2)
        temp2[i,t] <- u[i] + xi[i,t]
        specific[i,t] <- cut(temp2[i,t])
    }
    # area-specific trends
    xi[i,1:T] ~ car.normal(adj.tm[],weights.tm[],num.tm[],prec.xi[i])
    # area-specific intercepts (no smoothing)
    u[i] ~ dnorm(0,0.001)
    # hierarchical modelling of the local temporal variability
    prec.xi[i] <- pow(var.xi[i],-1)
    var.xi[i] <- exp(log.var.xi[i])
    log.var.xi[i] ~ dnorm(mean.log.var.xi,prec.log.var.xi)
    sigma.xi[i] <- pow(var.xi[i],0.5)
    }
    # hyper priors
    tau2 <- pow(sigma2,-2)
    sigma2 ~ dunif(0, 100)
    
    mean.log.var.xi ~ dnorm(0,0.001)
    prec.log.var.xi <- pow(var.log.var.xi,-1)
    var.log.var.xi <- pow(sd.log.var.xi,2)
    sd.log.var.xi ~ dunif(0,5)
    #sd.log.var.xi ~ dnorm(0,prec.sd.log.var.xi)I(0,)
    # sd.sd.log.var.xi <- 2.5
    #prec.sd.log.var.xi <- pow(sd.sd.log.var.xi,-2)
    
    
    # Specify weight matrix and adjacency matrix corresponding to RW(1) prior 
    # (Note - this could be given in the data file instead)
    
    for(t in 1:1) {
    weights.tm[t] <- 1;
    adj.tm[t] <- t+1;
    num.tm[t] <- 1
    }
    for(t in 2:(T-1)) {
    weights.tm[2+(t-2)*2] <- 1;
    adj.tm[2+(t-2)*2] <- t-1
    weights.tm[3+(t-2)*2] <- 1;
    adj.tm[3+(t-2)*2] <- t+1;
    num.tm[t] <- 2
    }
    for(t in T:T) {
    weights.tm[(T-2)*2 + 2] <- 1;
    adj.tm[(T-2)*2 + 2] <- t-1;
    num.tm[t] <- 1
    }
    
    
    
    } # end model
    ",fill = TRUE)
sink()


# Number of lakes
N <- length(unique(dat$lagoslakeid))
N

T <- length(unique(dat$Year))
T
head(dat)

N*T

# Response variable
y <- dat$median_secchi
length(y)


# Page 103 - Bayesian modeling using WinBUGS
Y <- matrix(y, c(N,T), byrow=T)
#write.csv(Y, 'Y.csv')

y.s <- structure(
  .Data=c(Y),
  .Dim=c(N,T)
)



# Load data
data <- list(y2=y.s, T = T, N = N)


# Initial values
inits <- function (){
  list (sigma2=runif(1) )
}


# Parameters monitored
parameters <- c('u','xi','sigma2','sigma.xi')


# MCMC settings
ni <- 15000
nt <- 3
nb <- 10000
nc <- 3



bugs.dir <- "C:/Program Files/WinBUGS14/"

start.time = Sys.time()         # Set timer 
# Call BUGS from R 

out <- bugs(data = data, inits = inits, parameters.to.save = parameters, 
            model.file = "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb,debug = F, bugs.directory=bugs.dir)

# 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Summarize posteriors
print(out, dig = 3)

# meanTrend <- array(NA,c(out[[1]]$n.sim,T) )
# ciTrend <- array(NA,c(out[[1]]$n.sim,T) )

SiteTrend <- array(NA,c(out$n.sim,T,N) )
dim(SiteTrend)
  
  
  # Summarize MCMC output for plotting
for(t in 1:T){
  for(n in 1:N){
    SiteTrend[,t,n] <- out$sims.list$xi[,n,t] + out$sims.list$u[,n]
  }
}	

SiteTrend2 <- matrix(NA, nrow=T, ncol=N)
for(t in 1:T){
  for(n in 1:N){
    SiteTrend2[t,n] <- mean(SiteTrend[,t,n])
  }
}	

head(SiteTrend2)
dim(SiteTrend2)

# Site=specific trend
SiteTrendLCI <- matrix(NA, nrow=T, ncol=N)
for(t in 1:T){
  for(n in 1:N){
    SiteTrendLCI[t,n] <- quantile(SiteTrend[,t,n],0.025)
  }
}	

SiteTrendUCI <- matrix(NA, nrow=T, ncol=N)
for(t in 1:T){
  for(n in 1:N){
    SiteTrendUCI[t,n] <- quantile(SiteTrend[,t,n],0.975)
  }
}	



#Setup the figure 173mm (W) x 120mm (H) 
dev.new(width=173/25.4,height=120/25.4)
m=cbind(c(1,2,3,4),c(5,5,5,5))
# par(oma=c(.8,.3,.2,.2),family="Arial",ps=10)
par(oma=c(.8,.3,.2,.2),ps=10)
layout(m,widths=c(2,2),heights=c(1,1,1,1))
par(mar=c(.5,2,0,0))

##### 4 trend figures...pull out and plot the raw data for each one (lake.data=subsetted by lake of interest from raw secchi data)
lake.data = secchi[which(secchi$lagoslakeid==802),]
plot(lake.data$Year,lake.data$median_secchi,type="p",xlab="",ylab="",xaxt="n",yaxt="n", pch=16)
points(c(1987:2011), SiteTrend2[,1], lty=1, lwd=2, type='l') 
points(c(1987:2011), SiteTrendLCI[,1], lty=2, lwd=2, type='l') 
points(c(1987:2011), SiteTrendUCI[,1], lty=2, lwd=2, type='l') 
title(main="Cluster 8",line=-1,adj=.2)
axis(side=1,labels = FALSE,tck=-0.03)
axis(side=2,labels = FALSE,tck=-0.03)
axis(side=2,line = -.75,lwd=0)
box(lwd=1)

lake.data = secchi[which(secchi$lagoslakeid==4733),]
plot(lake.data$Year,lake.data$median_secchi,type="p",xlab="",ylab="",xaxt="n",yaxt="n", pch=16)
title(main="Cluster 1",line=-1,adj=.8)
points(c(1987:2011), SiteTrend2[,3], lty=1, lwd=2, type='l') 
points(c(1987:2011), SiteTrendLCI[,3], lty=2, lwd=2, type='l') 
points(c(1987:2011), SiteTrendUCI[,3], lty=2, lwd=2, type='l') 
axis(side=1,labels = FALSE,tck=-0.03)
axis(side=2,labels = FALSE,tck=-0.03)
axis(side=2,line = -.75,lwd=0)
box(lwd=1)
mtext(side=2,"Secchi Depth (m)",line=1.2,at=1.8)

lake.data = secchi[which(secchi$lagoslakeid==5501),]
plot(lake.data$Year,lake.data$median_secchi,type="p",xlab="",ylab="",xaxt="n",yaxt="n",pch=16)
points(c(1987:2011), SiteTrend2[,4], lty=1, lwd=2, type='l') 
points(c(1987:2011), SiteTrendLCI[,4], lty=2, lwd=2, type='l') 
points(c(1987:2011), SiteTrendUCI[,4], lty=2, lwd=2, type='l') 
title(main="Cluster 2",line=-1,adj=.2)
axis(side=1,labels = FALSE,tck=-0.03)
axis(side=2,labels = FALSE,tck=-0.03)
axis(side=2,line = -.75,lwd=0)
box(lwd=1)

lake.data = secchi[which(secchi$lagoslakeid==858),]
plot(lake.data$Year,lake.data$median_secchi,type="p",xlab="",ylab="",xaxt="n",yaxt="n",pch=16)
points(c(1987:2011), SiteTrend2[,2], lty=1, lwd=2, type='l') 
points(c(1987:2011), SiteTrendLCI[,2], lty=2, lwd=2, type='l') 
points(c(1987:2011), SiteTrendUCI[,2], lty=2, lwd=2, type='l') 
title(main="Cluster 7",line=-1,adj=.2)
axis(side=1,labels = FALSE,tck=-0.03)
axis(side=2,labels = FALSE,tck=-0.03)
axis(side=1,line = -.75,lwd=0)
axis(side=2,line = -.75,lwd=0)
box(lwd=1)


#Generate left panel of plot with all lakes in WI plotted by cluster ID
map('county',region="Wisconsin",col=c('grey90'),fill=TRUE,resolution = 0,mar=c(0,0,0,0),border="grey50")
clust=1
points(WI_data$nhd_long[which(WI_data$cluster==clust)],WI_data$nhd_lat[which(WI_data$cluster==clust)],pch=15,col=rgb(1,0,0,.65))
clust=2
points(WI_data$nhd_long[which(WI_data$cluster==clust)],WI_data$nhd_lat[which(WI_data$cluster==clust)],pch=16,col=rgb(0,1,0,.65))
clust=3
points(WI_data$nhd_long[which(WI_data$cluster==clust)],WI_data$nhd_lat[which(WI_data$cluster==clust)],pch=17,col=rgb((255/255),(128/255),(0/255),.65))
clust=4
points(WI_data$nhd_long[which(WI_data$cluster==clust)],WI_data$nhd_lat[which(WI_data$cluster==clust)],pch=18,col=rgb((102/255),(0/255),(204/255),.65),cex=1.2)
clust=5
points(WI_data$nhd_long[which(WI_data$cluster==clust)],WI_data$nhd_lat[which(WI_data$cluster==clust)],pch=15,col=rgb(0,0,(204/255),.65))
clust=6
points(WI_data$nhd_long[which(WI_data$cluster==clust)],WI_data$nhd_lat[which(WI_data$cluster==clust)],pch=16,col=rgb(0,0,0,.65))
clust=7
points(WI_data$nhd_long[which(WI_data$cluster==clust)],WI_data$nhd_lat[which(WI_data$cluster==clust)],pch=17,col=rgb((0/255),(255/255),(255/255),.65))
clust=8
points(WI_data$nhd_long[which(WI_data$cluster==clust)],WI_data$nhd_lat[which(WI_data$cluster==clust)],pch=18,col=rgb((255/255),(0/255),(127/255),.65),cex=1.2)

legend('bottomleft',legend = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6","Cluster 7","Cluster 8"),cex=0.75,bty="n",pch=c(15,16,17,18,15,16,17,18),
       col = c(rgb(1,0,0,1),rgb(0,1,0,1),rgb((255/255),(128/255),(0/255),1),rgb((102/255),(0/255),(204/255),1),rgb(0,0,(204/255),1),
               rgb(0,0,0,1),rgb((0/255),(255/255),(255/255),1),rgb((255/255),(0/255),(127/255),1)),pt.cex = c(.75,.75,.75,1,.75,.75,.75,1))
map.scale(grconvertX(0.20,"npc"), grconvertY(.08, "npc"),col="black", metric = TRUE, ratio=FALSE, relwidth=0.075,cex=0.5)
northarrow(loc = c(-87.6,46.3),size = .25,cex = 0.75)

#######Add Arrows to plot##########
arrows(-95,41,-89.71,43.71,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05)
arrows(-95,43.5,-88.15,45.10,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05)
arrows(-95,45.4,-88.75,45.72,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05)
arrows(-95,47.4,-89.2,45.92,xpd=TRUE,col=rgb(0,0,0,.5),code=2,length=0.05)

# save.image()


