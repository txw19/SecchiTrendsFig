
library(reshape2)
library(R2WinBUGS)
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

# Convert to matrix
testwm <- as.matrix(testw)


#################################################################
########## BUGS CODE ############################################
#################################################################

# Define the model in the BUGS language and write a text file
sink("model.txt")
cat("
    model {
    # Common-trend Model
    for (i in 1:N) {
       for (t in 1:T) {
          y1[i,t] ~ dnorm(temp1[i,t], tau1)
            # temp1[i,t] <- alpha0 + eta[i] + gamma[t]
            temp1[i,t] <- alpha0 + gamma[t]
    }
      # eta[i] ~ dnorm(0,prec.eta)
    }
    
    # prior specifications for Model 1
    tau1 <- pow(sigma1,-2)
    sigma1 ~ dunif(0, 100)
    alpha0 ~ dnorm(0, 0.001)
    gamma[1:T] ~ car.normal(adj.tm[],weights.tm[],num.tm[],prec.gamma)
    prec.gamma <- pow(sigma.gamma,-2)
    sigma.gamma ~ dunif(0,10)
    # prec.eta <- pow(sigma.eta,-2)
    # sigma.eta ~ dunif(0,10)
    
    
    
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



# Loop through cluster and fit model
out <- list()

for(i in 1:8){
  # Select out a cluster
bayes <- testw[which(testw$cluster==i),3:27]
bayes <- as.matrix(bayes)
# Remove col names
colnames(bayes) <- NULL

# Number of lakes
N <- dim(bayes)[1]
N

# Number of years
T <- dim(bayes)[2]
T

# Load data
data <- list(y1 = bayes, T = T, N = N)


# Initial values
inits <- function (){
  list (alpha0 = rnorm(1), gamma=rnorm(T) )
}


# Parameters monitored
parameters <- c("alpha0",'gamma','sigma1')


# MCMC settings
ni <- 8000
nt <- 3
nb <- 5000
nc <- 3


bugs.dir <- "C:/Program Files/WinBUGS14/"

# start.time = Sys.time()         # Set timer 
# # Call BUGS from R 

out[[i]] <- bugs(data = data, inits = inits, parameters.to.save = parameters, 
            model.file = "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb,debug = F, bugs.directory=bugs.dir)

# # 
# end.time = Sys.time()
# elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
# cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# # Calculate computation time


} # End model loop

# Summarize posteriors
# print(out[[1]], dig = 3)

# Save BUGS output
saveRDS(out, file="secchi.bugs.rds")

meanTrend <- list()
ciTrend <- list()

for(i in 1:8){
FixedSiteTrend <- array(NA,c(out[[1]]$n.sim,T) )
dim(FixedSiteTrend)


# Summarize MCMC output for plotting
for(t in 1:T){
  FixedSiteTrend[,t] <- out[[i]]$sims.list$alpha0 + out[[i]]$sims.list$gamma[,t] 
}

dim(FixedSiteTrend)

# Posterior mean
meanTrend[[i]] <- apply(FixedSiteTrend, 2, mean)

# CIs
ciTrend[[i]] <- apply(FixedSiteTrend, 2, quantile, c(0.975, 0.025))
} # End MCMC processing loop

##############
######### PLOT 
##############
#Plot individual lake trends with common trend over the top
yrange = range(testw[,c(3:27)],na.rm=TRUE)

# dev.new(width=173/25.4,height=173/25.4)
# par(mfrow=c(4,2),oma=c(.5,1.5,.2,.2),mar=c(1,1,0,0) )
res <- 6
sitePlot <- c(1:8)
name_figure <- 'secchi_trend.png'
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

nf <- layout(matrix( c(1:8),nrow=4,ncol=2,byrow=T),  TRUE) 
layout.show(nf)
par(mar=c(0.5,0.5,0.5,0.5),oma=c(3,3,0,1),mai=c(0.1,0.1,0.1,0) )	


for(i in 1:length(clusters)){
  #Code to run DFA model for each cluster
   temp = data.matrix(testw[which(testw$cluster==clusters[i]),3:27])
  #Plotting all the data
  for (n in 1:nrow(temp)) {
    if(i == 1 | i == 3 | i==5) {
      if (n==1) {
        plot(x=c(1987:2011),y = temp[n,],ylim=yrange,type="l",col="grey",xlab="",xaxt="n",yaxt="n")
        title(main=paste0("Cluster ",clusters[i]),line=-1)
        axis(side=1,labels = FALSE,tck=-0.02)
        axis(side=2,labels = FALSE,tck=-0.02)
        axis(side=2,line = -.75,lwd=0)
        box(lwd=1)
        #if (i==3) mtext(side=2,"Standardized Secchi Depth",line=1.2,adj=4.5)
      }
    } else if (i == 7 & n == 1) {
        plot(x=c(1987:2011),y = temp[n,],ylim=yrange,type="l",col="grey",xlab="",xaxt="n",yaxt="n")
        title(main=paste0("Cluster ",clusters[i]),line=-1)
        axis(side=1,labels = FALSE,tck=-0.02)
        axis(side=2,labels = FALSE,tck=-0.02)
        axis(side=1,line = -.75,lwd=0)
        axis(side=2,line = -.75,lwd=0)
        box(lwd=1)
    } else if (i == 8 & n == 1)  {
      plot(x=c(1987:2011),y = temp[n,],ylim=yrange,type="l",col="grey",xlab="",xaxt="n",yaxt="n")
      title(main=paste0("Cluster ",clusters[i]),line=-1)
      axis(side=1,labels = FALSE,tck=-0.02)
      axis(side=2,labels = FALSE,tck=-0.02)
      axis(side=1,line = -.75,lwd=0)
      box(lwd=1)  
      } else if (n==1) {
      plot(x=c(1987:2011),y = temp[n,],ylim=yrange,type="l",col="grey",xlab="",xaxt="n",yaxt="n")
      title(main=paste0("Cluster ",clusters[i]),line=-1)
      axis(side=1,labels = FALSE,tck=-0.02)
      axis(side=2,labels = FALSE,tck=-0.02)
      }
    if(n > 1) lines(x=c(1987:2011),y = temp[n,],col="grey")
  }
   # Add random walk model fit
  points(c(1987:2011), meanTrend[[i]], lty=1, lwd=2, type='l') 
  points(c(1987:2011), ciTrend[[i]][1,], lty=2, lwd=2, type='l') 
  points(c(1987:2011), ciTrend[[i]][2,], lty=2, lwd=2, type='l') 
  
  mtext("Year", line = 1, side = 1, cex = 1, outer=T)
  mtext("Standardized Secchi Depth", line = 1, side = 2, cex = 1, outer=T)

}

par(def.par)
dev.off()
