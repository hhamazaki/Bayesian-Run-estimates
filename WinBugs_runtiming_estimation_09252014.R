################################################################
#  WinBUGS missing passage estimator program   04/09/2014      #
#  Main change is this assumes normal distribution error #
#  This program is designed to estimate missing run passage
#  using Hierarchical Bayesian method
#  The model specifies daily variance 
#  This is a R-code with R2WinBUGS Package and call for WinBUGS
#  You need to: 
#   1) Install R, 
#   2) Install R2WinBUGS and SDMTools Packages 
#   3) Install WinBUGS 
#   4) Prepare csv file with csv format
################################################################

################################################################
#  Model assumption:  Normal distribution of Escapement Counts 
################################################################

################################################################
# 1.0  Set working environment                                 #
################################################################
rm(list=ls(all=TRUE))
library(coda)
library(R2WinBUGS)
library(SDMTools)


# Set working directory: Where your file is located 
working_dir='C:/Projects/Kuskokwim_River/ESC/Test/'
bug_working_dir ='C:/Projects/Kuskokwim_River/ESC/Test/'
# Input file is the one you created 
input_file = '2016_MFG_Chum.csv'
# results file is the one you created 
results_file = '2016_MFG_Chum_bug_results.csv'
# Specify directory where WinBUGS program is located
bugs_dir = 'c:/WinBUGS14/'  

# Specify x-axis of graph and title.
start_day = 'Days from 6/15'
gtitle = 'MFG Sockeye'


################################################################
# 1.1   Import csv data to R                                   #
################################################################
# read table and put into data escape
escape <- read.csv(paste(working_dir,input_file,sep=''),header = TRUE)
# Find dimention of the dataset
dim(escape)
# x is days vector
x <- as.vector(escape$Day)
x2 <- as.Date(escape$Date,"%d-%b")
# nyrs is the number of years data have 
nyrs <- dim(escape)[2] -2
# ndays is the number of days  
ndays <- dim(escape)[1]
# Set an empty matrix y 


################################################################
# 1.2   Plot run timing                                        #
################################################################
#windows(h=6,w=12,record=TRUE)
#par(mfrow=c(3,4),mai=c(.6,.5,.1,.1))
#maxdat<-max(escape[,3:nyrs],na.rm=T)
#for(i in 1:nyrs){
#	plot(x2,escape[,i+2],col='blue', xlab= start_day, ylab='escapement')
#	legend('topright',bty ='n', c(gtitle,paste(substr(names(escape)[i+2],2,5))))
#}


################################################################
# 1.3   Convert data  to log scale                             #
################################################################

y <- matrix(0,nyrs,ndays)
# Convert escapement data to log(x+1) scale 
for (i in 1:ndays){
  for (j in 1:nyrs){
#     y[j,i]<-  log(escape[i,j+2]+0.01)
     y[j,i]<-  escape[i,j+2]+0.01  
     }
    }

################################################################
#  2.0: Create WinBUGS/OpenBUGS Model code					   #
################################################################
mod<-function() {
for(j in 1:nyrs) {
    for(i in 1:ndays){
     y[j,i] ~ dnorm(theta[j,i], tausq[j])
#    y[j,i] ~ dlnorm(log.theta[j,i], tausq[j])
#    y[j,i] ~ dpois(theta[j,i])	
#	log.theta[j,i] <- log(theta[j,i])
# Assume that run timing distribution takes log normal distribution 
    theta[j,i] <- a[j]*exp(-0.5*pow(log(x[i]/mu[j])/b[j],2))
# Assume that run timing distribution takes Extreme value distribution 
#   theta[j,i] <- a[j]*exp(-exp(-(x[i]-mu[j])/b[j])-(x[i]-mu[j])/b[j]+1)
# Assume that run timing distribution takes log-logistic distribution 
#   theta[j,i] <- (a[j]*(b[j]/mu[j])*pow((x[i]/mu[j]),b[j]-1))/pow(1+pow((x[i]/mu[j]),b[j]),2)   
 }
}
# a[] indicates the maximum height (amplitude) of the function a>0
# mu[] indicates the function peaks when x = mu mu>0 : Peak timing
# b[] indicates peak width of the function b>0 standard deviation

# Priors
for(i in 1:nyrs) {
# Normal distribution Positive only 
#  a: is independent not hierarhical 
   a[i] ~ dnorm(0,0.00001)%_%I(0,)
   b[i] ~ dnorm(b0,b0.prec)%_%I(0,)
   r[i] ~ dnorm(0,0.00001)%_%I(0,)
   mu[i] ~ dnorm(mu0,mu0.prec)%_%I(0,)   
     }  
		b0 ~ dnorm(0.5,0.001)%_%I(0,)
		mu0 ~ dnorm(50,0.001)%_%I(0,)
		b0.prec <-1/b0.ssq 
    b0.ssq <- b0.sigma*b0.sigma
    b0.sigma ~ dunif(0,100)  
		mu0.prec <-1/mu0.ssq 
    mu0.ssq <- mu0.sigma*mu0.sigma
    mu0.sigma ~ dunif(0,10)  
    
## This assumes that variance of each year is independent.     
 for(i in 1:nyrs) {    
	   tausq[i] <- pow(sigma[i],-2)
#		sigma[i] ~ dunif(0,100) 
		sigma[i] ~ dunif(0,100) 		
      }

# Backestimate escapement 
for(j in 1:nyrs){
    for(i in 1:ndays){ 
#    y2[j,i] <- y[j,i]
   }
  }
}

#write the model to a text file to be called by WinBUGS
bugfile<-paste(bug_working_dir,'/model.txt',sep='')
write.model(mod,bugfile)

################################################################
#  3.0: Create WinBUGS/OpenBUGS data file                      #  
#  Initial value an change based on your data and model       #   
################################################################

#  Create WinBUGS/OpenBUGS data file (datnew)
datnew<-list(nyrs=nyrs, ndays=ndays, x=x, y=y)

################################################################
# 3.1   Calculate initial value used for Bayesian method       #
################################################################
# zs is empirical standard deviation for each day
zs <- rep(0,nyrs)
for (i in 1:nyrs){
     zs[i]<-  sd(y[,i],na.rm=TRUE)
    }
# add 0 when zs is NA 
zs[is.na(zs)] <- 0  

# zs is empirical b for each year  
zb <- rep(0,nyrs)
for (i in 1:nyrs){
     zb[i]<-  wt.sd((x),y[i,])
    }

# z is bempirical mu for each year 
# za is empirical a for each year 
z <- rep(0,nyrs)
za <- rep(0,nyrs)
  for (j in 1:nyrs){
     za[j]<-  (max(escape[,j+2],na.rm=TRUE))    
     z[j]<-  sum(x*escape[,j+2],na.rm=TRUE)/sum(escape[,j+2],na.rm=TRUE)
    }
    
# Calculate Initial Values 
a <- za
mu <- z
b <- zb
sigma <- ifelse(zs>10,4,zs)
a0 <- median(a)
b0 <- mean(b)
mu0 <- median(mu)
a0.sigma <- 5
b0.sigma <- 1
mu0.sigma <- 5

################################################################
# 3.2   Save Initial value as list file                        #
################################################################
inits <- list(
a <- (a), 
mu <- mu,
b <- b,
sigma <- ifelse(sigma> 100, 50,sigma),
#a0 <- a0,
b0 <- b0,
mu0 <- mu0,
#a0.sigma <- a0.sigma,
b0.sigma <- b0.sigma ,
mu0.sigma <- mu0.sigma
)


################################################################
#  4.0: Run WinBUGS                                            #  
################################################################

# pass the initials to WinBUGS
inits<-list(inits)
#Define the parameters (nodes) of interest 
parameters <- c('a','b','mu','y','theta') 
	
#Run WinBUGS using the bugs() function
starttime=Sys.time()
sim <- bugs(data=datnew, inits=inits, parameters.to.save=parameters, model.file='model.txt',n.chains=1, 
	n.iter=5000,n.burnin=1000,n.thin=2,debug = TRUE, codaPkg=FALSE,  DIC=TRUE, bugs.directory = bugs_dir, working.directory=bug_working_dir)
   Sys.time()-starttime
print(sim$summary[,c(1:3,5,7)])
#write the data frame to a text file and save
bug_summary <- sim$summary
#write.csv(bug_summary,paste(working_dir,results_file,sep=''),row.names=T)    

#The bugs object is automatically imported back into R from bugs once the model is finished running.
post.samp <- sim$sims.array 
#create a data frame of posterior samples
post.samp <- as.data.frame(apply(post.samp,3,function(x) as.numeric(x) ))  
#create a data frame of posterior samples
nvars<-dim(post.samp)[2]
nsamps<-dim(post.samp)[1]


###############################################################
#  Graphics 
###############################################################

###############################################################
#  Plot estimates of each parameter: Optional 
###############################################################
#int<-25
#for(j in seq(1,nvars,int)){
#windows(h=6,w=12)
#par(mfrow=c(5,10),mai=c(0.2,0.2,0.2,0.2))

# Trace plost for Chain1
#for(i in 0:(int-1)){
#	mindat<-min(post.samp[,i+j])
#	maxdat<-max(post.samp[,i+j])
# plot density 
#	plot(density(post.samp[1:(nsamps),i+j]),col='blue',main=names(post.samp)[i+j],xlim=c(mindat,maxdat))
#	lines(density(post.samp[1:nsamps,i+j]),col='red')
# plot trace plot
#	plot(post.samp[1:(nsamps),i+j],col='blue',main=names(post.samp)[i+j],ylim=c(mindat,maxdat),type='l')
#	lines(post.samp[1:nsamps,i+j],col='red')
#}
#}


###############################################################
#  5.0 Extract data											  #
###############################################################

###############################################################
#  5.1 Extract predicted data 								
###############################################################
#create data with expected value
t1 <-(post.samp[,substr(names(post.samp),1,5)=='theta'])
#extract names:  this extracts only bracket part of the theta
t1.name <- substr(names(t1),6,15)
#change names:   this changes column name only bracket part of the theta
colnames(t1) <- t1.name
#y2 produces only data with missing passge dates:   
y2 <-(post.samp[,substr(names(post.samp),1,1)=='y'])
#extract names:  this extracts only bracket part of column name 
y2.name <- substr(names(y2),2,15)
#Extract predicted data based on y2's blacket 
t2 <-t1[,y2.name]
#extract names: this extracts only First part of bracket (year id)
tyear <- substr(names(y2),3,4)

###############################################################
#  5.2 Calculate median, 95% CI of missing dates           
###############################################################
y2med <- apply(t2[,1:dim(t2)[2]],2,median)
y2low <- apply(t2[,1:dim(t2)[2]],2,function(x) quantile(x, 0.025))
y2up <- apply(t2[,1:dim(t2)[2]],2,function(x) quantile(x, 0.975))


###############################################################
#  5.3 Calculate missing value for each year total:              
###############################################################
# Change names to year id 
colnames(t2) <- tyear
# Combine columns based on year id 
t3 <- as.data.frame(sapply(unique(colnames(t2)), function(x) rowSums(t2[, colnames(t2) == x, drop = FALSE])))
# Calculate median, 95% CI of missing dates  
# ym, ylow, and yup are median, 95% CI of annual passage of missing dates
ym <- round(apply(t3[,1:dim(t3)[2]],2,median),0)
ylow <- round(apply(t3[,1:dim(t3)[2]],2,function(x) quantile(x, 0.025)),0)
yup <- round(apply(t3[,1:dim(t3)[2]],2,function(x) quantile(x, 0.975)),0)


###############################################################
#  5.4 Create              
###############################################################

# 1.0: Extract years from ym
tname <- names(ym)
# 1.1: change 1, to 1, 2, to 2.... for all years
for(i in 1:length(tname)){
    if(substr(tname[i],2,2)== ','){
    tname[i] <- substr(tname[i],1,1)
    }
    }
# Chage name to numberic      
tname <- as.numeric(tname)
# Create vector that will include missing years
ym2 <- vector('numeric',nyrs)
yl2 <- vector('numeric',nyrs)
yu2 <- vector('numeric',nyrs)
t <- 1    
for(i in 1:nyrs){
# if year entry is the same as colum years, then put numbers to new vector
    if(i == tname[t]){ 
        ym2[i] <- ym[t] 
        yu2[i] <- yup[t]
        yl2[i] <- ylow[t]
# and move to next entry        
        t <- t+1
# if year entry is not the same as colum years, then put NA: Missing years         
        } else {  
        ym2[i] <- NA
        yu2[i] <- NA
        yl2[i] <- NA
        }
    }


###############################################################
#  Plot Run timing:  
###############################################################

Modelesc <- matrix(0,ndays,nyrs)
am <- bug_summary[substr(row.names(bug_summary),1,2)=='a[',5]
bm <- bug_summary[substr(row.names(bug_summary),1,2)=='b[',5]
mum <- bug_summary[substr(row.names(bug_summary),1,2)=='mu',5]

for (i in 1:ndays){
  for (j in 1:nyrs){
#     Expected log normal run timing   
      Modelesc[i,j]<- (am[j]*exp(-0.5*(log(x[i]/mum[j])/bm[j])^2))
#     Expected Extreme value normal run timing   
#      Modelesc[i,j] <-am[j]*exp(-exp(-(x[i]-mum[j])/bm[j])-(x[i]-mum[j])/bm[j]+1)
#     Expected log logistic run timing
#     Modelesc[i,j] <- (am[j]*(bm[j]/mum[j])*((x[i]/mum[j])^(bm[j]-1))/(1+(x[i]/mum[j])^bm[j])^2)-1   
     }
    }

est.esc <- matrix(0,ndays,nyrs)
for (i in 1:ndays){
  for (j in 1:nyrs){
     est.esc[i,j]<- ifelse(is.na(escape[i,j+2]), Modelesc[i,j],escape[i,j+2]) 
     }
    }
colSums(est.esc)

# y2m is a matrix of median passage estimates of all missing passage
y2m <- matrix(0,ndays,nyrs)
# y2u is a matrix of upper 95% CI passage estimates of all missing passage    
y2u <- matrix(0,ndays,nyrs)
# y2l is a matrix of upper 95% CI passage estimates of all missing passage    
y2l <- matrix(0,ndays,nyrs)

# Enter data to matrix
t <- 1
for (i in 1:nyrs){
  for (j in 1:ndays){
    if (is.na(escape[j,i+2])){
        y2m[j,i] <- y2med[t]
        y2u[j,i] <- y2up[t]
        y2l[j,i] <- y2low[t]
        t <- t+1
        } else {
        y2m[j,i] <- NA
        y2u[j,i] <- NA
        y2l[j,i] <- NA
        }
    }
  }

# plot graph
windows(h=6,w=12,record=TRUE)
par(mfrow=c(3,3),mai=c(.6,.5,.1,.1))
for(i in 1:nyrs){
	maxdat<-max(max(escape[,i+2],na.rm=T),max(y2u[,i],na.rm=T))
# plot observed passage
	plot(x2,escape[,i+2],col='blue', ylim = c(0,maxdat), xlab= start_day, ylab='escapement')
	legend('topright',bty ='n', c(gtitle,paste(substr(names(escape)[i+2],2,5)),paste('missing counts',ym2[i]),paste('95%CI ',yl2[i],' - ',yu2[i])))
# plot modeled run timing
	lines(x2,Modelesc[,i],col='red')
# plot 95% CI lines    
    arrows(x2,y0=y2u[,i],y1=y2l[,i],code=0)
# plot median passage stimate    
    points(x2,y2m[,i], pch=21, col='black',bg='white')
}

###############################################################
# 6.0 Data outputs        
###############################################################
###############################################################
# 6.1 Annual Observed and Esimated Missing Count         
###############################################################
# Extract Total observed counts 
esc.ob <- colSums(escape[,3:(nyrs+2)],na.rm = TRUE)
# Extract Year column 
year <- substr(names(esc.ob),2,5)
# Create data.frame with: year, observed, estimated, lower 95%CI, and upper 95%CI 
esc.sum <- data.frame(year,esc.ob, ym2,yl2,yu2)
# Calculate fraction of missing passage
esc.sum$p.est <- ym2/(esc.ob+ym2)
# Rename column name
names(esc.sum) <- c('year','observed','estimated','low.95%CI','upper.95%CI','% esimated')
# Write CSV file to working directory
write.csv(esc.sum,paste(working_dir,results_file,sep=''),row.names = FALSE)    

###############################################################
# 6.2  Daily Observed and Esimated Missing Count by year     
###############################################################
for(i in 1:nyrs){
# Create a data.frame with observed escapement, estimated, lower 95% CI, and upper 95% CI
esc.daily <- data.frame(escape$Date,escape[,2+i],round(y2m[,i],0),round(y2l[,i],0),round(y2u[,i],))
# Rename the column name 
names(esc.daily) <- c('Date','observed','estimated','low.95%CI','upper.95%CI')
# Write CSV file to working directory
write.csv(esc.daily,paste(working_dir,paste(gtitle,'_',year[i],'.csv',sep=''),sep=''),row.names = FALSE)    
}

