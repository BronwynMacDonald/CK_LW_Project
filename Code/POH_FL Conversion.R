# Bronwyn MacDonald
# May 9 2017

# Convert POH to Fork Length
# For project that is fitting von Bertalanffy Model to Chinook data
# For evaluation of size limits in ocean fisheries 


# Code that installs packages if need be - taken from Forecast.R
usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dep = TRUE)
    require(p, character.only = TRUE)
}

# Install packages
usePackage("snow")
usePackage("dclone")
usePackage("R2jags")
usePackage("lattice")
usePackage("superdiag")


# ======== Conversion of POH (spawning ground data) to Fork Length (required for applicability to ocean fisheries) ====== #

# 2008 data from Albion is primarily 4sub1's
POH_FL_data_2008_4_1 <- read.csv("Albion POH_FL Data_2008.csv")
plot(POH_FL_data_2008_4_1, xlim=c(300,1000), ylim=c(400,1200))

# 1981 data from Albion is primarily 5sub2's
POH_FL_data_1981_5_2 <- read.csv("Albion POH_FL Data_1981.csv")
points(POH_FL_data_1981_5_2, col=2)

# frequentist regression
freq.FL_POH_2008 <- lm(Fork.Length~POH.Length, data=POH_FL_data_2008_4_1)
summary(freq.FL_POH_2008)

# frequentist regression
freq.FL_POH_1981 <- lm(Fork.Length~POH.Length, data=POH_FL_data_1981_5_2)
summary(freq.FL_POH_1981)

# Bayesian regression
Regression <- function(POH_FL_data, new_data){
    n=length(POH_FL_data$POH.Length)
    attach(POH_FL_data)
    
 
    # Add prediction data 
    pred_POH <- new_data
    pred_FL <- rep(NA, length=length(pred_POH))
    
    x <- as.numeric(c(POH.Length, pred_POH))
    y <- as.numeric(c(Fork.Length,pred_FL))
    N <- length(x)
    
 
    # Model Specs
    n.chains=2
    #n.burnin=5000
    n.adapt=250
    n.update=250 
    n.iter=2000
    n.thin=4
    
    # Model
    regression.model <- function(){
      for(i in 1:N){
        # Likelihood function for each data point
        y[i] ~ dnorm(pred_FL[i], tau)
        # Linear regression equation
        pred_FL[i] <- alpha + beta*x[i] 
        
        # Track predicted values
        Rep[i] ~ dnorm(pred_FL[i],tau)
      }
      Track <- Rep[(n+1):N]
      alpha ~ dnorm(0.0,1.0E-4) # Prior for intercept
      beta ~ dnorm(0.0,1.0E-4) # Prior for slope 
      tau ~ dgamma(0.001,0.001) # Prior for uncertainty
    }
    
    # Inputs
    jags.data = list(N=N,x=x,y=y,n=n)
    #jags.data = list("N","x","y","n")
    # Initial Values
    Inits<-function(){
      list(alpha=100,beta=100)
      list(alpha=10,beta=0)
    }
    
    parameters = c("alpha","beta","tau","Track")  
    #parameters = c("alpha","beta","tau","Rep")  
    
    cl <- makeCluster(n.chains, "SOCK")
     
    # Run Model
    jags.regression = jags.parfit(cl,data=jags.data,inits=Inits,params=parameters,model=regression.model,
                                  n.chains=n.chains,n.iter=n.iter,n.adapt=n.adapt, n.update=n.update,n.thin=n.thin,DIC=T)		
    
    stopCluster(cl)
    detach(POH_FL_data)
    return(jags.regression)
    
} # end function that runs jags model



## end working May 11 2017

# Diagnostic Plots
Run_Plots <- function(jags.regression){
    regression.mcmc <- as.mcmc(jags.regression)
    summary(regression.mcmc)
    traceplot(jags.regression)
    xyplot(regression.mcmc, layout=c(2,2), aspect="fill")
    densityplot(regression.mcmc, layout=c(2,2), aspect="fill")
    gelman.plot(regression.mcmc)
    geweke.diag(regression.mcmc)
    raftery.diag(regression.mcmc)
    heidel.diag(regression.mcmc)
    
    # Run them all
    superdiag(regression.mcmc, burnin=5000)
}



#Run_Plots(regression_2008)

# ----MCMC Object Results ----



# 
# 
# coda.out_1 = read.bugs(bugs.model_1_coda)
# results_1 = summary(coda.out_1)$statistics
# quant_1 = summary(coda.out_1)$quantiles
# 
# plot(Julian_Day,Upper_Fraser)

